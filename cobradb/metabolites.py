from typing import Tuple, Dict, Optional

from sqlalchemy import select
from sqlalchemy.orm import Bundle
from cobradb.data_sources import DATA_SOURCE_NAMES
from cobradb.models import (
    Annotation,
    AnnotationLink,
    AnnotationProperty,
    Component,
    ComponentAnnotationMapping,
    ComponentIDMapping,
    ComponentReferenceMapping,
    DataSource,
    InChI,
    Model,
    ReferenceCompound,
    ReferenceCompoundAnnotationMapping,
    UniversalComponent,
    UniversalComponentReferenceMapping,
)
from libchebipy import ChebiEntity
import re
from pprint import pprint
import pandas as pd

FORMULA_PATTERN = re.compile(r"(([A-Z][a-z]?)([0-9])*)+")
FORMULA_PATTERN_SINGLE = re.compile(r"([A-Z][a-z]?)([0-9])*")


def fix_explicit_formula(formula):
    if not isinstance(formula, str):
        return False, None
    m = FORMULA_PATTERN.fullmatch(formula)
    if m is None:
        return False, None

    new_formula = ""
    is_original_formula = False
    for m in FORMULA_PATTERN_SINGLE.finditer(formula):
        atom = m[1]
        mult = m[2]
        if atom == "R":
            return False, None
        if mult is not None and int(mult) == "1":
            new_formula = new_formula + atom
            is_original_formula = False
        else:
            new_formula = new_formula + m[0]
    return is_original_formula, new_formula


def _formula_to_dict(formula):
    d = {}
    for m in FORMULA_PATTERN_SINGLE.finditer(formula):
        atom = m[1]
        mult = m[2]
        if mult is None:
            mult = 1
        mult = int(mult)
        d[atom] = mult
    return d


def are_explicit_formulae_equivalent(formula1, formula2):
    if formula1 is None or formula2 is None:
        return False
    d1 = _formula_to_dict(formula1)
    d2 = _formula_to_dict(formula2)
    return d1 == d2


def load_default_chebi_mapping():
    try:
        df = pd.read_csv("chebi_pH7_3_mapping.tsv", sep="\t", index_col=0, header=0)
        df.index = "CHEBI:" + df.index.astype(str)
        df["CHEBI_PH7_3"] = "CHEBI:" + df["CHEBI_PH7_3"].astype(str)
        return df["CHEBI_PH7_3"].to_dict()
    except:
        return {}


DEFAULT_CHEBI_MAPPING = load_default_chebi_mapping()


# class InChI:
#     def __init__(self, inchi_str):
#         self.string = inchi_str
#         self.formula, self.layers = InChI.parse_into_layers(inchi_str)
#
#     @staticmethod
#     def parse_into_layers(inchi_str):
#         l = inchi_str.split("/")
#         if (prefix := l.pop(0)) != "InChI=1S":
#             raise ValueError(f"InChI prefix not valid: '{prefix}'")
#         formula = l[0]
#         layers = {x[0]: x[1:] for x in l[1:]}
#         return formula, layers
#
#     def has_layer(self, layer_id):
#         return layer_id in self.layers
#

MAIN_RELATIONS = [
    "is_conjugate_acid_of",
    "is_conjugate_base_of",
    "is_tautomer_of",
]


def get_related_chebis(
    chebi: str, related_chebis: Optional[Dict[str, Optional[Tuple[str, str]]]] = None
) -> Dict[str, Optional[Tuple[str, str]]]:
    chebi_entity = ChebiEntity(chebi)

    if related_chebis is None:
        related_chebis = {chebi: None}

    new_chebis = set()

    outgoing_rel = chebi_entity.get_outgoings()
    for rel in outgoing_rel:
        if rel._Relation__typ in MAIN_RELATIONS:
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            if rel_chebi in related_chebis:
                continue
            related_chebis[rel_chebi] = (rel._Relation__typ, chebi)
            new_chebis.add(rel_chebi)

    for new_chebi in new_chebis:
        get_related_chebis(new_chebi, related_chebis)

    return related_chebis


def get_or_create_small_molecule_reference(
    chebi: str, session, cpd_cls=ReferenceCompound
) -> Tuple[bool, ReferenceCompound]:
    chebi_db = session.scalars(
        select(ReferenceCompound).filter(ReferenceCompound.bigg_id == chebi).limit(1)
    ).first()
    if chebi_db:
        return True, chebi_db
    chebi_entity = ChebiEntity(chebi)
    inchi_db = None
    if chebi_inchi := chebi_entity.get_inchi():
        inchi_obj = InChI.from_string(chebi_inchi)
        if inchi_obj is not None:
            inchi_db = session.scalars(
                select(InChI)
                .filter(
                    (InChI.key_major == inchi_obj.key_major)
                    & (InChI.key_minor == inchi_obj.key_minor)
                    & (InChI.key_proton == inchi_obj.key_proton)
                )
                .limit(1)
            ).first()
            if inchi_obj != inchi_db:
                inchi_db = None
            if inchi_db is None:
                inchi_db = inchi_obj
    cpd_dict = dict(
        bigg_id=chebi,
        name=chebi_entity.get_name(),
        html_name=chebi_entity.get_name(),
        charge=str(chebi_entity.get_charge()),
        formula=chebi_entity.get_formula(),
        inchi=inchi_db,
    )
    if cpd_cls is ReferenceCompound:
        cpd_dict["compound_type"] = "small_molecule"

    chebi_db = cpd_cls(**cpd_dict)
    session.add(chebi_db)

    if cpd_cls is ReferenceCompound:
        add_chebi_annotations(chebi_db, chebi_entity, session)

    return False, chebi_db


CHEBI_METABOLITE_PROPERTIES = {
    ("CHEBI", "CHEBI"): "CHEBI",
    ("KEGG COMPOUND accession", "KEGG COMPOUND"): "kegg.compound",
    ("DrugBank accession", "DrugBank"): "drugbank",
    ("KEGG DRUG accession", "KEGG DRUG"): "kegg.drug",
    ("Wikipedia accession", "Wikipedia"): "wikipedia.en",
    ("MetaCyc accession", "MetaCyc"): "metacyc.compound",
    ("HMDB accession", "HMDB"): "hmdb",
}


def add_chebi_annotations(reference_compound_db, chebi_entity, session):
    data_source_dbs = {}
    for (
        data_source_type,
        data_source_source,
    ), data_source_bigg_id in CHEBI_METABOLITE_PROPERTIES.items():
        data_source_db = session.scalars(
            select(DataSource)
            .filter(DataSource.bigg_id == data_source_bigg_id)
            .limit(1)
        ).first()
        if data_source_db is None:
            data_source_db = DataSource(
                bigg_id=data_source_bigg_id,
                name=DATA_SOURCE_NAMES[data_source_bigg_id],
                url_prefix=f"https://identifiers.org/",
            )
            session.add(data_source_db)
        data_source_dbs[data_source_bigg_id] = data_source_db

    chebi = reference_compound_db.bigg_id
    annotation_db = session.scalars(
        select(Annotation).filter(Annotation.bigg_id == chebi).limit(1)
    ).first()
    if not annotation_db:
        annotation_db = Annotation(
            bigg_id=chebi,
            type="chebi",
            default_data_source=data_source_dbs["CHEBI"],
        )
        session.add(annotation_db)

    for property, func in {
        "definition": chebi_entity.get_definition,
        "star": chebi_entity.get_star,
        "mass": chebi_entity.get_mass,
        "smiles": chebi_entity.get_smiles,
    }.items():
        prop_val = func()
        if isinstance(prop_val, str):
            prop_val = prop_val.strip()
            if len(prop_val) == 0:
                prop_val = None
        if prop_val is None:
            continue
        if isinstance(prop_val, int) and property.startswith("is_"):
            prop_val = bool(prop_val)
        prop_db = AnnotationProperty(key=property)
        prop_db.value = prop_val
        annotation_db.properties.append(prop_db)

    for name in chebi_entity.get_names():
        prop_db = AnnotationProperty(key="name")
        prop_db.value = name.get_name()
        annotation_db.properties.append(prop_db)

    alias_db = AnnotationLink(
        identifier=chebi,
        data_source=data_source_dbs["CHEBI"],
    )
    annotation_db.links.append(alias_db)

    for database_accession in chebi_entity.get_database_accessions():
        data_source_bigg_id = CHEBI_METABOLITE_PROPERTIES.get(
            (database_accession.get_type(), database_accession.get_source())
        )
        if data_source_bigg_id is None:
            continue
        prop_val = database_accession.get_accession_number()
        if isinstance(prop_val, str):
            prop_val = prop_val.strip()
            if len(prop_val) == 0:
                prop_val = None
        if prop_val is None:
            continue
        if isinstance(prop_val, str):
            prop_val = {prop_val}

        for val in prop_val:
            alias_db = AnnotationLink(
                identifier=val, data_source=data_source_dbs[data_source_bigg_id]
            )
            annotation_db.links.append(alias_db)

    annotation_mapping = ReferenceCompoundAnnotationMapping(
        reference_compound=reference_compound_db,
    )
    annotation_db.reference_compound_mappings.append(annotation_mapping)


BIGG_ID_PATTERN = re.compile(f"[a-zA-Z0-9][a-zA-Z0-9_]*")


def is_bigg_id_valid(proposed_bigg_id):
    return bool(BIGG_ID_PATTERN.fullmatch(proposed_bigg_id) is not None)


def create_metabolite(
    proposed_bigg_id: str, input_chebi: str, session, assure_present=None
):
    if not is_bigg_id_valid(proposed_bigg_id):
        return {"status": "error", "message": "invalid bigg id"}

    old_id_db = session.scalars(
        select(ComponentIDMapping)
        .filter(ComponentIDMapping.old_bigg_id == proposed_bigg_id)
        .limit(1)
    ).first()
    if old_id_db:
        return {
            "status": "error",
            "message": "bigg id is a deprecated id",
            "new_id": old_id_db.new_id,
        }

    universal_component_db = session.scalars(
        select(UniversalComponent)
        .filter(UniversalComponent.bigg_id == proposed_bigg_id)
        .limit(1)
    ).first()
    if universal_component_db:
        return {"status": "error", "message": "bigg id already exists"}

    all_chebis = get_related_chebis(input_chebi)
    all_chebis = list(all_chebis.keys())

    component_ref_mapping_db = session.execute(
        select(
            ComponentReferenceMapping,
            Bundle("ref_comp", ReferenceCompound.bigg_id),
            Bundle("uni_comp", UniversalComponent.bigg_id),
        )
        .join(ComponentReferenceMapping.reference_compound)
        .join(ComponentReferenceMapping.component)
        .join(Component.universal_component)
        .filter(ReferenceCompound.bigg_id.in_(all_chebis))
        .limit(1)
    ).one_or_none()
    if component_ref_mapping_db:
        return {
            "status": "error",
            "message": "chebi already associated",
            "chebi": component_ref_mapping_db.ref_comp.bigg_id,
            "bigg_id": component_ref_mapping_db.uni_comp.bigg_id,
        }

    references_db = {
        x: get_or_create_small_molecule_reference(x, session)[1] for x in all_chebis
    }

    if assure_present:
        assure_charge, assure_formula = assure_present
        if not any(
            (float(str(ref.charge)) == assure_charge)
            and are_explicit_formulae_equivalent(ref.formula, assure_formula)
            for ref in references_db.values()
        ):
            return {
                "status": "error",
                "message": "charge + formula combination not present",
            }

    default_chebi = None
    for ch in [input_chebi] + all_chebis:
        default_chebi = DEFAULT_CHEBI_MAPPING.get(ch)
        if default_chebi is not None:
            break
    if default_chebi not in references_db:
        default_chebi = None
    if default_chebi is None:
        default_chebi = input_chebi
    default_chebi_entity = ChebiEntity(default_chebi)

    universal_component_db = UniversalComponent(
        bigg_id=proposed_bigg_id,
        name=default_chebi_entity.get_name(),
    )
    session.add(universal_component_db)
    session.commit()

    successfully_added = []

    full_bids_created = {}
    for chebi, reference_db in references_db.items():
        full_bid = f"{proposed_bigg_id}:{reference_db.charge}"
        if full_bid not in full_bids_created:
            charge = reference_db.charge
            try:
                int_charge = int(str(charge))
            except:
                continue
            formula = reference_db.formula
            if not formula or formula == "nan":
                continue

            component_db = Component(
                bigg_id=full_bid,
                name=reference_db.name,
                formula=formula,
                charge=int_charge,
            )
            universal_component_db.components.append(component_db)
            full_bids_created[full_bid] = component_db
        component_db = full_bids_created[full_bid]

        component_ref_mapping_db = ComponentReferenceMapping(
            universal_component=universal_component_db,
            reference_compound=reference_db,
        )
        component_db.reference_mappings.append(component_ref_mapping_db)

        if chebi == default_chebi or default_chebi is None:
            universal_component_ref_mapping_db = session.scalars(
                select(UniversalComponentReferenceMapping)
                .filter(
                    UniversalComponentReferenceMapping.id == universal_component_db.id
                )
                .limit(1)
            ).first()
            if not universal_component_ref_mapping_db:
                universal_component_ref_mapping_db = UniversalComponentReferenceMapping(
                    id=universal_component_db.id, mapping=component_ref_mapping_db
                )
                session.add(universal_component_ref_mapping_db)

        session.commit()
        successfully_added.append((full_bid, reference_db.bigg_id))

    if not full_bids_created:
        session.commit()
        session.delete(universal_component_db)
    session.commit()

    return {"status": "success", "components_added": successfully_added}


def create_model_specific_metabolite(
    bigg_id, model_db_id, charge, name, formula, session
):
    if charge is None:
        charge = 0
    model_db = session.get(Model, model_db_id)
    new_universal_id = f"__{model_db.bigg_id}__{bigg_id}"
    new_bigg_id = f"__{model_db.bigg_id}__{bigg_id}:{charge}"

    universal_metabolite_db = session.scalars(
        select(UniversalComponent)
        .filter(UniversalComponent.bigg_id == new_universal_id)
        .limit(1)
    ).first()
    if not universal_metabolite_db:
        universal_metabolite_db = UniversalComponent(
            bigg_id=new_universal_id, name=name, model=model_db
        )
        session.add(universal_metabolite_db)

    metabolite_db = session.scalars(
        select(Component).filter(Component.bigg_id == new_bigg_id).limit(1)
    ).first()
    if not metabolite_db:
        metabolite_db = Component(
            bigg_id=new_bigg_id,
            universal_component=universal_metabolite_db,
            name=name,
            formula=formula,
            charge=charge,
            model=model_db,
        )
        session.add(metabolite_db)
    return universal_metabolite_db, metabolite_db
