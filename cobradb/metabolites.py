from typing import Tuple, Dict, Optional
from cobradb.models import (
    Component,
    ComponentIDMapping,
    ComponentReferenceMapping,
    ReferenceCompound,
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


class InChI:
    def __init__(self, inchi_str):
        self.string = inchi_str
        self.formula, self.layers = InChI.parse_into_layers(inchi_str)

    @staticmethod
    def parse_into_layers(inchi_str):
        l = inchi_str.split("/")
        if (prefix := l.pop(0)) != "InChI=1S":
            raise ValueError(f"InChI prefix not valid: '{prefix}'")
        formula = l[0]
        layers = {x[0]: x[1:] for x in l[1:]}
        return formula, layers

    def has_layer(self, layer_id):
        return layer_id in self.layers


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
    chebi: str, session
) -> Tuple[bool, ReferenceCompound]:
    chebi_db = (
        session.query(ReferenceCompound).filter(ReferenceCompound.id == chebi).first()
    )
    if chebi_db:
        return True, chebi_db
    chebi_entity = ChebiEntity(chebi)
    chebi_db = ReferenceCompound(
        id=chebi,
        name=chebi_entity.get_name(),
        html_name=chebi_entity.get_name(),
        compound_type="small_molecule",
        charge=str(chebi_entity.get_charge()),
        formula=chebi_entity.get_formula(),
    )
    session.add(chebi_db)
    return False, chebi_db


BIGG_ID_PATTERN = re.compile(f"[a-zA-Z0-9][a-zA-Z0-9_]*")


def is_bigg_id_valid(proposed_bigg_id):
    return bool(BIGG_ID_PATTERN.fullmatch(proposed_bigg_id) is not None)


def create_metabolite(proposed_bigg_id: str, input_chebi: str, session):
    if not is_bigg_id_valid(proposed_bigg_id):
        return {"status": "invalid bigg id"}

    old_id_db = (
        session.query(ComponentIDMapping)
        .filter(ComponentIDMapping.old_id == proposed_bigg_id)
        .first()
    )
    if old_id_db:
        return {"status": "bigg id is a deprecated id", "new_id": old_id_db.new_id}

    universal_component_db = (
        session.query(UniversalComponent)
        .filter(UniversalComponent.id == proposed_bigg_id)
        .first()
    )
    if universal_component_db:
        return {"status": "bigg id already exists"}

    all_chebis = get_related_chebis(input_chebi)
    all_chebis = list(all_chebis.keys())

    component_ref_mapping_db = (
        session.query(ComponentReferenceMapping)
        .filter(ComponentReferenceMapping.universal_id in all_chebis)
        .first()
    )
    if component_ref_mapping_db:
        return {
            "status": "chebi already associated",
            "chebi": component_ref_mapping_db.reference_id,
            "bigg_id": component_ref_mapping_db.universal_id,
        }

    references_db = {
        x: get_or_create_small_molecule_reference(x, session)[1] for x in all_chebis
    }

    default_chebi = DEFAULT_CHEBI_MAPPING.get(input_chebi)
    if default_chebi not in references_db:
        default_chebi = None
    if default_chebi is None:
        default_chebi = input_chebi
    default_chebi_entity = ChebiEntity(default_chebi)

    universal_component_db = UniversalComponent(
        id=proposed_bigg_id,
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
                id=full_bid,
                universal_id=universal_component_db.id,
                name=reference_db.name,
                formula=formula,
                charge=int_charge,
            )
            session.add(component_db)
            full_bids_created[full_bid] = component_db
        component_db = full_bids_created[full_bid]

        component_ref_mapping_db = ComponentReferenceMapping(
            component_id=component_db.id,
            universal_id=universal_component_db.id,
            reference_id=reference_db.id,
        )
        session.add(component_ref_mapping_db)

        if chebi == default_chebi or default_chebi is None:
            universal_component_ref_mapping_db = (
                session.query(UniversalComponentReferenceMapping)
                .filter(UniversalComponentReferenceMapping.id == proposed_bigg_id)
                .first()
            )
            if not universal_component_ref_mapping_db:
                session.commit()
                universal_component_ref_mapping_db = UniversalComponentReferenceMapping(
                    id=proposed_bigg_id, mapping_id=component_ref_mapping_db.id
                )
                session.add(universal_component_ref_mapping_db)

        successfully_added.append((full_bid, reference_db.id))

    if not full_bids_created:
        session.commit()
        session.delete(universal_component_db)
    session.commit()

    return {"status": "success", "components_added": successfully_added}


def create_model_specific_metabolite(bigg_id, model_id, charge, name, formula, session):
    if charge is None:
        charge = 0
    new_universal_id = f"__{model_id}__{bigg_id}"
    new_bigg_id = f"__{model_id}__{bigg_id}:{charge}"

    universal_metabolite_db = (
        session.query(UniversalComponent)
        .filter(UniversalComponent.id == new_universal_id)
        .first()
    )
    if not universal_metabolite_db:
        universal_metabolite_db = UniversalComponent(
            id=new_universal_id, name=name, model_specific=True
        )
        session.add(universal_metabolite_db)

    metabolite_db = session.query(Component).filter(Component.id == new_bigg_id).first()
    if not metabolite_db:
        metabolite_db = Component(
            id=new_bigg_id,
            universal_id=new_universal_id,
            name=name,
            formula=formula,
            charge=charge,
            model_specific=True,
        )
        session.add(metabolite_db)
    return universal_metabolite_db, metabolite_db
