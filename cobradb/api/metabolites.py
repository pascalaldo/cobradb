from typing import Any, List, Optional, Tuple, Type, TypeVar, Union

from sqlalchemy import select
from sqlalchemy.orm import Bundle, Session, subqueryload

from cobradb.data_sources import get_data_source_id
from cobradb.models import (
    Annotation,
    AnnotationLink,
    AnnotationProperty,
    BiGGBase,
    Component,
    ComponentIDMapping,
    ComponentReferenceMapping,
    InChI,
    Model,
    ReferenceCompound,
    ReferenceCompoundAnnotationMapping,
    ReferenceReactivePart,
    ReferenceReactivePartMatrix,
    UniversalComponent,
    UniversalComponentReferenceMapping,
)
from cobradb import chebi
from cobradb.chebi import DEFAULT_CHEBI_MAPPING, ChebiEntity, get_related_chebis
from cobradb.api import utils, reactions
import cobradb.api.bigg_ids as bigg_ids_api

CHEBI_METABOLITE_PROPERTIES = {
    ("CHEBI", "CHEBI"): "CHEBI",
    ("KEGG COMPOUND accession", "KEGG COMPOUND"): "kegg.compound",
    ("DrugBank accession", "DrugBank"): "drugbank",
    ("KEGG DRUG accession", "KEGG DRUG"): "kegg.drug",
    ("Wikipedia accession", "Wikipedia"): "wikipedia.en",
    ("MetaCyc accession", "MetaCyc"): "metacyc.compound",
    ("HMDB accession", "HMDB"): "hmdb",
}


def get_universal_component_by_bigg_id(
    session: Session, universal_bigg_id: str, model_id: Optional[int] = None
) -> Optional[UniversalComponent]:
    """Get the univeral component object from a universal BiGG ID. This
    function maps deprecated IDs to the correct new BiGG IDs.

    Arguments
    ---------

    session: SQLAlchemy session

    universal_bigg_id: The BiGG ID to look up
    """

    universal_component_db = session.scalars(
        select(UniversalComponent)
        .join(UniversalComponent.old_bigg_ids)
        .filter(ComponentIDMapping.old_bigg_id == universal_bigg_id)
        .limit(1)
    ).first()
    if universal_component_db is None:
        if model_id is None:
            model_sel = UniversalComponent.model_id == None
        else:
            model_sel = (UniversalComponent.model_id == None) | (
                UniversalComponent.model_id == model_id
            )
        universal_component_db = session.scalars(
            select(UniversalComponent)
            .filter(UniversalComponent.bigg_id == universal_bigg_id)
            .filter(model_sel)
            .limit(1)
        ).first()
    return universal_component_db


def _create_component(
    session: Session,
    bigg_id: str,
    universal_component_db: UniversalComponent,
    reference_compound_db: Optional[ReferenceCompound] = None,
    charge: Optional[Union[int, str]] = None,
    formula: Optional[str] = None,
    name: Optional[str] = None,
    model_db: Optional[Model] = None,
    require_formula: bool = True,
) -> Optional[Component]:
    if reference_compound_db is not None:
        charge = reference_compound_db.charge
        formula = reference_compound_db.formula
        if name is None:
            name = reference_compound_db.name
    try:
        int_charge = int(str(charge))
    except:
        return None

    if not formula or formula == "nan":
        if require_formula:
            return None
        else:
            formula = ""

    if not ":" in bigg_id:
        bigg_id = bigg_ids_api.create_component_bigg_id(bigg_id, charge=int_charge)

    component_db = Component(
        bigg_id=bigg_id,
        name=name,
        formula=formula,
        charge=int_charge,
        model=model_db,
    )
    universal_component_db.components.append(component_db)
    return component_db


def _create_universal_component(
    session: Session,
    bigg_id: str,
    name: Optional[str] = None,
    model_db: Optional[Model] = None,
):
    universal_component_db = UniversalComponent(
        bigg_id=bigg_id,
        name=name,
        model=model_db,
    )
    session.add(universal_component_db)
    return universal_component_db


def _create_inchi_object(
    session: Session,
    inchi: str,
) -> Optional[InChI]:
    """Create an InChI object from an InChI string if possible, or retrieve the corresponding object from the database.

    Arguments
    ---------

    session: SQLAlchemy session

    inchi: InChI string
    """
    if not inchi:
        return None

    inchi_db = None
    inchi_obj = InChI.from_string(inchi)
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
    return inchi_db


SMRT = TypeVar("SMRT", bound=BiGGBase)


def get_or_create_small_molecule_reference(
    session: Session,
    chebi: str,
    cpd_cls: Type[SMRT] = ReferenceCompound,
) -> Tuple[bool, Optional[SMRT]]:
    """Create a small molecule reference based on a ChEBI ID. Can either create a ReferenceCompound or ReferenceReactivePart.

    Arguments
    ---------
    session: SQLAlchemy session

    chebi: ChEBI ID as a string (CHEBI:00000 format)

    cpd_cls: Class of reference object to create. Default: ReferenceCompound."""
    chebi_db = utils.get_object_by_bigg_id(session, chebi, cpd_cls)
    if chebi_db:
        return True, chebi_db

    chebi_entity = ChebiEntity(chebi)
    if not chebi_entity:
        return False, None

    inchi_db = _create_inchi_object(session, chebi_entity.get_inchi())

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

    if isinstance(
        chebi_db, ReferenceCompound
    ):  # Typechecker is confused when using cpd_cls.
        add_chebi_annotations(session, chebi_db, chebi_entity)

    return False, chebi_db


def add_chebi_annotations(
    session: Session,
    reference_compound_db: ReferenceCompound,
    chebi_entity: ChebiEntity,
) -> None:
    """Add annotations from ChEBI to the reference compound.

    Arguments
    ---------
    session: SQLAlchemy session

    reference_compound_db: ReferenceCompound object to annotate.

    chebi_entity: ChebiEntity object to retrieve annotations from."""
    chebi = reference_compound_db.bigg_id
    annotation_db = utils.get_object_by_bigg_id(
        session,
        chebi,
        Annotation,
        opts=(subqueryload(Annotation.links), subqueryload(Annotation.properties)),
    )
    if not annotation_db:
        annotation_db = Annotation(
            bigg_id=chebi,
            type="chebi",
            default_data_source_id=get_data_source_id("CHEBI", session),
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
        data_source_id=get_data_source_id("CHEBI", session),
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
                identifier=val,
                data_source_id=get_data_source_id(data_source_bigg_id, session),
            )
            annotation_db.links.append(alias_db)

    annotation_mapping = ReferenceCompoundAnnotationMapping(
        reference_compound=reference_compound_db,
    )
    annotation_db.reference_compound_mappings.append(annotation_mapping)


def _create_reactive_part(
    session: Session,
    bigg_id: str,
    name: str,
    formula: str,
    charge: str,
    html_name: Optional[str] = None,
) -> Optional[ReferenceReactivePart]:
    if bigg_id.startswith("CHEBI:"):
        _existed, reactive_part_db = get_or_create_small_molecule_reference(
            session, bigg_id, cpd_cls=ReferenceReactivePart
        )
    else:
        reactive_part_db = ReferenceReactivePart(
            bigg_id=bigg_id,
            name=name,
            html_name=html_name,
            formula=formula,
            charge=charge,
        )
        session.add(reactive_part_db)
    return reactive_part_db


def _create_reference_compound(
    session: Session,
    bigg_id: str,
    name: str,
    compound_type: str,
    formula: Optional[str] = None,
    charge: Optional[str] = None,
    html_name: Optional[str] = None,
    reactive_parts: Optional[List[str]] = None,
):
    if bigg_id.startswith("CHEBI:"):
        _existed, compound_db = get_or_create_small_molecule_reference(session, bigg_id)
    else:
        compound_db = ReferenceCompound(
            bigg_id=bigg_id,
            name=name,
            html_name=html_name,
            formula=formula,
            charge=charge,
            compound_type=compound_type,
        )
        if reactive_parts is not None:
            for reactive_part in reactive_parts:
                reactive_part_db = utils.get_object_by_bigg_id(
                    session, reactive_part, ReferenceReactivePart
                )
                reactive_part_matrix_db = ReferenceReactivePartMatrix(
                    reactive_part=reactive_part_db,
                )
                compound_db.reactive_part_matrix.append(reactive_part_matrix_db)
        session.add(compound_db)
    session.commit()
    if compound_db is not None:
        chebi.add_reference_conversion_reactions(session, compound_db)
        reactions.add_reference_exchange_reactions(session, compound_db)
    return compound_db


def create_metabolite(
    session: Session,
    proposed_bigg_id: str,
    input_chebi: str,
    assure_present: Optional[Tuple[Any, Any]] = None,
):
    if not bigg_ids_api.is_bigg_id_valid(
        proposed_bigg_id, bigg_ids_api.BiGGIDType.UNIVERSAL_COMPONENT
    ):
        return {"status": "error", "message": "invalid bigg id"}

    universal_component_db = get_universal_component_by_bigg_id(
        session, proposed_bigg_id
    )
    if universal_component_db is not None:
        if universal_component_db.bigg_id == proposed_bigg_id:
            return {"status": "error", "message": "bigg id already exists"}
        else:
            return {
                "status": "error",
                "message": "bigg id is a deprecated id",
                "new_id": universal_component_db.bigg_id,
            }

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
        x: ref
        for x in all_chebis
        if (ref := get_or_create_small_molecule_reference(session, x)[1]) is not None
    }

    if assure_present:
        assure_charge, assure_formula = assure_present
        if not any(
            (float(str(ref.charge)) == assure_charge)
            and utils.are_explicit_formulae_equivalent(ref.formula, assure_formula)
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

    universal_component_db = _create_universal_component(
        session,
        bigg_id=proposed_bigg_id,
        name=default_chebi_entity.get_name(),
    )
    session.commit()

    successfully_added = []

    full_bids_created = {}
    for chebi, reference_db in references_db.items():
        full_bid = f"{proposed_bigg_id}:{reference_db.charge}"
        if full_bid not in full_bids_created:
            if (
                component_db := _create_component(
                    session,
                    full_bid,
                    universal_component_db,
                    reference_db,
                )
            ) is None:
                continue
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
    session: Session,
    bigg_id: str,
    model_db: Model,
    charge: Optional[Any],
    name: str,
    formula: str,
) -> Tuple[Optional[UniversalComponent], Optional[Component]]:
    if charge is None:
        charge = 0
    new_universal_id = bigg_ids_api.create_component_bigg_id(
        bigg_id,
        model_bigg_id=model_db.bigg_id,
    )
    new_bigg_id = bigg_ids_api.create_component_bigg_id(
        bigg_id,
        charge=charge,
        model_bigg_id=model_db.bigg_id,
    )

    if (
        universal_component_db := utils.get_object_by_bigg_id(
            session, new_universal_id, UniversalComponent
        )
    ) is None:
        universal_component_db = _create_universal_component(
            session,
            bigg_id=new_universal_id,
            name=name,
            model_db=model_db,
        )

    if (
        component_db := utils.get_object_by_bigg_id(session, new_bigg_id, Component)
    ) is None:
        component_db = _create_component(
            session,
            bigg_id=new_bigg_id,
            universal_component_db=universal_component_db,
            name=name,
            formula=formula,
            charge=charge,
            model_db=model_db,
            require_formula=False,
        )
    return universal_component_db, component_db


def create_component_id_mapping(
    session: Session, old_bigg_id: str, new_bigg_id: str
) -> Optional[ComponentIDMapping]:
    id_mapping_db = session.scalars(
        select(ComponentIDMapping)
        .filter(ComponentIDMapping.old_bigg_id == old_bigg_id)
        .limit(1)
    ).first()
    if id_mapping_db is not None:
        return None

    # Check that we're not making an entry inaccessible.
    if (
        utils.get_object_by_bigg_id(session, old_bigg_id, UniversalComponent)
        is not None
    ):
        return None

    if (
        new_component_db := utils.get_object_by_bigg_id(
            session, new_bigg_id, UniversalComponent
        )
    ) is None:
        return None
    id_mapping_db = ComponentIDMapping(
        old_bigg_id=old_bigg_id,
        new_universal_component=new_component_db,
    )
    session.add(id_mapping_db)
    return id_mapping_db
