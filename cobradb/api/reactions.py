from typing import Dict, List, Optional, Tuple
from sqlalchemy import select
from sqlalchemy.orm import Session

from cobradb.data_sources import get_data_source_id
from cobradb.models import (
    Annotation,
    AnnotationLink,
    Compartment,
    ReferenceCompound,
    ReferenceReaction,
    ReferenceReactionAnnotationMapping,
    ReferenceReactionParticipant,
)
from cobradb.api import utils


def _create_reference_reaction(
    session: Session,
    bigg_id: str,
    equation: str,
    participants: List[Tuple[int, float, str, Compartment]],
    annotations: Optional[Dict[str, List[str]]],
) -> Optional[ReferenceReaction]:
    reaction_db = ReferenceReaction(
        bigg_id=bigg_id,
        equation=equation,
    )
    for side_n, coefficient, compound_id, compartment in participants:
        compound_db = utils.get_object_by_bigg_id(
            session, compound_id, ReferenceCompound
        )
        participant_db = ReferenceReactionParticipant(
            compound=compound_db,
            side="L" if side_n < 0 else "R",
            coefficient=coefficient,
            compartment=compartment,
        )
        reaction_db.reaction_participants.append(participant_db)
    if len(reaction_db.reaction_participants) > 0:
        reaction_db.update_hash()
        session.add(reaction_db)

        add_reaction_annotations(session, reaction_db, annotations)
        return reaction_db
    return None


def add_reference_exchange_reactions(
    session: Session, compound_db: ReferenceCompound
) -> Optional[ReferenceReaction]:
    reaction_bigg_id = f"EX:{compound_db.bigg_id}"
    reference_reaction_db = session.scalars(
        select(ReferenceReaction)
        .filter(ReferenceReaction.bigg_id == reaction_bigg_id)
        .limit(1)
    ).first()
    if reference_reaction_db is not None:
        return
    compound_name = compound_db.name if compound_db.name else compound_db.id
    reference_reaction_db = ReferenceReaction(
        bigg_id=reaction_bigg_id,
        name=f"Exchange of {compound_name}.",
        equation=f"{compound_name} = ∅",
    )
    lhs_part = ReferenceReactionParticipant(
        compound=compound_db,
        side="L",
        coefficient="1",
        compartment="0",
    )
    reference_reaction_db.reaction_participants.append(lhs_part)
    reference_reaction_db.update_hash()
    session.add(reference_reaction_db)
    session.commit()
    return reference_reaction_db


def add_reaction_annotations(
    session: Session,
    reaction_db: ReferenceReaction,
    annotations: Optional[Dict[str, List[str]]],
):
    if annotations is None:
        return
    default_data_source_id = get_data_source_id("rhea", session)
    if default_data_source_id is None:
        print("Could not find RHEA data source.")
        return
    annotation_db = Annotation(
        bigg_id=reaction_db.bigg_id,
        default_data_source_id=default_data_source_id,
        type="rhea",
    )
    mapping = ReferenceReactionAnnotationMapping(
        reference_reaction=reaction_db,
    )
    annotation_db.reference_reaction_mappings.append(mapping)

    for namespace, identifiers in annotations.items():
        data_source_id = get_data_source_id(namespace, session)
        if not data_source_id:
            print(f"Unknown data source: {namespace}")
            continue
        for identifier in identifiers:
            link = AnnotationLink(
                data_source_id=data_source_id,
                identifier=identifier,
            )
            annotation_db.links.append(link)
    session.add(annotation_db)
