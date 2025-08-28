# -*- coding: utf-8 -*-

from cobradb.models import *
from cobradb import settings
from cobradb.util import timing

from sqlalchemy import func
from sqlalchemy import select
import re
import logging
import json

import time

from pprint import pprint
import re
from copy import deepcopy

COMP_IN_OUT = [
    ("c", "p"),
    ("c", "e"),
    ("p", "e"),
    ("n", "c"),
    ("n", "p"),
    ("n", "e"),
]


def load_bigg_id_data(filename):
    with open(filename, "r") as f:
        return json.load(f)


def match_reaction_data_with_db_entry(parsed_participants, db_entry, session):
    stmt = (
        select(
            ReferenceReactionParticipant,
            ReferenceCompound,
        )
        .join(
            ReferenceCompound,
            ReferenceCompound.id == ReferenceReactionParticipant.compound_id,
        )
        .where(ReferenceReactionParticipant.reaction_id == db_entry.id)
    )
    # print(stmt)
    # print(db_entry.id)
    reaction_mappings = [({}, {}, parsed_participants[:], [])]
    # participants_db = session.execute(stmt)
    participants_db = (
        session.query(
            ReferenceReactionParticipant,
            ReferenceCompound,
        )
        .join(
            ReferenceCompound,
            ReferenceCompound.id == ReferenceReactionParticipant.compound_id,
        )
        .filter(ReferenceReactionParticipant.reaction_id == db_entry.id)
    )
    if len(parsed_participants) != participants_db.count():
        return None
    for (
        participant,
        compound,
    ) in participants_db:
        requires_n = "n" in str(participant.coefficient).lower()

        # print(f"Participant: {participant} -> {compound}")
        stmt = (
            select(
                ComponentReferenceMapping,
                UniversalComponent,
                Component,
            )
            .join(
                UniversalComponent,
                UniversalComponent.id == ComponentReferenceMapping.universal_id,
            )
            .where(
                compound.id == ComponentReferenceMapping.reference_id,
            )
            .join(Component, Component.id == ComponentReferenceMapping.component_id)
        )
        # print(stmt)
        # print(compound.id)
        participant_mappings = []
        crmappings_db = session.execute(stmt)
        for crmapping, universal_component, component in crmappings_db:
            if requires_n and crmapping.reference_n is None:
                continue
            coefficient = participant.coefficient
            if requires_n:
                coefficient = coefficient.lower().replace("n", crmapping.reference_n)
            coefficient = float(coefficient)
            # print(f" {crmapping}: {universal_component}")
            ms = [
                p
                for p in parsed_participants
                if p[2] == universal_component.id
                and float(p[0]) == coefficient
                and ((p[4] == component.charge) or p[4] is None)
            ]
            if not ms:
                continue
            for m in ms:
                d_coef, d_side, d_universal_id, d_comp, d_charge = m
                participant_mappings.append(
                    (
                        (participant.side, d_side),
                        (
                            participant.compartment,
                            d_comp,
                        ),
                        m,
                        participant,
                        component,
                        crmapping,
                    ),
                )
        if not participant_mappings:
            reaction_mappings = []
            break
        old_front = reaction_mappings
        reaction_mappings = []
        for front in old_front:
            for participant_mapping in participant_mappings:
                if participant_mapping[2] not in front[2]:
                    continue
                new_item = deepcopy(front)
                for i, att in enumerate(["side", "comp"]):
                    if new_item is None:
                        break
                    if (att_m := participant_mapping[i][0]) in front[i]:
                        if participant_mapping[i][1] != front[i][att_m]:
                            print(
                                f"Inconsistent {att} mapping: prev: {att_m}:{front[i][att_m]}, new: {att_m}{participant_mapping[i][1]}"
                            )
                            new_item = None
                    elif participant_mapping[i][1] in front[i].values():
                        print(
                            f"Inconsistent {att} mapping: value:{participant_mapping[i][1]} already in {front[i]}"
                        )
                        new_item = None
                    else:
                        new_item[i][att_m] = participant_mapping[i][1]
                if new_item is not None:
                    new_item[3].append(
                        (
                            participant_mapping[3],
                            participant_mapping[4],
                            participant_mapping[5],
                        )
                    )
                    new_item[2].remove(participant_mapping[2])
                    reaction_mappings.append(new_item)
    successful_mappings = []
    for m in reaction_mappings:
        # TODO: reference_n should also be consistent
        if m[2]:
            continue
        if "in" in m[1] and "out" in m[1]:
            if not (m[1]["in"], m[1]["out"]) in COMP_IN_OUT:
                print(
                    f"mapping in: {m[1]['in']}, out: {m[1]['out']} not a known in/out pair."
                )
                continue
        successful_mappings.append(m)
    if not successful_mappings:
        print(f"No successful mappings: {reaction_mappings}")
        return None
    elif len(successful_mappings) > 1:
        print(f"Multiple successful mappings: {reaction_mappings}")
        return None
    else:
        return successful_mappings[0]


def parse_reaction_participants(participants, session):
    parsed_participants = []
    for i, side in enumerate("LR"):
        for p in participants[i]:
            full_id = p[1]
            if full_id == "":
                continue
            if ":" in full_id:
                universal_comp_comp_id, charge = full_id.rsplit(":")
            else:
                universal_comp_comp_id, charge = full_id, None
            universal_id, compartment = universal_comp_comp_id.rsplit("_", maxsplit=1)
            id_mapping_db = (
                session.query(ComponentIDMapping)
                .filter(ComponentIDMapping.old_id == universal_id)
                .first()
            )
            if id_mapping_db:
                print(f"Mapping: {id_mapping_db.old_id} -> {id_mapping_db.new_id}")
                universal_id = id_mapping_db.new_id
            parsed_participants.append((p[0], side, universal_id, compartment, charge))
    return parsed_participants


def find_reference_reaction(proposed_ids, parsed_participants, session):
    if proposed_ids is None:
        proposed_ids = []

    db_entries = []
    for rhea_id in proposed_ids:
        db_entry = (
            session.query(ReferenceReaction)
            .filter(ReferenceReaction.id == rhea_id)
            .first()
        )
        if db_entry:
            db_entries.append(db_entry)

    reference_db, mapping = None, None
    if db_entries:
        mappings = []
        for db_entry in db_entries:
            m = match_reaction_data_with_db_entry(
                parsed_participants, db_entry, session
            )
            if m is not None:
                mappings.append((db_entry, m))
        if len(mappings) > 1:
            print(f"No single unique mapping.")
            return None, None
        elif len(mappings) == 0:
            print("No mapping for proposed reference IDs")
        else:
            reference_db, mapping = mappings[0]

    if reference_db is None:
        db_entries = None
        for parsed_participant in parsed_participants:
            new_db_entries = (
                session.query(ReferenceReaction)
                .join(
                    ReferenceReactionParticipant,
                    ReferenceReactionParticipant.reaction_id == ReferenceReaction.id,
                )
                .join(
                    ReferenceCompound,
                    ReferenceCompound.id == ReferenceReactionParticipant.compound_id,
                )
                .join(
                    ComponentReferenceMapping,
                    ComponentReferenceMapping.reference_id == ReferenceCompound.id,
                )
                .filter(ComponentReferenceMapping.universal_id == parsed_participant[2])
                .distinct()
            )
            if db_entries is None:
                db_entries = set(new_db_entries)
            else:
                db_entries &= set(new_db_entries)
            if not db_entries:
                break
        if db_entries is None:
            return None, None
        db_entries = list(db_entries)
        # print(db_entries)
        for db_entry in db_entries:
            # print(f"Matching {db_entry.id}: {db_entry.equation}")
            m = match_reaction_data_with_db_entry(
                parsed_participants, db_entry, session
            )
            # print(m)
            if m is not None:
                reference_db, mapping = db_entry, m
                break

    return reference_db, mapping


@timing
def push_reactions(data, session):
    for n, (bigg_id, reaction_data) in enumerate(data.items()):
        # universal_reaction_db = (
        #     session.query(UniversalReaction)
        #     .filter(UniversalReaction.id == bigg_id)
        #     .first()
        # )
        # if universal_reaction_db is not None:
        #     continue
        print("###")
        print(f"{bigg_id}: {reaction_data['name']}")
        print(f" RHEA: {reaction_data.get('rhea')}")
        pprint(reaction_data["participants"])

        t0 = time.time()
        reaction_data["parsed_participants"] = parse_reaction_participants(
            reaction_data["participants"], session
        )
        print(reaction_data["parsed_participants"])
        print(f"TIME: PARSE: {(t1 := time.time()) - t0}")
        t0 = t1

        reference_db, m = find_reference_reaction(
            proposed_ids=reaction_data.get("rhea", [])
            + reaction_data.get("reference_ids", []),
            parsed_participants=reaction_data["parsed_participants"],
            session=session,
        )

        universal_reaction_db = (
            session.query(UniversalReaction)
            .filter(UniversalReaction.id == bigg_id)
            .first()
        )
        # if not universal_reaction_db:
        #     universal_reaction_db = UniversalReaction(
        #         id=bigg_id, name=db_entry.name, reference_id=db_entry.id
        #     )
        #     session.add(universal_reaction_db)
        # else:
        # if universal_reaction_db:
        #     logging.warn(f"Universal reaction already exists: {bigg_id}.")
        #     continue
        #
        universal_reaction_matrix_info = []
        reaction_matrix_info = []

        print(f"TIME: FINDREF: {(t1 := time.time()) - t0}")
        t0 = t1

        if reference_db is not None and m is not None:
            print("Using found reference")
            for participant, compound, crmapping in m[3]:
                compartment = m[1][participant.compartment]
                universal_id = f"{compound.universal_id}_{compartment}"
                comp_comp_id = f"{universal_id}:{compound.charge}"

                compartment_db = (
                    session.query(Compartment)
                    .filter(Compartment.id == compartment)
                    .first()
                )
                if not compartment_db:
                    compartment_db = Compartment(id=compartment, name=compartment)
                    session.add(compartment_db)

                universal_compartmentalized_component_db = (
                    session.query(UniversalCompartmentalizedComponent)
                    .filter(UniversalCompartmentalizedComponent.id == universal_id)
                    .first()
                )
                if not universal_compartmentalized_component_db:
                    universal_compartmentalized_component_db = (
                        UniversalCompartmentalizedComponent(
                            id=universal_id,
                            universal_component_id=compound.universal_id,
                            compartment_id=compartment_db.id,
                        )
                    )
                    session.add(universal_compartmentalized_component_db)

                compartmentalized_component_db = (
                    session.query(CompartmentalizedComponent)
                    .filter(CompartmentalizedComponent.id == comp_comp_id)
                    .first()
                )
                if not compartmentalized_component_db:
                    compartmentalized_component_db = CompartmentalizedComponent(
                        id=comp_comp_id,
                        component_id=compound.id,
                        universal_id=universal_compartmentalized_component_db.id,
                        compartment_id=compartment_db.id,
                    )
                    session.add(compartmentalized_component_db)

                part_coeff = participant.coefficient
                try:
                    part_coeff = float(part_coeff)
                except:
                    part_coeff = part_coeff.lower().replace("n", crmapping.reference_n)
                    part_coeff = float(part_coeff)
                coefficient = (-1 if m[0][participant.side] == "L" else 1) * float(
                    participant.coefficient
                )
                # TODO: Check if exists
                universal_reaction_matrix_info_item = dict(
                    universal_compartmentalized_component_id=universal_compartmentalized_component_db.id,
                    reference_reaction_participant_id=participant.id,
                    coefficient=coefficient,
                )
                universal_reaction_matrix_info.append(
                    universal_reaction_matrix_info_item
                )
                # universal_reaction_matrix_db = UniversalReactionMatrix(
                #     universal_id=universal_reaction_db.id,
                #     universal_compartmentalized_component_id=universal_compartmentalized_component_db.id,
                #     coefficient=coefficient,
                # )
                # session.add(universal_reaction_matrix_db)
                # session.commit()
                reaction_matrix_info.append(
                    dict(
                        reaction_matrix=universal_reaction_matrix_info_item,
                        compartmentalized_component_id=compartmentalized_component_db.id,
                        coefficient=coefficient,
                    )
                )
        else:
            for coeff, side, universal_id, compartment_id, charge in reaction_data[
                "parsed_participants"
            ]:
                print("Proceeding without reference")
                if charge is None:
                    print("Charge is None")
                    comp_db = (
                        session.query(Component)
                        .join(
                            ComponentReferenceMapping,
                            ComponentReferenceMapping.component_id == Component.id,
                        )
                        .join(
                            UniversalComponentReferenceMapping,
                            UniversalComponentReferenceMapping.mapping_id
                            == ComponentReferenceMapping.id,
                        )
                        .filter(UniversalComponentReferenceMapping.id == universal_id)
                        .first()
                    )
                    if comp_db is not None:
                        charge = comp_db.charge
                else:
                    print(f"Charge is {charge}")
                    charge = float(charge)
                    comp_db = (
                        session.query(Component)
                        .filter(
                            (Component.universal_id == universal_id)
                            & (Component.charge == charge)
                        )
                        .first()
                    )
                if charge is None and comp_db is None:
                    charge = 0
                    comp_db = (
                        session.query(Component)
                        .filter(
                            (Component.universal_id == universal_id)
                            & (Component.charge == charge)
                        )
                        .first()
                    )
                if charge is None or comp_db is None:
                    print(
                        f"No proper charge + component combination for {charge}, {universal_id}"
                    )
                    universal_reaction_matrix_info = None
                    break

                if isinstance(charge, float):
                    if charge.is_integer():
                        charge = int(charge)

                universal_comp_comp_id = f"{universal_id}_{compartment_id}"
                comp_comp_id = f"{universal_comp_comp_id}:{charge}"

                compartment_db = (
                    session.query(Compartment)
                    .filter(Compartment.id == compartment_id)
                    .first()
                )
                if not compartment_db:
                    compartment_db = Compartment(id=compartment_id, name=compartment_id)
                    session.add(compartment_db)

                universal_compartmentalized_component_db = (
                    session.query(UniversalCompartmentalizedComponent)
                    .filter(
                        UniversalCompartmentalizedComponent.id == universal_comp_comp_id
                    )
                    .first()
                )
                if not universal_compartmentalized_component_db:
                    universal_compartmentalized_component_db = (
                        UniversalCompartmentalizedComponent(
                            id=universal_comp_comp_id,
                            universal_component_id=universal_id,
                            compartment_id=compartment_db.id,
                        )
                    )
                    session.add(universal_compartmentalized_component_db)

                compartmentalized_component_db = (
                    session.query(CompartmentalizedComponent)
                    .filter(CompartmentalizedComponent.id == comp_comp_id)
                    .first()
                )
                if not compartmentalized_component_db:
                    compartmentalized_component_db = CompartmentalizedComponent(
                        id=comp_comp_id,
                        component_id=comp_db.id,
                        universal_id=universal_compartmentalized_component_db.id,
                        compartment_id=compartment_db.id,
                    )
                    session.add(compartmentalized_component_db)

                coefficient = (-1 if side == "L" else 1) * float(coeff)
                # TODO: Check if exists
                universal_reaction_matrix_info_item = dict(
                    universal_compartmentalized_component_id=universal_compartmentalized_component_db.id,
                    reference_reaction_participant_id=None,
                    coefficient=coefficient,
                )
                universal_reaction_matrix_info.append(
                    universal_reaction_matrix_info_item
                )
                reaction_matrix_info.append(
                    dict(
                        reaction_matrix=universal_reaction_matrix_info_item,
                        compartmentalized_component_id=compartmentalized_component_db.id,
                        coefficient=coefficient,
                    )
                )

        if not universal_reaction_matrix_info:
            print("No universal reaction info")
            continue
        universal_reaction_hash = UniversalReaction.generate_hash(
            universal_reaction_matrix_info
        )
        universal_reaction_bigg_id = None
        if universal_reaction_db is not None:
            universal_reaction_bigg_id = universal_reaction_db.id
        universal_reaction_db = (
            session.query(UniversalReaction)
            .filter(UniversalReaction.hash == universal_reaction_hash)
            .first()
        )
        print(f"universal hash 1: {universal_reaction_hash}")
        if universal_reaction_bigg_id is not None and universal_reaction_db is not None:
            if universal_reaction_bigg_id != universal_reaction_db.id:
                print(
                    "ERROR: No match between already existing universal reaction and newly proposed reaction."
                )
                # TODO: Probably rename
                continue
        if universal_reaction_db is None:
            print("Creating new universal reaction")
            reaction_name = reaction_data.get("name")
            if reaction_name is None and reference_db is not None:
                reaction_name = reference_db.name
            universal_reaction_db = UniversalReaction(
                id=bigg_id,
                name=reaction_name,
                reference_id=(None if reference_db is None else reference_db.id),
                hash=universal_reaction_hash,
            )
            session.add(universal_reaction_db)

            for urm in universal_reaction_matrix_info:
                urm["universal_id"] = universal_reaction_db.id
                urm_db = UniversalReactionMatrix(
                    universal_id=urm["universal_id"],
                    universal_compartmentalized_component_id=urm[
                        "universal_compartmentalized_component_id"
                    ],
                    reference_reaction_participant_id=urm[
                        "reference_reaction_participant_id"
                    ],
                    coefficient=urm["coefficient"],
                )
                session.add(urm_db)
                session.commit()
                urm["id"] = urm_db.id
        else:
            for urm in universal_reaction_matrix_info:
                urm_db = None
                for direction in [1, -1]:
                    urm_db = (
                        session.query(UniversalReactionMatrix)
                        .filter(
                            (
                                UniversalReactionMatrix.universal_id
                                == universal_reaction_db.id
                            )
                            & (
                                UniversalReactionMatrix.universal_compartmentalized_component_id
                                == urm["universal_compartmentalized_component_id"]
                            )
                            & (
                                UniversalReactionMatrix.coefficient
                                == direction * urm["coefficient"]
                            )
                        )
                        .first()
                    )
                    if urm_db is not None:
                        break
                if urm_db is None:
                    print("ERROR: Cannot find correct universal reaction matrix info.")
                urm["id"] = urm_db.id

        reaction_hash = Reaction.generate_hash(reaction_matrix_info)
        print(f"reaction hash 2: {reaction_hash}")
        reaction_db = (
            session.query(Reaction).filter(Reaction.id == reaction_hash).first()
        )
        if not reaction_db:
            print("Creating new reaction")
            reaction_db = Reaction(
                id=reaction_hash, universal_id=universal_reaction_db.id
            )
            session.add(reaction_db)
            for rm in reaction_matrix_info:
                reaction_matrix_db = ReactionMatrix(
                    reaction_id=reaction_db.id,
                    reaction_matrix_id=rm["reaction_matrix"]["id"],
                    compartmentalized_component_id=rm["compartmentalized_component_id"],
                )
                session.add(reaction_matrix_db)
            session.commit()
        else:
            print(reaction_db)
        session.commit()
        print(f"TIME: BIGLOOP: {(t1 := time.time()) - t0}")
        t0 = t1


@timing
def load_reactions(curated_reactions_filepath, session):
    logging.debug("Loading Curated Reactions reference data")

    data = load_bigg_id_data(curated_reactions_filepath)

    push_reactions(data, session)
