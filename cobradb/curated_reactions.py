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


def match_reaction_data_with_db_entry(reaction_data, db_entry, session):
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
    reaction_mappings = [({}, {}, reaction_data["parsed_participants"][:], [])]
    participants_db = session.execute(stmt)
    for (
        participant,
        compound,
    ) in participants_db:
        print(f"Participant: {participant} -> {compound}")
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
            print(f" {crmapping}: {universal_component}")
            ms = [
                p
                for p in reaction_data["parsed_participants"]
                if p[2] == universal_component.id
            ]
            if not ms:
                continue
            for m in ms:
                d_coef, d_side, d_universal_id, d_comp = m
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
                    ),
                )
        if participant_mappings:
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
                            (participant_mapping[3], participant_mapping[4])
                        )
                        new_item[2].remove(participant_mapping[2])
                        reaction_mappings.append(new_item)
    successful_mappings = []
    for m in reaction_mappings:
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


@timing
def push_reactions(data, session):
    for n, (bigg_id, reaction_data) in enumerate(data.items()):
        # if bigg_id != "UPP3MT":
        #     continue
        if reaction_data["rhea"]:
            print("###")
            print(f"{bigg_id}: {reaction_data['name']}")
            print(f" RHEA: {reaction_data['rhea']}")
            pprint(reaction_data["participants"])
            parsed_participants = []
            for i, side in enumerate("LR"):
                for p in reaction_data["participants"][i]:
                    universal_id, compartment = p[1].rsplit("_", maxsplit=1)
                    id_mapping_db = (
                        session.query(ComponentIDMapping)
                        .filter(ComponentIDMapping.old_id == universal_id)
                        .first()
                    )
                    if id_mapping_db:
                        print(
                            f"Mapping: {id_mapping_db.old_id} -> {id_mapping_db.new_id}"
                        )
                        universal_id = id_mapping_db.new_id
                    parsed_participants.append((p[0], side, universal_id, compartment))
            reaction_data["parsed_participants"] = parsed_participants

            db_entries = []
            for rhea_id in reaction_data["rhea"]:
                db_entry = (
                    session.query(ReferenceReaction)
                    .filter(ReferenceReaction.id == rhea_id)
                    .first()
                )
                if db_entry:
                    db_entries.append(db_entry)
            if not db_entries:
                print("!No DB entry for RHEA")
            else:
                mappings = []
                for db_entry in db_entries:
                    print(f"$Ref: {db_entry.equation}")
                    m = match_reaction_data_with_db_entry(
                        reaction_data, db_entry, session
                    )
                    if m is not None:
                        mappings.append((db_entry, m))
                if len(mappings) != 1:
                    print(f"No single unique mapping.")
                else:
                    db_entry, m = mappings[0]

                    universal_reaction_db = (
                        session.query(UniversalReaction)
                        .filter(Reaction.id == bigg_id)
                        .first()
                    )
                    if not universal_reaction_db:
                        universal_reaction_db = UniversalReaction(
                            id=bigg_id, name=db_entry.name, reference_id=db_entry.id
                        )
                        session.add(universal_reaction_db)
                    else:
                        logging.warn(f"Universal reaction already exists: {bigg_id}.")
                        continue

                    reaction_matrix_info = []

                    for participant, compound in m[3]:
                        compartment = m[1][participant.compartment]
                        flip = (m[0]["L"] == "R") if "L" in m[0] else (m["R"] == "L")
                        universal_id = f"{compound.universal_id}_{compartment}"
                        comp_comp_id = f"{universal_id}:{compound.charge}"

                        compartment_db = (
                            session.query(Compartment)
                            .filter(Compartment.id == compartment)
                            .first()
                        )
                        if not compartment_db:
                            compartment_db = Compartment(
                                id=compartment, name=compartment
                            )
                            session.add(compartment_db)

                        universal_compartmentalized_component_db = (
                            session.query(UniversalCompartmentalizedComponent)
                            .filter(
                                UniversalCompartmentalizedComponent.id == universal_id
                            )
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

                        coefficient = (
                            -1 if m[0][participant.side] == "L" else 1
                        ) * float(participant.coefficient)
                        # TODO: Check if exists
                        universal_reaction_matrix_db = UniversalReactionMatrix(
                            universal_id=universal_reaction_db.id,
                            universal_compartmentalized_component_id=universal_compartmentalized_component_db.id,
                            coefficient=coefficient,
                        )
                        session.add(universal_reaction_matrix_db)
                        session.commit()
                        reaction_matrix_info.append(
                            dict(
                                reaction_matrix_id=universal_reaction_matrix_db.id,
                                compartmentalized_component_id=compartmentalized_component_db.id,
                                coefficient=coefficient,
                            )
                        )
                    reaction_hash = Reaction.generate_hash(reaction_matrix_info)
                    reaction_db = (
                        session.query(Reaction)
                        .filter(Reaction.id == reaction_hash)
                        .first()
                    )
                    if not reaction_db:
                        reaction_db = Reaction(
                            id=reaction_hash, universal_id=universal_reaction_db.id
                        )
                        session.add(reaction_db)
                        for rm in reaction_matrix_info:
                            reaction_matrix_db = ReactionMatrix(
                                reaction_id=reaction_db.id,
                                reaction_matrix_id=rm["reaction_matrix_id"],
                                compartmentalized_component_id=rm[
                                    "compartmentalized_component_id"
                                ],
                            )
                            session.add(reaction_matrix_db)
                        session.commit()
        session.commit()


@timing
def load_reactions(curated_reactions_filepath, session):
    logging.debug("Loading Curated Reactions reference data")

    data = load_bigg_id_data(curated_reactions_filepath)

    push_reactions(data, session)
