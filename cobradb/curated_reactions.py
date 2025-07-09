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


def load_bigg_id_data(filename):
    with open(filename, "r") as f:
        return json.load(f)


def match_reaction_data_with_db_entry(reaction_data, db_entry, session):
    stmt = (
        select(
            ReferenceReactionParticipant,
            ReferenceCompound,
        )
        .where(ReferenceReactionParticipant.reaction_id == db_entry.id)
        .join(
            ReferenceCompound,
            ReferenceCompound.id == ReferenceReactionParticipant.compound_id,
        )
    )
    reaction_mapping = ({}, {})
    participants_db = session.execute(stmt)
    for (
        participant,
        compound,
    ) in participants_db:
        print(f"Participant: {participant} -> {compound}")
        stmt = (
            select(ComponentReferenceMapping, UniversalComponent, Component)
            .where(
                compound.id == ComponentReferenceMapping.reference_id,
            )
            .join(
                UniversalComponent,
                UniversalComponent.id == ComponentReferenceMapping.universal_id,
            )
            .join(Component, Component.id == ComponentReferenceMapping.component_id)
        )
        participant_mapping = None
        crmappings_db = session.execute(stmt)
        for crmapping, universal_component, component in crmappings_db:
            print(f" {crmapping}: {universal_component}, {component}")
            m = [
                p
                for p in reaction_data["parsed_participants"]
                if p[2] == universal_component.id
            ]
            if not m:
                continue
            d_coef, d_side, d_universal_id, d_comp = m[0]
            participant_mapping = (participant.side, d_side), (
                participant.compartment,
                d_comp,
            )
            break
        if participant_mapping:
            if (side_m := participant_mapping[0][0]) in reaction_mapping[0]:
                if participant_mapping[0][1] != reaction_mapping[0][side_m]:
                    print(
                        f"Inconsistent mapping: prev: {side_m}:{reaction_mapping[0][side_m]}, new: {side_m}{participant_mapping[0][1]}"
                    )
                    return False
            else:
                reaction_mapping[0][participant_mapping[0][0]] = participant_mapping[0][
                    1
                ]
            if (comp_m := participant_mapping[1][0]) in reaction_mapping[1]:
                if participant_mapping[1][1] != reaction_mapping[1][comp_m]:
                    print(
                        f"Inconsistent mapping: prev: {comp_m}:{reaction_mapping[1][side_m]}, new: {comp_m}{participant_mapping[1][1]}"
                    )
                    return False
            else:
                reaction_mapping[1][participant_mapping[1][0]] = participant_mapping[1][
                    1
                ]
    print(
        f"Successful mapping: sides: {reaction_mapping[0]}, comp: {reaction_mapping[1]}"
    )


@timing
def push_reactions(data, session):
    for bigg_id, reaction_data in data.items():
        if reaction_data["rhea"]:
            print("###")
            print(f"{bigg_id}: {reaction_data['name']}")
            print(f" RHEA: {reaction_data['rhea']}")
            pprint(reaction_data["participants"])
            parsed_participants = []
            for i, side in enumerate("LR"):
                for p in reaction_data["participants"][i]:
                    universal_id, compartment = p[1].rsplit("_", maxsplit=1)
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
                for db_entry in db_entries:
                    print(f"$Ref: {db_entry.equation}")
                    match_reaction_data_with_db_entry(reaction_data, db_entry, session)
                    # return
    # for ch, chebi_info in data["chebis"].items():
    #     chebi_db = (
    #         session.query(ReferenceCompound).filter(ReferenceCompound.id == ch).first()
    #     )
    #     if chebi_db:
    #         print(f"Chebi {ch} already exists: ({chebi_info}) ({chebi_db})")
    #         # if chebi_info["charge"] is not None and hasattr(chebi_db, "charge"):
    #         #     if int(chebi_db.charge) != chebi_info["charge"]:
    #         #         raise ValueError("Charge mismatch.")
    #         #     else:
    #         #         print("Matching charges")
    #     else:
    #         if chebi_info.get("formula"):
    #             print(f"Creating new entry for {ch}")
    #             chebi_db = ReferenceCompound(
    #                 id=ch,
    #                 name=chebi_info["name"],
    #                 formula=chebi_info["formula"],
    #                 charge=str(chebi_info.get("charge", 0)),
    #                 compound_type="small_molecule",
    #             )
    #             session.add(chebi_db)
    #         else:
    #             print(f"Skipping {ch}, no formula")
    # session.commit()
    #
    # for bid, bid_info in data["bigg_ids"].items():
    #     bigg_ids_handled = set()
    #     universal_component_db = UniversalComponent(
    #         id=bid,
    #     )
    #     session.add(universal_component_db)
    #     for ch in bid_info["chebis"]:
    #         chebi_db = (
    #             session.query(ReferenceCompound)
    #             .filter(ReferenceCompound.id == ch)
    #             .first()
    #         )
    #         if not chebi_db:
    #             continue
    #         if chebi_db.charge is None or chebi_db.formula is None:
    #             continue
    #         try:
    #             charge_int = int(chebi_db.charge)
    #         except:
    #             continue
    #         full_bid = f"{bid}:{charge_int}"
    #         if full_bid not in bigg_ids_handled:
    #             bigg_ids_handled.add(full_bid)
    #             component_db = Component(
    #                 id=full_bid,
    #                 universal_id=bid,
    #                 name=chebi_db.name,
    #                 formula=chebi_db.formula,
    #                 charge=charge_int,
    #             )
    #             session.add(component_db)
    #
    #             component_reference_db = ComponentReferenceMapping(
    #                 component_id=full_bid,
    #                 universal_id=bid,
    #                 reference_id=ch,
    #             )
    #             session.add(component_reference_db)
    # session.commit()
    #


@timing
def load_reactions(curated_reactions_filepath, session):
    logging.debug("Loading Curated Reactions reference data")

    data = load_bigg_id_data(curated_reactions_filepath)

    push_reactions(data, session)
