from sqlalchemy.orm import Session
from cobradb.models import (
    Compartment,
    CompartmentalizedComponent,
    Component,
    ComponentIDMapping,
    ComponentReferenceMapping,
    Model,
    Reaction,
    ReactionMatrix,
    ReferenceCompound,
    ReferenceReaction,
    ReferenceReactionParticipant,
    UniversalCompartmentalizedComponent,
    UniversalComponent,
    UniversalComponentReferenceMapping,
    UniversalReaction,
    UniversalReactionMatrix,
)
from cobradb.util import timing

from sqlalchemy import func
from sqlalchemy import select
import logging
import json

import time

from pprint import pprint
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
    reaction_mappings = [({}, {}, parsed_participants[:], [])]
    participants_query = select(
        ReferenceReactionParticipant,
    ).filter(ReferenceReactionParticipant.reaction_id == db_entry.id)

    if (
        len(parsed_participants)
        != session.execute(
            select(func.count()).select_from(participants_query.subquery())
        ).scalar_one()
    ):
        return None
    participants_db = session.scalars(participants_query)

    for participant in participants_db:
        requires_n = "n" in str(participant.coefficient).lower()

        # print(f"Participant: {participant} -> {compound}")
        crmappings_db = session.scalars(
            select(ComponentReferenceMapping)
            .join(
                ComponentReferenceMapping.universal_component,
            )
            .join(ComponentReferenceMapping.component)
            .where(
                participant.compound_id
                == ComponentReferenceMapping.reference_compound_id,
            )
        ).all()
        participant_mappings = []
        for crmapping in crmappings_db:
            if requires_n and crmapping.reference_n is None:
                continue
            coefficient = participant.coefficient
            if requires_n:
                coefficient = coefficient.lower().replace("n", crmapping.reference_n)
            coefficient = float(coefficient)
            ms = [
                p
                for p in parsed_participants
                if p[2] == crmapping.universal_component.bigg_id
                and float(p[0]) == coefficient
                and ((p[4] is None) or (float(p[4]) == crmapping.component.charge))
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
                        crmapping.component,
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
                            # participant_mapping[4],
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


def parse_reaction_participants(session: Session, participants):
    parsed_participants = []
    for i, side in enumerate("LR"):
        for p in participants[i]:
            full_id = p[1]
            if full_id == "":
                continue
            if ":" in full_id:
                universal_comp_comp_id, charge = full_id.rsplit(":", maxsplit=1)
            else:
                universal_comp_comp_id, charge = full_id, None
            universal_id, compartment = universal_comp_comp_id.rsplit("_", maxsplit=1)
            id_mapping_db = session.scalars(
                select(ComponentIDMapping)
                .join(ComponentIDMapping.new_universal_component)
                .filter(ComponentIDMapping.old_bigg_id == universal_id)
                .limit(1)
            ).first()
            if id_mapping_db:
                print(
                    f"Mapping: {id_mapping_db.old_bigg_id} -> {id_mapping_db.new_universal_component.bigg_id}"
                )
                universal_id = id_mapping_db.new_universal_component.bigg_id
            parsed_participants.append((p[0], side, universal_id, compartment, charge))
    return parsed_participants


def find_reference_reaction(proposed_ids, parsed_participants, session):
    if proposed_ids is None:
        proposed_ids = []

    db_entries = []
    for rhea_id in proposed_ids:
        db_entry = session.scalars(
            select(ReferenceReaction)
            .filter(ReferenceReaction.bigg_id == rhea_id)
            .limit(1)
        ).first()
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

        p_ucs = [p[2] for p in parsed_participants]
        db_rows = session.execute(
            select(
                ReferenceCompound.bigg_id,
                UniversalComponent.bigg_id,
            )
            .join(
                ComponentReferenceMapping,
                ComponentReferenceMapping.reference_compound_id == ReferenceCompound.id,
            )
            .join(ComponentReferenceMapping.universal_component)
            .filter(UniversalComponent.bigg_id.in_(p_ucs))
        ).all()
        db_rows_mapping = {k: v for v, k in db_rows}
        if any(x not in db_rows_mapping for x in p_ucs):
            return None, None

        p_info = [
            dict(
                reference_compound_bigg_id=db_rows_mapping[p[2]], coefficient=abs(p[0])
            )
            for p in parsed_participants
        ]
        hash_pattern = ReferenceReaction.generate_hash(p_info, pattern=True)

        db_entries = session.scalars(
            select(ReferenceReaction).where(
                ReferenceReaction.hash.regexp_match(hash_pattern)
            )
        ).all()
        # stmt = session.query(ReferenceReaction)
        #
        # for parsed_participant in parsed_participants:
        #     rrp_alias = aliased(ReferenceReactionParticipant)
        #     crm_alias = aliased(ComponentReferenceMapping)
        #     uc_alias = aliased(UniversalComponent)
        #     stmt = (
        #         stmt.join(
        #             rrp_alias,
        #             rrp_alias.reaction_id == ReferenceReaction.id,
        #         )
        #         .join(
        #             crm_alias,
        #             crm_alias.reference_compound_id == rrp_alias.compound_id,
        #         )
        #         .join(uc_alias, uc_alias.id == crm_alias.universal_component_id)
        #         .filter(uc_alias.bigg_id == parsed_participant[2])
        #     )
        # rrp_count_alias = aliased(ReferenceReactionParticipant)
        # subq = (
        #     session.query(
        #         rrp_count_alias.reaction_id,
        #         func.count(rrp_count_alias.id).label("count"),
        #     )
        #     .group_by(rrp_count_alias.reaction_id)
        #     .subquery()
        # )
        # stmt = stmt.join(subq, subq.c.reaction_id == ReferenceReaction.id).filter(
        #     subq.c.count == len(parsed_participants)
        # )
        #
        # db_entries = set(stmt.all())

        if not db_entries:
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
def push_reactions(session: Session, data):
    for n, (bigg_id, reaction_data) in enumerate(data.items()):
        # universal_reaction_db = (
        #     session.query(UniversalReaction)
        #     .filter(UniversalReaction.id == bigg_id)
        #     .first()
        # )
        # if universal_reaction_db is not None:
        #     continue
        # if bigg_id != "LIPOt2pp":
        #     continue
        print("###")
        print(f"{bigg_id}: {reaction_data['name']}")
        print(f" RHEA: {reaction_data.get('rhea')}")
        pprint(reaction_data["participants"])

        t0 = time.time()
        reaction_data["parsed_participants"] = parse_reaction_participants(
            session, reaction_data["participants"]
        )

        # Make sure there are no participants shared between the lhs and rhs
        cons_d = {}
        is_consistent = True
        for p in reaction_data["parsed_participants"]:
            k = (p[2], p[3], p[4])
            if k in cons_d:
                if cons_d[k] != p[1]:
                    is_consistent = False
                    break
            cons_d[k] = p[1]
        if not is_consistent:
            continue

        print(reaction_data["parsed_participants"])
        print(f"TIME: PARSE: {(t1 := time.time()) - t0}")
        t0 = t1

        reference_db, m = find_reference_reaction(
            proposed_ids=reaction_data.get("rhea", [])
            + reaction_data.get("reference_ids", []),
            parsed_participants=reaction_data["parsed_participants"],
            session=session,
        )

        reaction_model_bigg_id = reaction_data.get("model_bigg_id")
        if reaction_model_bigg_id is None:
            reaction_model = None
        else:
            reaction_model = session.scalars(
                select(Model).filter(Model.bigg_id == reaction_model_bigg_id).limit(1)
            ).first()

        if reaction_model is None:
            universal_reaction_model_cond = UniversalReaction.model_id == None
        else:
            universal_reaction_model_cond = (UniversalReaction.model_id == None) | (
                UniversalReaction.model_id == reaction_model.id
            )

        universal_reaction_db = session.scalars(
            select(UniversalReaction)
            .filter(UniversalReaction.bigg_id == bigg_id)
            .filter(universal_reaction_model_cond)
            .limit(1)
        ).first()

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

        is_biomass_reaction = reference_db is None and "BIOMASS" in bigg_id.upper()
        if is_biomass_reaction and reaction_model is not None:
            bigg_id = bigg_id.removeprefix(f"__{reaction_model.bigg_id}__")

        universal_reaction_matrix_info = []
        reaction_matrix_info = []

        print(f"TIME: FINDREF: {(t1 := time.time()) - t0}")
        t0 = t1

        if reference_db is not None and m is not None:
            print("Using found reference")
            for participant, crmapping in m[3]:
                component = session.scalars(
                    select(Component)
                    .where(Component.id == crmapping.component_id)
                    .join(Component.universal_component)
                ).first()
                compartment = m[1][participant.compartment]
                universal_id = f"{component.universal_component.bigg_id}_{compartment}"
                comp_comp_id = f"{universal_id}:{component.charge}"

                compartment_db = session.scalars(
                    select(Compartment)
                    .filter(Compartment.bigg_id == compartment)
                    .limit(1)
                ).first()
                if not compartment_db:
                    compartment_db = Compartment(bigg_id=compartment, name=compartment)
                    session.add(compartment_db)

                universal_compartmentalized_component_db = session.scalars(
                    select(UniversalCompartmentalizedComponent)
                    .filter(UniversalCompartmentalizedComponent.bigg_id == universal_id)
                    .limit(1)
                ).first()
                if not universal_compartmentalized_component_db:
                    universal_compartmentalized_component_db = (
                        UniversalCompartmentalizedComponent(
                            bigg_id=universal_id,
                            universal_component_id=component.universal_component_id,
                            compartment=compartment_db,
                        )
                    )
                    session.add(universal_compartmentalized_component_db)

                compartmentalized_component_db = session.scalars(
                    select(CompartmentalizedComponent)
                    .filter(CompartmentalizedComponent.bigg_id == comp_comp_id)
                    .limit(1)
                ).first()
                if not compartmentalized_component_db:
                    compartmentalized_component_db = CompartmentalizedComponent(
                        bigg_id=comp_comp_id,
                        universal_compartmentalized_component=universal_compartmentalized_component_db,
                        compartment=compartment_db,
                    )
                    component.compartmentalized_components.append(
                        compartmentalized_component_db
                    )

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
                    universal_compartmentalized_component=universal_compartmentalized_component_db,
                    reference_reaction_participant=participant,
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
                        compartmentalized_component=compartmentalized_component_db,
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
                    comp_db = session.scalars(
                        select(Component)
                        .join(
                            ComponentReferenceMapping,
                            ComponentReferenceMapping.component_id == Component.id,
                        )
                        .join(
                            UniversalComponentReferenceMapping,
                            UniversalComponentReferenceMapping.mapping_id
                            == ComponentReferenceMapping.id,
                        )
                        .join(Component.universal_component)
                        .filter(UniversalComponent.bigg_id == universal_id)
                        .limit(1)
                    ).first()
                    if comp_db is not None:
                        charge = comp_db.charge
                else:
                    print(f"Charge is {charge}")
                    charge = float(charge)
                    comp_db = session.scalars(
                        select(Component)
                        .join(Component.universal_component)
                        .filter(
                            (UniversalComponent.bigg_id == universal_id)
                            & (Component.charge == charge)
                        )
                        .limit(1)
                    ).first()
                if charge is None and comp_db is None:
                    charge = 0
                    comp_db = session.scalars(
                        select(Component)
                        .join(Component.universal_component)
                        .filter(
                            (UniversalComponent.bigg_id == universal_id)
                            & (Component.charge == charge)
                        )
                        .limit(1)
                    ).first()
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

                compartment_db = session.scalars(
                    select(Compartment)
                    .filter(Compartment.bigg_id == compartment_id)
                    .limit(1)
                ).first()
                if not compartment_db:
                    compartment_db = Compartment(
                        bigg_id=compartment_id, name=compartment_id
                    )
                    session.add(compartment_db)

                universal_compartmentalized_component_db = session.scalars(
                    select(UniversalCompartmentalizedComponent)
                    .filter(
                        UniversalCompartmentalizedComponent.bigg_id
                        == universal_comp_comp_id
                    )
                    .limit(1)
                ).first()
                if not universal_compartmentalized_component_db:
                    universal_compartmentalized_component_db = (
                        UniversalCompartmentalizedComponent(
                            bigg_id=universal_comp_comp_id,
                            universal_component=comp_db.universal_component,
                            compartment=compartment_db,
                        )
                    )
                    session.add(universal_compartmentalized_component_db)

                compartmentalized_component_db = session.scalars(
                    select(CompartmentalizedComponent)
                    .filter(CompartmentalizedComponent.bigg_id == comp_comp_id)
                    .limit(1)
                ).first()
                if not compartmentalized_component_db:
                    compartmentalized_component_db = CompartmentalizedComponent(
                        bigg_id=comp_comp_id,
                        component=comp_db,
                        universal_compartmentalized_component=universal_compartmentalized_component_db,
                        compartment=compartment_db,
                    )
                    session.add(compartmentalized_component_db)

                coefficient = (-1 if side == "L" else 1) * float(coeff)
                # TODO: Check if exists
                universal_reaction_matrix_info_item = dict(
                    universal_compartmentalized_component=universal_compartmentalized_component_db,
                    reference_reaction_participant=None,
                    coefficient=coefficient,
                )
                universal_reaction_matrix_info.append(
                    universal_reaction_matrix_info_item
                )
                reaction_matrix_info.append(
                    dict(
                        reaction_matrix=universal_reaction_matrix_info_item,
                        compartmentalized_component=compartmentalized_component_db,
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
            universal_reaction_bigg_id = universal_reaction_db.bigg_id
        universal_reaction_db = session.scalars(
            select(UniversalReaction)
            .filter(UniversalReaction.hash == universal_reaction_hash)
            .filter(universal_reaction_model_cond)
            .limit(1)
        ).first()
        print(f"universal hash 1: {universal_reaction_hash}")
        if universal_reaction_bigg_id is not None and universal_reaction_db is not None:
            if universal_reaction_bigg_id != universal_reaction_db.bigg_id:
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
                bigg_id=bigg_id,
                name=reaction_name,
                reference=reference_db,
                hash=universal_reaction_hash,
                model=reaction_model,
            )
            session.add(universal_reaction_db)

            for urm in universal_reaction_matrix_info:
                urm["universal_reaction"] = universal_reaction_db
                urm_db = UniversalReactionMatrix(
                    universal_reaction=urm["universal_reaction"],
                    universal_compartmentalized_component=urm[
                        "universal_compartmentalized_component"
                    ],
                    reference_reaction_participant=urm[
                        "reference_reaction_participant"
                    ],
                    coefficient=urm["coefficient"],
                )
                session.add(urm_db)
                # session.commit()
                urm["urm"] = urm_db
        else:
            for urm in universal_reaction_matrix_info:
                urm_db = None
                for direction in [1, -1]:
                    urm_db = session.scalars(
                        select(UniversalReactionMatrix)
                        .filter(
                            (
                                UniversalReactionMatrix.universal_reaction
                                == universal_reaction_db
                            )
                            & (
                                UniversalReactionMatrix.universal_compartmentalized_component
                                == urm["universal_compartmentalized_component"]
                            )
                            & (
                                UniversalReactionMatrix.coefficient
                                == direction * urm["coefficient"]
                            )
                        )
                        .limit(1)
                    ).first()
                    if urm_db is not None:
                        break
                if urm_db is None:
                    print("ERROR: Cannot find correct universal reaction matrix info.")
                urm["urm"] = urm_db

        reaction_hash = Reaction.generate_hash(reaction_matrix_info)
        print(f"reaction hash 2: {reaction_hash}")
        reaction_db = session.scalars(
            select(Reaction)
            .join(Reaction.model)
            .filter(
                (Reaction.hash == reaction_hash) & (Reaction.model == reaction_model)
            )
            .limit(1)
        ).first()

        if not reaction_db:
            print("Creating new reaction")
            copy_number = (
                session.scalar(
                    select(func.count(Reaction.id)).filter(
                        Reaction.universal_reaction == universal_reaction_db
                    )
                )
                + 1
            )
            reaction_id = Reaction.create_id(universal_reaction_db.bigg_id, copy_number)
            reaction_db = Reaction(
                bigg_id=reaction_id,
                hash=reaction_hash,
                model=reaction_model,
                copy_number=copy_number,
                universal_reaction=universal_reaction_db,
            )
            session.add(reaction_db)
            # session.commit()
            for rm in reaction_matrix_info:
                reaction_matrix_db = ReactionMatrix(
                    reaction=reaction_db,
                    universal_reaction_matrix=rm["reaction_matrix"]["urm"],
                    compartmentalized_component=rm["compartmentalized_component"],
                )
                session.add(reaction_matrix_db)
            # session.commit()
        else:
            print(reaction_db)

        if is_biomass_reaction:
            biomass_reference_db = session.scalars(
                select(ReferenceReaction)
                .filter(ReferenceReaction.bigg_id == "BiGGr:BIOMASS")
                .limit(1)
            ).first()
            if biomass_reference_db is None:
                biomass_reference_db = ReferenceReaction(
                    bigg_id="BiGGr:BIOMASS",
                    name="Biomass reaction",
                    hash="#BIOMASS",
                )
            universal_reaction_db.reference = biomass_reference_db
        session.commit()
        print(f"TIME: BIGLOOP: {(t1 := time.time()) - t0}")
        t0 = t1


@timing
def load_reactions(session: Session, curated_reactions_filepath):
    logging.debug("Loading Curated Reactions reference data")

    data = load_bigg_id_data(curated_reactions_filepath)

    push_reactions(session, data)
