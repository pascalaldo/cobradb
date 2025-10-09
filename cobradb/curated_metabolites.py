from typing import Any, Dict
from sqlalchemy import select
from sqlalchemy.orm import Session
from cobradb.api.bigg_ids import create_component_bigg_id
from cobradb.api import metabolites, utils
from cobradb.models import (
    Component,
    ComponentReferenceMapping,
    ReferenceCompound,
    UniversalComponent,
)
from cobradb.util import timing

import logging
import json


def load_bigg_id_data(filename):
    with open(filename, "r") as f:
        return json.load(f)


@timing
def push_metabolites(session: Session, data: Dict[str, Any]):
    for ch, chebi_info in data["chebis"].items():
        if not chebi_info.get("formula"):
            print(f"Skipping {ch}, no formula.")

        existed, chebi_db = metabolites.get_or_create_small_molecule_reference(
            session, ch
        )
        if existed:
            print(f"ChEBI entry '{ch}' already exists.")
        elif chebi_db is None:
            print(f"ChEBI entry '{ch}' could not be created.")
    session.commit()

    for bid, bid_info in data["bigg_ids"].items():
        if "##" in bid:
            continue

        metabolites.create_metabolite(session, bid, bid_info["chebis"][0])

    for old_id, new_id in data.get("bigg_id_mapping", {}).items():
        metabolites.create_component_id_mapping(session, old_id, new_id)
    session.commit()

    for bid, bid_info in data["bigg_ids"].items():
        if "##" in bid:
            continue
        bigg_ids_handled = set()
        if (
            universal_component_db := utils.get_object_by_bigg_id(
                session, bid, UniversalComponent
            )
        ) is None:
            universal_component_db = metabolites._create_universal_component(
                session,
                bigg_id=bid,
            )
        session.commit()
        # for old_id, new_id in data["bigg_id_mapping"].items():
        #     if new_id != bid:
        #         continue
        #     id_mapping_db = ComponentIDMapping(old_id=old_id, new_id=new_id)
        #     session.add(id_mapping_db)
        for ch in bid_info["chebis"]:
            if (
                chebi_db := utils.get_object_by_bigg_id(session, ch, ReferenceCompound)
            ) is None:
                continue
            if chebi_db.charge is None or chebi_db.formula is None:
                continue
            try:
                int_charge = int(str(chebi_db.charge))
            except:
                continue
            full_bid = create_component_bigg_id(bid, charge=int_charge)

            if full_bid not in bigg_ids_handled:
                bigg_ids_handled.add(full_bid)
                if (
                    component_db := utils.get_object_by_bigg_id(
                        session, full_bid, Component
                    )
                ) is None:
                    component_db = metabolites._create_component(
                        session,
                        full_bid,
                        universal_component_db=universal_component_db,
                        reference_compound_db=chebi_db,
                    )
            else:
                component_db = utils.get_object_by_bigg_id(session, full_bid, Component)

            if (
                component_reference_db := session.scalars(
                    select(ComponentReferenceMapping)
                    .filter(ComponentReferenceMapping.component == component_db)
                    .filter(ComponentReferenceMapping.reference_compound == chebi_db)
                ).first()
            ) is None:
                component_reference_db = ComponentReferenceMapping(
                    component=component_db,
                    universal_component=universal_component_db,
                    reference_compound=chebi_db,
                )
                session.add(component_reference_db)
    session.commit()


@timing
def load_bigg_ids(session: Session, curated_metabolites_filepath):
    logging.debug("Loading Curated Metabolites reference data")

    data = load_bigg_id_data(curated_metabolites_filepath)

    push_metabolites(session, data)
