# -*- coding: utf-8 -*-

from cobradb.models import *
from cobradb import metabolites, settings
from cobradb.util import timing

from sqlalchemy import func
import re
import logging
import json
from libchebipy import ChebiEntity

import time

from pprint import pprint
import re


def load_bigg_id_data(filename):
    with open(filename, "r") as f:
        return json.load(f)


@timing
def push_metabolites(data, session):
    for ch, chebi_info in data["chebis"].items():
        chebi_db = (
            session.query(ReferenceCompound).filter(ReferenceCompound.id == ch).first()
        )
        if chebi_db:
            print(f"Chebi {ch} already exists: ({chebi_info}) ({chebi_db})")
            # if chebi_info["charge"] is not None and hasattr(chebi_db, "charge"):
            #     if int(chebi_db.charge) != chebi_info["charge"]:
            #         raise ValueError("Charge mismatch.")
            #     else:
            #         print("Matching charges")
        else:
            if chebi_info.get("formula"):
                print(f"Creating new entry for {ch}")
                chebi_db = ReferenceCompound(
                    id=ch,
                    name=chebi_info["name"],
                    formula=chebi_info["formula"],
                    charge=str(chebi_info.get("charge", 0)),
                    compound_type="small_molecule",
                )
                session.add(chebi_db)
            else:
                print(f"Skipping {ch}, no formula")
    session.commit()

    for bid, bid_info in data["bigg_ids"].items():
        if "##" in bid:
            continue

        metabolites.create_metabolite(bid, bid_info["chebis"][0], session)

    for old_id, new_id in data.get("bigg_id_mapping", {}).items():
        id_mapping_db = (
            session.query(ComponentIDMapping)
            .filter(ComponentIDMapping.old_id == old_id)
            .first()
        )
        if id_mapping_db is not None:
            continue
        new_component_db = (
            session.query(UniversalComponent.id)
            .filter(UniversalComponent.id == new_id)
            .first()
        )
        if new_component_db is None:
            continue
        id_mapping_db = ComponentIDMapping(
            old_id=old_id,
            new_id=new_id,
        )
        session.add(id_mapping_db)
    session.commit()

    # for bid, bid_info in data["bigg_ids"].items():
    #     if "##" in bid:
    #         continue
    #     bigg_ids_handled = set()
    #     universal_component_db = UniversalComponent(
    #         id=bid,
    #     )
    #     session.add(universal_component_db)
    #     session.commit()
    #     for old_id, new_id in data["bigg_id_mapping"].items():
    #         if new_id != bid:
    #             continue
    #         id_mapping_db = ComponentIDMapping(old_id=old_id, new_id=new_id)
    #         session.add(id_mapping_db)
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
    #         component_reference_db = ComponentReferenceMapping(
    #             component_id=full_bid,
    #             universal_id=bid,
    #             reference_id=ch,
    #         )
    #         session.add(component_reference_db)
    session.commit()


@timing
def load_bigg_ids(curated_metabolites_filepath, session):
    logging.debug("Loading Curated Metabolites reference data")

    data = load_bigg_id_data(curated_metabolites_filepath)

    push_metabolites(data, session)
