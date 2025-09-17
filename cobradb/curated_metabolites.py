# -*- coding: utf-8 -*-

from sqlalchemy import select
from cobradb.models import *
from cobradb import metabolites
from cobradb.util import timing

import logging
import json


def load_bigg_id_data(filename):
    with open(filename, "r") as f:
        return json.load(f)


@timing
def push_metabolites(data, session):
    for ch, chebi_info in data["chebis"].items():
        chebi_db = session.scalars(
            select(ReferenceCompound).filter(ReferenceCompound.bigg_id == ch).limit(1)
        ).first()
        if chebi_db:
            print(f"Chebi {ch} already exists: ({chebi_info}) ({chebi_db})")
            # if chebi_info["charge"] is not None and hasattr(chebi_db, "charge"):
            #     if int(chebi_db.charge) != chebi_info["charge"]:
            #         raise ValueError("Charge mismatch.")
            #     else:
            #         print("Matching charges")
        else:
            if chebi_info.get("formula"):
                inchi_db = None
                if chebi_info.get("inchi"):
                    inchi_obj = InChI.from_string(chebi_info["inchi"])
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
                            session.add(inchi_db)
                            session.commit()
                print(f"Creating new entry for {ch}")
                chebi_db = ReferenceCompound(
                    bigg_id=ch,
                    name=chebi_info["name"],
                    formula=chebi_info["formula"],
                    charge=str(chebi_info.get("charge", 0)),
                    compound_type="small_molecule",
                    inchi=inchi_db,
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
        id_mapping_db = session.scalars(
            select(ComponentIDMapping)
            .filter(ComponentIDMapping.old_bigg_id == old_id)
            .limit(1)
        ).first()
        if id_mapping_db is not None:
            continue
        new_component_db = session.scalars(
            select(UniversalComponent.id)
            .filter(UniversalComponent.bigg_id == new_id)
            .limit(1)
        ).first()
        if new_component_db is None:
            continue
        id_mapping_db = ComponentIDMapping(
            old_bigg_id=old_id,
            new_id=new_component_db,
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
