from typing import Any, Dict
from sqlalchemy.orm import Session
from cobradb.api import metabolites
from cobradb.util import timing

import logging
import json


def load_bigg_id_data(filename):
    with open(filename, "r") as f:
        return json.load(f)


@timing
def push_metabolites(session: Session, data: Dict[str, Any]):
    for ch, chebi_info in data.get("chebis", {}).items():
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

    for bid, bid_info in data.get("bigg_ids", {}).items():
        if "##" in bid:
            continue

        if "complex" in bid_info:
            metabolites.create_complex_metabolite(session, bid, bid_info["complex"])
        else:
            if "chebi" in bid_info or "chebis" in bid_info:
                base_chebi = (
                    bid_info["chebi"] if "chebi" in bid_info else bid_info["chebis"][0]
                )
                ref_n = bid_info.get("reference_n")
                print(f"Creating metabolite {bid} for {base_chebi}")
                metabolites.create_metabolite(
                    session, bid, base_chebi, reference_n=ref_n
                )

            if "inchi" in bid_info:
                print("Creating entry for InChI")
                inchis = bid_info["inchi"]
                if isinstance(inchis, str):
                    inchis = [{"string": inchis, "name": bid_info.get("name")}]
                metabolites.create_metabolites_for_inchis(session, bid, inchis)

    for old_id, new_id in data.get("bigg_id_mapping", {}).items():
        metabolites.create_component_id_mapping(session, old_id, new_id)
    session.commit()


@timing
def load_bigg_ids(session: Session, curated_metabolites_filepath):
    logging.debug("Loading Curated Metabolites reference data")

    data = load_bigg_id_data(curated_metabolites_filepath)

    push_metabolites(session, data)
