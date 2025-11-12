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


@timing
def load_bigg_ids(session: Session, curated_metabolites_filepath):
    logging.debug("Loading Curated Metabolites reference data")

    data = load_bigg_id_data(curated_metabolites_filepath)

    push_metabolites(session, data)
