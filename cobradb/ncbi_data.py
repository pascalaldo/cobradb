from typing import Optional, Tuple, Union

from urllib.parse import quote
import requests

import time


def rate_limit(calls_per_second: float):
    seconds_per_call = 1 / calls_per_second

    def decorator(func):
        last_call = time.time()

        def wrapper(*args, **kwargs):
            nonlocal last_call
            elapsed = time.time() - last_call
            if (wait_time := (seconds_per_call - elapsed)) > 0:
                time.sleep(wait_time)

            last_call = time.time()
            return func(*args, **kwargs)

        return wrapper

    return decorator


@rate_limit(calls_per_second=5)
def get_organism_for_ncbi_assembly_accession(
    accession: str,
) -> Optional[Tuple[str, int, Optional[str]]]:
    r = requests.get(
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{quote(accession)}/dataset_report",
        headers={"accept": "application/json"},
    )
    if r.status_code != 200:
        return None

    response = r.json()
    if response.get("total_count") != 1:
        return None

    report = response.get("reports", [None])[0]
    if report is None:
        return None

    organism = report.get("organism")
    if organism is None:
        return None

    organism_name = organism.get("organism_name")
    tax_id = organism.get("tax_id")

    if organism_name is None or tax_id is None:
        return None

    strain = organism.get("strain")
    if (
        strain is None
        and (infraspecific := organism.get("infraspecific_names")) is not None
    ):
        strain = infraspecific.get("strain")

    return organism_name, tax_id, strain
