import re
from enum import Enum
from typing import Any, Optional

from cobradb.models import Component


class BiGGIDType(Enum):
    UNIVERSAL_COMPONENT = 1
    COMPONENT = 2
    UNIVERSAL_REACTION = 3
    REACTION = 4


BIGG_ID_PATTERNS = {
    BiGGIDType.UNIVERSAL_COMPONENT: re.compile(f"[a-zA-Z0-9][a-zA-Z0-9_]*"),
    BiGGIDType.COMPONENT: re.compile(f"[a-zA-Z0-9][a-zA-Z0-9_]*:\\-?[0-9]+"),
    BiGGIDType.UNIVERSAL_REACTION: re.compile(f"[A-Z0-9][a-zA-Z0-9_]*"),
    BiGGIDType.REACTION: re.compile(f"[A-Z0-9][a-zA-Z0-9_]*:[0-9]+"),
}


def is_bigg_id_valid(proposed_bigg_id: str, bigg_id_type: BiGGIDType) -> bool:
    return bool(BIGG_ID_PATTERNS[bigg_id_type].fullmatch(proposed_bigg_id) is not None)


def create_component_bigg_id(
    base_bigg_id: str,
    compartment_bigg_id: Optional[str] = None,
    charge: Optional[Any] = None,
    collection_bigg_id: Optional[str] = None,
) -> str:
    bigg_id = base_bigg_id
    if compartment_bigg_id is not None:
        bigg_id = f"{bigg_id}_{compartment_bigg_id}"
    if charge is not None:
        charge_str = Component.charge_to_string(charge)
        bigg_id = f"{bigg_id}:{charge_str}"
    if collection_bigg_id is not None:
        bigg_id = f"__{collection_bigg_id}__{bigg_id}"
    return bigg_id
