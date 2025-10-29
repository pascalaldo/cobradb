from typing import Any, Iterable, Optional, Tuple, Type, TypeVar
from sqlalchemy import select
from sqlalchemy.orm import Session
import re

from cobradb.models import BiGGBase

FORMULA_PATTERN = re.compile(r"(([A-Z][a-z]?)([0-9])*)+")
FORMULA_PATTERN_SINGLE = re.compile(r"([A-Z][a-z]?)([0-9])*")


T = TypeVar("T", bound=BiGGBase)


def get_object_by_bigg_id(
    session: Session,
    bigg_id: str,
    obj_class: Type[T],
    opts: Optional[Iterable[Any]] = None,
) -> Optional[T]:
    if opts is None:
        return session.scalars(
            select(obj_class).filter(obj_class.bigg_id == bigg_id).limit(1)
        ).first()
    else:
        return session.scalars(
            select(obj_class)
            .options(*opts)
            .filter(obj_class.bigg_id == bigg_id)
            .limit(1)
        ).first()


def fix_explicit_formula(formula, allow_R=False) -> Tuple[bool, Optional[str]]:
    if not isinstance(formula, str):
        return False, None
    m = FORMULA_PATTERN.fullmatch(formula)
    if m is None:
        return False, None

    new_formula = ""
    is_original_formula = False
    for m in FORMULA_PATTERN_SINGLE.finditer(formula):
        atom = m[1]
        mult = m[2]
        if atom == "R" and not allow_R:
            return False, None
        if mult is not None and int(mult) == "1":
            new_formula = new_formula + atom
            is_original_formula = False
        else:
            new_formula = new_formula + m[0]
    return is_original_formula, new_formula


def _formula_to_dict(formula):
    d = {}
    for m in FORMULA_PATTERN_SINGLE.finditer(formula):
        atom = m[1]
        mult = m[2]
        if mult is None:
            mult = 1
        mult = int(mult)
        d[atom] = mult
    return d


def are_explicit_formulae_equivalent(formula1, formula2):
    if formula1 is None or formula2 is None:
        return False
    d1 = _formula_to_dict(formula1)
    d2 = _formula_to_dict(formula2)
    return d1 == d2
