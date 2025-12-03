from typing import Any, Dict, Iterable, Optional, Tuple, Type, TypeVar, Union
from sqlalchemy import select
from sqlalchemy.orm import Session
import re

from cobradb.models import BiGGBase

FORMULA_PATTERN = re.compile(r"(([A-Z][a-z]?)([0-9])*)+")
FORMULA_PATTERN_SINGLE = re.compile(r"([A-Z][a-z]?)([0-9]*)((?=[A-Z])|$)")
FORMULA_DELTA_PATTERN_SINGLE = re.compile(r"([A-Z][a-z]?)([\+\-]?[0-9]*)((?=[A-Z])|$)")

NFORMULA_PATTERN = re.compile(
    r"^(?P<static1>([A-Z][a-z]?[0-9]*)*)(\((?P<variable>([A-Z][a-z]?[0-9]*)+)\)n)?\.?(?P<static2>([A-Z][a-z]?[0-9]*)*)$"
)
NCHARGE_PATTERN = re.compile(
    r"^(\(?(?P<static>-?[0-9]+)\)?)?(\((?P<variable>-?[0-9]+)\)n)?$"
)

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
        if mult == "":
            mult = None
        if atom == "R" and not allow_R:
            return False, None
        if mult is not None and int(mult) == "1":
            new_formula = new_formula + atom
            is_original_formula = False
        else:
            new_formula = new_formula + m[0]
    return is_original_formula, new_formula


class Formula:
    def __init__(
        self, formula: Optional[Union[str, Dict[str, Union[int, float]]]] = None
    ):
        if formula is None:
            self._d = {}
        elif isinstance(formula, str):
            self._d = Formula.interpret_formula_string(formula)
        elif isinstance(formula, dict):
            self._d = formula
        else:
            raise ValueError()

    @staticmethod
    def interpret_formula_string(formula: str) -> Dict[str, Union[int, float]]:
        d = {}
        for m in FORMULA_DELTA_PATTERN_SINGLE.finditer(formula):
            atom = m[1]
            mult = m[2]
            if mult is None or mult == "":
                mult = 1
            mult = int(mult)
            d[atom] = mult
        return d

    def __eq__(self, other: Any) -> bool:
        if other is None:
            return False
        if isinstance(other, str):
            other = Formula(other)
        if isinstance(other, Formula):
            return self._d == other._d
        return False

    def __str__(self) -> str:
        l = [
            (k if v == 1 else f"{k}{v}")
            for k in sorted(self._d.keys())
            if (v := self._d[k]) != 0
        ]
        return "".join(l)

    def grouped_str(self):
        plus_l = [
            (k if v == 1 else f"{k}{v}")
            for k in sorted(self._d.keys())
            if (v := self._d[k]) > 0
        ]
        minus_l = [
            (k if v == 1 else f"{k}{v}")
            for k in sorted(self._d.keys())
            if (v := -self._d[k]) > 0
        ]
        s = ""
        if plus_l:
            s = s + f"+({''.join(plus_l)}) "
        if minus_l:
            s = s + f"-({''.join(minus_l)})"
        s = s.strip()
        return s

    def __getitem__(self, atom: str):
        return self._d[atom]

    def copy(self) -> "Formula":
        return Formula(self._d.copy())

    @staticmethod
    def _add_sub(lhs_obj, rhs_obj, sub=False):
        if rhs_obj is None:
            return lhs_obj
        if isinstance(rhs_obj, (str, dict)):
            rhs_obj = Formula(rhs_obj)
        if not isinstance(rhs_obj, Formula):
            raise TypeError()
        for k, v in rhs_obj._d.items():
            if sub:
                lhs_obj._d[k] = lhs_obj._d.get(k, 0) - v
            else:
                lhs_obj._d[k] = lhs_obj._d.get(k, 0) + v
        lhs_keys = list(lhs_obj._d.keys())
        for k in lhs_keys:
            if lhs_obj._d[k] == 0:
                del lhs_obj._d[k]
        return lhs_obj

    def __iadd__(self, other):
        return Formula._add_sub(self, other)

    def __isub__(self, other):
        return Formula._add_sub(self, other, sub=True)

    def __add__(self, other) -> "Formula":
        result = self.copy()
        return Formula._add_sub(result, other)

    def __sub__(self, other) -> "Formula":
        result = self.copy()
        return Formula._add_sub(result, other, sub=True)

    @staticmethod
    def _mul(lhs_obj, rhs_obj):
        if not isinstance(rhs_obj, int):
            raise TypeError()
        for k in lhs_obj._d:
            lhs_obj._d[k] *= rhs_obj
        return lhs_obj

    def __mul__(self, other) -> "Formula":
        result = self.copy()
        return Formula._mul(result, other)

    def __rmul__(self, other) -> "Formula":
        result = self.copy()
        return Formula._mul(result, other)

    def __imul__(self, other):
        return Formula._mul(self, other)


def formula_to_dict(formula: str):
    d = {}
    for m in FORMULA_PATTERN_SINGLE.finditer(formula):
        atom = m[1]
        mult = m[2]
        if mult is None or mult == "":
            mult = 1
        mult = int(mult)
        d[atom] = mult
    return d


def are_explicit_formulae_equivalent(formula1, formula2):
    if formula1 is None or formula2 is None:
        return False
    d1 = formula_to_dict(formula1)
    d2 = formula_to_dict(formula2)
    return d1 == d2


class NFormula:
    def __init__(
        self,
        formula: Optional[Union[str, Tuple[Formula, Formula]]] = None,
    ):
        if formula is None:
            self._static = Formula()
            self._variable = Formula()
        elif isinstance(formula, str):
            self._static, self._variable = NFormula.interpret_formula_string(formula)
        elif isinstance(formula, tuple):
            self._static = formula[0]
            self._variable = formula[1]
        else:
            raise ValueError()

    def fill(self, n: int):
        return self._static + (n * self._variable)

    @staticmethod
    def interpret_formula_string(
        formula: str,
    ) -> Tuple[Formula, Formula]:
        m = NFORMULA_PATTERN.match(formula)
        if m is None:
            raise ValueError()
        f_static_1 = Formula(m.group("static1"))
        f_variable = Formula(m.group("variable"))
        f_static_2 = Formula(m.group("static2"))

        f_static = f_static_1 + f_static_2
        return f_static, f_variable

    def __eq__(self, other: Any) -> bool:
        if other is None:
            return False
        if isinstance(other, str):
            other = NFormula(other)
        if isinstance(other, NFormula):
            return (self._static == other._static) and (
                self._variable == other._variable
            )
        return False

    def __str__(self) -> str:
        s = str(self._static)
        variable_str = str(self._variable)
        if len(variable_str) != 0:
            s = f"{s}({variable_str})n"
        return s

    # def grouped_str(self):
    #     plus_l = [
    #         (k if v == 1 else f"{k}{v}")
    #         for k in sorted(self._d.keys())
    #         if (v := self._d[k]) > 0
    #     ]
    #     minus_l = [
    #         (k if v == 1 else f"{k}{v}")
    #         for k in sorted(self._d.keys())
    #         if (v := -self._d[k]) > 0
    #     ]
    #     s = ""
    #     if plus_l:
    #         s = s + f"+({''.join(plus_l)}) "
    #     if minus_l:
    #         s = s + f"-({''.join(minus_l)})"
    #     s = s.strip()
    #     return s

    def __getitem__(self, atom: str):
        return (self._static[atom], self._variable[atom])

    def copy(self) -> "NFormula":
        return NFormula((self._static.copy(), self._variable.copy()))

    @staticmethod
    def _add_sub(lhs_obj, rhs_obj, sub=False):
        if rhs_obj is None:
            return lhs_obj
        if isinstance(rhs_obj, str):
            rhs_obj = NFormula(rhs_obj)
        if isinstance(rhs_obj, Formula):
            rhs_obj = NFormula((rhs_obj, Formula()))
        if not isinstance(rhs_obj, NFormula):
            raise TypeError()
        if sub:
            lhs_obj._static -= rhs_obj._static
            lhs_obj._variable -= rhs_obj._variable
        else:
            lhs_obj._static += rhs_obj._static
            lhs_obj._variable += rhs_obj._variable
        return lhs_obj

    def __iadd__(self, other):
        return NFormula._add_sub(self, other)

    def __isub__(self, other):
        return NFormula._add_sub(self, other, sub=True)

    def __add__(self, other) -> "NFormula":
        result = self.copy()
        return NFormula._add_sub(result, other)

    def __sub__(self, other) -> "NFormula":
        result = self.copy()
        return NFormula._add_sub(result, other, sub=True)


class NCharge:
    def __init__(self, charge: Optional[Union[str, Tuple[int, int]]] = None):
        if charge is None:
            self._static = 0
            self._variable = 0
        elif isinstance(charge, tuple):
            self._static = charge[0]
            self._variable = charge[1]
        elif isinstance(charge, str):
            self._static, self._variable = NCharge.interpret_charge_string(charge)
        else:
            raise TypeError()

    @staticmethod
    def interpret_charge_string(charge: str):
        m = NCHARGE_PATTERN.match(charge)
        if m is None:
            raise ValueError()
        static_part = m.group("static")
        if static_part is None or static_part == "":
            static_part = 0
        variable_part = m.group("variable")
        if variable_part is None or variable_part == "":
            variable_part = 0
        return int(static_part), int(variable_part)

    def fill(self, n: int):
        return self._static + (n * self._variable)
