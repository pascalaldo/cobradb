import hashlib
from itertools import product

DUBLETS = [
    "".join(x) for x in product(*(([chr(ord("A") + i) for i in range(26)],) * 2))
]
TRIPLETS = [
    "".join(x)
    for x in product(*(([chr(ord("A") + i) for i in range(26)],) * 3))
    if not (
        x[0] == "E" or (x[0] == "T" and (x[1] < "T" or (x[1] == "T" and x[2] < "W")))
    )
]

PROTON_FLAGS = "AMLKJIHGFEDCBNOPQRSTUVWXYZA"


def _inchikeypart(inchi, form):
    s = ""
    x = list(hashlib.sha256(inchi).digest())
    for k in form:
        if k == 3:
            b0 = x[0]
            b1 = x[1]
            b = b0 | ((b1 & 0x3F) << 8)
            s = s + TRIPLETS[b]
            rem = (b1 & 0xC0) >> 6
            x = x[2:]
            for i in range(len(x)):
                xval = x[i]
                x[i] = ((xval & 0x3F) << 2) | rem
                rem = (xval & 0xC0) >> 6
        elif k == 2:
            b0 = x[0]
            b1 = x[1]
            b = b0 | ((b1 & 0x01) << 8)
            s = s + DUBLETS[b]
            rem = (b1 & 0xFE) >> 1
            x = x[2:]
            for i in range(len(x)):
                xval = x[i]
                x[i] = ((xval & 0x01) << 7) | rem
                rem = (xval & 0xFE) >> 1
        else:
            raise ValueError()
    return s


def inchi_object_to_inchikey(inchi):
    major_part = f"{inchi.formula}"
    for x in "chq":
        if (val := getattr(inchi, x)) is not None:
            major_part = f"{major_part}/{x}{val}"
    major_part = major_part.encode("ascii")
    k_major = _inchikeypart(major_part, form=[3, 3, 3, 3, 2])

    minor_part = ""
    for x in "btms":
        if (val := getattr(inchi, x)) is not None:
            minor_part = f"{minor_part}/{x}{val}"
    n = len(minor_part)
    if n > 0 and n < 255:
        minor_part = minor_part + minor_part
    minor_part = minor_part.encode("ascii")
    k_minor = _inchikeypart(minor_part, form=[3, 3, 2])

    n_protons = getattr(inchi, "p", 0)
    n_protons = 0 if n_protons is None else int(n_protons)
    proton_flag = PROTON_FLAGS[max(min(n_protons, 13), -13) + 13]

    standard_flag = "S"
    version_flag = "A"

    return (k_major, f"{k_minor}{standard_flag}{version_flag}", proton_flag)


def inchi_object_to_string(inchi):
    s = f"InChI=1S/{inchi.formula}"
    for x in "chqpbtms":
        if (val := getattr(inchi, x)) is not None:
            s = f"{s}/{x}{val}"
    return s


def string_to_inchi_dict(s, accept_unknown_layers=False):
    if not s.startswith("InChI=1S/"):
        return None
    s = s[9:]
    layers = s.split("/")
    d = {"formula": layers[0]}
    for l in layers[1:]:
        x = l[0]
        if x not in "chqpbtms":
            if not accept_unknown_layers:
                return None
        else:
            d[x] = l[1:]
    return d
