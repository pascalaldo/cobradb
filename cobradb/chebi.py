from libchebipy import ChebiEntity as lcpChebiEntity

from libchebipy._chebi_entity import parsers
import libchebipy._parsers as chebi_parsers

import pandas as pd

parsers.set_auto_update(False)
parsers.set_download_cache_path("/chebi/libChEBI")


def parse_structures():
    """COMMENT"""
    filename = chebi_parsers.get_file("structures.tsv.gz")

    df = pd.read_csv(filename, sep="\t", index_col=0, header=0)
    for _id, row in df.iterrows():
        cpd_id = int(row["compound_id"])
        inchi_key = row["standard_inchi_key"]
        if not pd.isna(inchi_key):
            chebi_parsers.__INCHI_KEYS[cpd_id] = chebi_parsers.Structure(
                inchi_key, chebi_parsers.Structure.InChIKey, 1
            )
        smiles = row["smiles"]
        if not pd.isna(smiles):
            chebi_parsers.__SMILES[cpd_id] = chebi_parsers.Structure(
                smiles, chebi_parsers.Structure.SMILES, 1
            )
    print("Successfully loaded all CHEBI InChI Keys and SMILES structures.")


class ChebiEntity(lcpChebiEntity):
    """Class changing some default behavior of ChebiEntity, since the original often fails."""

    # def __init__(self, chebi_id):
    #     super().__init__(
    #         chebi_id,
    #         parser="filesystem",
    #         auto_update=False,
    #         download_dir="/chebi/libChEBI",
    #     )

    def __init__(self, chebi_id):
        super().__init__(chebi_id)
