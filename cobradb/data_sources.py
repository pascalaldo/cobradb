from sqlalchemy import select
from cobradb.models import DataSource


DATA_SOURCE_NAMES = {
    "seed.compound": "SEED Compound",
    "kegg.compound": "KEGG Compound",
    "metacyc.compound": "MetaCyc Compound",
    "metanetx.chemical": "MetaNetX Chemical",
    "CHEBI": "CHEBI",
    "drugbank": "DrugBank",
    "kegg.drug": "KEGG Drug",
    "wikipedia.en": "Wikipedia",
    "hmdb": "HMDB",
    "rhea": "RHEA",
    "GO": "GO",
    "kegg.reaction": "KEGG Reaction",
    "metacyc.reaction": "MetaCyc Reaction",
    "ec-code": "EC",
}

DATA_SOURCE_ID_MAP = {}


def get_data_source_id(bigg_id, session):
    ds_id = DATA_SOURCE_ID_MAP.get(bigg_id)
    if ds_id is not None:
        return ds_id

    ds_db = session.scalars(
        select(DataSource).filter(DataSource.bigg_id == bigg_id).limit(1)
    ).first()
    if ds_db is None:
        ds_name = DATA_SOURCE_NAMES.get(bigg_id)
        if ds_name is None:
            return None
        ds_db = DataSource(
            bigg_id=bigg_id,
            name=ds_name,
            url_prefix=f"https://identifiers.org/{bigg_id}:",
        )
        session.add(ds_db)
        session.commit()
        ds_id = ds_db.id
        DATA_SOURCE_ID_MAP[bigg_id] = ds_id
        return ds_id

    ds_id = ds_db.id
    DATA_SOURCE_ID_MAP[bigg_id] = ds_id
    return ds_id
