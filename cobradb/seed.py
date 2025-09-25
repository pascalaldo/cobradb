from pathlib import Path

from sqlalchemy import select
from cobradb.models import (
    Annotation,
    AnnotationLink,
    AnnotationProperty,
    Component,
    ComponentAnnotationMapping,
    ComponentIDMapping,
    ComponentReferenceMapping,
    DataSource,
    InChI,
    ReferenceCompound,
    UniversalComponent,
)
from cobradb.data_sources import DATA_SOURCE_NAMES
from modelseedpy.biochem import from_local
import subprocess
import pandas as pd

SEED_METABOLITE_PROPERTIES = [
    "name",
    "mass",
    "is_core",
    "is_obsolete",
    "is_cofactor",
    "delta_g",
    "delta_g_err",
    "pka",
    "pkb",
    "is_abstract",
    "smiles",
]
SEED_METABOLITE_ALIASES = {
    "seed.compound": "seed.compound",
    "KEGG": "kegg.compound",
    "MetaCyc": "metacyc.compound",
    "metanetx.chemical": "metanetx.chemical",
}


def download_or_update(directory):
    directory = Path(directory)
    git_directory = directory / "ModelSEEDDatabase"

    if git_directory.exists():
        process = subprocess.run(
            ["git", "pull"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            cwd=git_directory,
        )
    else:
        process = subprocess.run(
            [
                "git",
                "clone",
                # "-b",
                # "dev",
                "https://github.com/ModelSEED/ModelSEEDDatabase.git",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            cwd=directory,
        )
    if process.returncode != 0:
        print("Could not download or update the ModelSEED Database.")


def load_modelseed_database(directory):
    directory = Path(directory)
    db = from_local(str(directory / "ModelSEEDDatabase"))
    return db


def push_model_seed_metabolites(modelseed_db, session):
    data_source_dbs = {}
    for data_source_bigg_id in SEED_METABOLITE_ALIASES.values():
        data_source_db = session.scalars(
            select(DataSource)
            .filter(DataSource.bigg_id == data_source_bigg_id)
            .limit(1)
        ).first()
        if data_source_db is None:
            data_source_db = DataSource(
                bigg_id=data_source_bigg_id,
                name=DATA_SOURCE_NAMES[data_source_bigg_id],
                url_prefix=f"https://identifiers.org/{data_source_bigg_id}:",
            )
            session.add(data_source_db)
        data_source_dbs[data_source_bigg_id] = data_source_db
    session.commit()

    for cpd_id in modelseed_db.compounds:
        cpd = modelseed_db.get_seed_compound(cpd_id)

        aliases = cpd.aliases
        cpd_bigg_ids = aliases.get("BiGG", [])
        if isinstance(cpd_bigg_ids, str):
            cpd_bigg_ids = [cpd_bigg_ids]

        inchi_key = cpd.inchikey
        inchi_key = None if pd.isna(inchi_key) or len(inchi_key) == 0 else inchi_key

        charge = cpd.data.get("charge")
        charge = None if pd.isna(charge) else int(charge)

        if charge is None:
            charge = 0

        bigg_component_db = None
        for cpd_bigg_id in cpd_bigg_ids:
            bigg_component_db = session.scalars(
                select(Component)
                .join(Component.universal_component)
                .join(UniversalComponent.old_bigg_ids)
                .filter(ComponentIDMapping.old_bigg_id == cpd_bigg_id)
                .filter(Component.charge == charge)
                .limit(1)
            ).first()
            if bigg_component_db is None:
                bigg_component_db = session.scalars(
                    select(Component)
                    .join(Component.universal_component)
                    .filter(UniversalComponent.bigg_id == cpd_bigg_id)
                    .filter(Component.charge == charge)
                    .limit(1)
                ).first()
            if bigg_component_db is not None:
                break

        inchi_key_component_db = None
        if inchi_key is not None:
            inchi_key_component_db = session.execute(
                select(Component, InChI)
                .join(Component.reference_mappings)
                .join(ComponentReferenceMapping.reference_compound)
                .join(ReferenceCompound.inchi)
                .filter(Component.charge == charge)
                .filter(InChI.key == inchi_key)
            ).all()

            inchi = cpd.inchi
            inchi = None if pd.isna(inchi) or len(inchi) == 0 else inchi

            # Select proper match based on full InChI, if possible
            if inchi is None:
                if len(inchi_key_component_db) == 1:
                    inchi_key_component_db = inchi_key_component_db[0][0]
                else:
                    inchi_key_component_db = None
            else:
                for cdb, ikdb in inchi_key_component_db:
                    if ikdb.to_string() == inchi:
                        inchi_key_component_db = cdb
                        break
                else:
                    inchi_key_component_db = None

        if bigg_component_db is None and inchi_key_component_db is None:
            continue

        cpd_names = cpd.data.get("aliases", "")
        cpd_names = cpd_names.split("|")
        cpd_names = [x[5:] for x in cpd_names if x.startswith("Name:")]
        cpd_names = cpd_names[0] if len(cpd_names) > 0 else ""
        cpd_names = [x.strip() for x in cpd_names.split(";")]

        annotation_bigg_id = f"seed.metabolite:{cpd_id}"
        annotation_db = session.scalars(
            select(Annotation).filter(Annotation.bigg_id == annotation_bigg_id).limit(1)
        ).first()
        if annotation_db is None:
            annotation_db = Annotation(
                bigg_id=annotation_bigg_id,
                type="seed",
                default_data_source=data_source_dbs["seed.compound"],
            )
            session.add(annotation_db)

            for property in SEED_METABOLITE_PROPERTIES:
                prop_val = getattr(cpd, property, None)
                if prop_val is None:
                    prop_val = cpd.data.get(property)
                if isinstance(prop_val, str):
                    prop_val = prop_val.strip()
                    if len(prop_val) == 0:
                        prop_val = None
                if prop_val is None:
                    continue
                if isinstance(prop_val, int) and property.startswith("is_"):
                    prop_val = bool(prop_val)
                prop_db = AnnotationProperty(key=property)
                prop_db.value = prop_val
                annotation_db.properties.append(prop_db)

            for name in cpd_names:
                if name == getattr(cpd, "name", None):
                    continue
                prop_db = AnnotationProperty(key="name")
                prop_db.value = name
                annotation_db.properties.append(prop_db)

            for property, namespace in SEED_METABOLITE_ALIASES.items():
                if namespace == "seed.compound":
                    prop_val = cpd_id
                else:
                    prop_val = aliases.get(property)
                    if isinstance(prop_val, str):
                        prop_val = prop_val.strip()
                        if len(prop_val) == 0:
                            prop_val = None
                    if prop_val is None:
                        continue
                if isinstance(prop_val, str):
                    prop_val = {prop_val}

                for val in prop_val:
                    alias_db = AnnotationLink(
                        identifier=val, data_source=data_source_dbs[namespace]
                    )
                    annotation_db.links.append(alias_db)

        if bigg_component_db is not None and inchi_key_component_db is not None:
            if bigg_component_db.bigg_id == inchi_key_component_db.bigg_id:
                annotation_mapping = ComponentAnnotationMapping(
                    component=bigg_component_db,
                    bigg_id_match=True,
                    inchi_match=True,
                )
                annotation_db.component_mappings.append(annotation_mapping)
            else:
                annotation_mapping_1 = ComponentAnnotationMapping(
                    component=bigg_component_db,
                    bigg_id_match=True,
                    inchi_match=False,
                )
                annotation_db.component_mappings.append(annotation_mapping_1)
                annotation_mapping_2 = ComponentAnnotationMapping(
                    component=bigg_component_db,
                    bigg_id_match=False,
                    inchi_match=True,
                )
                annotation_db.component_mappings.append(annotation_mapping_2)
        elif bigg_component_db is not None:
            annotation_mapping = ComponentAnnotationMapping(
                component=bigg_component_db,
                bigg_id_match=True,
                inchi_match=(None if inchi_key is None else False),
            )
            annotation_db.component_mappings.append(annotation_mapping)
        else:
            annotation_mapping = ComponentAnnotationMapping(
                component=inchi_key_component_db,
                bigg_id_match=(None if not cpd_bigg_ids is None else False),
                inchi_match=True,
            )
            annotation_db.component_mappings.append(annotation_mapping)
    session.commit()
