# -*- coding: utf-8 -*-

from pathlib import Path
from cobra.io import save_json_model, save_yaml_model, write_sbml_model
from cobradb.api.bigg_ids import create_component_bigg_id
from cobradb.api.metabolites import (
    create_model_specific_metabolite,
    get_universal_component_by_bigg_id,
)
from cobradb.curated_reactions import push_reactions

from cobradb.models import (
    AlreadyLoadedError,
    Chromosome,
    Compartment,
    CompartmentalizedComponent,
    Component,
    ComponentReferenceMapping,
    Gene,
    GeneReactionMatrix,
    Genome,
    Model,
    ModelCompartmentalizedComponent,
    ModelCount,
    ModelGene,
    ModelReaction,
    Publication,
    PublicationModel,
    Reaction,
    ReactionMatrix,
    ReferenceReaction,
    Synonym,
    UniversalCompartmentalizedComponent,
    UniversalReaction,
)
from cobradb import settings
from cobradb import parse
from cobradb.api import utils
from cobradb.util import (
    get_or_create_data_source,
    format_formula,
    scrub_name,
    check_none,
    timing,
)

from sqlalchemy import func, select, update
import re
import logging
from collections import defaultdict
import os
from difflib import SequenceMatcher


class GenbankNotFound(Exception):
    pass


def get_model_list():
    """Get the models that are available, as SBML, in ome_data/models"""
    return [
        x.replace(".xml", "").replace(".mat", "")
        for x in os.listdir(settings.model_directory)
        if ".xml" in x or ".mat" in x
    ]


def check_for_model(name):
    """Check for model, case insensitive, and ignore periods and underscores"""

    def min_name(n):
        return n.lower().replace(".", "").replace(" ", "").replace("_", "")

    for x in get_model_list():
        if min_name(name) == min_name(x):
            return x
    return None


def _sim(name1, name2):
    """Return true if the names are similar"""
    clean1, clean2 = [n.lower().replace(" ", "") for n in [name1, name2]]
    return SequenceMatcher(None, clean1, clean2).ratio() > 0.7


def improve_name(session, db, new_name):
    """If the new_name is a better descriptive name for the reaction or metabolite,
    then update.

    """
    cur_name = db.name
    cur_id = db.id
    # New name is not None and not similar to bigg_id
    if (
        new_name is not None
        and not _sim(new_name, cur_id)
        and (cur_name is None or _sim(cur_name, cur_id))
    ):
        logging.debug("Replacing name %s with %s" % (cur_name, new_name))
        db.name = new_name
    # session.commit()


@timing
def load_model(model_filepath, pub_ref, genome_ref, session):
    """Load a model into the database. Returns the bigg_id for the new model.

    Arguments
    ---------

    model_filepath: the path to the file where model is stored.

    pub_ref: a publication PMID or doi for the model, as a tuple like this:

        ('doi', '10.1128/ecosalplus.10.2.1')

        ('pmid', '21988831')

        Can be None

    genome_ref: A tuple specifying the genome accession type and value. The
    first element can be ncbi_accession, ncbi_assembly, or organism.

    session: An instance of Session.

    """
    # apply id normalization
    logging.debug("Parsing SBML")
    model, old_parsed_ids = parse.load_and_normalize(model_filepath)
    model_bigg_id = model.id

    # check that the model doesn't already exist
    if (
        model_db := utils.get_object_by_bigg_id(session, model_bigg_id, Model)
    ) is not None:
        raise AlreadyLoadedError(
            f"Model {model_bigg_id} already loaded",
            bigg_id=model_db.bigg_id,
            db_id=model_db.id,
        )

    # check for a genome annotation for this model
    if genome_ref is not None and genome_ref[0] == "organism":
        genome_id = None
        organism = genome_ref[1]
    elif genome_ref is not None and genome_ref[0] in [
        "ncbi_accession",
        "ncbi_assembly",
    ]:
        genome_db = session.scalars(
            select(Genome)
            .filter(Genome.accession_type == genome_ref[0])
            .filter(Genome.accession_value == genome_ref[1])
            .limit(1)
        ).first()
        if genome_db is None:
            raise GenbankNotFound(
                "Genome for model {} not found with genome_ref {}".format(
                    model_bigg_id, genome_ref
                )
            )
        genome_id = genome_db.id
        organism = genome_db.organism
    else:
        logging.info(
            "No Genome reference or organism provided for model {}".format(
                model_bigg_id
            )
        )
        genome_id = None
        organism = None

    # Load the model objects. Remember: ORDER MATTERS! So don't mess around.
    logging.debug("Loading objects for model {}".format(model.id))
    published_filename = os.path.basename(model_filepath)
    model_db_id = load_new_model(
        session, model, genome_id, pub_ref, published_filename, organism
    )

    # metabolites/components and linkouts
    # get compartment names
    if os.path.exists(settings.compartment_names):
        with open(settings.compartment_names, "r") as f:
            compartment_names = {}
            for line in f.readlines():
                sp = [x.strip() for x in line.split("\t")]
                try:
                    compartment_names[sp[0]] = sp[1]
                except IndexError:
                    continue
    else:
        logging.warning("No compartment names file")
        compartment_names = {}
    final_metabolite_ids = load_metabolites(
        session,
        model_db_id,
        model,
        compartment_names,
        old_parsed_ids["metabolites"],
    )

    # # reactions
    model_db_rxn_ids = load_reactions(
        session,
        model_db_id,
        model,
        old_parsed_ids["reactions"],
        final_metabolite_ids,
    )
    #
    # genes
    model_db_rxn_ids = {}
    load_genes(session, model_db_id, model, model_db_rxn_ids, old_parsed_ids["genes"])

    # count model objects for the model summary web page
    load_model_count(session, model_db_id)

    session.commit()

    # TODO: Some potential of security issues here. Not a problem as long as the inputs are controlled by maintainers.
    model_output_path = Path("/models/models/") / model_bigg_id
    logging.warning(f"Writing corrected model to: {model_output_path}")
    for gz in ["", ".gz"]:
        write_sbml_model(model, model_output_path.with_suffix(f".biggr.sbml{gz}"))
        save_json_model(model, model_output_path.with_suffix(f".biggr.json{gz}"))
        save_yaml_model(model, model_output_path.with_suffix(f".biggr.yaml{gz}"))

    return model_bigg_id, model_db_id


@timing
def load_new_model(session, model, genome_db_id, pub_ref, published_filename, organism):
    """Load the model.

    Arguments:
    ---------

    session: A SQLAlchemy session.

    model: A COBRApy model.

    genome_db_id: The database ID of the genome. Can be None.

    pub_ref: a publication PMID or doi for the model, as a string like this:

        doi:10.1128/ecosalplus.10.2.1

        pmid:21988831

        Can be None

    organism: The organism. Can be None.

    Returns:
    -------

    The database ID of the new model row.

    """
    model_db = Model(
        bigg_id=model.id,
        genome_id=genome_db_id,
        published_filename=published_filename,
        organism=organism,
    )
    session.add(model_db)
    if pub_ref is not None:
        session.commit()
        ref_type, ref_id = pub_ref
        publication_db = session.scalars(
            select(Publication)
            .filter(Publication.reference_type == ref_type)
            .filter(Publication.reference_id == ref_id)
            .limit(1)
        ).first()
        if publication_db is None:
            publication_db = Publication(reference_type=ref_type, reference_id=ref_id)
        publication_model_db = None
        if publication_db.id is not None and model_db.id is not None:
            publication_model_db = session.scalars(
                select(PublicationModel)
                .filter(PublicationModel.publication_id == publication_db.id)
                .filter(PublicationModel.model_id == model_db.id)
                .limit(1)
            ).first()
        if publication_model_db is None:
            publication_model_db = PublicationModel(publication=publication_db)
            model_db.publication_models.append(publication_model_db)
    session.commit()
    return model_db.id


@timing
def load_metabolites(
    session, model_db_id, model, compartment_names, old_metabolite_ids
):
    """Load the metabolites as components and model components.

    Arguments:
    ---------

    session: An SQLAlchemy session.

    model_id: The database ID for the model.

    model: The COBRApy model.

    old_metabolite_ids: A dictionary where keys are new IDs and values are old
    IDs for compartmentalized metabolites.
    Returns
    -------

    comp_comp_db_ids: A dictionary where keys are the original compartmentalized
    metabolite ids and the values are the database IDs for the compartmentalized
    components.

    final_metabolite_ids: A new dictionary where keys are original
    compartmentalized metabolite IDs from the model and values are the new
    compartmentalized metabolite IDs.

    """
    final_metabolite_ids = {}

    model_db = session.get(Model, model_db_id)

    # for each metabolite in the model
    for metabolite in model.metabolites:
        model_specific_id = metabolite.id
        metabolite_id = model_specific_id.removeprefix(f"__{model_db.bigg_id}__")
        metabolite_id = parse.remove_duplicate_tag(metabolite_id)
        print(f"Metabolite: {model_specific_id} -> {metabolite_id}")

        try:
            component_bigg_id, compartment_bigg_id = parse.split_compartment(
                metabolite_id
            )
        except Exception:
            logging.error(
                (
                    "Could not find compartment for metabolite %s in"
                    "model %s" % (metabolite_id, model_db.bigg_id)
                )
            )
            continue

        if (
            universal_component_db := get_universal_component_by_bigg_id(
                session, component_bigg_id
            )
        ) is not None:
            new_universal_bigg_id = universal_component_db.bigg_id
        else:
            new_universal_bigg_id = component_bigg_id

        print(f"# {new_universal_bigg_id}: {component_bigg_id}")
        print(f"## <> {universal_component_db}")

        # look for the formula in these places
        formula_fns = [
            lambda m: getattr(m, "formula", None),  # support cobra v0.3 and 0.4
            lambda m: m.notes.get("FORMULA", None),
            lambda m: m.notes.get("FORMULA1", None),
        ]
        # Cast to string, but not for None
        strip_str_or_none = lambda v: str(v).strip() if v is not None else None
        # Ignore the empty string
        ignore_empty_str = lambda s: s if s != "" else None
        # Use a generator for lazy evaluation
        values = (
            ignore_empty_str(strip_str_or_none(formula_fn(metabolite)))
            for formula_fn in formula_fns
        )
        # Get the first non-null result. Otherwise _formula = None.
        _formula = format_formula(next(filter(None, values), None))
        # Check for non-valid formulas
        if parse.invalid_formula(_formula):
            logging.warning(
                "Invalid formula %s for metabolite %s in model %s"
                % (_formula, metabolite_id, model.id)
            )
            _formula = None
        _is_orig, _formula = utils.fix_explicit_formula(_formula)

        # get charge
        try:
            charge = int(metabolite.charge)
            # check for float charge
            if charge != metabolite.charge:
                logging.warning(
                    "Could not load charge {} for {} in model {}".format(
                        metabolite.charge, metabolite_id, model.id
                    )
                )
                charge = None
        except Exception:
            if hasattr(metabolite, "charge") and metabolite.charge is not None:
                logging.debug(
                    "Could not convert charge to integer for metabolite {} in model {}: {}".format(
                        metabolite_id, model.id, metabolite.charge
                    )
                )
            charge = None

        if charge is None:
            charge = 0

        new_biggr_id = create_component_bigg_id(new_universal_bigg_id, charge=charge)

        print(f"% {new_biggr_id}")

        if (
            metabolite_db := utils.get_object_by_bigg_id(
                session, new_biggr_id, Component
            )
        ) is not None:
            charge_zero = 0 if charge is None else charge
            if (
                metabolite_db.charge != charge_zero
                or not utils.are_explicit_formulae_equivalent(
                    metabolite_db.formula, _formula
                )
            ):
                logging.warning(
                    f"Found component, but charge or formula did not match: {metabolite_db} (charge: {metabolite_db.charge}, formula: {metabolite_db.formula}) != (charge: {charge}, formula: {_formula})"
                )
                metabolite_db = None

        # if necessary, add the new metabolite, and keep track of the ID
        if metabolite_db is None:
            print("Creating a model-specific metabolite entry.")
            # make the new metabolite
            universal_component_db, metabolite_db = create_model_specific_metabolite(
                bigg_id=new_universal_bigg_id,
                model_db=model_db,
                charge=charge,
                formula=_formula,
                name="",
                session=session,
            )

        print(f"%% <> {metabolite_db} ({universal_component_db})")
        if metabolite_db is None or universal_component_db is None:
            raise ValueError("Error parsing metabolites")

        final_metabolite_ids[metabolite_id] = metabolite_db.bigg_id

        compartment_db = utils.get_object_by_bigg_id(
            session, compartment_bigg_id, Compartment
        )
        if compartment_db is None:
            raise ValueError(f"Could not find compartment {compartment_bigg_id}")

        universal_comp_comp_id = create_component_bigg_id(
            universal_component_db.bigg_id,
            compartment_bigg_id=compartment_db.bigg_id,
        )
        comp_comp_id = create_component_bigg_id(
            universal_component_db.bigg_id,
            compartment_bigg_id=compartment_db.bigg_id,
            charge=metabolite_db.charge,
        )

        # if there is no compartmentalized component, add a new one
        if (
            universal_comp_component_db := utils.get_object_by_bigg_id(
                session, universal_comp_comp_id, UniversalCompartmentalizedComponent
            )
        ) is None:
            universal_comp_component_db = UniversalCompartmentalizedComponent(
                bigg_id=universal_comp_comp_id,
                universal_component=universal_component_db,
                compartment=compartment_db,
            )
            session.add(universal_comp_component_db)

        # if there is no compartmentalized component, add a new one
        if (
            comp_component_db := utils.get_object_by_bigg_id(
                session, comp_comp_id, CompartmentalizedComponent
            )
        ) is None:
            comp_component_db = CompartmentalizedComponent(
                bigg_id=comp_comp_id,
                component=metabolite_db,
                universal_compartmentalized_component=universal_comp_component_db,
                compartment=compartment_db,
            )
            session.add(comp_component_db)

        # if there is no model compartmentalized component, add a new one
        model_comp_comp_db = session.scalars(
            select(ModelCompartmentalizedComponent)
            .filter(
                ModelCompartmentalizedComponent.compartmentalized_component
                == comp_component_db
            )
            .filter(ModelCompartmentalizedComponent.model == model_db)
            .limit(1)
        ).first()
        if model_comp_comp_db is None:
            model_comp_comp_db = ModelCompartmentalizedComponent(
                model=model_db,
                compartmentalized_component=comp_component_db,
            )
            session.add(model_comp_comp_db)

        if universal_comp_component_db is not None:
            new_id = universal_comp_component_db.bigg_id
            if new_id != metabolite.id:
                if model.metabolites.has_id(new_id):
                    other_metabolite = model.metabolites.get_by_id(new_id)
                    if other_metabolite.charge == metabolite.charge:
                        raise Exception("Two identical metabolites.")
                    other_metabolite.id = create_component_bigg_id(
                        other_metabolite.id, charge=other_metabolite.charge
                    )
                    metabolite.id = comp_component_db.bigg_id
                else:
                    metabolite.id = new_id
        if comp_component_db is not None and metabolite_db is not None:
            metabolite.annotation["biggr"] = comp_component_db.bigg_id
            metabolite.annotation["sbo"] = "SBO:0000247"

            comp_ref_map = session.scalars(
                select(ComponentReferenceMapping)
                .join(ComponentReferenceMapping.reference_compound)
                .filter(ComponentReferenceMapping.component == metabolite_db)
            ).all()
            if comp_ref_map:
                chebis = []
                for crm in comp_ref_map:
                    if crm.reference_compound.bigg_id.startswith("CHEBI:"):
                        chebis.append(crm.reference_compound.bigg_id)
                    metabolite.annotation["sbo"] = crm.reference_compound.get_sbo()

                if chebis:
                    metabolite.annotation["chebi"] = chebis

    return final_metabolite_ids


def compartmentalized_id_to_universal_compartmentalized_id(comp_comp_id):
    try:
        universal_comp_comp_id, _charge = comp_comp_id.rsplit(":", maxsplit=1)
        return universal_comp_comp_id
    except:
        return comp_comp_id


@timing
def load_reactions(
    session,
    model_db_id,
    model,
    old_reaction_ids,
    final_metabolite_ids,
):
    """Load the reactions and stoichiometries into the model.

    TODO if the reaction is already loaded, we need to check the stoichometry
    has. If that doesn't match, then add a new reaction with an incremented ID
    (e.g. ACALD_1)

    Arguments
    ---------

    session: An SQLAlchemy session.

    model_db_id: The database ID for the model.

    model: The COBRApy model.

    old_reaction_ids: A dictionary where keys are new IDs and values are old IDs
    for reactions.

    comp_comp_db_ids: A dictionary where keys are the original compartmentalized
    metabolite ids and the values are the database IDs for the compartmentalized
    components.

    final_metabolite_ids: A new dictionary where keys are original
    compartmentalized metabolite IDs from the model and values are the new
    compartmentalized metabolite IDs.

    Returns
    -------

    A dictionary with keys for reaction BiGG IDs in the model and values for the
    associated ModelReaction.id in the database.

    """

    model_db = session.get(Model, model_db_id)
    model_db_rxn_ids = {}
    for reaction in model.reactions:
        # Drop duplicates label
        model_specific_id = reaction.id
        reaction_id = model_specific_id.removeprefix(f"__{model_db.bigg_id}__")
        reaction_id = parse.remove_duplicate_tag(reaction_id)
        reaction_id, model_reaction_copy_number = parse.split_id_and_copy_tag(
            reaction_id
        )
        # TODO: Make sure to keep copy number equal model and biggr

        # if reaction_id != "EX_glc__D_e" and reaction_id != "EX___iML1515__glc__D_e":
        #     continue

        participants = [
            dict(
                compartmentalized_component_bigg_id=(
                    m.id if ":" in m.id else f"{m.id}:{int(m.charge)}"
                ),
                coefficient=coeff,
            )
            for m, coeff in reaction.metabolites.items()
        ]
        reaction_hash = Reaction.generate_hash(participants)
        # print(f"reaction hash 0: {reaction_hash}")

        # Get the reaction
        reaction_db = session.scalars(
            select(Reaction)
            .filter(
                (Reaction.hash == reaction_hash)
                & ((Reaction.model_id == None) | (Reaction.model == model_db))
            )
            .join(Reaction.universal_reaction)
            .limit(1)
        ).first()

        if reaction_db is None:
            # universal_participants_d = {}
            # for m, coeff in reaction.metabolites.items():
            #     comp_comp_bigg_id = m.id if ":" in m.id else f"{m.id}:{int(m.charge)}"
            #     default_cc_alias = aliased(CompartmentalizedComponent)
            #     default_c_alias = aliased(Component)
            #     default_db = session.execute(
            #         select(default_cc_alias, default_c_alias, UniversalCompartmentalizedComponent, Compartment, CompartmentalizedComponent, Component).join(
            #             default_c_alias,
            #             default_cc_alias.component_id == default_c_alias.id,
            #         )
            #         .join(default_cc_alias.compartment)
            #         .join(default_cc_alias.universal_compartmentalized_component)
            #         .join(default_c_alias.reference_mappings)
            #         .join(ComponentReferenceMapping.universal_component_reference_mapping)
            #         .join(ComponentReferenceMapping.component)
            #         .join(CompartmentalizedComponent, (CompartmentalizedComponent.component_id == Component.id) & (CompartmentalizedComponent.compartment_id == default_cc_alias.compartment_id))
            #         .filter(CompartmentalizedComponent.bigg_id == comp_comp_bigg_id)
            #         .limit(1)
            #     ).first()
            #
            #     if default_db is None:
            #         ucc_bigg_id = compartmentalized_id_to_universal_compartmentalized_id(
            #                 comp_comp_db_ids.get(m.id, m.id)
            #             )
            #         u_part = [dict(
            #             universal_compartmentalized_component_bigg_id=ucc_bigg_id,
            #             coefficient=coeff,
            #         )]
            #     else:
            #         default_cc_db, default_c_db, default_ucc_db, compartment_db, cc_db, c_db = default_db
            #         u_part = [dict(
            #             universal_compartmentalized_component_bigg_id=default_ucc_db,
            #             coefficient=coeff,
            #         )]
            #         if default_c_db.charge != c_db.charge:
            #             u_part.append(
            #                 dict(
            #                     universal_compartmentalized_component_bigg_id=f"h_{compartment_db.bigg_id}",
            #                     coefficient=(c_db.charge - default_c_db.charge)*coeff,
            #                 )
            #             )
            #     for p in u_part:
            #         ucc_bigg_id = p["univeral_compartmentalized_component_bigg_id"]
            #         if ucc_bigg_id in universal_participants_d:
            #             universal_participants_d[ucc_bigg_id]["coefficient"] += p["coefficient"]
            #         else:
            #             universal_participants_d[ucc_bigg_id] = p
            # universal_participants = [v for v in universal_participants_d.values() if v["coefficient"] != 0]
            universal_participants = [
                dict(
                    universal_compartmentalized_component_bigg_id=compartmentalized_id_to_universal_compartmentalized_id(
                        m.id
                    ),
                    coefficient=coeff,
                )
                for m, coeff in reaction.metabolites.items()
            ]
            universal_reaction_hash = UniversalReaction.generate_hash(
                universal_participants
            )
            print(f"universal hash 0: {universal_reaction_hash}")
            universal_reaction_db = session.scalars(
                select(UniversalReaction)
                .filter(UniversalReaction.hash == universal_reaction_hash)
                .limit(1)
            ).first()
            if universal_reaction_db is not None:
                logging.warn(
                    f"\t UniversalReaction {reaction_id}: {universal_reaction_db}"
                )
                # TODO: Generate reaction with alternative component variants
                reaction_data = {
                    universal_reaction_db.bigg_id: {
                        "name": reaction.name,
                        "participants": [
                            [
                                (abs(coeff), x["compartmentalized_component_bigg_id"])
                                for x in participants
                                if (coeff := float(x["coefficient"])) < 0
                            ],
                            [
                                (abs(coeff), x["compartmentalized_component_bigg_id"])
                                for x in participants
                                if (coeff := float(x["coefficient"])) >= 0
                            ],
                        ],
                    }
                }
                print("! Creating new reaction variant.")
                push_reactions(session, reaction_data)
                session.commit()
                print(f"reaction hash 3: {reaction_hash}")
                reaction_db = session.scalars(
                    select(Reaction)
                    .filter(
                        (Reaction.hash == reaction_hash)
                        & ((Reaction.model_id == None) | (Reaction.model == model_db))
                    )
                    .join(Reaction.universal_reaction)
                    .limit(1)
                ).first()
            else:
                # Check for exchange reactions
                is_exchange = False
                if (
                    len(participants) == 1
                    and abs(float(participants[0]["coefficient"])) == 1
                ):
                    universal_compartmentalized_component_bigg_id = (
                        universal_participants[0][
                            "universal_compartmentalized_component_bigg_id"
                        ]
                    )
                    if (
                        reaction_id
                        != f"EX_{universal_compartmentalized_component_bigg_id}"
                    ):
                        print(f"Wrong name for exchange reaction: {reaction_id}")
                        reaction_id = (
                            f"EX_{universal_compartmentalized_component_bigg_id}"
                        )
                    is_exchange = True

                reaction_model_bigg_id = None
                if not is_exchange:
                    print("! Creating model-specific reaction.")
                    reaction_id = model_specific_id
                    reaction_model_bigg_id = model_db.bigg_id
                else:
                    print("! Creating exchange reaction.")

                reaction_data = {
                    reaction_id: {
                        "name": reaction.name,
                        "participants": [
                            [
                                (abs(coeff), x["compartmentalized_component_bigg_id"])
                                for x in participants
                                if (coeff := float(x["coefficient"])) < 0
                            ],
                            [
                                (abs(coeff), x["compartmentalized_component_bigg_id"])
                                for x in participants
                                if (coeff := float(x["coefficient"])) >= 0
                            ],
                        ],
                        "model_bigg_id": reaction_model_bigg_id,
                    }
                }
                push_reactions(session, reaction_data)
                session.commit()
                reaction_db = session.scalars(
                    select(Reaction)
                    .filter(
                        (Reaction.hash == reaction_hash)
                        & ((Reaction.model_id == None) | (Reaction.model == model_db))
                    )
                    .join(Reaction.universal_reaction)
                    .limit(1)
                ).first()

        if reaction_db is None:
            print("ERROR: Reaction was not correctly created.")

        logging.warn(f"Reaction {reaction_id}: {reaction_db}")

        reaction_matrix_db = session.scalars(
            select(ReactionMatrix)
            .join(ReactionMatrix.universal_reaction_matrix)
            .join(ReactionMatrix.compartmentalized_component)
            .filter(ReactionMatrix.reaction == reaction_db)
        ).all()
        model_reaction_is_reversed = False
        for rm_db in reaction_matrix_db:
            # TODO: May fail in some cases where metabolites occur at both sides
            for m, coeff in reaction.metabolites.items():
                comp_comp_id = f"{m.id}:{m.charge}"
                if comp_comp_id == rm_db.compartmentalized_component.bigg_id:
                    if rm_db.universal_reaction_matrix.coefficient == coeff:
                        pass
                    elif rm_db.universal_reaction_matrix.coefficient == -1 * coeff:
                        model_reaction_is_reversed = True
                    else:
                        logging.error(
                            f"Coefficients do not match for {reaction_id}, {comp_comp_id}"
                        )
            break

        # If the reaction is reversed, then switch upper and lower bound
        lower_bound = (
            -reaction.upper_bound
            if model_reaction_is_reversed
            else reaction.lower_bound
        )
        upper_bound = (
            -reaction.lower_bound
            if model_reaction_is_reversed
            else reaction.upper_bound
        )

        # TODO: Flip reaction in model if necessary

        # subsystem
        subsystem = check_none(reaction.subsystem.strip())

        copy_number = (
            session.scalars(
                select(func.count(ModelReaction.id))
                .join(ModelReaction.reaction)
                .filter(
                    Reaction.universal_reaction_id == reaction_db.universal_reaction_id
                )
                .filter(ModelReaction.model == model_db)
            ).first()
            + 1
        )

        reaction.id = (
            reaction_db.universal_reaction.bigg_id
            if copy_number == 1
            else f"{reaction_db.universal_reaction.bigg_id}:{copy_number}"
        )

        uni_ref_db = session.execute(
            select(UniversalReaction, ReferenceReaction)
            .outerjoin(
                ReferenceReaction,
                UniversalReaction.reference_id == ReferenceReaction.id,
            )
            .filter(UniversalReaction.id == reaction_db.universal_reaction_id)
            .limit(1)
        ).first()
        reaction.annotation["sbo"] = "SBO:0000176"
        if uni_ref_db is not None:
            universal_db, reference_db = uni_ref_db
            reaction.annotation["sbo"] = universal_db.get_sbo(reference_db)
            if reference_db is not None:
                if reference_db.bigg_id.startswith("RHEA:"):
                    reaction.annotation["rhea"] = reference_db.bigg_id

        # make a new reaction
        model_reaction_id = ModelReaction.create_id(
            model_db.bigg_id, reaction_db.universal_reaction.bigg_id, copy_number
        )
        model_reaction_db = ModelReaction(
            bigg_id=model_reaction_id,
            model=model_db,
            reaction=reaction_db,
            gene_reaction_rule=reaction.gene_reaction_rule,
            original_gene_reaction_rule=reaction.gene_reaction_rule,
            upper_bound=upper_bound,
            lower_bound=lower_bound,
            objective_coefficient=reaction.objective_coefficient,
            copy_number=copy_number,
            subsystem=subsystem,
        )
        session.add(model_reaction_db)
        session.commit()
    return model_db_rxn_ids


# find gene functions
def _match_gene_by_fns(fn_list, session, gene_id, chromosome_ids):
    """Go through each funciton and look for a match."""
    for fn in fn_list:
        match, is_alternative_transcript = fn(session, gene_id, chromosome_ids)
        if len(match) > 0:
            if len(match) > 1:
                logging.warning(
                    "Multiple matches for gene {} with function {}. Using the first match.".format(
                        gene_id, fn.__name__
                    )
                )
            return match[0], is_alternative_transcript
    return None, False


def _by_bigg_id(session, gene_id, chromosome_ids):
    # look for a matching model gene
    gene_db = session.scalars(
        select(Gene)
        .filter(func.lower(Gene.bigg_id) == func.lower(gene_id))
        .filter(Gene.chromosome_id.in_(chromosome_ids))
    ).all()
    return gene_db, False


def _by_name(session, gene_id, chromosome_ids):
    gene_db = session.scalars(
        select(Gene)
        .filter(func.lower(Gene.name) == func.lower(gene_id))
        .filter(Gene.chromosome_id.in_(chromosome_ids))
    ).all()
    return gene_db, False


def _by_synonym(session, gene_id, chromosome_ids):
    gene_db = session.scalars(
        select(Gene)
        .join(Synonym, Synonym.ome_id == Gene.id)
        .filter(Gene.chromosome_id.in_(chromosome_ids))
        .filter(func.lower(Synonym.synonym) == func.lower(gene_id))
    ).all()
    return gene_db, False


def _by_alternative_transcript(session, gene_id, chromosome_ids):
    """Function to check for the alternative transcript match."""
    check = re.match(r"(.*)_AT[0-9]{1,2}$", gene_id)
    if not check:
        gene_db = []
    else:
        # find the old gene
        gene_db = session.scalars(
            select(Gene)
            .filter(Gene.chromosome_id.in_(chromosome_ids))
            .filter(func.lower(Gene.bigg_id) == func.lower(check.group(1)))
            .filter(Gene.alternative_transcript_of.is_(None))
        ).all()
    return gene_db, True


def _by_alternative_transcript_name(session, gene_id, chromosome_ids):
    """Function to check for the alternative transcript match."""
    check = re.match(r"(.*)_AT[0-9]{1,2}$", gene_id)
    if not check:
        gene_db = []
    else:
        # find the old gene
        gene_db = session.scalars(
            select(Gene)
            .filter(Gene.chromosome_id.in_(chromosome_ids))
            .filter(func.lower(Gene.name) == func.lower(check.group(1)))
            .filter(Gene.alternative_transcript_of.is_(None))
        ).all()
    return gene_db, True


def _by_alternative_transcript_synonym(session, gene_id, chromosome_ids):
    """Function to check for the alternative transcript match."""
    check = re.match(r"(.*)_AT[0-9]{1,2}$", gene_id)
    if not check:
        gene_db = []
    else:
        # find the old gene
        gene_db = session.scalars(
            select(Gene)
            .join(Synonym, Synonym.ome_id == Gene.id)
            .filter(Gene.chromosome_id.in_(chromosome_ids))
            .filter(func.lower(Synonym.synonym) == func.lower(check.group(1)))
            .filter(Gene.alternative_transcript_of.is_(None))
        ).all()
    return gene_db, True


def _by_bigg_id_no_underscore(session, gene_id, chromosome_ids):
    """Matches for T maritima genes"""
    # look for a matching model gene
    gene_db = session.scalars(
        select(Gene)
        .filter(func.lower(Gene.bigg_id) == func.lower(gene_id.replace("_", "")))
        .filter(Gene.chromosome_id.in_(chromosome_ids))
    ).all()
    return gene_db, False


def _replace_gene_str(rule, old_gene, new_gene):
    return re.sub(r"\b" + old_gene + r"\b", new_gene, rule)


@timing
def load_genes(session, model_db_id, model, model_db_rxn_ids, old_gene_ids):
    """Load the genes for this model.

    Arguments:
    ---------

    session: An SQLAlchemy session.

    model_db_id: The database ID for the model.

    model: The COBRApy model.

    model_db_rxn_ids: A dictionary with keys for reactions in the model and
    values for the associated ModelReaction.id in the database.

    old_gene_ids: A dictionary where keys are new IDs and values are old IDs for
    genes.

    """
    # only grab this once
    data_source_id = get_or_create_data_source(session, "old_bigg_id")

    # find the model in the db
    model_db = session.get(Model, model_db_id)

    # find the chromosomes in the db
    chromosome_ids = session.scalars(
        select(Chromosome.id).filter(Chromosome.genome_id == model_db.genome_id)
    ).all()
    chromosome_ids = [c for c in chromosome_ids]
    if len(chromosome_ids) == 0:
        logging.warning("No chromosomes for model %s" % model_db.bigg_id)

    context = {}

    # keep track of the gene-reaction associations
    gene_bigg_id_to_model_reaction_db_ids = defaultdict(set)
    for reaction in model.reactions:
        # find the ModelReaction that corresponds to this particular reaction in
        # the model
        # if reaction.id not in model_db_rxn_ids:
        #     continue
        universal_reaction_id, copy_number = ModelReaction.interpret_id(reaction.id)
        model_reaction_db = session.scalars(
            select(ModelReaction)
            .join(ModelReaction.reaction)
            .join(Reaction.universal_reaction)
            .filter(
                (UniversalReaction.bigg_id == universal_reaction_id)
                & (ModelReaction.copy_number == copy_number)
                & (ModelReaction.model == model_db)
            )
            .limit(1)
        ).first()

        if model_reaction_db is None:
            logging.error(
                "Could not find ModelReaction {} for {} in model {}. Cannot load GeneReactionMatrix entries".format(
                    model_db_rxn_ids[reaction.id], reaction.id, model.id
                )
            )
            continue
        for gene in reaction.genes:
            gene_bigg_id_to_model_reaction_db_ids[gene.id].add(model_reaction_db.id)

    # load the genes
    for gene in model.genes:
        if len(chromosome_ids) == 0:
            gene_db = None
            is_alternative_transcript = False
        else:
            # find a matching gene
            fns = [
                _by_bigg_id,
                _by_name,
                # _by_synonym,
                # _by_alternative_transcript,
                # _by_alternative_transcript_name,
                # _by_alternative_transcript_synonym,
                _by_bigg_id_no_underscore,
            ]
            gene_db, is_alternative_transcript = _match_gene_by_fns(
                fns, session, gene.id, chromosome_ids
            )

        dup_gene_alternative_transcript = (not gene_db) and is_alternative_transcript
        if not gene_db:
            # add
            if len(chromosome_ids) > 0:
                logging.warning(
                    "Gene not in genbank file: {} from model {}".format(
                        gene.id, model.id
                    )
                )
            gene_db = session.scalars(
                select(Gene).filter(Gene.bigg_id == gene.id).limit(1)
            ).first()
            if gene_db is None:
                gene_db = Gene(
                    bigg_id=gene.id,
                    # name is optional in cobra 0.4b2. This will probably change back.
                    name=scrub_name(getattr(gene, "name", None)),
                    mapped_to_genbank=False,
                )
                session.add(gene_db)
            # session.commit()

        context[gene.id] = {
            "gene": gene,
            "gene_db": gene_db,
            "dup_gene_alternative_transcript": dup_gene_alternative_transcript,
        }
        if dup_gene_alternative_transcript:
            # duplicate gene for the alternative transcript
            old_gene_db = gene_db
            ome_gene = {}
            ome_gene["bigg_id"] = gene.bigg_id
            ome_gene["name"] = old_gene_db.name
            ome_gene["leftpos"] = old_gene_db.leftpos
            ome_gene["rightpos"] = old_gene_db.rightpos
            ome_gene["chromosome_id"] = old_gene_db.chromosome_id
            ome_gene["strand"] = old_gene_db.strand
            ome_gene["mapped_to_genbank"] = True
            ome_gene["alternative_transcript_of"] = old_gene_db.id
            gene_db = Gene(**ome_gene)
            session.add(gene_db)
            # session.commit()
            context[gene.id]["old_gene_db"] = old_gene_db
    session.commit()

    for gene_id, ctx in context.items():
        # add the deprecated id if necessary
        gene_db, old_gene_db, dup_gene_alternative_transcript = (
            ctx["gene_db"],
            ctx.get("old_gene_db"),
            ctx["dup_gene_alternative_transcript"],
        )

        # if dup_gene_alternative_transcript:
        #     # duplicate all the synonyms
        #     synonyms_db = (
        #         session.query(Synonym).filter(Synonym.ome_id == old_gene_db.id).all()
        #     )
        #     for syn_db in synonyms_db:
        #         # add a new synonym
        #         ome_synonym = {}
        #         ome_synonym["type"] = syn_db.type
        #         ome_synonym["ome_id"] = gene_db.id
        #         ome_synonym["synonym"] = syn_db.synonym
        #         ome_synonym["data_source_id"] = syn_db.data_source_id
        #         synonym_object = Synonym(**ome_synonym)
        #         session.add(synonym_object)
        #
        # add model gene
        model_gene_db = session.scalars(
            select(ModelGene)
            .filter(ModelGene.gene == gene_db)
            .filter(ModelGene.model == model_db)
            .limit(1)
        ).first()
        if model_gene_db is None:
            model_gene_db = ModelGene(gene_id=gene_db.id, model=model_db)
            session.add(model_gene_db)
            # session.commit()
        ctx["model_gene_db"] = model_gene_db

    session.commit()
    for gene_id, ctx in context.items():
        gene_db, model_gene_db = ctx["gene_db"], ctx["model_gene_db"]
        # add old gene synonym
        old_bigg_synonyms = {}
        # for old_bigg_id in old_gene_ids[gene_id]:
        #     synonym_db = (
        #         session.query(Synonym)
        #         .filter(Synonym.type == "gene")
        #         .filter(Synonym.ome_id == gene_db.id)
        #         .filter(Synonym.synonym == old_bigg_id)
        #         .filter(Synonym.data_source_id == data_source_id)
        #         .first()
        #     )
        #     if synonym_db is None:
        #         synonym_db = Synonym(
        #             type="gene",
        #             ome_id=gene_db.id,
        #             synonym=old_bigg_id,
        #             data_source_id=data_source_id,
        #         )
        #         session.add(synonym_db)
        #         # session.commit()
        #     old_bigg_synonyms[old_bigg_id] = synonym_db
        ctx["old_bigg_synonyms"] = old_bigg_synonyms
    session.commit()
    for gene_id, ctx in context.items():
        gene_db, model_gene_db, old_bigg_synonyms = (
            ctx["gene_db"],
            ctx["model_gene_db"],
            ctx["old_bigg_synonyms"],
        )
        # for old_bigg_id, synonym_db in old_bigg_synonyms.items():
        #     # add OldIDSynonym
        #     old_id_db = (
        #         session.query(OldIDSynonym)
        #         .filter(OldIDSynonym.type == "model_gene")
        #         .filter(OldIDSynonym.ome_id == model_gene_db.id)
        #         .filter(OldIDSynonym.synonym_id == synonym_db.id)
        #         .first()
        #     )
        #     if old_id_db is None:
        #         old_id_db = OldIDSynonym(
        #             type="model_gene", ome_id=model_gene_db.id, synonym_id=synonym_db.id
        #         )
        #         session.add(old_id_db)
        #         # session.commit()

        # find model reaction
        try:
            model_reaction_db_ids = gene_bigg_id_to_model_reaction_db_ids[gene_id]
        except KeyError:
            # error message above
            continue

        for mr_db_id in model_reaction_db_ids:
            # add to the GeneReactionMatrix, if not already present
            found_gene_reaction_row = (
                session.scalars(
                    select(func.count(GeneReactionMatrix.id))
                    .filter(GeneReactionMatrix.model_gene == model_gene_db)
                    .filter(GeneReactionMatrix.model_reaction_id == mr_db_id)
                ).first()
                > 0
            )
            if not found_gene_reaction_row:
                new_object = GeneReactionMatrix(
                    model_gene=model_gene_db, model_reaction_id=mr_db_id
                )
                session.add(new_object)

            # update the gene_reaction_rule if the gene id has changed
            if gene_id != gene_db.bigg_id:
                mr = session.get(ModelReaction, mr_db_id)
                new_rule = _replace_gene_str(
                    mr.gene_reaction_rule, gene_id, gene_db.bigg_id
                )
                session.execute(
                    update(ModelReaction),
                    [{"id": mr_db_id, "gene_reaction_rule": new_rule}],
                )
    session.commit()


@timing
def load_model_count(session, model_db_id):
    metabolite_count = session.scalars(
        select(func.count(ModelCompartmentalizedComponent.id)).filter(
            ModelCompartmentalizedComponent.model_id == model_db_id
        )
    ).first()
    reaction_count = session.scalars(
        select(func.count(ModelReaction.id)).filter(
            ModelReaction.model_id == model_db_id
        )
    ).first()
    gene_count = session.scalars(
        select(func.count(ModelGene.id)).filter(ModelGene.model_id == model_db_id)
    ).first()
    mc = ModelCount(
        model_id=model_db_id,
        gene_count=gene_count,
        metabolite_count=metabolite_count,
        reaction_count=reaction_count,
    )
    session.add(mc)
