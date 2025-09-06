# -*- coding: utf-8 -*-

from cobradb.metabolites import get_or_create_small_molecule_reference
from cobradb.models import *
from cobradb import settings
from cobradb.util import timing

from sqlalchemy import func
import re
import logging

import time

from rdflib import Graph, Namespace
from rdflib.namespace import RDF
from libchebipy import ChebiEntity
from pprint import pprint
import re

RHEA = Namespace("http://rdf.rhea-db.org/")
RDF = Namespace("http://www.w3.org/2000/01/rdf-schema#")
CHEBI = Namespace("http://purl.obolibrary.org/obo/CHEBI_")
BIOPAX = Namespace("http://www.biopax.org/release/biopax-level3.owl#")
EC = Namespace("http://purl.uniprot.org/enzyme/")

CHEBI_URI_PATTERN = re.compile(r"https?://purl\.obolibrary\.org/obo/CHEBI_([0-9]+)")
EC_URI_PATTERN = re.compile(r"https?://purl\.uniprot\.org/enzyme/([0-9\.\-n]+)")


def create_hierarchical_conversion_reaction(lhs_chebi, rhs_chebi, session):
    if lhs_chebi.get_charge() != rhs_chebi.get_charge():
        return
    reaction_id = f"HIER:{lhs_chebi.get_id()}_{rhs_chebi.get_id()}"
    reference_reaction_db = (
        session.query(ReferenceReaction)
        .filter(ReferenceReaction.id == reaction_id)
        .first()
    )
    if reference_reaction_db is not None:
        return
    reference_reaction_db = ReferenceReaction(
        id=reaction_id,
        name=f"Conversion of {lhs_chebi.get_id()} to {rhs_chebi.get_id()}, because one is an instance of the other.",
        equation=f"{lhs_chebi.get_name()} = {rhs_chebi.get_name()}",
    )
    session.add(reference_reaction_db)
    lhs_part = ReferenceReactionParticipant(
        reaction_id=reaction_id,
        compound_id=lhs_chebi.get_id(),
        side="L",
        coefficient="1",
        compartment="0",
    )
    session.add(lhs_part)
    rhs_part = ReferenceReactionParticipant(
        reaction_id=reaction_id,
        compound_id=rhs_chebi.get_id(),
        side="R",
        coefficient="1",
        compartment="0",
    )
    session.add(rhs_part)


CHEBI_PROTONATION_RELATIONS = ["is_conjugate_acid_of", "is_conjugate_base_of"]


def add_reference_exchange_reactions(compound_db, session):
    reaction_id = f"EX:{compound_db.id}"
    reference_reaction_db = (
        session.query(ReferenceReaction)
        .filter(ReferenceReaction.id == reaction_id)
        .first()
    )
    if reference_reaction_db is not None:
        return
    compound_name = compound_db.name if compound_db.name else compound_db.id
    reference_reaction_db = ReferenceReaction(
        id=reaction_id,
        name=f"Exchange of {compound_name}.",
        equation=f"{compound_name} = ∅",
    )
    session.add(reference_reaction_db)
    lhs_part = ReferenceReactionParticipant(
        reaction_id=reaction_id,
        compound_id=compound_db.id,
        side="L",
        coefficient="1",
        compartment="0",
    )
    session.add(lhs_part)
    session.commit()


def add_reference_conversion_reactions(compound_db, session):
    if not compound_db.id.startswith("CHEBI:"):
        return
    chebi_entity = ChebiEntity(compound_db.id)
    main_charge = chebi_entity.get_charge()

    outgoing_rel = chebi_entity.get_outgoings()
    for rel in outgoing_rel:
        if rel._Relation__typ in CHEBI_PROTONATION_RELATIONS:
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            rel_chebi_db = (
                session.query(ReferenceCompound)
                .filter(ReferenceCompound.id == rel_chebi)
                .first()
            )
            if rel_chebi_db is None:
                continue
            rel_chebi_entity = ChebiEntity(rel_chebi)
            rel_charge = rel_chebi_entity.get_charge()
            if main_charge == rel_charge:
                print("Same charge")
                continue
            elif main_charge < rel_charge:
                lhs_chebi = chebi_entity
                rhs_chebi = rel_chebi_entity
            else:
                lhs_chebi = rel_chebi_entity
                rhs_chebi = chebi_entity
            n_h_plus = rhs_chebi.get_charge() - lhs_chebi.get_charge()
            reaction_id = f"PROT:{lhs_chebi.get_id()}_{rhs_chebi.get_id()}"
            reference_reaction_db = (
                session.query(ReferenceReaction)
                .filter(ReferenceReaction.id == reaction_id)
                .first()
            )
            if reference_reaction_db is not None:
                continue
            reference_reaction_db = ReferenceReaction(
                id=reaction_id,
                name=f"Protonation of {lhs_chebi.get_id()} to {rhs_chebi.get_id()}.",
                equation=f"{lhs_chebi.get_name()} + {n_h_plus} H+ = {rhs_chebi.get_name()}",
            )
            session.add(reference_reaction_db)
            lhs_part = ReferenceReactionParticipant(
                reaction_id=reaction_id,
                compound_id=lhs_chebi.get_id(),
                side="L",
                coefficient="1",
                compartment="0",
            )
            session.add(lhs_part)
            rhs_part = ReferenceReactionParticipant(
                reaction_id=reaction_id,
                compound_id=rhs_chebi.get_id(),
                side="R",
                coefficient="1",
                compartment="0",
            )
            session.add(rhs_part)
            proton_part = ReferenceReactionParticipant(
                reaction_id=reaction_id,
                compound_id="CHEBI:15378",
                side="L",
                coefficient="1",
                compartment="0",
            )
            session.add(proton_part)
            session.commit()
        elif rel._Relation__typ == "is_a":
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            rel_chebi_db = (
                session.query(ReferenceCompound)
                .filter(ReferenceCompound.id == rel_chebi)
                .first()
            )
            rel_chebi_entity = ChebiEntity(rel_chebi)
            if rel_chebi_db is None:
                # Look for relations where one step is skipped in our DB.
                for rel_2 in rel_chebi_entity.get_outgoings():
                    if rel_2._Relation__typ == "is_a":
                        rel_2_chebi = f"CHEBI:{rel_2._Relation__target_chebi_id}"
                        rel_chebi_db = (
                            session.query(ReferenceCompound)
                            .filter(ReferenceCompound.id == rel_2_chebi)
                            .first()
                        )
                        if rel_chebi_db is None:
                            continue
                        rel_chebi_entity = ChebiEntity(rel_2_chebi)
                        create_hierarchical_conversion_reaction(
                            chebi_entity, rel_chebi_entity, session
                        )
                session.commit()
                continue
            create_hierarchical_conversion_reaction(
                chebi_entity, rel_chebi_entity, session
            )
            session.commit()

    incoming_rel = chebi_entity.get_incomings()
    for rel in incoming_rel:
        if rel._Relation__typ == "is_a":
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            rel_chebi_db = (
                session.query(ReferenceCompound)
                .filter(ReferenceCompound.id == rel_chebi)
                .first()
            )
            rel_chebi_entity = ChebiEntity(rel_chebi)
            if rel_chebi_db is None:
                # Look for relations where one step is skipped in our DB.
                for rel_2 in rel_chebi_entity.get_incomings():
                    if rel_2._Relation__typ == "is_a":
                        rel_2_chebi = f"CHEBI:{rel_2._Relation__target_chebi_id}"
                        rel_chebi_db = (
                            session.query(ReferenceCompound)
                            .filter(ReferenceCompound.id == rel_2_chebi)
                            .first()
                        )
                        if rel_chebi_db is None:
                            continue
                        rel_chebi_entity = ChebiEntity(rel_2_chebi)
                        create_hierarchical_conversion_reaction(
                            rel_chebi_entity, chebi_entity, session
                        )
                session.commit()
                continue
            create_hierarchical_conversion_reaction(
                rel_chebi_entity, chebi_entity, session
            )
            session.commit()


@timing
def push_rhea_reference(rhea_db, session):
    for rp_id, reactive_part in rhea_db["reactive_parts"].items():
        if rp_id.startswith("CHEBI:"):
            get_or_create_small_molecule_reference(
                rp_id, session, cpd_cls=ReferenceReactivePart
            )
        else:
            reactive_part_db = ReferenceReactivePart(
                id=rp_id,
                name=reactive_part["name"],
                html_name=reactive_part["html_name"],
                formula=reactive_part["formula"],
                charge=reactive_part["charge"],
            )
            session.add(reactive_part_db)
    session.commit()
    for cp_id, compound in rhea_db["compounds"].items():
        if cp_id.startswith("CHEBI:"):
            _existed, compound_db = get_or_create_small_molecule_reference(
                cp_id, session
            )
        else:
            compound_db = ReferenceCompound(
                id=cp_id,
                name=compound["name"],
                html_name=compound["html_name"],
                formula=compound.get("formula"),
                charge=compound.get("charge"),
                compound_type=compound["type"],
            )
            session.add(compound_db)
            for reactive_part in compound.get("reactive_parts", []):
                reactive_part_matrix_db = ReferenceReactivePartMatrix(
                    compound_id=cp_id,
                    reactive_part_id=reactive_part,
                )
                session.add(reactive_part_matrix_db)
        session.commit()
        add_reference_conversion_reactions(compound_db, session)
        add_reference_exchange_reactions(compound_db, session)
    session.commit()
    for rx_id, reaction in rhea_db["reactions"].items():
        reaction_db = ReferenceReaction(
            id=rx_id,
            equation=reaction["equation"],
        )
        session.add(reaction_db)
        for side_n, coefficient, compound_id, compartment in reaction["participants"]:
            participant_db = ReferenceReactionParticipant(
                reaction_id=rx_id,
                compound_id=compound_id,
                side="L" if side_n < 0 else "R",
                coefficient=coefficient,
                compartment=compartment,
            )
            session.add(participant_db)
    session.commit()


@timing
def load_rhea(rhea_filepath, session):
    logging.debug("Loading RHEA reference data")

    graph = load_rhea_rdf(rhea_filepath)
    rhea_db = extract_reactions(graph)

    push_rhea_reference(rhea_db, session)


# Load RDF file
def load_rhea_rdf(file_path):
    g = Graph()
    g.parse(file_path, format="xml")  # RDF/XML format
    return g


def get_single_object(graph, subject, predicate):
    o = graph.value(subject=subject, predicate=predicate, default=None, any=False)
    if o is None:
        raise ValueError(f"No object found for {subject}, {predicate}")
    return o


def get_single_object_or_default(graph, subject, predicate, default=None):
    o = graph.value(subject=subject, predicate=predicate, default=None, any=False)
    if o is None:
        return default
    else:
        return o.toPython()


def get_chebi_from_uri(uri):
    m = CHEBI_URI_PATTERN.fullmatch(uri)
    if not m:
        raise ValueError(f"Could not parse CHEBI URI: '{uri}'")
    return f"CHEBI:{m.group(1)}"


def get_ec_from_uri(uri):
    m = EC_URI_PATTERN.fullmatch(uri)
    if not m:
        raise ValueError(f"Could not parse EC URI: '{uri}'")
    return m.group(1)


def get_small_molecule_info(graph, compound):
    chebi = get_chebi_from_uri(
        get_single_object(graph, compound, RHEA.chebi).toPython()
    )
    rhea_accession = get_single_object(graph, compound, RHEA.accession).toPython()
    name = get_single_object(graph, compound, RHEA.name).toPython()
    html_name = get_single_object(graph, compound, RHEA.htmlName).toPython()
    formula = get_single_object_or_default(graph, compound, RHEA.formula)
    charge = get_single_object_or_default(graph, compound, RHEA.charge)
    return {
        "type": "small_molecule",
        "rhea_accession": rhea_accession,
        "name": name,
        "chebi": chebi,
        "charge": charge,
        "formula": formula,
        "html_name": html_name,
    }


def get_reactive_part_info(graph, reactive_part):
    chebi = get_chebi_from_uri(
        get_single_object(graph, reactive_part, RHEA.chebi).toPython()
    )
    name = get_single_object(graph, reactive_part, RHEA.name).toPython()
    html_name = get_single_object(graph, reactive_part, RHEA.htmlName).toPython()
    formula = get_single_object_or_default(graph, reactive_part, RHEA.formula)
    charge = get_single_object_or_default(graph, reactive_part, RHEA.charge)

    return {
        "type": "reactive_part",
        "name": name,
        "chebi": chebi,
        "charge": charge,
        "formula": formula,
        "html_name": html_name,
    }


def get_generic_polypeptide_info(graph, compound):
    reactive_parts = [
        get_reactive_part_info(graph, reactive_part)
        for reactive_part in graph.objects(compound, RHEA.reactivePart)
    ]
    reactive_parts = {
        reactive_part["chebi"]: reactive_part for reactive_part in reactive_parts
    }

    rhea_accession = get_single_object(graph, compound, RHEA.accession).toPython()
    name = get_single_object(graph, compound, RHEA.name).toPython()
    html_name = get_single_object(graph, compound, RHEA.htmlName).toPython()
    charge = get_single_object_or_default(graph, compound, RHEA.charge)

    return {
        "type": "generic_polypeptide",
        "rhea_accession": rhea_accession,
        "name": name,
        "charge": charge,
        "html_name": html_name,
        "reactive_parts": list(reactive_parts.keys()),
    }, reactive_parts


def get_generic_polynucleotide_info(graph, compound):
    reactive_parts = [
        get_reactive_part_info(graph, reactive_part)
        for reactive_part in graph.objects(compound, RHEA.reactivePart)
    ]
    reactive_parts = {
        reactive_part["chebi"]: reactive_part for reactive_part in reactive_parts
    }

    rhea_accession = get_single_object(graph, compound, RHEA.accession).toPython()
    name = get_single_object(graph, compound, RHEA.name).toPython()
    html_name = get_single_object(graph, compound, RHEA.htmlName).toPython()
    charge = get_single_object_or_default(graph, compound, RHEA.charge)

    return {
        "type": "generic_polynucleotide",
        "rhea_accession": rhea_accession,
        "name": name,
        "charge": charge,
        "html_name": html_name,
        "reactive_parts": list(reactive_parts.keys()),
    }, reactive_parts


def get_polymer_info(graph, compound):
    chebi = get_chebi_from_uri(
        get_single_object(graph, compound, RHEA.underlyingChebi).toPython()
    )
    rhea_accession = get_single_object(graph, compound, RHEA.accession).toPython()
    name = get_single_object(graph, compound, RHEA.name).toPython()
    html_name = get_single_object(graph, compound, RHEA.htmlName).toPython()
    name = get_single_object(graph, compound, RHEA.polymerizationIndex).toPython()
    formula = get_single_object_or_default(graph, compound, RHEA.formula)
    charge = get_single_object_or_default(graph, compound, RHEA.charge)
    return {
        "type": "polymer",
        "rhea_accession": rhea_accession,
        "name": name,
        "chebi": chebi,
        "charge": charge,
        "formula": formula,
        "html_name": html_name,
    }


def get_coefficients(graph):
    contains_variants = graph.subjects(RDF.subPropertyOf, RHEA.contains)
    return {
        variant: get_single_object(graph, variant, RHEA.coefficient).toPython()
        for variant in contains_variants
    }


def determine_participant_stoichiometry(graph, coefficients, side, participant):
    predicates = graph.predicates(side, participant, unique=True)
    for predicate in predicates:
        if predicate == RHEA.contains:
            continue
        return coefficients[predicate]
    raise ValueError("No coefficient found")


def determine_compartment(graph, participant):
    c = get_single_object_or_default(graph, participant, RHEA.location)
    if c is None:
        return "0"
    if c == "http://rdf.rhea-db.org/In":
        return "in"
    if c == "http://rdf.rhea-db.org/Out":
        return "out"
    raise ValueError("Unknown comparment")


# Extract and print reaction info
def extract_reactions(graph: Graph):
    coefficients = get_coefficients(graph)
    rhea_db = {"reactions": {}, "compounds": {}, "reactive_parts": {}}
    for reaction in graph.subjects(
        predicate=RDF.subClassOf,
        object=RHEA.Reaction,
    ):
        try:
            equation = get_single_object(graph, reaction, RHEA.equation).toPython()
        except ValueError:
            equation = None
        # direction = get_single_object(graph, reaction, RHEA.direction).split("#")[-1]
        # status = get_single_object(graph, reaction, RHEA.status)
        accession = get_single_object(graph, reaction, RHEA.accession).toPython()
        ec = [get_ec_from_uri(x) for x in graph.objects(reaction, RHEA.ec)]

        reaction_info = {
            "accession": accession,
            "equation": equation,
            # "direction": direction,
            # "status": status,
            "ec": ec,
        }
        participants = []

        for side in graph.objects(reaction, RHEA.side):
            side_letter = side.rsplit("_", maxsplit=1)[-1].upper()
            if side_letter == "R":
                side_n = 1
            elif side_letter == "L":
                side_n = -1
            else:
                raise ValueError(f"Unknown reaction side (not L/R): {side}")
            for participant in graph.objects(side, RHEA.contains):
                coefficient = determine_participant_stoichiometry(
                    graph, coefficients, side, participant
                )
                compartment = determine_compartment(graph, participant)
                compound = get_single_object(graph, participant, RHEA.compound)
                info, rp_info = {}, {}
                for compound_subclass in graph.objects(compound, RDF.subClassOf):
                    if compound_subclass == RHEA.SmallMolecule:
                        info = get_small_molecule_info(graph, compound)
                        break
                    elif compound_subclass == RHEA.GenericPolypeptide:
                        info, rp_info = get_generic_polypeptide_info(graph, compound)
                        break
                    elif compound_subclass == RHEA.Polymer:
                        info = get_polymer_info(graph, compound)
                        break
                    elif compound_subclass == RHEA.GenericPolynucleotide:
                        info, rp_info = get_generic_polynucleotide_info(graph, compound)
                        break
                    else:
                        if not str(compound_subclass).startswith(
                            "http://purl.obolibrary.org/obo/CHEBI"
                        ):
                            raise ValueError(
                                f"#### OTHER CLASS {compound_subclass} ###"
                            )
                rhea_db["compounds"][info["rhea_accession"]] = info
                rhea_db["reactive_parts"].update(rp_info)
                participants.append(
                    (side_n, coefficient, info["rhea_accession"], compartment)
                )

        reaction_info["participants"] = participants
        rhea_db["reactions"][accession] = reaction_info
    return rhea_db
