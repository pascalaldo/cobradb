# -*- coding: utf-8 -*-

from typing import Any, Dict
from sqlalchemy.orm import Session
from cobradb.api.metabolites import (
    _create_reactive_part,
    _create_reference_compound,
)
from cobradb.api.reactions import _create_reference_reaction


from cobradb.util import timing

import re
import logging

from rdflib import Graph, Namespace
import re

RHEA = Namespace("http://rdf.rhea-db.org/")
RDF = Namespace("http://www.w3.org/2000/01/rdf-schema#")
CHEBI = Namespace("http://purl.obolibrary.org/obo/CHEBI_")
BIOPAX = Namespace("http://www.biopax.org/release/biopax-level3.owl#")
EC = Namespace("http://purl.uniprot.org/enzyme/")

CHEBI_URI_PATTERN = re.compile(r"https?://purl\.obolibrary\.org/obo/CHEBI_([0-9]+)")

EC_URI_PATTERN = re.compile(r"https?://purl\.uniprot\.org/enzyme/([0-9\.\-n]+)")
GO_URI_PATTERN = re.compile(r"https?://purl\.obolibrary\.org/obo/GO_([0-9]+)")
KEGG_URI_PATTERN = re.compile(r"https?://identifiers\.org/kegg\.reaction/(R[0-9]+)")
METACYC_URI_PATTERN = re.compile(r"https?://identifiers\.org/biocyc/METACYC:(.+)")

ANNOTATION_PATTERNS = {
    "ec-code": EC_URI_PATTERN,
    "GO": GO_URI_PATTERN,
    "kegg.reaction": KEGG_URI_PATTERN,
    "metacyc.reaction": METACYC_URI_PATTERN,
}


@timing
def push_rhea_reference(session: Session, rhea_db: Dict[str, Any]) -> None:
    n_reactive_parts = len(rhea_db["reactive_parts"])
    for i, (rp_id, reactive_part) in enumerate(rhea_db["reactive_parts"].items()):
        print(f"RHEA: Reactive Part {i+1}/{n_reactive_parts}")
        reactive_part_db = _create_reactive_part(
            session,
            rp_id,
            name=reactive_part["name"],
            html_name=reactive_part["html_name"],
            formula=reactive_part["formula"],
            charge=reactive_part["charge"],
        )
    session.commit()
    session.close()

    n_compounds = len(rhea_db["compounds"])
    for i, (cp_id, compound) in enumerate(rhea_db["compounds"].items()):
        print(f"RHEA: Compound {i+1}/{n_compounds}")
        compound_db = _create_reference_compound(
            session,
            cp_id,
            name=compound["name"],
            compound_type=compound["type"],
            formula=compound.get("formula"),
            charge=compound.get("charge"),
            html_name=compound.get("html_name"),
            reactive_parts=compound.get("reactive_parts"),
        )
    session.commit()
    session.close()

    n_reactions = len(rhea_db["reactions"])
    for i, (rx_id, reaction) in enumerate(rhea_db["reactions"].items()):
        print(f"RHEA: Reaction {i+1}/{n_reactions}")
        reaction_db = _create_reference_reaction(
            session,
            rx_id,
            equation=reaction["equation"],
            participants=reaction["participants"],
            annotations=reaction.get("annotations"),
        )
    session.commit()
    session.close()


@timing
def load_rhea(session: Session, rhea_filepath):
    logging.warning("Loading RHEA reference data")

    graph = load_rhea_rdf(rhea_filepath)

    # graph.serialize(destination=f"/chebi/rhea.ttl")
    # graph.serialize(destination=f"/chebi/rhea.nt", format="nt")

    logging.warning("Extracting reactions")
    rhea_db = extract_reactions(graph)

    logging.warning("Pushing reactions to DB")
    push_rhea_reference(session, rhea_db)


# Load RDF file
def load_rhea_rdf(file_path):
    g = Graph()
    # g.parse(file_path, format="xml")  # RDF/XML format
    g.parse(file_path)
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


def add_annotation(uri, annotations):
    for namespace, pattern in ANNOTATION_PATTERNS.items():
        m = pattern.fullmatch(uri)
        if m:
            identifier = m.group(1)
            if namespace in annotations:
                annotations[namespace].append(identifier)
            else:
                annotations[namespace] = [identifier]
            break


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
        # ec = [get_ec_from_uri(x) for x in graph.objects(reaction, RHEA.ec)]

        all_reaction_variants = [reaction]
        all_reaction_variants.extend(graph.objects(reaction, RHEA.directionalReaction))
        all_reaction_variants.extend(
            graph.objects(reaction, RHEA.bidirectionalReaction)
        )

        annotations = {"rhea": []}
        for r in all_reaction_variants:
            annotations["rhea"].append(
                str(
                    get_single_object(graph, r, RHEA.accession).toPython()
                ).removeprefix("RHEA:")
            )
            for x in graph.objects(r, RDF.seeAlso):
                add_annotation(x, annotations)
            for x in graph.objects(r, RHEA.ec):
                add_annotation(x, annotations)

        reaction_info = {
            "accession": accession,
            "equation": equation,
            # "direction": direction,
            # "status": status,
            "annotations": annotations,
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
