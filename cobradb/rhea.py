# -*- coding: utf-8 -*-

from sqlalchemy import select
from cobradb.chebi import ChebiEntity
from cobradb.data_sources import get_data_source_id
from cobradb.metabolites import get_or_create_small_molecule_reference
from cobradb.models import *
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


def create_hierarchical_conversion_reaction(lhs_chebi, rhs_chebi, session):
    if lhs_chebi.get_charge() != rhs_chebi.get_charge():
        return
    reaction_bigg_id = f"HIER:{lhs_chebi.get_id()}_{rhs_chebi.get_id()}"
    reference_reaction_db = session.scalars(
        select(ReferenceReaction)
        .filter(ReferenceReaction.bigg_id == reaction_bigg_id)
        .limit(1)
    ).first()
    if reference_reaction_db is not None:
        return
    reference_reaction_db = ReferenceReaction(
        bigg_id=reaction_bigg_id,
        name=f"Conversion of {lhs_chebi.get_id()} to {rhs_chebi.get_id()}, because one is an instance of the other.",
        equation=f"{lhs_chebi.get_name()} = {rhs_chebi.get_name()}",
    )
    lhs_compound = session.scalars(
        select(ReferenceCompound)
        .filter(ReferenceCompound.bigg_id == lhs_chebi.get_id())
        .limit(1)
    ).first()
    lhs_part = ReferenceReactionParticipant(
        compound=lhs_compound,
        side="L",
        coefficient="1",
        compartment="0",
    )
    reference_reaction_db.reaction_participants.append(lhs_part)
    rhs_compound = session.scalars(
        select(ReferenceCompound)
        .filter(ReferenceCompound.bigg_id == rhs_chebi.get_id())
        .limit(1)
    ).first()
    rhs_part = ReferenceReactionParticipant(
        compound=rhs_compound,
        side="R",
        coefficient="1",
        compartment="0",
    )
    reference_reaction_db.reaction_participants.append(rhs_part)
    reference_reaction_db.update_hash()
    session.add(reference_reaction_db)


CHEBI_PROTONATION_RELATIONS = ["is_conjugate_acid_of", "is_conjugate_base_of"]


def add_reference_exchange_reactions(compound_db, session):
    reaction_bigg_id = f"EX:{compound_db.bigg_id}"
    reference_reaction_db = session.scalars(
        select(ReferenceReaction)
        .filter(ReferenceReaction.bigg_id == reaction_bigg_id)
        .limit(1)
    ).first()
    if reference_reaction_db is not None:
        return
    compound_name = compound_db.name if compound_db.name else compound_db.id
    reference_reaction_db = ReferenceReaction(
        bigg_id=reaction_bigg_id,
        name=f"Exchange of {compound_name}.",
        equation=f"{compound_name} = ∅",
    )
    lhs_part = ReferenceReactionParticipant(
        compound=compound_db,
        side="L",
        coefficient="1",
        compartment="0",
    )
    reference_reaction_db.reaction_participants.append(lhs_part)
    reference_reaction_db.update_hash()
    session.add(reference_reaction_db)
    session.commit()


def add_reference_conversion_reactions(compound_db, session):
    if not compound_db.bigg_id.startswith("CHEBI:"):
        return
    chebi_entity = ChebiEntity(compound_db.bigg_id)
    main_charge = chebi_entity.get_charge()

    outgoing_rel = chebi_entity.get_outgoings()
    for rel in outgoing_rel:
        if rel._Relation__typ in CHEBI_PROTONATION_RELATIONS:
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            rel_chebi_db = session.scalars(
                select(ReferenceCompound)
                .filter(ReferenceCompound.bigg_id == rel_chebi)
                .limit(1)
            ).first()
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
            reaction_bigg_id = f"PROT:{lhs_chebi.get_id()}_{rhs_chebi.get_id()}"
            reference_reaction_db = session.scalars(
                select(ReferenceReaction)
                .filter(ReferenceReaction.bigg_id == reaction_bigg_id)
                .limit(1)
            ).first()
            if reference_reaction_db is not None:
                continue
            reference_reaction_db = ReferenceReaction(
                bigg_id=reaction_bigg_id,
                name=f"Protonation of {lhs_chebi.get_id()} to {rhs_chebi.get_id()}.",
                equation=f"{lhs_chebi.get_name()} + {n_h_plus} H+ = {rhs_chebi.get_name()}",
            )
            lhs_compound = session.scalars(
                select(ReferenceCompound)
                .filter(ReferenceCompound.bigg_id == lhs_chebi.get_id())
                .limit(1)
            ).first()
            lhs_part = ReferenceReactionParticipant(
                compound=lhs_compound,
                side="L",
                coefficient="1",
                compartment="0",
            )
            reference_reaction_db.reaction_participants.append(lhs_part)
            rhs_compound = session.scalars(
                select(ReferenceCompound)
                .filter(ReferenceCompound.bigg_id == rhs_chebi.get_id())
                .limit(1)
            ).first()
            rhs_part = ReferenceReactionParticipant(
                compound=rhs_compound,
                side="R",
                coefficient="1",
                compartment="0",
            )
            reference_reaction_db.reaction_participants.append(rhs_part)
            proton_compound = session.scalars(
                select(ReferenceCompound)
                .filter(ReferenceCompound.bigg_id == "CHEBI:15378")
                .limit(1)
            ).first()
            proton_part = ReferenceReactionParticipant(
                compound=proton_compound,
                side="L",
                coefficient="1",
                compartment="0",
            )
            reference_reaction_db.reaction_participants.append(proton_part)
            reference_reaction_db.update_hash()
            session.add(reference_reaction_db)
            session.commit()
        elif rel._Relation__typ == "is_a":
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            rel_chebi_db = session.scalars(
                select(ReferenceCompound)
                .filter(ReferenceCompound.bigg_id == rel_chebi)
                .limit(1)
            ).first()
            rel_chebi_entity = ChebiEntity(rel_chebi)
            if rel_chebi_db is None:
                # Look for relations where one step is skipped in our DB.
                for rel_2 in rel_chebi_entity.get_outgoings():
                    if rel_2._Relation__typ == "is_a":
                        rel_2_chebi = f"CHEBI:{rel_2._Relation__target_chebi_id}"
                        rel_chebi_db = session.scalars(
                            select(ReferenceCompound)
                            .filter(ReferenceCompound.bigg_id == rel_2_chebi)
                            .limit(1)
                        ).first()
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
            rel_chebi_db = session.scalars(
                select(ReferenceCompound)
                .filter(ReferenceCompound.bigg_id == rel_chebi)
                .limit(1)
            ).first()
            rel_chebi_entity = ChebiEntity(rel_chebi)
            if rel_chebi_db is None:
                # Look for relations where one step is skipped in our DB.
                for rel_2 in rel_chebi_entity.get_incomings():
                    if rel_2._Relation__typ == "is_a":
                        rel_2_chebi = f"CHEBI:{rel_2._Relation__target_chebi_id}"
                        rel_chebi_db = session.scalars(
                            select(ReferenceCompound)
                            .filter(ReferenceCompound.bigg_id == rel_2_chebi)
                            .limit(1)
                        ).first()
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


def add_reaction_annotations(reaction_db, reaction, session):
    default_data_source_id = get_data_source_id("rhea", session)
    if default_data_source_id is None:
        print("Could not find RHEA data source.")
        return
    annotation_db = Annotation(
        bigg_id=reaction["accession"],
        default_data_source_id=default_data_source_id,
        type="rhea",
    )
    mapping = ReferenceReactionAnnotationMapping(
        reference_reaction=reaction_db,
    )
    annotation_db.reference_reaction_mappings.append(mapping)
    annotations = reaction.get("annotations", {})
    for namespace, identifiers in annotations.items():
        data_source_id = get_data_source_id(namespace, session)
        if not data_source_id:
            print(f"Unknown data source: {namespace}")
            continue
        for identifier in identifiers:
            link = AnnotationLink(
                data_source_id=data_source_id,
                identifier=identifier,
            )
            annotation_db.links.append(link)
    session.add(annotation_db)


@timing
def push_rhea_reference(rhea_db, session):
    n_reactive_parts = len(rhea_db["reactive_parts"])
    for i, (rp_id, reactive_part) in enumerate(rhea_db["reactive_parts"].items()):
        print(f"RHEA: Reactive Part {i+1}/{n_reactive_parts}")
        if rp_id.startswith("CHEBI:"):
            _existed, reactive_part_db = get_or_create_small_molecule_reference(
                rp_id, session, cpd_cls=ReferenceReactivePart
            )
        else:
            reactive_part_db = ReferenceReactivePart(
                bigg_id=rp_id,
                name=reactive_part["name"],
                html_name=reactive_part["html_name"],
                formula=reactive_part["formula"],
                charge=reactive_part["charge"],
            )
            session.add(reactive_part_db)
    session.commit()
    session.close()

    n_compounds = len(rhea_db["compounds"])
    for i, (cp_id, compound) in enumerate(rhea_db["compounds"].items()):
        print(f"RHEA: Compound {i+1}/{n_compounds}")
        if cp_id.startswith("CHEBI:"):
            _existed, compound_db = get_or_create_small_molecule_reference(
                cp_id, session
            )
        else:
            compound_db = ReferenceCompound(
                bigg_id=cp_id,
                name=compound["name"],
                html_name=compound["html_name"],
                formula=compound.get("formula"),
                charge=compound.get("charge"),
                compound_type=compound["type"],
            )
            for reactive_part in compound.get("reactive_parts", []):
                reactive_part_db = session.scalars(
                    select(ReferenceReactivePart)
                    .filter(ReferenceReactivePart.bigg_id == reactive_part)
                    .limit(1)
                ).first()
                reactive_part_matrix_db = ReferenceReactivePartMatrix(
                    reactive_part=reactive_part_db,
                )
                compound_db.reactive_part_matrix.append(reactive_part_matrix_db)
            session.add(compound_db)
        session.commit()
        add_reference_conversion_reactions(compound_db, session)
        add_reference_exchange_reactions(compound_db, session)
    session.commit()
    session.close()

    n_reactions = len(rhea_db["reactions"])
    for i, (rx_id, reaction) in enumerate(rhea_db["reactions"].items()):
        print(f"RHEA: Reaction {i+1}/{n_reactions}")
        reaction_db = ReferenceReaction(
            bigg_id=rx_id,
            equation=reaction["equation"],
        )
        for side_n, coefficient, compound_id, compartment in reaction["participants"]:
            compound_db = session.scalars(
                select(ReferenceCompound)
                .filter(ReferenceCompound.bigg_id == compound_id)
                .limit(1)
            ).first()
            participant_db = ReferenceReactionParticipant(
                compound=compound_db,
                side="L" if side_n < 0 else "R",
                coefficient=coefficient,
                compartment=compartment,
            )
            reaction_db.reaction_participants.append(participant_db)
        if len(reaction_db.reaction_participants) > 0:
            reaction_db.update_hash()
            session.add(reaction_db)

            add_reaction_annotations(reaction_db, reaction, session)

    session.commit()
    session.close()


@timing
def load_rhea(rhea_filepath, session):
    logging.warning("Loading RHEA reference data")

    graph = load_rhea_rdf(rhea_filepath)

    # graph.serialize(destination=f"/chebi/rhea.ttl")
    # graph.serialize(destination=f"/chebi/rhea.nt", format="nt")

    logging.warning("Extracting reactions")
    rhea_db = extract_reactions(graph)

    logging.warning("Pushing reactions to DB")
    push_rhea_reference(rhea_db, session)


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
