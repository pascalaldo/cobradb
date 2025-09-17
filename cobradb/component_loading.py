# -*- coding: utf-8 -*-

from cobradb.models import *
from cobradb import settings
from cobradb.util import scrub_gene_id, get_or_create_data_source, get_or_create, timing

import sys, os, math, re
from os.path import basename
from warnings import warn
from sqlalchemy import select, text, or_, and_, func
import logging
import six
import itertools as it
import gzip
import time


class BadGenomeError(Exception):
    pass


def _load_gb_file(genbank_file_handle):
    """Load the Genbank file.

    Arguments
    ---------

    genbank_file_handle: The handle to the genbank file.

    """
    # imports
    from Bio import SeqIO

    # load the genbank file
    logging.debug("Loading file: %s" % genbank_file_handle.name)
    try:
        gb_file = SeqIO.parse(genbank_file_handle, "gb")
    except IOError:
        raise BadGenomeError("File '%s' not found" % genbank_file_handle.name)
    except Exception as e:
        raise BadGenomeError(
            'BioPython failed to parse %s with error "%s"'
            % (genbank_file_handle.name, e.message)
        )
    return gb_file


def get_genbank_accessions(genbank_filepath, fast=False):
    """Load the file and return the NCBI Accession and Assembly IDs (if available).

    Returns a dictionary of accessions with keys: 'ncbi_accession',
    'ncbi_assembly', 'ncbi_bioproject'.

    Arguments
    ---------

    genbank_filepath: The path to the genbank file.

    fast: If True, then only look in the first 100 lines. Faster because we do
    not load the whole file.

    """
    out = {"ncbi_assembly": None, "ncbi_accession": None, "ncbi_bioproject": None}

    if fast:
        # try to find the BioProject ID in the first 100 lines. Otherwise, use
        # the full SeqIO.read
        line_limit = 100
        regex_dict = {
            k: re.compile(v)
            for k, v in six.iteritems(
                {
                    "ncbi_accession": r"VERSION\s+([\w.-]+)[^\w.-]",
                    "ncbi_assembly": r"Assembly:\s*([\w.-]+)[^\w.-]",
                    "ncbi_bioproject": r"BioProject:\s*([\w.-]+)[^\w.-]",
                }
            )
        }
        with open(genbank_filepath, "r") as f:
            for i, line in enumerate(f.readlines()):
                for key, regex in six.iteritems(regex_dict):
                    match = regex.search(line)
                    if match is not None:
                        out[key] = match.group(1)
                if i > line_limit:
                    break
    else:
        # load the genbank file
        with open(genbank_filepath, "r") as f:
            gb_file = _load_gb_file(f)
            out["ncbi_accession"] = gb_file.id
            for value in it.chain.from_iterable(x.split() for x in gb_file.dbxrefs):
                if "Assembly" in value:
                    out["ncbi_assembly"] = value.split(":")[1]
                if "BioProject" in value:
                    out["ncbi_bioproject"] = value.split(":")[1]

    return out


def load_gene_synonym(session, gene_db, synonym, data_source_id):
    """Load the synonym for this gene from the given genome."""
    data_source_id = get_or_create_data_source(session, data_source_id)
    synonym_db, _ = get_or_create(
        session,
        Synonym,
        type="gene",
        ome_id=gene_db.id,
        synonym=synonym,
        data_source_id=data_source_id,
    )
    return synonym_db.id


def collect_gene_synonym(synonym_collection, gene_db, synonym, data_source_id):
    if not data_source_id in synonym_collection:
        synonym_collection[data_source_id] = set()
    synonym_collection[data_source_id].add((synonym, gene_db))


def insert_synonyms(session, synonym_collection):
    from sqlalchemy.dialects.postgresql import insert

    # syns = []
    for data_source_id, synonyms in synonym_collection.items():
        data_source_db = get_or_create_data_source(session, data_source_id)
        for syn, gene_db in synonyms:
            synonym = Synonym(
                type="gene",
                ome_id=gene_db.id,
                synonym=syn,
            )
            data_source_db.synonyms.append(synonym)
    # stmt = insert(Synonym).values(syns)
    # stmt = stmt.on_conflict_do_nothing(
    #     index_elements=[
    #         Synonym.type,
    #         Synonym.ome_id,
    #         Synonym.synonym,
    #         Synonym.data_source,
    #     ]
    # )
    # session.execute(stmt)
    # session.commit()


def _get_qual(feat, name, get_first=False):
    """Get a non-null attribute from the feature."""
    try:
        qual = feat.qualifiers[name]
    except KeyError:
        if get_first:
            return None
        else:
            return []

    def nonempty_str(s):
        s = s.strip()
        return None if s == "" else s

    if get_first:
        return nonempty_str(qual[0])
    else:
        return [y for y in (nonempty_str(x) for x in qual) if y is not None]


def _get_geneid(feature):
    """Get the value of GeneID from db_xref, or else return None."""
    db_xref = _get_qual(feature, "db_xref")
    if not db_xref:
        return None
    for ref in db_xref:
        splitrefs = [x.strip() for x in ref.split(":")]
        if len(splitrefs) == 2 and splitrefs[0].lower() == "geneid":
            return splitrefs[1]
    return None


@timing
def load_assembly(assembly_id, assembly_path, chromosome_accessions, session):
    """Load the genome and chromosomes."""

    # check that the genome doesn't already exist
    if (
        session.scalar(
            select(func.count(Genome.id))
            .filter(Genome.accession_type == "ncbi_assembly")
            .filter(Genome.accession_value == assembly_id)
        )
        > 0
    ):
        raise AlreadyLoadedError(f"Assembly {assembly_id} already loaded")

    logging.debug(f"Adding new genome: {assembly_id}")
    genome_db = Genome(accession_type="ncbi_assembly", accession_value=assembly_id)
    session.add(genome_db)
    # session.commit()

    if assembly_path is None:
        logging.warning(f"No assembly file found for {assembly_id}")
        return
    if chromosome_accessions is None:
        chromosome_accessions = []
    with gzip.open(assembly_path, "rt") as f:
        gb_file = _load_gb_file(f)
        for i, record in enumerate(gb_file):
            if not record.id in chromosome_accessions:
                logging.warning(f"Skipping chromosome [{i+1} of ?] {record.id}")
                continue
            logging.info(f"Loading chromosome [{i+1} of ?] {record.id}")
            load_chromosome(record, genome_db, session)
    session.commit()


def first_uppercase(val):
    if val is None:
        return val
    return str(val[0]).upper()


@timing
def load_chromosome(record, genome_db, session):
    chromosome_db = None
    if genome_db.id is not None:
        chromosome_db = session.scalars(
            select(Chromosome)
            .filter(Chromosome.ncbi_accession == record.id)
            .filter(Chromosome.genome_id == genome_db.id)
            .limit(1)
        ).first()
    if not chromosome_db:
        logging.debug("Loading new chromosome: {}".format(record.id))
        chromosome_db = Chromosome(ncbi_accession=record.id)
        genome_db.chromosomes.append(chromosome_db)
    else:
        logging.debug("Chromosome already loaded: %s" % record.id)

    # update genome
    if genome_db.organism is None:
        logging.warning(f"Organism: {record.annotations['organism']}")
        genome_db.organism = record.annotations["organism"]

    bigg_id_warnings = 0
    duplicate_genes_warnings = 0
    warning_num = 5
    synonym_collection = {}
    added_gene_bigg_ids = {}
    for i, feature in enumerate(record.features):
        # update genome with the source information
        if genome_db.taxon_id is None and feature.type == "source":
            for ref in _get_qual(feature, "db_xref"):
                if "taxon" == ref.split(":")[0]:
                    genome_db.taxon_id = ref.split(":")[1]
                    break
            continue

        # only read in CDSs
        if feature.type != "CDS":
            continue

        # bigg_id required
        bigg_id = None
        gene_name = None
        refseq_name = None
        locus_tag = None

        # get bigg_id if possible from locus_tag or and GeneID
        found_tag = _get_qual(feature, "locus_tag", True)
        found_gene_id = _get_geneid(feature)
        if found_tag is not None:
            locus_tag = found_tag
            bigg_id = scrub_gene_id(found_tag)
        elif found_gene_id is not None:
            bigg_id = scrub_gene_id(found_gene_id)

        # get name
        found_name = _get_qual(feature, "gene", True)
        if found_name is not None:
            gene_name = found_name
            refseq_name = found_name

        # warn about no locus_tag / bigg_id
        if gene_name is not None and bigg_id is None:
            if bigg_id_warnings <= warning_num:
                msg = (
                    "No locus_tag for gene. Using Gene name as bigg_id: %s" % gene_name
                )
                if bigg_id_warnings == warning_num:
                    msg += " (Warnings limited to %d)" % warning_num
                logging.warning(msg)
                bigg_id_warnings += 1
            bigg_id = scrub_gene_id(gene_name)
            gene_name = bigg_id
        elif bigg_id is None:
            logging.warning(
                (
                    "No locus_tag or gene name for gene %d in chromosome "
                    "%s" % (i, chromosome_db.ncbi_accession)
                )
            )
            continue

        gene_db = None
        if bigg_id in added_gene_bigg_ids:
            gene_db = added_gene_bigg_ids[bigg_id]
        elif chromosome_db.id is not None:
            gene_db = session.scalars(
                select(Gene)
                .filter(Gene.bigg_id == bigg_id)
                .filter(Gene.chromosome_id == chromosome_db.id)
                .limit(1)
            ).first()

        if gene_db is None:
            # get the strand and positions
            strand = None
            if feature.location.strand == 1:
                strand = "+"
            elif feature.location.strand == -1:
                strand = "-"
            leftpos = int(feature.location.start) + 1
            rightpos = int(feature.location.end)

            dna_sequence = str(feature.extract(record.seq)).upper()
            protein_sequence = first_uppercase(
                feature.qualifiers.get("translation", None)
            )

            # finally, create the gene
            gene_db = Gene(
                bigg_id=bigg_id,
                locus_tag=locus_tag,
                name=gene_name,
                leftpos=leftpos,
                rightpos=rightpos,
                strand=strand,
                dna_sequence=dna_sequence,
                protein_sequence=protein_sequence,
                mapped_to_genbank=True,
            )
            chromosome_db.genome_regions.append(gene_db)
            added_gene_bigg_ids[bigg_id] = gene_db

        else:
            # warn about duplicate genes.
            #
            # TODO The only downside to loading CDS's this way is that the
            # leftpos and rightpos correspond to a CDS, not the whole gene. So
            # these need to be fixed eventually.
            if duplicate_genes_warnings <= warning_num:
                msg = "Duplicate genes %s on chromosome %s" % (
                    bigg_id,
                    chromosome_db.ncbi_accession,
                )
                if duplicate_genes_warnings == warning_num:
                    msg += " (Warnings limited to %d)" % warning_num
                logging.warning(msg)
                duplicate_genes_warnings += 1

        # synonym_collection = {}
        # load the synonyms for the gene
        if locus_tag is not None:
            collect_gene_synonym(
                synonym_collection, gene_db, locus_tag, "refseq_locus_tag"
            )

        if refseq_name is not None:
            collect_gene_synonym(
                synonym_collection, gene_db, refseq_name, "refseq_name"
            )

        for ref in _get_qual(feature, "gene_synonym"):
            synonyms = [x.strip() for x in ref.split(";")]
            for syn in synonyms:
                collect_gene_synonym(synonym_collection, gene_db, syn, "refseq_synonym")

        for ref in _get_qual(feature, "db_xref"):
            splitrefs = [x.strip() for x in ref.split(":")]
            if len(splitrefs) == 2:
                collect_gene_synonym(
                    synonym_collection, gene_db, splitrefs[1], splitrefs[0]
                )

        for ref in _get_qual(feature, "old_locus_tag"):
            for syn in [x.strip() for x in ref.split(";")]:
                collect_gene_synonym(
                    synonym_collection, gene_db, syn, "refseq_old_locus_tag"
                )

        for ref in _get_qual(feature, "note"):
            for value in [x.strip() for x in ref.split(";")]:
                sp = value.split(":")
                if len(sp) == 2 and sp[0] == "ORF_ID":
                    collect_gene_synonym(
                        synonym_collection, gene_db, sp[1], "refseq_orf_id"
                    )
    session.commit()
    insert_synonyms(session, synonym_collection)
    session.commit()
