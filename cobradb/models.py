# -*- coding: utf-8 -*-

"""Module to implement ORM to the ome database"""

from cobradb.settings import db_connection_string

from sqlalchemy import (
    ForeignKey,
    Column,
    Integer,
    String,
    Float,
    Table,
    LargeBinary,
    Boolean,
    create_engine,
    MetaData,
    Enum,
    DateTime,
    UniqueConstraint,
)
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm.session import close_all_sessions

# Connect to postgres
engine = create_engine(db_connection_string)
Base = declarative_base()
metadata = MetaData()
Session = sessionmaker(bind=engine)

COL_ID_STR = String(100)
COL_NAME_STR = String(400)
COL_HTML_NAME_STR = String(1000)
COL_FORMULA_STR = String(1000)
COL_CHARGE_STR = String(100)
COL_COEFFICIENT_STR = String(100)
COL_COMPARTMENT_STR = String(10)
COL_EQUATION_STR = String(5000)

# Make the enums
_enum_l = [
    Enum(
        "component",
        "reaction",
        "gene",
        "compartmentalized_component",
        name="synonym_type",
    ),
    Enum("pmid", "doi", name="reference_type"),
    Enum(
        "model_reaction",
        "model_compartmentalized_component",
        "model_gene",
        name="old_id_synonym_type",
    ),
    Enum("is_version", name="is_version"),
    Enum("component", "reaction", name="deprecated_id_types"),
    Enum(
        "model_compartmentalized_component",
        "model_reaction",
        name="escher_map_matrix_type",
    ),
    Enum(
        "small_molecule",
        "generic_polypeptide",
        "generic_polynucleotide",
        "polymer",
        name="compound_type",
    ),
    Enum("L", "R", name="reaction_side"),
]
custom_enums = {x.name: x for x in _enum_l}

# ------------
# Exceptions
# ------------


class NotFoundError(Exception):
    pass


class AlreadyLoadedError(Exception):
    pass


# --------
# Tables
# --------


class DatabaseVersion(Base):
    __tablename__ = "database_version"

    is_version = Column(custom_enums["is_version"], primary_key=True)
    date_time = Column(DateTime, nullable=False)

    __table_args__ = (UniqueConstraint("is_version"),)

    def __init__(self, date_time):
        self.is_version = "is_version"
        self.date_time = date_time


class Genome(Base):
    __tablename__ = "genome"

    id = Column(Integer, primary_key=True)
    accession_type = Column(String(200), nullable=False)
    accession_value = Column(String(200), nullable=False)
    organism = Column(String(200), nullable=True)
    taxon_id = Column(String(200), nullable=True)
    ncbi_assembly_id = Column(String(200), nullable=True)

    __table_args__ = (UniqueConstraint("accession_type", "accession_value"),)

    def __repr__(self):
        return (
            "<cobradb Genome(id={self.id}, accession_type={self.accession_type}, "
            "accession_value={self.accession_value})>".format(self=self)
        )


class Chromosome(Base):
    __tablename__ = "chromosome"

    id = Column(Integer, primary_key=True)
    ncbi_accession = Column(String(200))
    genome_id = Column(Integer, ForeignKey("genome.id"))

    __table_args__ = (UniqueConstraint("ncbi_accession", "genome_id"),)

    def __repr__(self):
        return "<cobradb Chromosome(id={self.id}, ncbi_accession={self.ncbi_accession}, genome_id={self.genome_id})>".format(
            self=self
        )


class GenomeRegion(Base):
    __tablename__ = "genome_region"
    id = Column(Integer, primary_key=True)
    chromosome_id = Column(Integer, ForeignKey("chromosome.id"))
    bigg_id = Column(String, nullable=False)
    leftpos = Column(Integer, nullable=True)
    rightpos = Column(Integer, nullable=True)
    strand = Column(String(1), nullable=True)
    type = Column(String(20))
    dna_sequence = Column(String, nullable=True)
    protein_sequence = Column(String, nullable=True)

    __table_args__ = (UniqueConstraint("bigg_id", "chromosome_id"),)

    __mapper_args__ = {"polymorphic_identity": "genome_region", "polymorphic_on": type}

    def __repr__(self):
        return "<cobradb GenomeRegion(id={self.id}, leftpos={self.leftpos}, rightpos={self.rightpos})>".format(
            self=self
        )


class ReferenceCompound(Base):
    __tablename__ = "reference_compound"

    id = Column(COL_ID_STR, primary_key=True)
    name = Column(COL_NAME_STR, nullable=False)
    html_name = Column(COL_HTML_NAME_STR, nullable=True)
    compound_type = Column(custom_enums["compound_type"], nullable=False)
    charge = Column(COL_CHARGE_STR, nullable=True)
    formula = Column(COL_FORMULA_STR, nullable=True)

    def __repr__(self):
        return (
            "<cobradb ReferenceCompound(id={self.id}, name={self.name}, "
            "comound_type={self.compound_type})>".format(self=self)
        )


class UniversalComponent(Base):
    __tablename__ = "universal_component"

    id = Column(COL_ID_STR, primary_key=True)
    name = Column(COL_NAME_STR, nullable=True)
    model_specific = Column(Boolean, nullable=False, default=False)


class Component(Base):
    __tablename__ = "component"

    id = Column(COL_ID_STR, primary_key=True)
    universal_id = Column(COL_ID_STR, ForeignKey(UniversalComponent.id), nullable=False)
    name = Column(String, nullable=True)
    # type = Column(String(20))
    formula = Column(String, nullable=False)
    charge = Column(Integer, nullable=False)
    model_specific = Column(Boolean, nullable=False, default=False)

    # __mapper_args__ = {"polymorphic_identity": "component", "polymorphic_on": type}

    def __repr__(self):
        return f"Component ({self.id}): {self.name}"


class ComponentReferenceMapping(Base):
    __tablename__ = "component_reference_mapping"

    id = Column(Integer, primary_key=True)
    component_id = Column(COL_ID_STR, ForeignKey(Component.id), nullable=False)
    universal_id = Column(COL_ID_STR, ForeignKey(UniversalComponent.id), nullable=False)
    reference_id = Column(COL_ID_STR, ForeignKey(ReferenceCompound.id), nullable=False)
    reference_n = Column(Integer, nullable=True)


class UniversalComponentReferenceMapping(Base):
    __tablename__ = "univeral_component_reference_mapping"
    id = Column(COL_ID_STR, ForeignKey(UniversalComponent.id), primary_key=True)
    mapping_id = Column(
        Integer, ForeignKey(ComponentReferenceMapping.id), nullable=False
    )


class DataSource(Base):
    __tablename__ = "data_source"

    id = Column(Integer, primary_key=True)
    bigg_id = Column(COL_ID_STR, nullable=False)
    name = Column(String(100))
    url_prefix = Column(String)

    __table_args__ = (UniqueConstraint("bigg_id"),)

    def __repr__(self):
        return (
            "<cobradb DataSource(id={self.id}, bigg_id={self.bigg_id}, "
            "name={self.name}, url_prefix={self.url_prefix})>"
        ).format(self=self)


class Synonym(Base):
    __tablename__ = "synonym"
    id = Column(Integer, primary_key=True)
    ome_id = Column(COL_ID_STR)
    synonym = Column(COL_ID_STR)
    type = Column(custom_enums["synonym_type"])
    data_source_id = Column(Integer, ForeignKey("data_source.id", ondelete="CASCADE"))

    __table_args__ = (UniqueConstraint("ome_id", "synonym", "type", "data_source_id"),)

    def __repr__(self):
        return (
            '<cobradb Synonym(id=%d, synonym="%s", type="%s", ome_id=%d, data_source_id=%d)>'
            % (self.id, self.synonym, self.type, self.ome_id, self.data_source_id)
        )


class Publication(Base):
    __tablename__ = "publication"
    id = Column(Integer, primary_key=True)
    reference_type = Column(custom_enums["reference_type"])
    reference_id = Column(String)

    __table_args__ = (UniqueConstraint("reference_type", "reference_id"),)


class PublicationModel(Base):
    __tablename__ = "publication_model"
    model_id = Column(
        COL_ID_STR, ForeignKey("model.id", ondelete="CASCADE"), primary_key=True
    )
    publication_id = Column(
        Integer, ForeignKey("publication.id", ondelete="CASCADE"), primary_key=True
    )

    __table_args__ = (UniqueConstraint("model_id", "publication_id"),)


class OldIDSynonym(Base):
    __tablename__ = "old_id_model_synonym"
    id = Column(Integer, primary_key=True)
    type = Column(custom_enums["old_id_synonym_type"])
    synonym_id = Column(
        Integer, ForeignKey("synonym.id", ondelete="CASCADE"), nullable=False
    )
    ome_id = Column(Integer, nullable=False)

    __table_args__ = (UniqueConstraint("synonym_id", "ome_id"),)

    def __repr__(self):
        return '<cobradb OldIDSynonym(id=%d, type="%s", ome_id=%d, synonym_id=%d)>' % (
            self.id,
            self.type,
            self.ome_id,
            self.synonym_id,
        )


class GenomeRegionMap(Base):
    __tablename__ = "genome_region_map"

    genome_region_id_1 = Column(
        Integer, ForeignKey("genome_region.id"), primary_key=True
    )
    genome_region_id_2 = Column(
        Integer, ForeignKey("genome_region.id"), primary_key=True
    )
    distance = Column(Integer)

    __table_args__ = (UniqueConstraint("genome_region_id_1", "genome_region_id_2"),)

    def __repr__(self):
        return "GenomeRegionMap (%d <--> %d) distance:%d" % (
            self.genome_region_id_1,
            self.genome_region_id_2,
            self.distance,
        )


class DeprecatedID(Base):
    __tablename__ = "deprecated_id"

    id = Column(Integer, primary_key=True)
    type = Column(custom_enums["deprecated_id_types"])
    deprecated_id = Column(COL_ID_STR)
    ome_id = Column(COL_ID_STR)

    __table_args__ = (UniqueConstraint("type", "deprecated_id", "ome_id"),)

    def __repr__(self):
        return '<cobradb DeprecatedID(type="%s", deprecated_id="%s", ome_id=%d)>' % (
            self.type,
            self.deprecated_id,
            self.ome_id,
        )


class Model(Base):
    __tablename__ = "model"

    id = Column(COL_ID_STR, primary_key=True)
    genome_id = Column(
        Integer, ForeignKey("genome.id", onupdate="CASCADE", ondelete="CASCADE")
    )
    organism = Column(String(200), nullable=True)
    published_filename = Column(String, nullable=True)

    def __repr__(self):
        return "<cobradb Model(id={self.id}, bigg_id={self.bigg_id})>".format(self=self)


class ModelGene(Base):
    __tablename__ = "model_gene"

    id = Column(Integer, primary_key=True)
    model_id = Column(
        COL_ID_STR,
        ForeignKey("model.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )
    gene_id = Column(
        Integer,
        ForeignKey("gene.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )

    __table_args__ = (UniqueConstraint("model_id", "gene_id"),)


class ModelReaction(Base):
    __tablename__ = "model_reaction"

    id = Column(Integer, primary_key=True)
    reaction_id = Column(
        COL_ID_STR,
        ForeignKey("reaction.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )
    model_id = Column(
        COL_ID_STR,
        ForeignKey("model.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )
    copy_number = Column(Integer, nullable=False)

    objective_coefficient = Column(Float, nullable=False)
    lower_bound = Column(Float, nullable=False)
    upper_bound = Column(Float, nullable=False)
    gene_reaction_rule = Column(String, nullable=False)
    original_gene_reaction_rule = Column(String, nullable=True)
    subsystem = Column(String, nullable=True)

    __table_args__ = (UniqueConstraint("reaction_id", "model_id", "copy_number"),)

    def __repr__(self):
        return "<cobradb ModelCeaction(id={self.id}, reaction_id={self.reaction_id}, model_id={self.model_id}, copy_number={self.copy_number})>".format(
            self=self
        )


class GeneReactionMatrix(Base):
    __tablename__ = "gene_reaction_matrix"

    id = Column(Integer, primary_key=True)
    model_gene_id = Column(
        Integer,
        ForeignKey("model_gene.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )
    model_reaction_id = Column(
        Integer,
        ForeignKey("model_reaction.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )

    __table_args__ = (UniqueConstraint("model_gene_id", "model_reaction_id"),)

    def __repr__(self):
        return "<cobradb GeneReactionMatrix(id={self.id}, model_gene_id={self.model_gene_id}, model_reaction_id={self.model_reaction_id})>".format(
            self=self
        )


class UniversalCompartmentalizedComponent(Base):
    __tablename__ = "universal_compartmentalized_component"
    id = Column(COL_ID_STR, primary_key=True)
    universal_component_id = Column(
        COL_ID_STR,
        ForeignKey(UniversalComponent.id, onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )
    compartment_id = Column(
        COL_ID_STR,
        ForeignKey("compartment.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )


class CompartmentalizedComponent(Base):
    __tablename__ = "compartmentalized_component"
    id = Column(COL_ID_STR, primary_key=True)
    component_id = Column(
        COL_ID_STR,
        ForeignKey("component.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )
    universal_id = Column(
        COL_ID_STR,
        ForeignKey(UniversalCompartmentalizedComponent.id),
        nullable=False,
    )
    compartment_id = Column(
        COL_ID_STR,
        ForeignKey("compartment.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )

    __table_args__ = (UniqueConstraint("compartment_id", "component_id"),)


class ModelCompartmentalizedComponent(Base):
    __tablename__ = "model_compartmentalized_component"
    id = Column(Integer, primary_key=True)
    model_id = Column(
        COL_ID_STR,
        ForeignKey("model.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )
    compartmentalized_component_id = Column(
        COL_ID_STR, ForeignKey("compartmentalized_component.id"), nullable=False
    )
    formula = Column(String, nullable=True)
    charge = Column(Integer, nullable=True)

    __table_args__ = (UniqueConstraint("compartmentalized_component_id", "model_id"),)


class Compartment(Base):
    __tablename__ = "compartment"
    id = Column(COL_ID_STR, primary_key=True)
    name = Column(String)

    def __repr__(self):
        return "<cobradb Compartment(id={self.id})>".format(self=self)


class EscherMap(Base):
    __tablename__ = "escher_map"
    id = Column(Integer, primary_key=True)
    map_name = Column(String, nullable=False)
    map_data = Column(LargeBinary, nullable=False)
    model_id = Column(COL_ID_STR, ForeignKey(Model.id), nullable=False)
    priority = Column(Integer, nullable=False)

    __table_args__ = (UniqueConstraint("map_name"),)


class EscherMapMatrix(Base):
    __tablename__ = "escher_map_matrix"
    id = Column(Integer, primary_key=True)
    ome_id = Column(Integer, nullable=False)
    type = Column(custom_enums["escher_map_matrix_type"], nullable=False)
    escher_map_id = Column(Integer, ForeignKey(EscherMap.id), nullable=False)
    # the reaction id or node id
    escher_map_element_id = Column(String(50))

    __table_args__ = (UniqueConstraint("ome_id", "type", "escher_map_id"),)


class ModelCount(Base):
    __tablename__ = "model_count"
    id = Column(Integer, primary_key=True)
    model_id = Column(
        COL_ID_STR,
        ForeignKey("model.id", onupdate="CASCADE", ondelete="CASCADE"),
        nullable=False,
    )
    reaction_count = Column(Integer)
    gene_count = Column(Integer)
    metabolite_count = Column(Integer)


class Gene(GenomeRegion):
    __tablename__ = "gene"

    id = Column(
        Integer,
        ForeignKey("genome_region.id", onupdate="CASCADE", ondelete="CASCADE"),
        primary_key=True,
    )
    name = Column(String, nullable=True)
    bigg_id = Column(String, nullable=False)
    locus_tag = Column(String, nullable=True)
    mapped_to_genbank = Column(Boolean, nullable=False)
    alternative_transcript_of = Column(Integer, ForeignKey("gene.id"), nullable=True)

    __mapper_args__ = {"polymorphic_identity": "gene"}

    def __repr__(self):
        return "<cobradb Gene(id=%d, bigg_id=%s, name=%s)>" % (
            self.id,
            self.bigg_id,
            self.name,
        )


class ReferenceReactivePart(Base):
    __tablename__ = "reference_reactive_part"

    id = Column(COL_ID_STR, primary_key=True)
    name = Column(COL_NAME_STR, nullable=False)
    html_name = Column(COL_HTML_NAME_STR, nullable=True)
    formula = Column(COL_FORMULA_STR, nullable=True)
    charge = Column(COL_CHARGE_STR, nullable=True)

    def __repr__(self):
        return (
            "<cobradb ReferenceReactivePart(id={self.id}, "
            "accession={self.accession}>"
        ).format(self=self)


class ReferenceReactivePartMatrix(Base):
    __tablename__ = "reference_reactive_part_matrix"

    id = Column(Integer, primary_key=True)
    compound_id = Column(COL_ID_STR, ForeignKey(ReferenceCompound.id), nullable=False)
    reactive_part_id = Column(
        COL_ID_STR, ForeignKey(ReferenceReactivePart.id), nullable=False
    )

    def __repr__(self):
        return (
            "<cobradb ReferenceReactivePartMatrix(id={self.id}, "
            "compound_id={self.compound_id}, "
            "reactive_part_id={self.reactive_part_id}>"
        ).format(self=self)


class ReferenceReaction(Base):
    __tablename__ = "reference_reaction"

    id = Column(COL_ID_STR, primary_key=True)
    name = Column(COL_NAME_STR, nullable=True)
    equation = Column(COL_EQUATION_STR, nullable=True)

    def __repr__(self):
        return (
            "<cobradb ReferenceReaction(id={self.id}, " "equation='{self.equation}'>"
        ).format(self=self)


class ReferenceReactionParticipant(Base):
    __tablename__ = "reference_reaction_participant"

    id = Column(Integer, primary_key=True)
    reaction_id = Column(COL_ID_STR, ForeignKey(ReferenceReaction.id), nullable=False)
    compound_id = Column(COL_ID_STR, ForeignKey(ReferenceCompound.id), nullable=False)
    side = Column(custom_enums["reaction_side"], nullable=False)
    coefficient = Column(COL_COEFFICIENT_STR, nullable=False)
    compartment = Column(COL_COMPARTMENT_STR, nullable=False)

    def __repr__(self):
        return (
            f"<cobradb ReferenceReactionParticipant(id={self.id}, "
            f"reaction_id={self.reaction_id}, "
            f"compound_id={self.compound_id}, side={self.side}, "
            f"coefficient={self.coefficient}, compartment={self.compartment})"
        )


class Reaction(Base):
    __tablename__ = "reaction"

    id = Column(COL_ID_STR, primary_key=True)
    # type = Column(String(20))
    name = Column(String, nullable=True)
    reference_id = Column(COL_ID_STR, ForeignKey(ReferenceReaction.id), nullable=False)
    # reaction_hash = Column(String, nullable=False)
    # pseudoreaction = Column(Boolean, default=False)

    # __mapper_args__ = {"polymorphic_identity": "reaction", "polymorphic_on": type}

    def __repr__(self):
        return "<cobradb Reaction(id=%d)>" % (
            self.id,
            ", pseudoreaction" if self.pseudoreaction else "",
        )


class ReactionMatrix(Base):
    __tablename__ = "reaction_matrix"
    id = Column(Integer, primary_key=True)
    reaction_id = Column(COL_ID_STR, ForeignKey(Reaction.id), nullable=False)
    compartmentalized_component_id = Column(
        COL_ID_STR,
        ForeignKey(
            CompartmentalizedComponent.id, onupdate="CASCADE", ondelete="CASCADE"
        ),
        nullable=False,
    )
    coefficient = Column(Float)

    __table_args__ = (
        UniqueConstraint("reaction_id", "compartmentalized_component_id"),
    )


class ComponentIDMapping(Base):
    __tablename__ = "component_id_mapping"

    old_id = Column(COL_ID_STR, primary_key=True)
    new_id = Column(COL_ID_STR, ForeignKey(UniversalComponent.id), nullable=False)
