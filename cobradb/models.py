# -*- coding: utf-8 -*-

"""Module to implement ORM to the ome database"""

import datetime
import math
from operator import itemgetter
from typing import List, Optional

from cobradb.inchi import (
    inchi_object_to_inchikey,
    inchi_object_to_string,
    string_to_inchi_dict,
)
from cobradb.settings import db_connection_string

from sqlalchemy import (
    ForeignKey,
    String,
    Boolean,
    create_engine,
    MetaData,
    Enum,
    DateTime,
    UniqueConstraint,
    func,
)
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.orm import DeclarativeBase
from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.ext.hybrid import hybrid_property

# Connect to postgres
engine = create_engine(db_connection_string)
metadata = MetaData()
Session = sessionmaker(bind=engine)

HASH_STR_MAX_LEN = 5000
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
    Enum("passed", "failed", "skipped", name="test_result"),
]
custom_enums = {x.name: x for x in _enum_l}

COMPOUND_TYPE_TO_SBO = {
    "small_molecule": "SBO:0000247",
    "generic_polypeptide": "SBO:0000252",
    "generic_polynucleotide": "SBO:0000246",
    "polymer": "SBO:0000248",
}

# ------------
# Exceptions
# ------------


class NotFoundError(Exception):
    pass


class AlreadyLoadedError(Exception):
    def __init__(self, message, bigg_id=None, db_id=None):
        super().__init__(message)
        self.bigg_id = bigg_id
        self.db_id = db_id


# --------
# Tables
# --------


class Base(DeclarativeBase):
    pass


class DatabaseVersion(Base):
    __tablename__ = "database_version"

    is_version = mapped_column(custom_enums["is_version"], primary_key=True)
    date_time: Mapped[datetime.datetime] = mapped_column(DateTime)

    __table_args__ = (UniqueConstraint("is_version"),)

    def __init__(self, date_time):
        self.is_version = "is_version"
        self.date_time = date_time


class InChI(Base):
    __tablename__ = "inchi"

    id: Mapped[int] = mapped_column(primary_key=True)
    formula: Mapped[str]
    c: Mapped[Optional[str]]
    h: Mapped[Optional[str]]
    q: Mapped[Optional[str]]
    p: Mapped[Optional[str]]
    b: Mapped[Optional[str]]
    t: Mapped[Optional[str]]
    m: Mapped[Optional[str]]
    s: Mapped[Optional[str]]

    key_major: Mapped[str] = mapped_column(String(14))
    key_minor: Mapped[str] = mapped_column(String(10))
    key_proton: Mapped[str] = mapped_column(String(1))

    reference_compounds: Mapped[List["ReferenceCompound"]] = relationship(
        back_populates="inchi"
    )
    reference_reactive_parts: Mapped[List["ReferenceReactivePart"]] = relationship(
        back_populates="inchi"
    )

    @hybrid_property
    def key(self):
        return self.key_major + "-" + self.key_minor + "-" + self.key_proton

    def to_string(self):
        return inchi_object_to_string(self)

    def calculate_key_parts(self):
        major, minor, proton = inchi_object_to_inchikey(self)
        self.key_major = major
        self.key_minor = minor
        self.key_proton = proton

    @classmethod
    def from_string(cls, inchi):
        d = string_to_inchi_dict(inchi)
        if d is None:
            return None
        o = cls(**d)
        o.calculate_key_parts()
        return o

    def __eq__(self, other):
        if other is None:
            return False
        for x in ["formula", "c", "h", "q", "p", "b", "t", "m", "s"]:
            if getattr(self, x) != getattr(other, x):
                return False
        return True


class Genome(Base):
    __tablename__ = "genome"

    id: Mapped[int] = mapped_column(primary_key=True)
    accession_type: Mapped[str]
    accession_value: Mapped[str]
    organism: Mapped[Optional[str]]
    taxon_id: Mapped[Optional[str]]
    ncbi_assembly_id: Mapped[Optional[str]]

    chromosomes: Mapped[List["Chromosome"]] = relationship(back_populates="genome")
    models: Mapped[List["Model"]] = relationship(back_populates="genome")

    __table_args__ = (UniqueConstraint("accession_type", "accession_value"),)

    def __repr__(self):
        return (
            "<cobradb Genome(id={self.id}, accession_type={self.accession_type}, "
            "accession_value={self.accession_value})>".format(self=self)
        )


class Chromosome(Base):
    __tablename__ = "chromosome"

    id: Mapped[int] = mapped_column(primary_key=True)
    ncbi_accession: Mapped[str]

    genome_id: Mapped[int] = mapped_column(ForeignKey("genome.id"))
    genome: Mapped["Genome"] = relationship(back_populates="chromosomes")

    genome_regions: Mapped[List["GenomeRegion"]] = relationship(
        back_populates="chromosome"
    )

    __table_args__ = (UniqueConstraint("ncbi_accession", "genome_id"),)

    def __repr__(self):
        return "<cobradb Chromosome(id={self.id}, ncbi_accession={self.ncbi_accession}, genome_id={self.genome_id})>".format(
            self=self
        )


class GenomeRegion(Base):
    __tablename__ = "genome_region"

    id: Mapped[int] = mapped_column(primary_key=True)

    chromosome_id: Mapped[Optional[int]] = mapped_column(ForeignKey("chromosome.id"))
    chromosome: Mapped[Optional["Chromosome"]] = relationship(
        back_populates="genome_regions"
    )

    bigg_id: Mapped[str]
    leftpos: Mapped[Optional[int]]
    rightpos: Mapped[Optional[int]]
    strand: Mapped[Optional[str]] = mapped_column(String(1))
    type: Mapped[str] = mapped_column(String(20))
    dna_sequence: Mapped[Optional[str]]
    protein_sequence: Mapped[Optional[str]]

    __table_args__ = (UniqueConstraint("bigg_id", "chromosome_id"),)

    __mapper_args__ = {"polymorphic_identity": "genome_region", "polymorphic_on": type}

    def __repr__(self):
        return "<cobradb GenomeRegion(id={self.id}, leftpos={self.leftpos}, rightpos={self.rightpos})>".format(
            self=self
        )


class ReferenceCompound(Base):
    __tablename__ = "reference_compound"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]
    name: Mapped[str]
    html_name: Mapped[Optional[str]]
    compound_type = mapped_column(custom_enums["compound_type"], nullable=False)
    charge: Mapped[Optional[str]]
    formula: Mapped[Optional[str]]

    inchi_id: Mapped[Optional[int]] = mapped_column(ForeignKey(InChI.id))
    inchi: Mapped[Optional[InChI]] = relationship(back_populates="reference_compounds")

    reactive_part_matrix: Mapped[List["ReferenceReactivePartMatrix"]] = relationship(
        back_populates="compound"
    )

    reaction_participants: Mapped[List["ReferenceReactionParticipant"]] = relationship(
        back_populates="compound"
    )

    reference_mappings: Mapped[List["ComponentReferenceMapping"]] = relationship(
        back_populates="reference_compound"
    )

    __table_args__ = (UniqueConstraint("bigg_id"),)

    def get_sbo(self):
        return COMPOUND_TYPE_TO_SBO.get(str(self.compound_type), "small_molecule")

    def __repr__(self):
        return (
            "<cobradb ReferenceCompound(id={self.id}, name={self.name}, "
            "comound_type={self.compound_type})>".format(self=self)
        )


class UniversalComponent(Base):
    __tablename__ = "universal_component"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]
    name: Mapped[Optional[str]]
    components: Mapped[List["Component"]] = relationship(
        back_populates="universal_component"
    )
    reference_mapping: Mapped["UniversalComponentReferenceMapping"] = relationship(
        back_populates="universal_component"
    )

    model_id: Mapped[Optional[int]] = mapped_column(ForeignKey("model.id"))
    model: Mapped[Optional["Model"]] = relationship(
        back_populates="model_namespace_universal_components"
    )

    universal_compartmentalized_components: Mapped[
        List["UniversalCompartmentalizedComponent"]
    ] = relationship(back_populates="universal_component")

    old_bigg_ids: Mapped[List["ComponentIDMapping"]] = relationship(
        back_populates="new_universal_component"
    )

    reference_mappings: Mapped[List["ComponentReferenceMapping"]] = relationship(
        back_populates="universal_component"
    )

    __table_args__ = (UniqueConstraint("bigg_id"),)


class Component(Base):
    __tablename__ = "component"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    universal_component_id: Mapped[int] = mapped_column(
        ForeignKey(UniversalComponent.id)
    )
    universal_component: Mapped[UniversalComponent] = relationship(
        back_populates="components"
    )

    name: Mapped[Optional[str]]
    # type = Column(String(20))
    formula: Mapped[Optional[str]]
    charge: Mapped[int]

    model_id: Mapped[Optional[int]] = mapped_column(ForeignKey("model.id"))
    model: Mapped[Optional["Model"]] = relationship(
        back_populates="model_namespace_components"
    )

    compartmentalized_components: Mapped[List["CompartmentalizedComponent"]] = (
        relationship(back_populates="component")
    )

    reference_mappings: Mapped[List["ComponentReferenceMapping"]] = relationship(
        back_populates="component"
    )

    # __mapper_args__ = {"polymorphic_identity": "component", "polymorphic_on": type}

    __table_args__ = (UniqueConstraint("bigg_id"),)

    @staticmethod
    def charge_to_string(coefficient):
        if coefficient is None:
            return 0
        if isinstance(coefficient, str):
            return coefficient
        if isinstance(coefficient, int):
            return str(coefficient)
        if isinstance(coefficient, float):
            if coefficient.is_integer():
                return str(int(coefficient))
            else:
                return str(coefficient)

    def __repr__(self):
        return f"Component ({self.id}): {self.name}"


class ComponentReferenceMapping(Base):
    __tablename__ = "component_reference_mapping"

    id: Mapped[int] = mapped_column(primary_key=True)

    component_id: Mapped[int] = mapped_column(ForeignKey(Component.id))
    component: Mapped[Component] = relationship(back_populates="reference_mappings")

    universal_component_id: Mapped[int] = mapped_column(
        ForeignKey(UniversalComponent.id)
    )
    universal_component: Mapped[UniversalComponent] = relationship(
        back_populates="reference_mappings"
    )

    reference_compound_id: Mapped[int] = mapped_column(ForeignKey(ReferenceCompound.id))
    reference_compound: Mapped[ReferenceCompound] = relationship(
        back_populates="reference_mappings"
    )

    reference_n: Mapped[Optional[int]]


class UniversalComponentReferenceMapping(Base):
    __tablename__ = "universal_component_reference_mapping"

    id: Mapped[int] = mapped_column(ForeignKey(UniversalComponent.id), primary_key=True)
    universal_component: Mapped[UniversalComponent] = relationship(
        back_populates="reference_mapping"
    )

    mapping_id: Mapped[int] = mapped_column(ForeignKey(ComponentReferenceMapping.id))
    mapping: Mapped[ComponentReferenceMapping] = relationship()


class DataSource(Base):
    __tablename__ = "data_source"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]
    name: Mapped[Optional[str]]
    url_prefix: Mapped[Optional[str]]

    synonyms: Mapped[List["Synonym"]] = relationship(back_populates="data_source")

    __table_args__ = (UniqueConstraint("bigg_id"),)

    def __repr__(self):
        return (
            "<cobradb DataSource(id={self.id}, bigg_id={self.bigg_id}, "
            "name={self.name}, url_prefix={self.url_prefix})>"
        ).format(self=self)


class Synonym(Base):
    __tablename__ = "synonym"

    id: Mapped[int] = mapped_column(primary_key=True)
    ome_id: Mapped[int]
    synonym: Mapped[str]
    type = mapped_column(custom_enums["synonym_type"])

    data_source_id: Mapped[int] = mapped_column(
        ForeignKey("data_source.id", ondelete="CASCADE")
    )
    data_source: Mapped[DataSource] = relationship(back_populates="synonyms")

    __table_args__ = (UniqueConstraint("ome_id", "synonym", "type", "data_source_id"),)

    def __repr__(self):
        return (
            '<cobradb Synonym(id=%d, synonym="%s", type="%s", ome_id=%d, data_source_id=%d)>'
            % (self.id, self.synonym, self.type, self.ome_id, self.data_source_id)
        )


class Publication(Base):
    __tablename__ = "publication"

    id: Mapped[int] = mapped_column(primary_key=True)
    reference_type = mapped_column(custom_enums["reference_type"])

    reference_id: Mapped[str]

    publication_models: Mapped[List["PublicationModel"]] = relationship(
        back_populates="publication"
    )

    __table_args__ = (UniqueConstraint("reference_type", "reference_id"),)


class PublicationModel(Base):
    __tablename__ = "publication_model"

    model_id: Mapped[int] = mapped_column(
        ForeignKey("model.id", ondelete="CASCADE"), primary_key=True
    )
    model: Mapped["Model"] = relationship(back_populates="publication_models")

    publication_id: Mapped[int] = mapped_column(
        ForeignKey("publication.id", ondelete="CASCADE"), primary_key=True
    )
    publication: Mapped["Publication"] = relationship(
        back_populates="publication_models"
    )

    __table_args__ = (UniqueConstraint("model_id", "publication_id"),)


class OldIDSynonym(Base):
    __tablename__ = "old_id_model_synonym"

    id: Mapped[int] = mapped_column(primary_key=True)
    type = mapped_column(custom_enums["old_id_synonym_type"])
    synonym_id: Mapped[int] = mapped_column(
        ForeignKey("synonym.id", ondelete="CASCADE"),
    )
    ome_id: Mapped[int]

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

    genome_region_id_1: Mapped[int] = mapped_column(
        ForeignKey("genome_region.id"), primary_key=True
    )
    genome_region_1: Mapped[GenomeRegion] = relationship(
        foreign_keys=[genome_region_id_1]
    )

    genome_region_id_2: Mapped[int] = mapped_column(
        ForeignKey("genome_region.id"), primary_key=True
    )
    genome_region_2: Mapped[GenomeRegion] = relationship(
        foreign_keys=[genome_region_id_2]
    )

    distance: Mapped[int]

    __table_args__ = (UniqueConstraint("genome_region_id_1", "genome_region_id_2"),)

    def __repr__(self):
        return "GenomeRegionMap (%d <--> %d) distance:%d" % (
            self.genome_region_id_1,
            self.genome_region_id_2,
            self.distance,
        )


class DeprecatedID(Base):
    __tablename__ = "deprecated_id"

    id: Mapped[int] = mapped_column(primary_key=True)
    type = mapped_column(custom_enums["deprecated_id_types"])
    deprecated_id: Mapped[str]
    ome_id: Mapped[int]

    __table_args__ = (UniqueConstraint("type", "deprecated_id", "ome_id"),)

    def __repr__(self):
        return '<cobradb DeprecatedID(type="%s", deprecated_id="%s", ome_id=%d)>' % (
            self.type,
            self.deprecated_id,
            self.ome_id,
        )


class Model(Base):
    __tablename__ = "model"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    genome_id: Mapped[int] = mapped_column(
        ForeignKey("genome.id", onupdate="CASCADE", ondelete="CASCADE")
    )
    genome: Mapped[Genome] = relationship(back_populates="models")

    organism: Mapped[Optional[str]]
    published_filename: Mapped[Optional[str]]
    date_modified: Mapped[datetime.datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now()
    )

    model_namespace_components: Mapped[List[Component]] = relationship(
        back_populates="model"
    )
    model_namespace_universal_components: Mapped[List[UniversalComponent]] = (
        relationship(back_populates="model")
    )
    model_namespace_reactions: Mapped[List["Reaction"]] = relationship(
        back_populates="model"
    )
    model_genes: Mapped[List["ModelGene"]] = relationship(back_populates="model")
    model_reactions: Mapped[List["ModelReaction"]] = relationship(
        back_populates="model"
    )
    model_compartmentalized_components: Mapped[
        List["ModelCompartmentalizedComponent"]
    ] = relationship(back_populates="model")

    model_count: Mapped[Optional["ModelCount"]] = relationship(back_populates="model")
    publication_models: Mapped[List["PublicationModel"]] = relationship(
        back_populates="model"
    )
    escher_maps: Mapped[List["EscherMap"]] = relationship(back_populates="model")
    memote_results: Mapped[List["MemoteResult"]] = relationship(back_populates="model")

    __table_args__ = (UniqueConstraint("bigg_id"),)

    def __repr__(self):
        return "<cobradb Model(id={self.id}, bigg_id={self.bigg_id})>".format(self=self)


class ModelGene(Base):
    __tablename__ = "model_gene"

    id: Mapped[int] = mapped_column(primary_key=True)
    model_id: Mapped[int] = mapped_column(
        ForeignKey("model.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    model: Mapped["Model"] = relationship(back_populates="model_genes")

    gene_id: Mapped[int] = mapped_column(
        ForeignKey("gene.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    gene: Mapped["Gene"] = relationship(back_populates="model_genes")

    reaction_matrix: Mapped[List["GeneReactionMatrix"]] = relationship(
        back_populates="model_gene"
    )

    memote_results: Mapped[List["MemoteResult"]] = relationship(
        back_populates="model_gene"
    )

    __table_args__ = (UniqueConstraint("model_id", "gene_id"),)


class ModelReaction(Base):
    __tablename__ = "model_reaction"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    reaction_id: Mapped[int] = mapped_column(
        ForeignKey("reaction.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    reaction: Mapped["Reaction"] = relationship(back_populates="model_reactions")

    model_id: Mapped[int] = mapped_column(
        ForeignKey("model.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    model: Mapped["Model"] = relationship(back_populates="model_reactions")

    copy_number: Mapped[int]

    objective_coefficient: Mapped[float]
    lower_bound: Mapped[float]
    upper_bound: Mapped[float]

    gene_reaction_rule: Mapped[str]
    original_gene_reaction_rule: Mapped[Optional[str]]
    subsystem: Mapped[str]

    reaction_matrix: Mapped[List["GeneReactionMatrix"]] = relationship(
        back_populates="model_reaction"
    )

    memote_results: Mapped[List["MemoteResult"]] = relationship(
        back_populates="model_reaction"
    )

    __table_args__ = (
        UniqueConstraint("reaction_id", "model_id", "copy_number"),
        UniqueConstraint("bigg_id"),
    )

    @staticmethod
    def create_id(model_id, universal_reaction_id, copy_number):
        if copy_number == 1:
            return f"{model_id}|{universal_reaction_id}"
        else:
            return f"{model_id}|{universal_reaction_id}:{copy_number}"

    def __repr__(self):
        return "<cobradb ModelCeaction(id={self.id}, reaction_id={self.reaction_id}, model_id={self.model_id}, copy_number={self.copy_number})>".format(
            self=self
        )

    @staticmethod
    def interpret_id(bigg_id):
        copy_number = 1
        if ":" in bigg_id:
            bigg_id, copy_number = bigg_id.rsplit(":", maxsplit=1)
            copy_number = int(copy_number)
        return bigg_id, copy_number


class GeneReactionMatrix(Base):
    __tablename__ = "gene_reaction_matrix"

    id: Mapped[int] = mapped_column(primary_key=True)

    model_gene_id: Mapped[int] = mapped_column(
        ForeignKey("model_gene.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    model_gene: Mapped[ModelGene] = relationship(back_populates="reaction_matrix")

    model_reaction_id: Mapped[int] = mapped_column(
        ForeignKey("model_reaction.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    model_reaction: Mapped[ModelReaction] = relationship(
        back_populates="reaction_matrix"
    )

    __table_args__ = (UniqueConstraint("model_gene_id", "model_reaction_id"),)

    def __repr__(self):
        return "<cobradb GeneReactionMatrix(id={self.id}, model_gene_id={self.model_gene_id}, model_reaction_id={self.model_reaction_id})>".format(
            self=self
        )


class UniversalCompartmentalizedComponent(Base):
    __tablename__ = "universal_compartmentalized_component"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    universal_component_id: Mapped[int] = mapped_column(
        ForeignKey(UniversalComponent.id, onupdate="CASCADE", ondelete="CASCADE"),
    )
    universal_component: Mapped[UniversalComponent] = relationship(
        back_populates="universal_compartmentalized_components"
    )

    compartment_id: Mapped[int] = mapped_column(
        ForeignKey("compartment.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    compartment: Mapped["Compartment"] = relationship(
        back_populates="universal_compartmentalized_components"
    )

    compartmentalized_components: Mapped[List["CompartmentalizedComponent"]] = (
        relationship(back_populates="universal_compartmentalized_component")
    )

    matrix: Mapped[List["UniversalReactionMatrix"]] = relationship(
        back_populates="universal_compartmentalized_component"
    )

    __table_args__ = (UniqueConstraint("bigg_id"),)


class CompartmentalizedComponent(Base):
    __tablename__ = "compartmentalized_component"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    component_id: Mapped[int] = mapped_column(
        ForeignKey("component.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    component: Mapped[Component] = relationship(
        back_populates="compartmentalized_components"
    )

    compartment_id: Mapped[int] = mapped_column(
        ForeignKey("compartment.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    compartment: Mapped["Compartment"] = relationship(
        back_populates="compartmentalized_components"
    )

    universal_compartmentalized_component_id: Mapped[int] = mapped_column(
        ForeignKey(UniversalCompartmentalizedComponent.id)
    )
    universal_compartmentalized_component: Mapped[
        UniversalCompartmentalizedComponent
    ] = relationship(back_populates="compartmentalized_components")

    model_compartmentalized_components: Mapped[
        List["ModelCompartmentalizedComponent"]
    ] = relationship(back_populates="compartmentalized_component")

    matrix: Mapped[List["ReactionMatrix"]] = relationship(
        back_populates="compartmentalized_component"
    )

    __table_args__ = (
        UniqueConstraint("compartment_id", "component_id"),
        UniqueConstraint("bigg_id"),
    )


class ModelCompartmentalizedComponent(Base):
    __tablename__ = "model_compartmentalized_component"

    id: Mapped[int] = mapped_column(primary_key=True)

    model_id: Mapped[int] = mapped_column(
        ForeignKey("model.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    model: Mapped["Model"] = relationship(
        back_populates="model_compartmentalized_components"
    )

    compartmentalized_component_id: Mapped[int] = mapped_column(
        ForeignKey("compartmentalized_component.id")
    )
    compartmentalized_component: Mapped[CompartmentalizedComponent] = relationship(
        back_populates="model_compartmentalized_components"
    )

    memote_results: Mapped[List["MemoteResult"]] = relationship(
        back_populates="model_compartmentalized_component"
    )

    __table_args__ = (UniqueConstraint("compartmentalized_component_id", "model_id"),)


class Compartment(Base):
    __tablename__ = "compartment"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]
    name: Mapped[str]

    universal_compartmentalized_components: Mapped[
        List["UniversalCompartmentalizedComponent"]
    ] = relationship(back_populates="compartment")
    compartmentalized_components: Mapped[List["CompartmentalizedComponent"]] = (
        relationship(back_populates="compartment")
    )

    __table_args__ = (UniqueConstraint("bigg_id"),)

    def __repr__(self):
        return "<cobradb Compartment(id={self.id})>".format(self=self)


class EscherMap(Base):
    __tablename__ = "escher_map"

    id: Mapped[int] = mapped_column(primary_key=True)

    map_name: Mapped[str]
    map_data: Mapped[str]

    model_id: Mapped[int] = mapped_column(ForeignKey(Model.id))
    model: Mapped["Model"] = relationship(back_populates="escher_maps")

    priority: Mapped[int]

    matrix: Mapped[List["EscherMapMatrix"]] = relationship(back_populates="escher_map")

    __table_args__ = (UniqueConstraint("map_name"),)


class EscherMapMatrix(Base):
    __tablename__ = "escher_map_matrix"

    id: Mapped[int] = mapped_column(primary_key=True)
    ome_id: Mapped[int]
    type = mapped_column(custom_enums["escher_map_matrix_type"], nullable=False)

    escher_map_id: Mapped[int] = mapped_column(ForeignKey(EscherMap.id))
    escher_map: Mapped[EscherMap] = relationship(back_populates="matrix")
    # the reaction id or node id
    escher_map_element_id: Mapped[str]

    __table_args__ = (UniqueConstraint("ome_id", "type", "escher_map_id"),)


class ModelCount(Base):
    __tablename__ = "model_count"

    id: Mapped[int] = mapped_column(primary_key=True)

    model_id: Mapped[int] = mapped_column(
        ForeignKey("model.id", onupdate="CASCADE", ondelete="CASCADE"),
    )
    model: Mapped["Model"] = relationship(back_populates="model_count")

    reaction_count: Mapped[int]
    gene_count: Mapped[int]
    metabolite_count: Mapped[int]


class Gene(GenomeRegion):
    __tablename__ = "gene"

    id: Mapped[int] = mapped_column(
        ForeignKey("genome_region.id", onupdate="CASCADE", ondelete="CASCADE"),
        primary_key=True,
    )
    bigg_id: Mapped[str]
    name: Mapped[Optional[str]]
    locus_tag: Mapped[Optional[str]]
    mapped_to_genbank: Mapped[bool]

    alternative_transcript_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey("gene.id")
    )
    alternative_transcript_of: Mapped[Optional["Gene"]] = relationship(
        foreign_keys=[alternative_transcript_id]
    )

    model_genes: Mapped[List[ModelGene]] = relationship(back_populates="gene")

    __table_args__ = (UniqueConstraint("bigg_id"),)
    __mapper_args__ = {"polymorphic_identity": "gene"}

    def __repr__(self):
        return "<cobradb Gene(id=%d, bigg_id=%s, name=%s)>" % (
            self.id,
            self.bigg_id,
            self.name,
        )


class ReferenceReactivePart(Base):
    __tablename__ = "reference_reactive_part"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    name: Mapped[str]
    html_name: Mapped[Optional[str]]
    formula: Mapped[Optional[str]]
    charge: Mapped[Optional[str]]

    inchi_id: Mapped[Optional[int]] = mapped_column(ForeignKey(InChI.id))
    inchi: Mapped[Optional[InChI]] = relationship(
        back_populates="reference_reactive_parts"
    )

    matrix: Mapped[List["ReferenceReactivePartMatrix"]] = relationship(
        back_populates="reactive_part"
    )

    __table_args__ = (UniqueConstraint("bigg_id"),)

    def __repr__(self):
        return (
            "<cobradb ReferenceReactivePart(id={self.id}, "
            "accession={self.accession}>"
        ).format(self=self)


class ReferenceReactivePartMatrix(Base):
    __tablename__ = "reference_reactive_part_matrix"

    id: Mapped[int] = mapped_column(primary_key=True)

    compound_id: Mapped[int] = mapped_column(ForeignKey(ReferenceCompound.id))
    compound: Mapped[ReferenceCompound] = relationship(
        back_populates="reactive_part_matrix"
    )

    reactive_part_id: Mapped[int] = mapped_column(ForeignKey(ReferenceReactivePart.id))
    reactive_part: Mapped[ReferenceReactivePart] = relationship(back_populates="matrix")

    def __repr__(self):
        return (
            "<cobradb ReferenceReactivePartMatrix(id={self.id}, "
            "compound_id={self.compound_id}, "
            "reactive_part_id={self.reactive_part_id}>"
        ).format(self=self)


class ReferenceReaction(Base):
    __tablename__ = "reference_reaction"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    name: Mapped[Optional[str]]
    equation: Mapped[Optional[str]]

    hash: Mapped[str]

    reaction_participants: Mapped[List["ReferenceReactionParticipant"]] = relationship(
        back_populates="reaction"
    )

    universal_reactions: Mapped[List["UniversalReaction"]] = relationship(
        back_populates="reference"
    )

    __table_args__ = (UniqueConstraint("bigg_id"),)

    @staticmethod
    def generate_hash(participants, pattern=False):
        part_dict = {}
        for p in participants:
            if isinstance(p, dict):
                rc = p.get("reference_compound")
                if rc is None:
                    rc_id = p["reference_compound_bigg_id"]
                else:
                    rc_id = rc.bigg_id
                coefficient = str(p["coefficient"]).lower()
            else:
                rc_id = p.compound.bigg_id
                coefficient = p.coefficient
            if "n" in coefficient:
                coefficient = math.inf
            else:
                try:
                    coefficient = abs(float(coefficient))
                except:
                    coefficient = math.inf
            part_dict[rc_id] = part_dict.get(rc_id, 0) + coefficient
        if pattern:
            hash_str = "/".join(
                f"(({Reaction.coefficient_to_string(coeff)})|N)\\${p_id}"
                for p_id, coeff in sorted(
                    part_dict.items(),
                    key=itemgetter(0),
                )
            )
            hash_str = f"^{hash_str}$"
        else:
            hash_str = "/".join(
                f"{Reaction.coefficient_to_string(coeff)}${p_id}"
                for p_id, coeff in sorted(
                    part_dict.items(),
                    key=itemgetter(0),
                )
            )
        return hash_str

    def update_hash(self):
        self.hash = ReferenceReaction.generate_hash(self.reaction_participants)

    def __repr__(self):
        return (
            "<cobradb ReferenceReaction(id={self.id}, " "equation='{self.equation}'>"
        ).format(self=self)


class ReferenceReactionParticipant(Base):
    __tablename__ = "reference_reaction_participant"

    id: Mapped[int] = mapped_column(primary_key=True)

    reaction_id: Mapped[int] = mapped_column(ForeignKey(ReferenceReaction.id))
    reaction: Mapped[ReferenceReaction] = relationship(
        back_populates="reaction_participants"
    )

    compound_id: Mapped[int] = mapped_column(ForeignKey(ReferenceCompound.id))
    compound: Mapped[ReferenceCompound] = relationship(
        back_populates="reaction_participants"
    )

    side = mapped_column(custom_enums["reaction_side"], nullable=False)
    coefficient: Mapped[str]
    compartment: Mapped[str]

    universal_reaction_matrix: Mapped[List["UniversalReactionMatrix"]] = relationship(
        back_populates="reference_reaction_participant"
    )

    def __repr__(self):
        return (
            f"<cobradb ReferenceReactionParticipant(id={self.id}, "
            f"reaction_id={self.reaction_id}, "
            f"compound_id={self.compound_id}, side={self.side}, "
            f"coefficient={self.coefficient}, compartment={self.compartment})"
        )


class UniversalReaction(Base):
    __tablename__ = "universal_reaction"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    hash: Mapped[str]
    # type = Column(String(20))
    name: Mapped[Optional[str]]

    reference_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey(ReferenceReaction.id)
    )
    reference: Mapped[Optional[ReferenceReaction]] = relationship(
        back_populates="universal_reactions"
    )
    # reaction_hash = Column(String, nullable=False)
    # pseudoreaction = Column(Boolean, default=False)
    reactions: Mapped[List["Reaction"]] = relationship(
        back_populates="universal_reaction"
    )
    matrix: Mapped[List["UniversalReactionMatrix"]] = relationship(
        back_populates="universal_reaction"
    )

    __table_args__ = (UniqueConstraint("bigg_id"),)
    # __mapper_args__ = {"polymorphic_identity": "reaction", "polymorphic_on": type}

    def get_sbo(self):
        bare_id = str(self.id)
        if bare_id.startswith("__"):
            bare_id = bare_id[2:]
            if "__" in bare_id:
                _model, bare_id = bare_id.split("__", maxsplit=1)
            else:
                bare_id = str(self.id)
        if bare_id.startswith("EX_"):
            return "SBO:0000627"
        if bare_id.startswith("SK_"):
            return "SBO:0000632"
        if bare_id.startswith("DM_"):
            return "SBO:0000628"
        if bare_id.startswith("BIOMASS"):
            return "SBO:0000629"

        return "SBO:0000176"

    def __repr__(self):
        return "<cobradb Reaction(id=%s)>" % (self.id,)

    @staticmethod
    def generate_hash(components):
        comp_dict = {}
        for c in components:
            cc = c.get("universal_compartmentalized_component")
            if cc is None:
                cc_id = c["universal_compartmentalized_component_bigg_id"]
            else:
                cc_id = cc.bigg_id
            comp_dict[cc_id] = comp_dict.get(cc_id, 0) + c["coefficient"]
            # if comp_dict[cc_id] == 0:
            #     del comp_dict[cc_id]
        sorting_1 = sorted(comp_dict.items(), key=lambda x: (x[0], abs(x[1])))
        first_item = next(iter(sorting_1))
        if first_item[1] > 0:
            comp_dict = {k: -1 * v for k, v in comp_dict.items()}
        hash_str = "/".join(
            f"{Reaction.coefficient_to_string(coeff)}${cc_id}"
            for cc_id, coeff in sorted(
                comp_dict.items(),
                key=itemgetter(0, 1),
            )
        )
        if len(hash_str) > HASH_STR_MAX_LEN:
            import hashlib

            hash_str = hashlib.sha256(hash_str.encode()).hexdigest()
        return hash_str


class Reaction(Base):
    __tablename__ = "reaction"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    hash: Mapped[str]

    model_id: Mapped[Optional[int]] = mapped_column(ForeignKey(Model.id))
    model: Mapped[Optional["Model"]] = relationship(
        back_populates="model_namespace_reactions"
    )

    copy_number: Mapped[int]

    universal_reaction_id: Mapped[int] = mapped_column(ForeignKey(UniversalReaction.id))
    universal_reaction: Mapped[UniversalReaction] = relationship(
        back_populates="reactions"
    )

    matrix: Mapped[List["ReactionMatrix"]] = relationship(back_populates="reaction")
    model_reactions: Mapped[List["ModelReaction"]] = relationship(
        back_populates="reaction"
    )

    __table_args__ = (UniqueConstraint("hash", "model_id"), UniqueConstraint("bigg_id"))

    def __repr__(self):
        return "<cobradb Reaction(id=%s)>" % (self.id,)

    @staticmethod
    def create_id(universal_id, copy_number):
        if copy_number != 1:
            return f"{universal_id}:{copy_number}"
        else:
            return universal_id

    @staticmethod
    def coefficient_to_string(coefficient):
        if coefficient is None:
            return 0
        if isinstance(coefficient, str):
            return coefficient
        if isinstance(coefficient, int):
            return str(coefficient)
        if isinstance(coefficient, float):
            if math.isinf(coefficient):
                return "N"
            elif coefficient.is_integer():
                return str(int(coefficient))
            else:
                return str(coefficient)

    @staticmethod
    def generate_hash(components):
        comp_dict = {}
        for c in components:
            cc = c.get("compartmentalized_component")
            if cc is None:
                cc_id = c["compartmentalized_component_bigg_id"]
            else:
                cc_id = cc.bigg_id
            comp_dict[cc_id] = comp_dict.get(cc_id, 0) + c["coefficient"]
            # if comp_dict[cc_id] == 0:
            #     del comp_dict[cc_id]
        sorting_1 = sorted(comp_dict.items(), key=lambda x: (x[0], abs(x[1])))
        first_item = next(iter(sorting_1))
        if first_item[1] > 0:
            comp_dict = {k: -1 * v for k, v in comp_dict.items()}
        hash_str = "/".join(
            f"{Reaction.coefficient_to_string(coeff)}${cc_id}"
            for cc_id, coeff in sorted(
                comp_dict.items(),
                key=itemgetter(0, 1),
            )
        )

        if len(hash_str) > HASH_STR_MAX_LEN:
            import hashlib

            hash_str = hashlib.sha256(hash_str.encode()).hexdigest()
        return hash_str


class UniversalReactionMatrix(Base):
    __tablename__ = "universal_reaction_matrix"

    id: Mapped[int] = mapped_column(primary_key=True)

    universal_reaction_id: Mapped[int] = mapped_column(ForeignKey(UniversalReaction.id))
    universal_reaction: Mapped[UniversalReaction] = relationship(
        back_populates="matrix"
    )

    reference_reaction_participant_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey(ReferenceReactionParticipant.id)
    )
    reference_reaction_participant: Mapped[Optional[ReferenceReactionParticipant]] = (
        relationship(back_populates="universal_reaction_matrix")
    )

    universal_compartmentalized_component_id: Mapped[int] = mapped_column(
        ForeignKey(
            UniversalCompartmentalizedComponent.id,
            onupdate="CASCADE",
            ondelete="CASCADE",
        ),
    )
    universal_compartmentalized_component: Mapped[
        UniversalCompartmentalizedComponent
    ] = relationship(back_populates="matrix")

    coefficient: Mapped[float]

    reaction_matrices: Mapped[List["ReactionMatrix"]] = relationship(
        back_populates="universal_reaction_matrix"
    )
    __table_args__ = (
        UniqueConstraint(
            "universal_reaction_id", "universal_compartmentalized_component_id"
        ),
    )


class ReactionMatrix(Base):
    __tablename__ = "reaction_matrix"

    id: Mapped[int] = mapped_column(primary_key=True)

    reaction_id: Mapped[int] = mapped_column(ForeignKey(Reaction.id))
    reaction: Mapped[Reaction] = relationship(back_populates="matrix")

    universal_reaction_matrix_id: Mapped[int] = mapped_column(
        ForeignKey(UniversalReactionMatrix.id)
    )
    universal_reaction_matrix: Mapped[UniversalReactionMatrix] = relationship(
        back_populates="reaction_matrices"
    )

    compartmentalized_component_id: Mapped[int] = mapped_column(
        ForeignKey(
            CompartmentalizedComponent.id, onupdate="CASCADE", ondelete="CASCADE"
        ),
    )
    compartmentalized_component: Mapped[CompartmentalizedComponent] = relationship(
        back_populates="matrix"
    )


class ComponentIDMapping(Base):
    __tablename__ = "component_id_mapping"

    old_bigg_id: Mapped[str] = mapped_column(primary_key=True)
    new_id: Mapped[int] = mapped_column(ForeignKey(UniversalComponent.id))
    new_universal_component: Mapped[UniversalComponent] = relationship(
        back_populates="old_bigg_ids"
    )


# MEMOTE-related models
class MemoteTest(Base):
    __tablename__ = "memote_test"

    id: Mapped[int] = mapped_column(primary_key=True)
    bigg_id: Mapped[str]

    name: Mapped[str]
    summary: Mapped[Optional[str]]
    format_type: Mapped[str] = mapped_column(String(20))

    invert: Mapped[Optional[bool]]
    invert_result: Mapped[Optional[bool]]
    field: Mapped[Optional[str]] = mapped_column(String(20))

    results: Mapped[List["MemoteResult"]] = relationship(back_populates="test")

    __table_args__ = (UniqueConstraint("bigg_id"),)


class MemoteResult(Base):
    __tablename__ = "memote_result"

    id: Mapped[int] = mapped_column(primary_key=True)

    model_id: Mapped[int] = mapped_column(ForeignKey(Model.id))
    model: Mapped["Model"] = relationship(back_populates="memote_results")

    test_id: Mapped[int] = mapped_column(ForeignKey(MemoteTest.id))
    test: Mapped[MemoteTest] = relationship(back_populates="results")

    model_reaction_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey(ModelReaction.id)
    )
    model_reaction: Mapped[Optional[ModelReaction]] = relationship(
        back_populates="memote_results"
    )

    model_compartmentalized_component_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey(ModelCompartmentalizedComponent.id)
    )
    model_compartmentalized_component: Mapped[
        Optional[ModelCompartmentalizedComponent]
    ] = relationship(back_populates="memote_results")

    model_gene_id: Mapped[Optional[int]] = mapped_column(ForeignKey(ModelGene.id))
    model_gene: Mapped[Optional[ModelGene]] = relationship(
        back_populates="memote_results"
    )

    message: Mapped[Optional[str]]
    data: Mapped[Optional[float]]
    duration: Mapped[Optional[float]]
    metric: Mapped[Optional[float]]
    result = mapped_column(custom_enums["test_result"], nullable=True)

    data_count: Mapped[Optional[int]]
