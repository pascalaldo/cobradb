from functools import partial
import logging
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Optional,
    Tuple,
    Union,
    # assert_type,
)
from cobradb import get_data, parse
from sqlalchemy import select
from sqlalchemy.orm import Session
from cobradb.models import (
    ComponentIDMapping,
    EscherModule,
    ModelReaction,
    ReactionMatrix,
    UniversalComponent,
)
from biggr_maps import map, pathway, template
from pathlib import Path
import os


class EscherModuleDefinition:
    def __init__(
        self, bigg_id: str, name: str, description: str, exclude_pseudo: bool = True
    ):
        self.bigg_id: str = bigg_id
        self.name: str = name
        self.description: str = description
        self.exclude_pseudo: bool = exclude_pseudo
        self._module_db: Optional[EscherModule] = None

    def get_or_create(self, session: Session):
        if self._module_db is not None:
            return self._module_db
        self._module_db = session.scalars(
            select(EscherModule).filter(EscherModule.bigg_id == self.bigg_id).limit(1)
        ).first()
        if self._module_db is None:
            self._module_db = EscherModule(
                bigg_id=self.bigg_id,
                name=self.name,
                description=self.description,
            )
            session.add(self._module_db)
            session.commit()
        return self._module_db

    def select_model_reactions(
        self, session: Session, model_reactions: Iterable[ModelReaction]
    ) -> List[ModelReaction]:
        raise NotImplementedError("Not implemented in base class.")

    def build_map(self, session: Session, model_reactions: Iterable[ModelReaction]):
        m = map.Map(
            name=self.name,
            description=self.description,
            canvas=(-1000, -1000, 2000, 2000),
        )
        return m


def _reactants_sorting(vals):
    node: map.MetaboliteNode = vals[1]
    if node.node_is_primary:
        return 0
    if node.bigg_id.startswith("h_"):
        return 2
    return 1


class EscherBackboneModuleDefinition(EscherModuleDefinition):
    def __init__(
        self,
        backbone: Iterable[map.MetaboliteNode],
        sequential: bool = False,
        placement_kwargs: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        self.backbone = {b.bigg_id: b for b in backbone}
        self.sequential = sequential
        self.placement_kwargs = {} if placement_kwargs is None else placement_kwargs
        super().__init__(**kwargs)

    def _sel_f(self, matched_d: Dict[int, bool]):
        if self.sequential:
            for pos, bool_coeff in matched_d.items():
                if (pos + 1) not in matched_d:
                    continue
                if matched_d[pos + 1] == bool_coeff:
                    continue
                return True
            return False
        else:
            return len(set(matched_d.values())) == 2

    def _sel_bb_f(self, matched_d: Dict[int, bool]):
        if self.sequential:
            for pos, bool_coeff in matched_d.items():
                if (pos + 1) not in matched_d:
                    continue
                if matched_d[pos + 1] == bool_coeff:
                    continue
                return [pos, pos + 1]
            raise ValueError()
        else:
            return list(matched_d.keys())

    def _match_with_bb(
        self,
        model_reaction: ModelReaction,
        backbone_bigg_ids: List[str],
        return_early: bool = False,
    ) -> Dict[int, bool]:
        matched_bb_pos = {}
        for reaction_matrix in model_reaction.reaction.matrix:
            cc_bigg_id = reaction_matrix.compartmentalized_component.bigg_id
            if cc_bigg_id in backbone_bigg_ids:
                bool_coeff = reaction_matrix.universal_reaction_matrix.coefficient > 0
                matched_bb_pos[backbone_bigg_ids.index(cc_bigg_id)] = bool_coeff
            if return_early and self._sel_f(matched_bb_pos):
                break
        return matched_bb_pos

    def select_model_reactions(
        self, session: Session, model_reactions: Iterable[ModelReaction]
    ) -> List[ModelReaction]:
        selected_model_reactions = []
        backbone_bigg_ids = list(self.backbone)
        for model_reaction in model_reactions:
            if (
                self.exclude_pseudo
                and model_reaction.reaction.universal_reaction.is_pseudo
            ):
                continue
            matched_bb_pos = self._match_with_bb(
                model_reaction, backbone_bigg_ids, return_early=True
            )
            if self._sel_f(matched_bb_pos):
                selected_model_reactions.append(model_reaction)
        return selected_model_reactions

    def _build_reaction(
        self, model_reaction: ModelReaction, backbone: Dict[str, map.MetaboliteNode]
    ) -> List[Tuple[Union[int, float], map.MetaboliteNode]]:
        backbone_bigg_ids = list(backbone)
        matched_bb_pos = self._match_with_bb(model_reaction, backbone_bigg_ids)
        pos = self._sel_bb_f(matched_bb_pos)
        backbone_bigg_ids = [backbone_bigg_ids[i] for i in pos]
        backbone = {x: backbone[x] for x in backbone_bigg_ids}

        reactants = []
        for reaction_matrix in model_reaction.reaction.matrix:
            cc_db = reaction_matrix.compartmentalized_component
            bigg_id = cc_db.model_compartmentalized_components[0].bigg_id
            if (node := backbone.get(cc_db.bigg_id)) is None:
                name = cc_db.component.name
                if name is None:
                    name = cc_db.component.universal_component.name
                if name is None:
                    name = "Unknown"
                node = map.MetaboliteNode(
                    bigg_id=bigg_id,
                    name=name,
                    node_is_primary=False,
                )
            else:
                node.bigg_id = bigg_id
            reactants.append(
                (reaction_matrix.universal_reaction_matrix.coefficient, node)
            )
        reactants = list(sorted(reactants, key=_reactants_sorting))
        return reactants

    def build_map(
        self, session: Session, model_reactions: Iterable[ModelReaction]
    ) -> map.Map:
        m = super().build_map(session, model_reactions)

        backbone = {bigg_id: node.copy() for bigg_id, node in self.backbone.items()}
        for node in backbone.values():
            m.add_node(node)

        for model_reaction in model_reactions:
            if (
                self.exclude_pseudo
                and model_reaction.reaction.universal_reaction.is_pseudo
            ):
                continue

            reaction_info = self._build_reaction(model_reaction, backbone)
            name = model_reaction.reaction.universal_reaction.name
            if name is None:
                name = "Unknown"

            reversibility = True
            if model_reaction.lower_bound >= 0:
                reversibility = False
            elif model_reaction.upper_bound <= 0:
                reversibility = False
                # Flip the reaction, not the best for data mapping, but otherwise the
                # reversibility is not correct.
                reaction_info = [(-coeff, node) for coeff, node in reaction_info]

            reaction = pathway.place_reaction_on_backbone(
                map=m,
                name=name,
                bigg_id=model_reaction.bigg_id,
                reaction_info=reaction_info,
                reversibility=reversibility,
                **self.placement_kwargs,
            )
            m.add_reaction(reaction)

        return m


class EscherTemplateModuleDefinition(EscherModuleDefinition):
    def __init__(
        self,
        bigg_id: str,
        template_filename: Union[str, os.PathLike],
        bigg_ids_in_template_are_universal: bool = False,
        required_fraction_matching_reactions: float = 0.5,
        **kwargs,
    ):
        self.bigg_id = bigg_id
        self.template_filename = Path(template_filename)
        self.bigg_ids_in_template_are_universal = bigg_ids_in_template_are_universal
        self.required_fraction_matching_reactions = required_fraction_matching_reactions
        super().__init__(bigg_id=bigg_id, **kwargs)

    def load_template(self, session: Session):
        full_path = get_data("escher_templates", self.template_filename)
        with open(full_path, "r") as f:
            m, reactions, nodes = template.load_as_template(f)

        if self.bigg_ids_in_template_are_universal:
            node_bigg_ids = [
                parse.split_compartment(node.bigg_id)[0]
                for node in nodes.values()
                if node.node_type == "metabolite"
            ]
            bigg_id_mapping_db = session.execute(
                select(ComponentIDMapping.old_bigg_id, UniversalComponent.bigg_id)
                .join(ComponentIDMapping.new_universal_component)
                .filter(ComponentIDMapping.old_bigg_id.in_(node_bigg_ids))
            ).all()
            bigg_id_mapping = {
                row.old_bigg_id: row.bigg_id for row in bigg_id_mapping_db
            }
            for node in nodes.values():
                if node.node_type != "metabolite":
                    continue
                old_bigg_id, compartment = parse.split_compartment(node.bigg_id)
                if old_bigg_id not in bigg_id_mapping:
                    continue
                new_bigg_id = f"{bigg_id_mapping[old_bigg_id]}_{compartment}"
                node.bigg_id = new_bigg_id

        return m, reactions, nodes

    def select_model_reactions(
        self, session: Session, model_reactions: Iterable[ModelReaction]
    ) -> List[ModelReaction]:
        m, reactions, _ = self.load_template(session)

        selected_model_reactions = []

        def _select_model_reactions(
            mr: ModelReaction,
            rx: map.Reaction,
            _mr_bigg_ids: Dict[str, str],
            _mr_matrix_dbs: Dict[str, ReactionMatrix],
        ):
            selected_model_reactions.append(mr)
            rx.finalized = True

        self.map_template_to_model_reactions(
            model_reactions, m, reactions, _select_model_reactions
        )
        if (
            len(selected_model_reactions) / len(reactions)
        ) < self.required_fraction_matching_reactions:
            return []
        return selected_model_reactions

    def _build_reaction(
        self,
        model_reaction: ModelReaction,
        reaction: map.AutoReactionWithOptionalMetabolites,
        model_reaction_bigg_ids: Dict[str, str],
        model_reaction_matrix_dbs: Dict[str, ReactionMatrix],
    ):

        for i in range(len(reaction.metabolites)):
            coefficient, node = reaction.metabolites[i]
            bigg_id = node.bigg_id

            reaction_matrix = model_reaction_matrix_dbs[bigg_id]

            if (
                new_coefficient := reaction_matrix.universal_reaction_matrix.coefficient
            ) != coefficient:
                reaction.metabolites[i] = (new_coefficient, node)

            name = reaction_matrix.compartmentalized_component.component.name
            if name is None:
                name = (
                    reaction_matrix.compartmentalized_component.component.universal_component.name
                )
            if name is not None:
                node.name = name
            del model_reaction_bigg_ids[bigg_id]
            del model_reaction_matrix_dbs[bigg_id]

        absolute_sides = []
        default_placement_opts = map.PlacementOptions()
        for bigg_id, metabolite_info in reaction.optional_metabolites.items():
            # Calculate which side of the reaction metabolites are shown to
            # provide smarter placement of additional metabolites.

            ref_node = reaction.mid_marker
            plus_minus = metabolite_info["coefficient"] > 0
            if reaction.multi_markers[plus_minus] is not None:
                ref_node = reaction.multi_markers[plus_minus]
            _, _, _, _, _, effictive_angle_delta = reaction.calculate_placement(
                ref_node,
                metabolite_info["node"],
                plus_minus,
                0,
                0,
                metabolite_info.get("b1_b2"),
                placement_opts=default_placement_opts,
            )
            absolute_side = (effictive_angle_delta >= 0) == bool(plus_minus)
            absolute_sides.append(absolute_side)

            if bigg_id not in model_reaction_bigg_ids:
                continue
            reaction_matrix = model_reaction_matrix_dbs[bigg_id]

            node = metabolite_info["node"]
            if (
                new_coefficient := reaction_matrix.universal_reaction_matrix.coefficient
            ) != metabolite_info["coefficient"]:
                metabolite_info["coefficient"] = new_coefficient

            name = reaction_matrix.compartmentalized_component.component.name
            if name is None:
                name = (
                    reaction_matrix.compartmentalized_component.component.universal_component.name
                )
            if name is not None:
                node.name = name

            reaction.add_metabolite(**metabolite_info)
            del model_reaction_bigg_ids[bigg_id]
            del model_reaction_matrix_dbs[bigg_id]

        placement_opts = None
        if len(set(absolute_sides)) == 1:
            placement_opts = map.PlacementOptions(
                placement_f=partial(
                    map.AutoReaction.same_side_placement,
                    absolute_side=absolute_sides[0],
                )
            )
        for bigg_id, reaction_matrix in model_reaction_matrix_dbs.items():
            cc_db = reaction_matrix.compartmentalized_component
            name = cc_db.component.name
            if name is None:
                name = cc_db.component.universal_component.name
            if name is None:
                name = "Unknown"
            node = map.MetaboliteNode(
                bigg_id=cc_db.bigg_id,
                name=name,
                node_is_primary=False,
            )
            reaction.add_metabolite(
                coefficient=reaction_matrix.universal_reaction_matrix.coefficient,
                node=node,
                placement_opts=placement_opts,
            )
        reaction.bigg_id = model_reaction.bigg_id
        reaction.finalized = True

    def map_template_to_model_reactions(
        self,
        model_reactions: Iterable[ModelReaction],
        m: map.Map,
        reactions: List[map.AutoReactionWithOptionalMetabolites],
        build_reaction_f: Callable[
            [
                ModelReaction,
                map.Reaction,
                Dict[str, str],
                Dict[str, ReactionMatrix],
            ],
            None,
        ],
    ) -> Dict[str, str]:
        bigg_id_mapping = {}

        for model_reaction in model_reactions:
            if (
                self.exclude_pseudo
                and model_reaction.reaction.universal_reaction.is_pseudo
            ):
                continue

            model_reaction_bigg_ids = {}
            model_reaction_matrix_dbs = {}
            for reaction_matrix in model_reaction.reaction.matrix:
                if self.bigg_ids_in_template_are_universal:
                    bigg_id = (
                        reaction_matrix.compartmentalized_component.universal_compartmentalized_component.bigg_id
                    )
                    model_cc_bigg_id = reaction_matrix.compartmentalized_component.model_compartmentalized_components[
                        0
                    ].bigg_id
                    model_reaction_bigg_ids[bigg_id] = model_cc_bigg_id
                else:
                    bigg_id = reaction_matrix.compartmentalized_component.model_compartmentalized_components[
                        0
                    ].bigg_id
                    model_reaction_bigg_ids[bigg_id] = bigg_id
                model_reaction_matrix_dbs[bigg_id] = reaction_matrix

            for bigg_id_1, bigg_id_2 in model_reaction_bigg_ids.items():
                if bigg_id_1 not in bigg_id_mapping:
                    continue
                if bigg_id_mapping[bigg_id_1] == bigg_id_2:
                    continue
                logging.error(
                    f"BiGG ID Clash: map BiGG ID: {bigg_id_1} -> new {bigg_id_2} (previous mapping: {bigg_id_mapping[bigg_id_1]})"
                )
                break
            else:
                for reaction in reactions:
                    if reaction.finalized:
                        continue
                    for _, metabolite in reaction.metabolites:
                        if metabolite.bigg_id not in model_reaction_bigg_ids:
                            break
                    else:
                        bigg_id_mapping.update(model_reaction_bigg_ids)
                        build_reaction_f(
                            model_reaction,
                            reaction,
                            model_reaction_bigg_ids,
                            model_reaction_matrix_dbs,
                        )
                        break
                # else:
                #     logging.error(f"Could not match reaction: {model_reaction.bigg_id}")
        return bigg_id_mapping

    def build_map(
        self, session: Session, model_reactions: Iterable[ModelReaction]
    ) -> map.Map:
        m, reactions, _ = self.load_template(session)

        bigg_id_mapping = self.map_template_to_model_reactions(
            model_reactions, m, reactions, self._build_reaction
        )

        for reaction in reactions:
            if not reaction.finalized:
                continue
            m.add_reaction(reaction)

        for node in m.nodes.values():
            if node.node_type == "metabolite":
                # assert_type(node, map.MetaboliteNode)
                node.bigg_id = bigg_id_mapping.get(node.bigg_id, node.bigg_id)

        return m


ESCHER_MODULE_DEFINITION_LIST = [
    EscherBackboneModuleDefinition(
        bigg_id="ubiquinone",
        name="Ubiquinone Reduction/Oxidation",
        description="Auto-generated map of all reactions involved in the reversible redox cycle of ubiquinone-8.",
        backbone=[
            map.MetaboliteNode(
                bigg_id="q8h2_c:0",
                name="ubiquinol-8",
                node_is_primary=True,
                x=0,
                y=500,
            ),
            map.MetaboliteNode(
                bigg_id="q8_c:0",
                name="ubiquinone-8",
                node_is_primary=True,
                x=0,
                y=-500,
            ),
        ],
        placement_kwargs=dict(
            placement_f=partial(pathway.alternating_pathways_sides, centered=False),
            add_metabolite_opts=dict(
                placement_opts=map.PlacementOptions(
                    placement_f=map.AutoReaction.same_side_placement
                )
            ),
        ),
    ),
    EscherBackboneModuleDefinition(
        bigg_id="menaquinone",
        name="Menaquinone Reduction/Oxidation",
        description="Auto-generated map of all reactions involved in the reversible redox cycle of menaquinone-8.",
        backbone=[
            map.MetaboliteNode(
                bigg_id="mql8_c:0",
                name="ubiquinol-8",
                node_is_primary=True,
                x=0,
                y=500,
            ),
            map.MetaboliteNode(
                bigg_id="mqn8_c:0",
                name="ubiquinone-8",
                node_is_primary=True,
                x=0,
                y=-500,
            ),
        ],
        placement_kwargs=dict(
            placement_f=partial(pathway.alternating_pathways_sides, centered=False),
            add_metabolite_opts=dict(
                placement_opts=map.PlacementOptions(
                    placement_f=map.AutoReaction.same_side_placement
                )
            ),
        ),
    ),
    EscherBackboneModuleDefinition(
        bigg_id="2demethylmenaquinone",
        name="2-Demethylmenaquinone Reduction/Oxidation",
        description="Auto-generated map of all reactions involved in the reversible redox cycle of 2-demethylmenaquinone-8.",
        backbone=[
            map.MetaboliteNode(
                bigg_id="2dmmql8_c:0",
                name="ubiquinol-8",
                node_is_primary=True,
                x=0,
                y=500,
            ),
            map.MetaboliteNode(
                bigg_id="2dmmq8_c:0",
                name="ubiquinone-8",
                node_is_primary=True,
                x=0,
                y=-500,
            ),
        ],
        placement_kwargs=dict(
            placement_f=partial(pathway.alternating_pathways_sides, centered=False),
            add_metabolite_opts=dict(
                placement_opts=map.PlacementOptions(
                    placement_f=map.AutoReaction.same_side_placement
                )
            ),
        ),
    ),
    EscherTemplateModuleDefinition(
        bigg_id="ecoli_central_metabolism",
        template_filename="simple_central_metabolism.json",
        bigg_ids_in_template_are_universal=True,
        required_fraction_matching_reactions=0.5,
        name="E. coli Central Metabolism",
        description="Model reactions mapped onto a map of the central metabolism of E. coli.",
    ),
    EscherTemplateModuleDefinition(
        bigg_id="nucleotide_histidine_biosynthesis",
        template_filename="nucleotide_histidine_biosynthesis.json",
        bigg_ids_in_template_are_universal=True,
        required_fraction_matching_reactions=0.6,
        name="Nucleotide and Histidine Biosynthesis",
        description="Model reactions mapped onto a map of nucleotide and histidine biosynthesis in E. coli.",
    ),
    EscherBackboneModuleDefinition(
        bigg_id="arginine_biosynthesis",
        name="L-Arginine Biosynthesis (via L-Ornithine)",
        description="Auto-generated map of all reactions interconverting metabolites in the L-arginine biosynthesis pathway, via L-ornithine and L-citrulline.",
        sequential=True,
        backbone=[
            map.MetaboliteNode(
                bigg_id="glu__L_c:-1",
                name="L-glutamate",
                node_is_primary=True,
                x=-1500,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="acglu_c:-2",
                name="N-acetylglutamate",
                node_is_primary=True,
                x=-1000,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="acg5p_c:-3",
                name="N-acetyl-L-gamma-glutamyl phosphate",
                node_is_primary=True,
                x=-500,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="acg5sa_c:-1",
                name="2-acetamido-5-oxopentanoate",
                node_is_primary=True,
                x=0,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="acorn_c:0",
                name="N(2)-acetyl-L-ornithine",
                node_is_primary=True,
                x=500,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="orn__L_c:1",
                name="L-ornithine",
                node_is_primary=True,
                x=500,
                y=600,
            ),
            map.MetaboliteNode(
                bigg_id="citr__L_c:0",
                name="L-citrulline",
                node_is_primary=True,
                x=0,
                y=600,
            ),
            map.MetaboliteNode(
                bigg_id="argsuc_c:-1",
                name="L-argininosuccinate",
                node_is_primary=True,
                x=-500,
                y=600,
            ),
            map.MetaboliteNode(
                bigg_id="arg__L_c:1",
                name="L-arginine",
                node_is_primary=True,
                x=-1000,
                y=600,
            ),
        ],
        placement_kwargs=dict(
            placement_f=partial(pathway.alternating_pathways_sides, centered=True),
            add_metabolite_opts=dict(
                placement_opts=map.PlacementOptions(
                    placement_f=map.AutoReaction.alternating_side_placement
                )
            ),
        ),
    ),
    EscherBackboneModuleDefinition(
        bigg_id="aspartate_asparagine_biosynthesis",
        name="L-Aspartate and L-Asparagine Biosynthesis",
        description="Auto-generated map of all reactions interconverting metabolites in the L-aspartate and L-asparagine biosynthesis pathway.",
        sequential=True,
        backbone=[
            map.MetaboliteNode(
                bigg_id="oaa_c:-2",
                name="oxaloacetate",
                node_is_primary=True,
                x=-1500,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="asp__L_c:-1",
                name="L-aspartate",
                node_is_primary=True,
                x=-1000,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="asn__L_c:0",
                name="L-asparagine",
                node_is_primary=True,
                x=-500,
                y=0,
            ),
        ],
        placement_kwargs=dict(
            placement_f=partial(pathway.alternating_pathways_sides, centered=True),
            add_metabolite_opts=dict(
                placement_opts=map.PlacementOptions(
                    placement_f=map.AutoReaction.alternating_side_placement
                )
            ),
        ),
    ),
    EscherBackboneModuleDefinition(
        bigg_id="cysteine_biosynthesis",
        name="L-Cysteine Biosynthesis",
        description="Auto-generated map of all reactions interconverting metabolites in the L-cysteine biosynthesis pathway from L-serine.",
        sequential=False,
        backbone=[
            map.MetaboliteNode(
                bigg_id="ser__L_c:0",
                name="L-serine",
                node_is_primary=True,
                x=-500,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="acser_c:0",
                name="O-acetyl-L-serine",
                node_is_primary=True,
                x=-1000,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="slfcys:-1",
                name="S-sulfo-L-cysteine",
                node_is_primary=True,
                x=-1500,
                y=0,
            ),
            map.MetaboliteNode(
                bigg_id="cys__L:0",
                name="L-cysteine",
                node_is_primary=True,
                x=-1500,
                y=0,
            ),
        ],
        placement_kwargs=dict(
            placement_f=partial(pathway.alternating_pathways_sides, centered=True),
            add_metabolite_opts=dict(
                placement_opts=map.PlacementOptions(
                    placement_f=map.AutoReaction.alternating_side_placement
                )
            ),
        ),
    ),
]
ESCHER_MODULE_DEFINITIONS = {m.bigg_id: m for m in ESCHER_MODULE_DEFINITION_LIST}
