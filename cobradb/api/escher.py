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
    def __init__(self, bigg_id: str, name: str, description: str):
        self.bigg_id: str = bigg_id
        self.name: str = name
        self.description: str = description
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
        placement_kwargs: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        self.backbone = {b.bigg_id: b for b in backbone}
        self.placement_kwargs = {} if placement_kwargs is None else placement_kwargs
        super().__init__(**kwargs)

    def select_model_reactions(
        self, session: Session, model_reactions: Iterable[ModelReaction]
    ) -> List[ModelReaction]:
        selected_model_reactions = []
        for model_reaction in model_reactions:
            matched_bigg_ids = set()
            for reaction_matrix in model_reaction.reaction.matrix:
                cc_bigg_id = reaction_matrix.compartmentalized_component.bigg_id
                if cc_bigg_id in self.backbone:
                    matched_bigg_ids.add(cc_bigg_id)
                if len(matched_bigg_ids) >= 2:
                    break
            if len(matched_bigg_ids) >= 2:
                selected_model_reactions.append(model_reaction)
        return selected_model_reactions

    def _build_reaction(
        self, model_reaction: ModelReaction, backbone: Dict[str, map.MetaboliteNode]
    ) -> List[Tuple[Union[int, float], map.MetaboliteNode]]:
        reactants = []
        for reaction_matrix in model_reaction.reaction.matrix:
            cc_db = reaction_matrix.compartmentalized_component
            if (node := backbone.get(cc_db.bigg_id)) is None:
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
            reaction_info = self._build_reaction(model_reaction, backbone)
            name = model_reaction.reaction.universal_reaction.name
            if name is None:
                name = "Unknown"
            pathway.place_reaction_on_backbone(
                map=m,
                name=name,
                bigg_id=model_reaction.bigg_id,
                reaction_info=reaction_info,
                **self.placement_kwargs,
            )

        return m


class EscherTemplateModuleDefinition(EscherModuleDefinition):
    def __init__(
        self,
        bigg_id: str,
        template_filename: Union[str, os.PathLike],
        bigg_ids_in_template_are_universal: bool = False,
        required_fraction_matching_reactions: float = 0.5,
    ):
        self.bigg_id = bigg_id
        self.template_filename = Path(template_filename)
        self.bigg_ids_in_template_are_universal = bigg_ids_in_template_are_universal
        self.required_fraction_matching_reactions = required_fraction_matching_reactions
        super().__init__(bigg_id=bigg_id, name="Template", description="Template")

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
                else:
                    logging.error(f"Could not match reaction: {model_reaction.bigg_id}")
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
        description="Ubiquinone Reduction/Oxidation",
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
        description="Menaquinone Reduction/Oxidation",
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
    EscherTemplateModuleDefinition(
        bigg_id="central_metabolism_1",
        template_filename="iJO1366.Central metabolism.json",
        bigg_ids_in_template_are_universal=True,
        required_fraction_matching_reactions=0.0,
    ),
]
ESCHER_MODULE_DEFINITIONS = {m.bigg_id: m for m in ESCHER_MODULE_DEFINITION_LIST}
