from functools import partial
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union
from sqlalchemy import select
from sqlalchemy.orm import Session
from cobradb.models import EscherModule, ModelReaction
from biggr_maps import map, pathway


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
        self, model_reactions: Iterable[ModelReaction]
    ) -> List[ModelReaction]:
        raise NotImplementedError("Not implemented in base class.")

    def build_map(self, model_reactions: Iterable[ModelReaction]):
        m = map.Map(
            name=self.name,
            description=self.description,
            canvas=(-1000, -1000, 2000, 2000),
        )
        return m


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
        self, model_reactions: Iterable[ModelReaction]
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
        self, model_reaction: ModelReaction
    ) -> List[Tuple[Union[int, float], map.MetaboliteNode]]:
        reactants = []
        for reaction_matrix in model_reaction.reaction.matrix:
            cc_db = reaction_matrix.compartmentalized_component
            if (node := self.backbone.get(cc_db.bigg_id)) is None:
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
        return reactants

    def build_map(self, model_reactions: Iterable[ModelReaction]) -> map.Map:
        m = super().build_map(model_reactions)

        for node in self.backbone.values():
            node.identifier = None
            m.add_node(node)

        for model_reaction in model_reactions:
            reaction_info = self._build_reaction(model_reaction)
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
            placement_f=partial(pathway.alternating_pathways_sides, centered=False)
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
            placement_f=partial(pathway.alternating_pathways_sides, centered=False)
        ),
    ),
]
ESCHER_MODULE_DEFINITIONS = {m.bigg_id: m for m in ESCHER_MODULE_DEFINITION_LIST}
