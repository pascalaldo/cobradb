import logging
from sqlalchemy import select
from sqlalchemy.orm import Session, joinedload, subqueryload

from cobradb.api.escher import ESCHER_MODULE_DEFINITIONS
from cobradb.models import (
    CompartmentalizedComponent,
    Model,
    ModelCompartmentalizedComponent,
    ModelReaction,
    ModelReactionEscherMapping,
    Reaction,
    ReactionMatrix,
)


def create_escher_modules(session: Session):
    models_db = session.scalars(select(Model)).all()
    for model_db in models_db:
        logging.warning(f"Creating escher modules for model {model_db.bigg_id}")

        model_reactions = session.scalars(
            select(ModelReaction)
            .options(
                joinedload(ModelReaction.reaction)
                .subqueryload(Reaction.matrix)
                .joinedload(ReactionMatrix.compartmentalized_component)
                .subqueryload(
                    CompartmentalizedComponent.model_compartmentalized_components.and_(
                        ModelCompartmentalizedComponent.model_id == model_db.id
                    )
                )
            )
            .filter(ModelReaction.model_id == model_db.id)
            .order_by(ModelReaction.bigg_id)
        ).all()

        for escher_module in ESCHER_MODULE_DEFINITIONS.values():
            logging.warning(f"Escher Module: {escher_module.bigg_id}")
            escher_module_db = escher_module.get_or_create(session)
            selected_reactions = escher_module.select_model_reactions(
                session, model_reactions
            )
            print(f"Selected reactions: {len(selected_reactions)}")

            for r in selected_reactions:
                r.escher_mappings.append(
                    ModelReactionEscherMapping(escher_module=escher_module_db)
                )
        session.commit()
