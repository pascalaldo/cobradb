from typing import Any, Dict, List
from sqlalchemy import and_, literal, or_, select
from sqlalchemy.orm import Session

from cobradb.models import ModelCollection, Taxon


def find_common_parent_for_taxons(session: Session, tax_id_1: int, tax_id_2: int):
    # This could maybe be done more efficiently, but I'm not that good with SQL
    topq_1 = select(Taxon.id, Taxon.parent_id)
    topq_1 = topq_1.filter(Taxon.id == tax_id_1)
    topq_1 = topq_1.cte("cte_1", recursive=True)

    topq_2 = select(Taxon.id, Taxon.parent_id)
    topq_2 = topq_2.filter(Taxon.id == tax_id_2)
    topq_2 = topq_2.cte("cte_2", recursive=True)

    bottomq_1 = select(Taxon.id, Taxon.parent_id)
    bottomq_1 = bottomq_1.join(topq_1, Taxon.id == topq_1.c.parent_id)

    bottomq_2 = select(Taxon.id, Taxon.parent_id)
    bottomq_2 = bottomq_2.join(topq_2, Taxon.id == topq_2.c.parent_id)

    recursive_q_1 = select(topq_1.union(bottomq_1)).cte("cte_3")
    recursive_q_2 = select(topq_2.union(bottomq_2)).cte("cte_4")

    recursive_q_1_id = select(recursive_q_1.c.id)
    recursive_q_2_id = select(recursive_q_2.c.id)
    q = (
        select(Taxon.parent_id)
        .filter(
            and_(
                Taxon.parent_id.in_(recursive_q_1_id),
                Taxon.parent_id.in_(recursive_q_2_id),
                or_(
                    and_(
                        Taxon.id.in_(recursive_q_1_id), ~Taxon.id.in_(recursive_q_2_id)
                    ),
                    and_(
                        ~Taxon.id.in_(recursive_q_1_id), Taxon.id.in_(recursive_q_2_id)
                    ),
                ),
            )
        )
        .limit(1)
    )
    return session.scalars(q).first()


def update_collection_data(session: Session, collection_list: List[Dict[str, Any]]):
    collections_db = session.scalars(select(ModelCollection))
    for collection in collections_db:
        collection_data = next(
            (x for x in collection_list if x.get("bigg_id") == collection.bigg_id), {}
        )

        tax_id = collection_data.get("tax_id")
        if tax_id is not None:
            if (
                session.scalars(
                    select(Taxon).filter(Taxon.id == tax_id).limit(1)
                ).first()
                is None
            ):
                tax_id = None

        if tax_id is not None:
            collection.taxon_id = tax_id
        else:
            if collection.taxon_id is None:
                # Find common parent
                parent_tax_id = None
                for model in collection.models:
                    if model.taxon_id is None:
                        continue
                    if parent_tax_id is None:
                        parent_tax_id = model.taxon_id
                    else:
                        parent_tax_id = find_common_parent_for_taxons(
                            session, parent_tax_id, model.taxon_id
                        )
                if parent_tax_id is not None:
                    collection.taxon_id = parent_tax_id

        collection.oneliner = collection_data.get("oneliner")
        collection.description = collection_data.get("description")

    session.commit()
