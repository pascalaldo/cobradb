from os import PathLike
from typing import Union
from sqlalchemy import delete, select, update
from sqlalchemy.orm import Session

from cobradb.models import Taxon, TaxonomicRank


def _parse_nodes_line(line: str):
    try:
        tax_id, parent_tax_id, rank, _ = line.split("|", maxsplit=3)
        tax_id = int(tax_id.strip())
        rank = rank.strip()
        parent_tax_id = int(parent_tax_id.strip())
    except:
        print(f"Failed to parse line: {line.strip()}")
        return None, None, None
    return tax_id, parent_tax_id, rank


def load_taxonomy(
    session: Session,
    taxonomy_nodes_path: Union[str, PathLike],
    taxonomy_names_path: Union[str, PathLike],
    batch_size: int = 10000,
):
    print("Loading taxonomies")
    session.execute(delete(Taxon))
    session.commit()

    ranks = {
        name: rank_id
        for rank_id, name in session.execute(
            select(TaxonomicRank.id, TaxonomicRank.name)
        ).all()
    }
    with open(taxonomy_nodes_path, "r") as f:
        for i, line in enumerate(f):
            tax_id, _parent_tax_id, rank = _parse_nodes_line(line)
            if tax_id is None:
                continue
            if rank not in ranks:
                rank_db = TaxonomicRank(name=rank)
                session.add(rank_db)
                session.commit()
                ranks[rank] = rank_db.id
            rank_id = ranks[rank]

            session.add(
                Taxon(
                    id=tax_id,
                    rank_id=rank_id,
                )
            )
            if (i + 1) % batch_size == 0:
                print(f"Creating entry: {i}")
                session.commit()
    session.commit()
    session.close()
    with open(taxonomy_nodes_path, "r") as f:
        updates = []
        for i, line in enumerate(f):
            tax_id, parent_tax_id, _rank = _parse_nodes_line(line)
            if tax_id is None:
                continue

            updates.append({"id": tax_id, "parent_id": parent_tax_id})
            if (i + 1) % batch_size == 0:
                print(f"Updating parent id: {i}")
                session.execute(update(Taxon), updates)
                updates = []
                session.commit()
        if updates:
            session.execute(update(Taxon), updates)

    session.commit()
    with open(taxonomy_names_path, "r") as f:
        updates = []
        i = 0
        for line in f:
            try:
                tax_id, name, _, name_class, _ = line.split("|", maxsplit=4)
                tax_id = int(tax_id.strip())
                name = name.strip()
                name_class = name_class.strip()
            except:
                print(f"Failed to parse line: {line.strip()}")
                continue
            if name_class != "scientific name":
                continue
            updates.append({"id": tax_id, "name": name})
            if (i + 1) % batch_size == 0:
                print(f"Updating name: {i}")
                session.execute(update(Taxon), updates)
                updates = []
                session.commit()
            i = i + 1
        if updates:
            session.execute(update(Taxon), updates)

    session.commit()
    session.close()
