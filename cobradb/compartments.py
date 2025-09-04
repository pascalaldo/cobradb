import os
import logging

from cobradb.models import Compartment


def load_compartments(compartment_names_file, session):
    if os.path.exists(compartment_names_file):
        with open(compartment_names_file, "r") as f:
            compartment_names = {}
            for line in f.readlines():
                sp = [x.strip() for x in line.split("\t")]
                try:
                    compartment_names[sp[0]] = sp[1]
                except IndexError:
                    continue
    else:
        logging.warning("No compartment names file")
        compartment_names = {}

    for k, v in compartment_names.items():
        compartment_db = session.query(Compartment).filter(Compartment.id == k).first()
        if compartment_db is None:
            compartment_db = Compartment(id=k, name=v)
            session.add(compartment_db)
        else:
            compartment_db.name = v
        session.commit()
