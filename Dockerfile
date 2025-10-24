FROM python:3.10
ENV PYTHONUNBUFFERED=1
RUN mkdir /app
WORKDIR /app

# Install dependencies

RUN apt-get update && apt-get install -y \
  libxml2-dev
COPY requirements.txt /app
RUN pip install -r requirements.txt

RUN git clone https://github.com/pascalaldo/bigg_models_data.git bigg_models_data

COPY ./ /app
RUN python setup.py install
# Install reqs  

#COPY settings.ini /app/settings.ini
#COPY settings.small.ini /app/settings.ini

# CMD ["bin/new_load_db", "--drop-all"]
# CMD ["bin/new_load_db", "--drop-all", "--skip-genomes", "--skip-curated-reactions", "--skip-models"]
# CMD ["bin/new_load_db", "--skip-genomes"]
# CMD ["bin/new_load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites"]
# CMD ["bin/new_load_db", "--skip-rhea", "--skip-curated-metabolites"]
CMD ["bin/new_load_db", "--drop-models", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-seed-metabolites", "--skip-seed-reactions"]
# CMD ["bin/new_load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-compartments", "--skip-models", "--skip-seed-metabolites", "--skip-seed-reactions"]
# CMD ["bin/new_load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-compartments", "--skip-models"]
# CMD ["bin/new_load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions"]
# CMD ["bin/new_load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-models", "--skip-memote", "--skip-compartments", "--skip-seed-metabolites"]
# CMD ["bin/new_load_db", "--drop-all", "--skip-models"]
# CMD ["bin/new_load_db", "--skip-rhea"]

# CMD ["bin/test_parts"]
