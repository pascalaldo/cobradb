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
#CMD ["bin/load_db", "--drop-all"]
CMD ["bin/load_db", "--drop-all"]
