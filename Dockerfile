FROM python:3.10
ENV PYTHONUNBUFFERED=1
RUN mkdir /app
WORKDIR /app

# Install dependencies

RUN apt-get update && apt-get install -y \
  libxml2-dev

COPY ./ /app
RUN python setup.py install
# Install reqs  

#COPY settings.ini /app/settings.ini

