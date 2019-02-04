FROM python:3.6

RUN apt-get update && apt-get install -y --no-install-recommends

RUN apt-get install -y \
        build-essential \
        gcc \
        gfortran \
        python3-dev \
        wget \
        libpng-dev \
        git

RUN pip install -U pip
RUN pip install poetry

WORKDIR /gcgc
COPY ./pyproject.toml ./poetry.lock ./

RUN poetry config settings.virtualenvs.create false
RUN poetry install
RUN poetry cache:clear pypi --all

WORKDIR /gcgc
COPY ./ ./
RUN poetry install \
        && find -name "*.pyc" -delete \
        && find -name "__pycache__" -delete

WORKDIR /
