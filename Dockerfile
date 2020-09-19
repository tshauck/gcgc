ARG PY_VERSION=3.7
FROM python:$PY_VERSION
RUN echo $PY_VERSION

RUN pip install -U pip poetry

WORKDIR /gcgc
COPY ./pyproject.toml ./poetry.lock ./

RUN POETRY_VIRTUALENVS_CREATE=false poetry install -E third_party

WORKDIR /gcgc
COPY ./ ./
RUN poetry install -E third_party

WORKDIR /workspace
