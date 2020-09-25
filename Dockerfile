ARG PY_VERSION=3.7
FROM python:$PY_VERSION

RUN pip install -U pip poetry

WORKDIR /gcgc
COPY ./ ./

# Install dataclasses if we're on 3.6
RUN python3.6 -m pip install dataclass || true
RUN POETRY_VIRTUALENVS_CREATE=false poetry install -E pytorch -E hf

WORKDIR /workspace
