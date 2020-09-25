ARG PY_VERSION=3.7
FROM python:$PY_VERSION

RUN pip install -U pip

WORKDIR /gcgc
COPY ./ ./

RUN pip install -e .[dev,hf,pytorch]

WORKDIR /workspace
