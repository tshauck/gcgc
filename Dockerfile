ARG PY_VERSION=3.7
FROM python:$PY_VERSION
RUN echo $PY_VERSION

RUN curl --proto '=https' --tlsv1.2 -sSf https://just.systems/install.sh | bash -s -- --to /usr/bin/

RUN pip install pip==19.0.1 poetry==1.0.5

WORKDIR /gcgc
COPY ./pyproject.toml ./poetry.lock ./

RUN POETRY_VIRTUALENVS_CREATE=false poetry install

WORKDIR /gcgc
COPY ./ ./
RUN poetry install

WORKDIR /workspace
ENTRYPOINT ["gcgc"]
