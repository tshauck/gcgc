ARG PY_VERSION=3.6
FROM python:$PY_VERSION
RUN echo $PY_VERSION

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
            build-essential \
            cmake \
            libgoogle-perftools-dev \
            pkg-config \
            python3-dev \
            wget \
        && rm -rf /var/lib/apt/lists/*

RUN wget -O /tmp/piece.tar.gz https://github.com/google/sentencepiece/archive/v0.1.85.tar.gz

WORKDIR /tmp/
RUN tar -xzf /tmp/piece.tar.gz

WORKDIR /tmp/sentencepiece-0.1.85
RUN mkdir build && cd build && cmake .. && make && make install && ldconfig -v

RUN pip install pip==19.0.1 poetry==0.12.17

WORKDIR /gcgc
COPY ./pyproject.toml ./poetry.lock ./

RUN poetry config settings.virtualenvs.create false
RUN poetry install
RUN poetry cache:clear pypi --all

WORKDIR /gcgc
COPY ./ ./
RUN poetry install

WORKDIR /workspace
ENTRYPOINT ["gcgc"]
