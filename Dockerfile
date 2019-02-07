FROM python:3.6

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
            build-essential \
            gcc \
            gfortran \
            python3-dev \
            wget \
            libpng-dev \
            git \
        && rm -rf /var/lib/apt/lists/*

RUN pip install pip==19.0.1
RUN pip install poetry==0.12.11

WORKDIR /gcgc
COPY ./pyproject.toml ./poetry.lock ./

RUN poetry config settings.virtualenvs.create false
RUN poetry install
RUN poetry cache:clear pypi --all

WORKDIR /gcgc
COPY ./ ./
RUN poetry install \
        && find . -name "*.pyc" -delete \
        && find . -name "__pycache__" -delete

ENTRYPOINT ["gcgc"]
