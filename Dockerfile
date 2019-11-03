FROM python:3.6

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
            build-essential \
            python3-dev \
            wget \
        && rm -rf /var/lib/apt/lists/*

RUN pip install pip==19.3.1

WORKDIR /gcgc

COPY ./dev-requirements.txt ./dev-requirements.txt
RUN pip install -r ./dev-requirements.txt

COPY ./ ./
RUN python setup.py install
