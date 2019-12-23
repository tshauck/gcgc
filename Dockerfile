FROM python:3.7

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
            build-essential \
            cmake \
            libgoogle-perftools-dev \
            pkg-config \
            python3-dev \
            wget \
        && rm -rf /var/lib/apt/lists/*

RUN pip install pip==19.3.1

WORKDIR /gcgc

COPY ./dev-requirements.txt ./dev-requirements.txt
RUN pip install -r ./dev-requirements.txt

RUN wget -O /tmp/piece.tar.gz https://github.com/google/sentencepiece/archive/v0.1.85.tar.gz

WORKDIR /tmp/
RUN tar -xzf /tmp/piece.tar.gz

WORKDIR /tmp/sentencepiece-0.1.85
RUN mkdir build && cd build && cmake .. && make && make install && ldconfig -v

WORKDIR /gcgc
COPY ./ ./
RUN pip install .[sentencepiece]

ENTRYPOINT ["gcgc"]
