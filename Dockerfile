ARG PY_VERSION=3.7
FROM python:$PY_VERSION

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install

RUN curl -sL https://taskfile.dev/install.sh | sh
RUN mv ./bin/task /usr/local/bin/task

RUN pip install -U pip torch==1.6

WORKDIR /gcgc

COPY ./ ./
RUN pip install -e .[dev,hf,pytorch]
