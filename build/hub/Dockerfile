FROM python

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-key C99B11DEB97541F0 && \
    apt-add-repository https://cli.github.com/packages && \
    apt update && \
    apt install gh && \
    python -m pip install --pre -U commitizen

ENTRYPOINT ["/bin/bash", "/src/build/hub/make-release.sh"]
