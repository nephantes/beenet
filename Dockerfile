FROM ubuntu:22.04
LABEL author="alper.kucukural@umassmed.edu"  description="Docker image containing all requirements for the dolphinnext/honeycomb pipeline"

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git
    
RUN wget --quiet https://bioinfo.umassmed.edu/pub/beenet && \
    chmod 755 beenet && mv beenet /usr/bin/. 