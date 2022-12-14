FROM ubuntu:22.04
LABEL author="alper.kucukural@umassmed.edu"  description="Docker image containing all requirements for the dolphinnext/honeycomb pipeline"

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git
    
RUN wget --cipher 'DEFAULT:!DH' https://bioinfo.umassmed.edu/pub/beenet 
RUN wget --cipher 'DEFAULT:!DH' https://bioinfo.umassmed.edu/pub/STAR
RUN chmod 755 beenet && mv beenet /usr/bin/. 
RUN chmod 755 beenet && mv STAR /usr/bin/. 