FROM ubuntu:xenial

RUN apt -y update

RUN apt -y install \
	git \
	wget \
	vim \
	zlib1g-dev \
	bzip2 \
	libbz2-dev \
	xz-utils \
    liblzma-dev \
	curl \
    libcurl4-openssl-dev \
	libssl-dev \
	ncurses-dev \
	graphviz \
    unzip \
    zip
    
RUN wget --quiet https://bioinfo.umassmed.edu/pub/beenet && \
    chmod 755 beenet && mv beenet /usr/bin/. 