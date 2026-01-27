
FROM ubuntu:20.04
LABEL Description="Telseq docker" Version="0.0.3"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        automake \
        autoconf \
        autotools-dev \
        build-essential \
        cmake \
        libhts-dev \
        libjemalloc-dev \
        libsparsehash-dev \
        libz-dev \
        wget \
        zlib1g-dev \
        pkg-config \
        binutils \
        autoconf autotools-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev


#htslib
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 -O htslib.tar.bz2 &&  \
    tar -xjvf htslib.tar.bz2 && \
    cd htslib-1.18 && \
    autoconf -i && \
    ./configure && \
    make && \
    make install

# build telseq
RUN apt-get install -y gdb
COPY src src
RUN cd src && \
    ./autogen.sh && \
    ./configure --prefix=/usr/local && \
    make && \
    make install
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
CMD ["/bin/bash"]