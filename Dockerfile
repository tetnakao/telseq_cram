FROM ubuntu:24.04 AS build

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        autotools-dev \
        build-essential \
        ca-certificates \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libsparsehash-dev \
        libssl-dev \
        pkg-config \
        wget \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# htslib
RUN mkdir -p /deps && \
    cd /deps && \
    wget -q https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 -O htslib.tar.bz2 && \
    tar -xjf htslib.tar.bz2 && \
    cd htslib-1.18 && \
    ./configure && \
    make -j"$(nproc)" && \
    make install

# telseq
COPY src /src
RUN cd /src && \
    ./autogen.sh && \
    ./configure --prefix=/usr/local && \
    make -j"$(nproc)" && \
    make install


FROM ubuntu:24.04

LABEL Description="Telseq docker" Version="0.0.5"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libbz2-1.0 \
        libcurl4t64 \
        liblzma5 \
        libssl3t64 \
        libstdc++6 \
        zlib1g && \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /usr/local/bin /usr/local/bin
COPY --from=build /usr/local/lib /usr/local/lib
RUN ldconfig && ! ldd /usr/local/bin/telseq | grep "not found"

ENV LD_LIBRARY_PATH=/usr/local/lib
CMD ["/bin/bash"]
