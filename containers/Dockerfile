# Container for the Nextflow version of the ATX epigenomic preprocessing
# pipeline. Build from the repository root:
#
#   docker build -f containers/Dockerfile -t atx-fastq2frags:$(cat version) .
#
# Scripts and reference assets live under /root, because
# scripts/atac_analysis_summary_prep.py hard-codes those paths.
ARG BASE_IMAGE=ubuntu:22.04
FROM ${BASE_IMAGE}

ARG BBMAP_VERSION=39.01
ARG CHROMAP_REF=li_dev4

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        default-jre-headless \
        git \
        libz-dev \
        pigz \
        procps \
        python3 \
        python3-dev \
        python3-pip \
        python-is-python3 \
        tabix && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /root

# BBMap (bbduk.sh)
RUN curl -fsSL "https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz/download" \
        -o bbmap.tar.gz && \
    tar -xzf bbmap.tar.gz && \
    rm bbmap.tar.gz
ENV PATH="/root/bbmap:${PATH}"

# Chromap, AtlasXomics branch supporting multi-segment barcodes (--read-format bc:..,bc:..)
RUN git clone --depth 1 -b "${CHROMAP_REF}" https://github.com/haowenz/chromap.git /root/chromap && \
    make -C /root/chromap && \
    ln -s /root/chromap/chromap /usr/local/bin/chromap

COPY requirements.txt /root/requirements.txt
RUN python3 -m pip install --no-cache-dir --upgrade pip && \
    python3 -m pip install --no-cache-dir -r /root/requirements.txt

COPY barcodes /root/barcodes
COPY blacklist /root/blacklist
COPY chrom_sizes /root/chrom_sizes
COPY references /root/references
COPY scripts /root/scripts
RUN gcc -O3 -o /root/scripts/substitute_filter /root/scripts/substitute_filter.c

# Fragment filtering and lane QC helpers used by bin/ scripts.
COPY lib /root/lib
ENV PYTHONPATH="/root/lib:/root/scripts"
ENV NUMEXPR_MAX_THREADS=42

COPY version /root/version

# Tasks may run as a non-root user (docker -u, Singularity): make /root
# readable and point library caches at a writable location.
RUN chmod -R a+rX /root
ENV MPLCONFIGDIR=/tmp/mplconfig \
    NUMBA_CACHE_DIR=/tmp/numba_cache
