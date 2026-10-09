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
        tabix \
        unzip && \
    rm -rf /var/lib/apt/lists/*

# AWS CLI v2. Nextflow's AWS Batch executor uses it inside the task container
# to stage files to and from S3 (not needed when running with Fusion).
RUN curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscliv2.zip && \
    unzip -q /tmp/awscliv2.zip -d /tmp && \
    /tmp/aws/install && \
    rm -rf /tmp/aws /tmp/awscliv2.zip

WORKDIR /root

# BBMap (bbduk.sh)
RUN curl -fsSL "https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz/download" \
        -o bbmap.tar.gz && \
    tar -xzf bbmap.tar.gz && \
    rm bbmap.tar.gz
ENV PATH="/root/bbmap:${PATH}"

# Chromap, AtlasXomics branch supporting multi-segment barcodes (--read-format bc:..,bc:..)
# CHROMAP_SORT_BUFFER_GIB sets the fixed buffer Chromap allocates for sorting
# and deduplicating BED output (upstream: 10 GiB). Lower it only for small test
# runs on laptops; see "Local testing" in the README.
ARG CHROMAP_SORT_BUFFER_GIB=10
RUN git clone --depth 1 -b "${CHROMAP_REF}" https://github.com/haowenz/chromap.git /root/chromap && \
    if [ "${CHROMAP_SORT_BUFFER_GIB}" != "10" ]; then \
        grep -q 'max_mem_size = 10 \* ((uint64_t)1 << 30);' /root/chromap/src/mapping_writer.h && \
        sed -i "s|max_mem_size = 10 \* ((uint64_t)1 << 30);|max_mem_size = ${CHROMAP_SORT_BUFFER_GIB} * ((uint64_t)1 << 30);|" \
            /root/chromap/src/mapping_writer.h; \
    fi && \
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
