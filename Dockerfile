FROM nfcore/base
LABEL authors="davismcc@gmail.com" \
    maintainer="Davis McCarthy <davismcc@gmail.com>" \
    description="Docker image containing all requirements for davismcc/nf-hipsci-endo pipeline"

# Install container-wide requrements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
        curl \
        g++ \
        gcc \
        gfortran \
        git \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgsl-dev \
        libgsl2 \
        liblzma-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssh2-1-dev \
        libssl-dev \
        libxml2-dev \
        libzmq3-dev \
        make \
        pandoc \
        pandoc-citeproc \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/davismcc-nf-hipsci-endo/bin:$PATH

# FROM rocker/verse:3.5.0
# RUN mkdir -p /usr/local/lib/R/site-library
# ADD install.R /tmp/
# RUN R -f /tmp/install.R
