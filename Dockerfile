FROM nfcore/base
LABEL authors="davismcc@gmail.com" \
    maintainer="Davis McCarthy <davismcc@gmail.com>" \
    description="Docker image containing all requirements for davismcc/nf-hipsci-endo pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml python=3.6 && conda clean -a
ENV PATH /opt/conda/envs/davismcc-nf-hipsci-endo/bin:$PATH

FROM broadinstitute/gatk

# FROM rocker/verse
# FROM bioconductor/release_core2
# RUN mkdir -p /usr/local/lib/R/site-library
# ADD install.R /tmp/
# RUN R -f /tmp/install.R
