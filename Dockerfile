FROM nfcore/base
MAINTAINER Davis McCarthy <davismcc@gmail.com>
LABEL authors="davismcc@gmail.com" \
    description="Docker image containing all requirements for davismcc/nf-hipsci-endo pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/davismcc-nf-hipsci-endo/bin:$PATH
