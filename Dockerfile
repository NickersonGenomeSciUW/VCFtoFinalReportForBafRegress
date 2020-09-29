FROM ubuntu:latest

ARG bcftools_version=1.11

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    r-base \
    r-cran-randomforest \
    python \
    libbz2-dev zlib1g-dev liblzma-dev \
    wget

WORKDIR /root/tools

# Install BafRegress
RUN wget 'https://genome.sph.umich.edu/w/images/d/d5/BafRegress.tar.gz' && \
  tar xf BafRegress.tar.gz && \
  rm BafRegress.tar.gz

# Install bcftools
RUN wget -q https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 && \
  tar xf bcftools-${bcftools_version}.tar.bz2 && \
  rm bcftools-${bcftools_version}.tar.bz2 && \
  cd bcftools-${bcftools_version} && \
  ./configure --prefix=/root/tools/bcftools && \
  make && \
  make install && \
  cd ../ && \
  rm -r bcftools-${bcftools_version}  # Don't need the source directory now that we've installed the binary.

RUN mkdir /root/data
COPY GDA.MAF.txt.gz /root/data
RUN gzip -d /root/data/GDA.MAF.txt.gz
COPY parseVcfToBAFRegress.py ./parseVcfToBAFRegress.py





