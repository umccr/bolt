FROM continuumio/miniconda3:24.5.0-0 as build

RUN \
  conda install conda-libmamba-solver

RUN \
  echo -e > ~/.condarc '\
solver: libmamba\n\
channels:\n\
  - conda-forge\n\
  - bioconda\n\
  - defaults'

COPY ./conda/env/bolt_env.yml /tmp/
RUN \
  conda env update --prefix /env/ --file /tmp/bolt_env.yml

COPY ./ /tmp/bolt/
RUN \
  conda run --prefix /env/ pip install --no-deps /tmp/bolt/

RUN \
  conda install --prefix /env/ bash

FROM quay.io/bioconda/base-glibc-busybox-bash:2.1.0

COPY --from=build /env/ /env/

ENV PATH="/env/bin:${PATH}"
