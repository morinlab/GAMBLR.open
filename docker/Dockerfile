# Mirrors envs/gamblr_collective.yaml plus the six GAMBLR packages
# (installed from GitHub on top), so the same environment is available as
# either a conda env (envs/gamblr_collective.yaml directly) or a container
# (this image, which Apptainer/Singularity can pull straight from GHCR via
# `apptainer pull docker://ghcr.io/morinlab/gamblr-collective:<tag>`) --
# see LCR-modules' two supported install modes.
#
# Built from envs/gamblr_collective.yaml itself (not a separately
# maintained package list), so this image can never drift out of sync with
# the conda environment it's meant to mirror.
#
# Build with docker/release_container.sh -- see docker/RELEASING.md.

FROM mambaorg/micromamba:2

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/gamblr_collective.yaml /tmp/gamblr_collective.yaml

RUN micromamba install -y -n base -f /tmp/gamblr_collective.yaml && \
    micromamba clean --all --yes

# Required by micromamba's own base image convention: activates the env
# for every RUN instruction from this point on.
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Which ref of each GAMBLR package to install -- default "master" for a
# normal release build; release_container.sh can override these for
# pre-merge testing of open PR branches.
ARG GAMBLR_DATA_REF=master
ARG GAMBLR_HELPERS_REF=master
ARG GAMBLR_UTILS_REF=master
ARG GAMBLR_VIZ_REF=master
ARG GAMBLR_PREDICT_REF=master
ARG GAMBLR_OPEN_REF=master
ENV GAMBLR_DATA_REF=${GAMBLR_DATA_REF} \
    GAMBLR_HELPERS_REF=${GAMBLR_HELPERS_REF} \
    GAMBLR_UTILS_REF=${GAMBLR_UTILS_REF} \
    GAMBLR_VIZ_REF=${GAMBLR_VIZ_REF} \
    GAMBLR_PREDICT_REF=${GAMBLR_PREDICT_REF} \
    GAMBLR_OPEN_REF=${GAMBLR_OPEN_REF}

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/install_gamblr.R /tmp/install_gamblr.R
RUN Rscript /tmp/install_gamblr.R

CMD ["R"]
