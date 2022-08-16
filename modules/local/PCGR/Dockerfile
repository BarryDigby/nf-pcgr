FROM nfcore/base:1.14

LABEL \
    author="Barry Digby" \
    description="PCGR image for nf-core pipelines" \
    maintainer="bdigby@fredhutch.org"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add mamba installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-mamba/bin/:$PATH

# Setup default ARG variables
ARG PCGR_VERSION="1.0.3"
ARG PCGR_REPO="https://raw.githubusercontent.com/sigven/pcgr/v${PCGR_VERSION}/conda/env/lock/"
ARG PLATFORM="linux"

# Follow install instructions from PCGR documentation
RUN mamba create --file ${PCGR_REPO}/pcgr-${PLATFORM}-64.lock --prefix ./pcgr
RUN mamba create --file ${PCGR_REPO}/pcgrr-${PLATFORM}-64.lock --prefix ./pcgrr

# Add PCGR installation dir to PATH (instead of doing 'mamba activate')
ENV PATH /opt/mambaforge/envs/pcgr/bin/:$PATH

# Data bundle vars
ARG GENOME="grch38"
ARG BUNDLE_VERSION="20220203"
ARG BUNDLE="pcgr.databundle.${GENOME}.${BUNDLE_VERSION}.tgz"

# Download and unpack to PCGR bin
RUN wget http://insilico.hpc.uio.no/pcgr/${BUNDLE}
RUN gzip -dc ${BUNDLE} | tar xvf -
RUN mv data/ /opt/mambaforge/envs/pcgr/bin

# Dump the details of the installed packages to a file for posterity
RUN mamba env export --name pcgr > nf-core-pcgr-${PCGR_VERSION}.yml
