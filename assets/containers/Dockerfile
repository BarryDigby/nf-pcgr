FROM nfcore/base:1.14

LABEL \
    author="Barry Digby" \
    description="Custom container to reformat VCF files for PCGR" \
    maintainer="bdigby@fredhutch.org"

# unzip
# install main packages:
RUN apt-get update; apt-get clean all;
RUN apt-get install --yes unzip

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add mamba installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/pysam/bin/:$PATH

# Download xcmds
RUN wget --no-check-certificate --tries=0 https://github.com/gudeqing/xcmds/archive/refs/heads/master.zip && \
    unzip master.zip && \
    pip install xcmds-master/dist/*

## N.B: must use 'python3.6' on cmdline not 'python'
