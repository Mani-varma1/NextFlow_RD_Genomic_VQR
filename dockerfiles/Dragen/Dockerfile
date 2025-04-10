FROM continuumio/miniconda3

# Set environment to non-interactive
ENV DEBIAN_FRONTEND=noninteractive

# Add necessary channels and install dragen-os
RUN conda update -n base -c defaults conda && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n dragen-env dragen-os && \
    conda clean -a -y

# Activate the environment by default
SHELL ["conda", "run", "-n", "dragen-env", "/bin/bash", "-c"]

# Set working directory
WORKDIR /data

# Default command
CMD ["dragen-os", "--help"]