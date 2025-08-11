# Base docker with miniconda
FROM docker.io/continuumio/miniconda3:latest

# Connect this docker image to the biobb_workflows repository
LABEL org.opencontainers.image.source https://github.com/mmb-irb/MDDB-workflow

# Define working dir
WORKDIR /app

# Install libgl1 for avoiding problems with VMD
RUN apt-get update && \
    apt-get install -y libgl1

# Clone repo
RUN git clone https://github.com/mmb-irb/MDDB-workflow.git

# Enable libmamba as solver
RUN conda config --set solver libmamba

# Create new environment
RUN conda env create --file=/app/MDDB-workflow/envs/environment.yml && conda clean -afy

# Install mwf SW
WORKDIR /app/MDDB-workflow
RUN conda run -n mwf_env python -m pip install -e .

# Define working dir
WORKDIR /app

# Export python path
ENV PYTHONPATH="${PYTHONPATH}:/app/MDDB-workflow"

# Code to run when container is started
RUN echo '#!/bin/sh' > entrypoint.sh && \
    echo 'conda run --no-capture-output -n mwf_env "$@"' >> entrypoint.sh && \
    chmod +x entrypoint.sh

# Serve the app
ENTRYPOINT ["/app/entrypoint.sh"]
