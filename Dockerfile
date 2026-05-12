# Dockerfile (thin image for workflows / notebooks)
FROM pgc-images.sbgenomics.com/rokita-lab/openpedcanverse-base:2025-12-01

LABEL maintainer="jrokita@childrensnational.org" \
      org.opencontainers.image.title="openpedcanverse" \
      org.opencontainers.image.description="Project-specific image layered on openpedcanverse-base"

# Recommended: set working dir for your project code
WORKDIR /rocker-build/

# If you have a small set of extra R/Python deps *for this repo only*,
# install them here. Keep this light so GH Actions doesn’t OOM.
# Example:
# RUN R -q -e "install.packages(c('here','argparse'), repos='https://cloud.r-project.org')"
# RUN python3 -m pip install --no-cache-dir some-extra-package==1.2.3

# Install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libxml2-dev \
        libutf8proc-dev \
        cmake \
        pkg-config \
        libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*
    
    
#Install Arrow 
RUN R -e "install.packages(c('devtools', 'arrow'), repos='https://cloud.r-project.org')"

#Install conumee
RUN R -e "devtools::install_github('hovestadtlab/conumee2', subdir='conumee2', dependencies=TRUE)"
