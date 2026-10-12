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



#Install conumee

RUN R -e "devtools::install_github( \
  'hovestadtlab/conumee2', \
  subdir='conumee2', \
  ref='0a21b9d508506d14ccbc3c8820d5974b26ff3824', \
  dependencies=TRUE)"


# Install duckplyr and nanoparquet for fast lookups against large lesion tables
# (SNV, fusion) and parquet files in the molecular subtyping modules
RUN R -e "install.packages('remotes', repos = 'https://cloud.r-project.org'); \
  remotes::install_version('nanoparquet', version = '0.4.3', repos = 'https://cloud.r-project.org'); \
  remotes::install_version('duckdb', version = '1.5.0', repos = 'https://cloud.r-project.org'); \
  remotes::install_version('duckplyr', version = '1.2.1', repos = 'https://cloud.r-project.org')"
