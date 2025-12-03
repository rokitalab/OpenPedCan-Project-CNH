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


