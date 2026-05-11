#!/usr/bin/env bash


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/mcr/v83/runtime/glnxa64:/opt/mcr/v83/bin/glnxa64:/opt/mcr/v83/sys/os/glnxa64
export XAPPLRESDIR=/opt/mcr/v83/X11/app-defaults


/home/rstudio/gistic_install/share/gistic2-2.0.23-0/gp_gistic2_from_seg \
  -b test-out/test-gistic \
  -seg test-out/-IlluminaHumanMethylationEPIC-gistic.seg \
  -refgene ~/gistic_install/share/gistic2-2.0.23-0/refgenefiles/hg19.mat \
  -genegistic 1 \
  -smallmem 1 \
  -broad 1 \
  -brlen 0.7 \
  -conf 0.99 \
  -armpeel 1 \
  -ta 0.2 \
  -td 0.2
