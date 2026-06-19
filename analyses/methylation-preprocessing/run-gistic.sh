#!/usr/bin/env bash

mkdir -p test-out/test-gistic-EPIC

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/mcr/v83/runtime/glnxa64:/opt/mcr/v83/bin/glnxa64:/opt/mcr/v83/sys/os/glnxa64
export XAPPLRESDIR=/opt/mcr/v83/X11/app-defaults


/home/rstudio/gistic_install/share/gistic2-2.0.23-0/gp_gistic2_from_seg \
  -b test-out/test-gistic-EPIC \
  -seg test-out/test-IlluminaHumanMethylationEPIC-gistic.seg \
  -refgene ~/gistic_install/share/gistic2-2.0.23-0/refgenefiles/hg19.mat \
  -genegistic 1 \
  -smallmem 1 \
  -broad 1 \
  -brlen 0.7 \
  -conf 0.99 \
  -armpeel 1 \
  -ta 0.2 \
  -td 0.2

mkdir -p test-out/test-gistic-EPICv2

/home/rstudio/gistic_install/share/gistic2-2.0.23-0/gp_gistic2_from_seg \
  -b test-out/test-gistic-EPICv2 \
  -seg test-out/test-IlluminaHumanMethylationEPICv2-gistic.seg \
  -refgene ~/gistic_install/share/gistic2-2.0.23-0/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat \
  -genegistic 1 \
  -smallmem 1 \
  -broad 1 \
  -brlen 0.7 \
  -conf 0.99 \
  -armpeel 1 \
  -ta 0.2 \
  -td 0.2
  
mkdir -p test-out/test-gistic-450k

/home/rstudio/gistic_install/share/gistic2-2.0.23-0/gp_gistic2_from_seg \
  -b test-out/test-gistic-450k \
  -seg test-out/test-IlluminaHumanMethylation450k-gistic.seg \
  -refgene ~/gistic_install/share/gistic2-2.0.23-0/refgenefiles/hg19.mat \
  -genegistic 1 \
  -smallmem 1 \
  -broad 1 \
  -brlen 0.7 \
  -conf 0.99 \
  -armpeel 1 \
  -ta 0.2 \
  -td 0.2
