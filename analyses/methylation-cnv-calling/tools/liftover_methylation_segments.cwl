class: CommandLineTool
cwlVersion: v1.2
id: liftover_methylation_segments
label: Liftover Methylation CNV Segments
doc: |-
  Lift legacy methylation-array GISTIC segments from hg19 to hg38 and combine
  them with EPICv2 segments, which already use hg38 coordinates.

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: pgc-images.sbgenomics.com/sicklera/openpedcanverse:latest
- class: ResourceRequirement
  ramMin: $(inputs.ram * 1000)
  coresMin: $(inputs.cores)
- class: InitialWorkDirRequirement
  listing:
  - entryname: 02-liftover-methylation-segments.R
    writable: false
    entry:
      $include: ../scripts/02-liftover-methylation-segments.R
  - $(inputs.gistic_segments)

baseCommand: [Rscript, --vanilla, 02-liftover-methylation-segments.R]
arguments:
- prefix: --seg_dir
  valueFrom: .

inputs:
  gistic_segments:
    type: File[]
    doc: "Per-array GISTIC SEG files produced by methylation CNV calling."
  chain_file:
    type: File?
    doc: "hg19-to-hg38 chain file; required for EPICv1 or 450k input."
    inputBinding:
      prefix: --chain_file

  ram:
    type: int
    default: 8
    doc: "Minimum RAM in GB."
  cores:
    type: int
    default: 1
    doc: "Minimum CPU cores."

outputs:
  combined_hg38_segments:
    type: File
    outputBinding:
      glob: combined_hg38.gistic.seg
