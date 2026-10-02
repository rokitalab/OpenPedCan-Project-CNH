class: CommandLineTool
cwlVersion: v1.2
id: liftover_methylation_segments
label: Liftover Methylation CNV Segments
doc: |-
  Lift legacy methylation-array GISTIC segments from hg19 to hg38 and combine
  them with EPICv2 segments, which already use hg38 coordinates.

requirements:
- class: InlineJavascriptRequirement
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

outputs:
  combined_hg38_segments:
    type: File
    outputBinding:
      glob: combined_hg38.gistic.seg
