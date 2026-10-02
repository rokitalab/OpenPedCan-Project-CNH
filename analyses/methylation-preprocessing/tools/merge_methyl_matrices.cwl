class: CommandLineTool
cwlVersion: v1.2
id: merge_methyl_matrices
doc: |-
  Intersect probe sets and merge per-array methylation matrices.

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: pgc-images.sbgenomics.com/sicklera/openpedcanverse:latest
- class: InitialWorkDirRequirement
  listing:
  - entryname: 02-merge-methyl-matrices.R
    writable: false
    entry:
      $include: ../scripts/02-merge-methyl-matrices.R
  - $(inputs.methylation_files)
baseCommand: [Rscript, --vanilla, 02-merge-methyl-matrices.R]
arguments:
- prefix: --output_dir
  valueFrom: .
inputs:
  run_merge:
    type: boolean
    default: true
    doc: "Workflow control input; the merge step runs only when this is true."
  methylation_files:
    type: 'File[]'
    doc: "Per-array parquet outputs from preprocess_illumina_arrays, including detection p-values."
  output_basename:
    type: string
    inputBinding:
      prefix: --output_prefix
outputs:
  merged_values:
    type: 'File[]'
    outputBinding:
      glob:
      - '*-IlluminaHumanMethylationEPICv1-EPICv2-methyl-*.parquet'
      - '*-IlluminaHumanMethylationEPICv1-EPICv2-450k-methyl-*.parquet'
      - '*-IlluminaHumanMethylationEPICv1-450k-methyl-*.parquet'
      - '*-IlluminaHumanMethylationEPICv2-450k-methyl-*.parquet'
