class: CommandLineTool
cwlVersion: v1.2
id: methylation_cnv_calls
label: Methylation CNV Calling
doc: |-
  Call copy-number segments from a QC-filtered methylation RGChannelSet using
  Noob normalization and conumee2.

requirements:
- class: DockerRequirement
  dockerPull: pgc-images.sbgenomics.com/sicklera/openpedcanverse:latest
- class: InitialWorkDirRequirement
  listing:
  - entryname: 01-methylation-cnv-calls.R
    writable: false
    entry:
      $include: ../scripts/01-methylation-cnv-calls.R

baseCommand: [Rscript, --vanilla, 01-methylation-cnv-calls.R]

inputs:
  rg_set_file:
    type: File
    doc: "QC-filtered RGChannelSet produced by methylation preprocessing (.qs2)."
    inputBinding:
      prefix: --rg_set_file
  manifest_file:
    type: File
    doc: "Manifest containing file_name, Bioassay_ID, sample_type, and platform columns."
    inputBinding:
      prefix: --manifest_file
  output_basename:
    type: string
    doc: "Prefix for the output segment files."
    inputBinding:
      prefix: --output_basename
  array_type:
    type: string
    default: EPIC
    doc: "Array type: EPIC, EPICv2, or 450k."
    inputBinding:
      prefix: --array_type

outputs:
  segments:
    type: File
    outputBinding:
      glob: '*-segments.seg'
  gistic_segments:
    type: File
    outputBinding:
      glob: '*-gistic.seg'
