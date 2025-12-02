class: CommandLineTool
cwlVersion: v1.2
id: preprocess_illumina_arrays
doc: |-
  Prepocess raw Illumina Infinium HumanMethylation BeadArrays (450K, and 850k)
  intensities using minfi into usable methylation measurements (Beta and M values)
  and copy number (cn-values) for OpenPedCan.

requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: pgc-images.sbgenomics.com/sicklera/openpedcanverse:latest
- class: ResourceRequirement
  ramMin: $(inputs.ram * 1000)
  coresMin: $(inputs.cores)
- class: InitialWorkDirRequirement
  listing:
  - entryname: 01-preprocess-illumina-arrays.R
    writable: false
    entry:
      $include: ../scripts/01-preprocess-illumina-arrays.R
baseCommand: [Rscript, --vanilla, 01-preprocess-illumina-arrays.R]
arguments:
- position: 99
  prefix: ''
  shellQuote: false
  valueFrom: |
    1>&2
inputs:
  input_idats_dir: { type: Directory, loadListing: shallow_listing, inputBinding: { prefix: "--base_dir", position: 1 }, doc: "Directory containing the IDATs to process." }
  output_basename: { type: 'string', doc: "Prefix string for output file name.", inputBinding: { position: 1, prefix: "--output_basename"} }
  manifest_file: {type: File, inputBinding: { prefix: "--manifest_file", position: 1 }, doc: "Manifest file containing 'file_name' and 'Bioassay_ID' columns"}
  funnorm: { type: 'boolean?', inputBinding: { prefix: "--funnorm", position: 1 }, doc: "If set, use funnorm for normalization" }
  snp_filter: { type: 'boolean?', inputBinding: { prefix: "--snp_filter", position: 1 }, doc: "If set, drops the probes that contain either a SNP at the CpG interrogation or at the single nucleotide extension." }
  ram: { type: 'int?', default: 32, doc: "GB of RAM to allocate to the task." }
  cores: { type: 'int?', default: 16, inputBinding: { prefix: "--n_cores", position: 1 }, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  beta_values:
    type: File
    outputBinding:
      glob: '*-methyl-beta-values-masked.qs2'
  m_values_unmasked:
    type: File
    outputBinding:
      glob: '*-methyl-m-values-unmasked.qs2'
  m_values_masked:
    type: File
    outputBinding:
      glob: '*-methyl-m-values-masked.qs2'
  cn_values:
    type: File
    outputBinding:
      glob: '*-methyl-cn-values.qs2'
  zero_mad:
    type: File
    outputBinding:
      glob: '*-zero-mad.txt'
