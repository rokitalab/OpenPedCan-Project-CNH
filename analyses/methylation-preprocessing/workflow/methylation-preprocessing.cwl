cwlVersion: v1.2
class: Workflow
id: methylation-preprocessing
label: Methylation Preprocessing
doc: |-
  # Methylation Preprocessing
  
  Prepocess raw Illumina Infinium
  intensities using minfi into usable methylation measurements (Beta and M values)
  and copy number (cn-values) for OpenPedCan.

  This workflow can preprocess multiple array types at one time and will generate separate
  outputs for each array type in the input folder.

  ### Required Inputs
  input_idats_dir: { type: Directory, doc: "Directory containing the IDAT files to process." }
   
  ### Optional Inputs
  controls_present: { type: 'boolean?', doc: "If set, preprocesses the Illumina methylation array dataset assuming presence of either normal and tumor samples or samples of mutiple cancer groups or both." }
  snp_filter: { type: 'boolean?', doc: "If set, drops the probes that contain either a SNP at the CpG interrogation or at the single nucleotide extension." }
  ram: { type: 'int?', default: 32, doc: "GB of RAM to allocate to the task." }
  cores: { type: 'int?', default: 16, doc: "Minimum reserved number of CPU cores for the task." }

requirements:
- class: SubworkflowFeatureRequirement
- class: ScatterFeatureRequirement
inputs:
  input_idats_dir: { type: Directory, doc: "Directory containing the IDAT files to process." }
  manifest_file: {type: File, doc: "Manifest file containing 'file_name' and 'Bioassay_ID' columns"}
  controls_present: { type: 'boolean?', doc: "If set, preprocesses the Illumina methylation array dataset assuming presence of either normal and tumor samples or samples of mutiple cancer groups or both." }
  snp_filter: { type: 'boolean?', doc: "If set, drops the probes that contain either a SNP at the CpG interrogation or at the single nucleotide extension." }
  ram: { type: 'int?', default: 32, doc: "GB of RAM to allocate to the task." }
  cores: { type: 'int?', default: 16, doc: "Minimum reserved number of CPU cores for the task." }

outputs:
  beta_values: {type: 'File[]', outputSource: preprocess_illumina_arrays/beta_values }
  m_values: {type: 'File[]', outputSource: preprocess_illumina_arrays/m_values }
  cn_values: {type: 'File[]', outputSource: preprocess_illumina_arrays/cn_values }
  
steps:
  unzip_and_sort_files:
    run: ../tools/unzip_and_sort_files.cwl
    in:
      input_idats_dir: input_idats_dir
      ram: ram
      cores: cores
    out: [array_dirs]
  
  preprocess_illumina_arrays:
    run: ../tools/preprocess_illumina_arrays.cwl
    scatter: [input_idats_dir]
    in:
      input_idats_dir: unzip_and_sort_files/array_dirs
      manifest_file: manifest_file
      controls_present: controls_present
      snp_filter: snp_filter
      ram: ram
      cores: cores
    out: [beta_values, m_values, cn_values]
  

