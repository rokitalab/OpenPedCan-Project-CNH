class: CommandLineTool
cwlVersion: v1.2
id: unzip_and_sort_files
doc: |-
  Unzip and sort idats files into separate directories based on methylation array type.

requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: pgc-images.sbgenomics.com/rokita-lab/openpedcanverse:latest
- class: ResourceRequirement
  ramMin: $(inputs.ram * 1000)
  coresMin: $(inputs.cores)
- class: InitialWorkDirRequirement
  listing:
  - entryname: 00-unzip-and-sort.R
    writable: false
    entry:
      $include: ../scripts/00-unzip-and-sort.R
baseCommand: [Rscript, --vanilla, 00-unzip-and-sort.R]
arguments:
- position: 99
  prefix: ''
  shellQuote: false
  valueFrom: |
    1>&2
inputs:
  input_idats_dir: { type: Directory, loadListing: shallow_listing, inputBinding: { prefix: "--base_dir", position: 1 }, doc: "Directory containing the IDATs to process." }
  ram: { type: 'int?', default: 32, doc: "GB of RAM to allocate to the task." }
  cores: { type: 'int?', default: 16, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  array_dirs:
    type: 'Directory[]'
    outputBinding:
      glob: 'output_dir/*'
