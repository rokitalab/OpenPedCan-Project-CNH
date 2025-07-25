cwltool  --user-space-docker-cmd="docker" --outdir outputs/current_Rscript workflow/methylation-preprocessing.cwl workflow/methylation-preprocessing-input.yml

 docker run -it --rm -v $PWD:$PWD dmiller15/minfi:4.2.0 bash