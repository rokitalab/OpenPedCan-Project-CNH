# RAM optimization on methylation preprocessing script 01-preprocess-illumina-arrays.R

# pull example samples from S3
# bash mount_s3.sh bti-private-us-east-1-prd-fonseca-lab
# 10 samples, 1 sample (205566000176_R05C01) not in manifest but in previous PBTA cohort
# path1=data/bti-private-us-east-1-prd-fonseca-lab/rare-brain-tumor-program/source/MethylationData/20250624-Methylation/


# only copy "IlluminaHumanMethylationEPICv2 types in one manifest
grep EPICv2  data-modeling/generate-bioassay-ids/output/* |less -S |grep itt334|cut -f6|less > epicv2_idats.txt
while read file; do 
    find data/bti-private-us-east-1-prd-fonseca-lab/rare-brain-tumor-program/source/MethylationData/ -name "$file" -type f -exec cp {} data/test_epicv2_idats/ \;; 
done < /home/ubuntu/epicv2_idats.txt


Rscript --vanilla scripts/01-preprocess-illumina-arrays.R --base_dir=data/test1_epicv2_idats/ \
  --manifest_file=itt334-rbt-methylation_IDs_assigned.tsv