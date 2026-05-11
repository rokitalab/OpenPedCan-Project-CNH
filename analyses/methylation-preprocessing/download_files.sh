#!/usr/bin/env bash

set -euo pipefail

mkdir -p input-test

filelist="idat_filesnames.txt"

# ---- FIRST BUCKET ----
aws s3 ls s3://bti-private-us-east-1-prd-rood-lab/ --recursive --profile cnh-sso > bucket_listing.txt
awk '{print $4}' bucket_listing.txt > bucket_files.txt

while IFS= read -r fname || [ -n "$fname" ]; do
    fname=$(echo "$fname" | tr -d '\r')

    echo "Searching (rood): $fname"

    match=$(grep -E "/$fname\$" bucket_files.txt | head -n 1 || true)

    if [ -n "$match" ]; then
        echo "Found: $match"
        aws s3 cp "s3://bti-private-us-east-1-prd-rood-lab/$match" input-test/ --profile cnh-sso
    else
        echo "Missing (rood): $fname" >> missing_files.txt
    fi
done < "$filelist"

rm bucket_listing.txt bucket_files.txt

# ---- SECOND BUCKET ----
aws s3 ls s3://bti-private-us-east-1-prd-fonseca-lab/ --recursive --profile cnh-sso > bucket_listing.txt
awk '{print $4}' bucket_listing.txt > bucket_files.txt

while IFS= read -r fname || [ -n "$fname" ]; do
    fname=$(echo "$fname" | tr -d '\r')

    echo "Searching (fonseca): $fname"

    match=$(grep -E "/$fname\$" bucket_files.txt | head -n 1 || true)

    if [ -n "$match" ]; then
        echo "Found: $match"
        aws s3 cp "s3://bti-private-us-east-1-prd-fonseca-lab/$match" input-test/ --profile cnh-sso
    else
        echo "Missing (fonseca): $fname" >> missing_files_final.txt
    fi
done < "$filelist"