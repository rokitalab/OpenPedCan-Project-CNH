#!/usr/bin/env bash

mkdir -p input-test

aws s3 ls s3://bti-private-us-east-1-prd-rood-lab/ --recursive > bucket_listing.txt --profile cnh-sso

awk '{print $4}' bucket_listing.txt > bucket_files.txt

fname=idat_filesnames.txt

while read -r fname; do
    match=$(grep -F "/$fname" bucket_files.txt || true)

    if [ -n "$match" ]; then
        aws s3 cp "s3://bti-private-us-east-1-prd-rood-lab/$match" input-test/ --profile cnh-sso
    else
        echo "$fname" >> missing_files.txt
    fi
done < "$fname"

rm bucket_listing.txt
rm bucket_files.txt

aws s3 ls s3://bti-private-us-east-1-prd-fonseca-lab/ --recursive > bucket_listing.txt --profile cnh-sso

awk '{print $4}' bucket_listing.txt > bucket_files.txt


while read -r fname; do
    match=$(grep -F "/$fname" bucket_files.txt || true)

    if [ -n "$match" ]; then
        aws s3 cp "s3://bti-private-us-east-1-prd-fonseca-lab/$match" input-test/ --profile cnh-sso
    else
        echo "$fname" >> missing_files_final.txt
    fi
done < "$fname"
