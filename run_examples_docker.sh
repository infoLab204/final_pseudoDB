#!/bin/bash
set -e

# Make output directory
mkdir -p example_out

# Run Docker-based examples
/usr/bin/time --verbose docker run \
  --user $(id -u):$(id -g) \
  -v ./example_data:/data \
  -v ./example_out:/output \
  pseudodb \
  -sp human_chr22 \
  -fa /data/chr22.fa \
  -s /data/docker_sample_list_chr22.txt \
  -o /output/docker_human_chr22_case1 \
  -t 16 -sl > ./example_out/docker_human_chr22_case1.log 2>&1

/usr/bin/time --verbose docker run \
  --user $(id -u):$(id -g) \
  -v ./example_data:/data \
  -v ./example_out:/output \
  pseudodb \
  -sp human_chr22 \
  -fa /data/chr22.fa \
  -s /data/docker_sample_list_chr22.txt \
  -o /output/docker_human_chr22_case2 \
  -db /data/chr22_dbSNP.vcf.gz \
  -dn dbSNP \
  -t 16 > ./example_out/docker_human_chr22_case2.log 2>&1

/usr/bin/time --verbose docker run \
  --user $(id -u):$(id -g) \
  -v ./example_data:/data \
  -v ./example_out:/output \
  pseudodb \
  -sp human_chr22 \
  -fa /data/chr22.fa \
  -s /data/docker_sample_list_chr22.txt \
  -o /output/docker_human_chr22_case3 \
  -db /data/chr22_pseudoDB.vcf.gz \
  -dn pseudoDB \
  -t 16 > ./example_out/docker_human_chr22_case3.log 2>&1
