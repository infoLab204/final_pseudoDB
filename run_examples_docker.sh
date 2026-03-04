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
	-fa /data/human_chr22.fa \
	-s /data/human_chr22_sample_list_docker.txt \
	-o /output/human_chr22_case1_docker \
	-t 16 -sl > ./example_out/human_chr22_case1_docker.log 2>&1

/usr/bin/time --verbose docker run \
	--user $(id -u):$(id -g) \
	-v ./example_data:/data \
	-v ./example_out:/output \
	pseudodb \
	-sp human_chr22 \
	-fa /data/human_chr22.fa \
	-s /data/human_chr22_sample_list_docker.txt \
	-o /output/human_chr22_case2_docker \
	-db /data/human_chr22_dbSNP.vcf.gz \
	-dn dbSNP \
	-t 16 > ./example_out/human_chr22_case2_docker.log 2>&1

/usr/bin/time --verbose docker run \
	--user $(id -u):$(id -g) \
	-v ./example_data:/data \
	-v ./example_out:/output \
	pseudodb \
	-sp human_chr22 \
	-fa /data/human_chr22.fa \
	-s /data/human_chr22_sample_list_docker.txt \
	-o /output/human_chr22_case3_docker \
	-db /data/human_chr22_pseudoDB.vcf.gz \
	-dn pseudoDB \
	-t 16 > ./example_out/human_chr22_case3_docker.log 2>&1
