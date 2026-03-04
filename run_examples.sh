#!/bin/bash
set -e

# Make output directory
mkdir -p example_out

# Run Conda-based examples
/usr/bin/time --verbose pseudoDB \
	-sp human_chr22 \
	-fa ./example_data/human_chr22.fa \
	-s ./example_data/human_chr22_sample_list.txt \
	-o ./example_out/human_chr22_case1 \
	-t 16 -sl > ./example_out/human_chr22_case1.log 2>&1

/usr/bin/time --verbose pseudoDB \
	-sp human_chr22 \
	-fa ./example_data/human_chr22.fa \
	-s ./example_data/human_chr22_sample_list.txt \
	-o ./example_out/human_chr22_case2 \
	-db ./example_data/human_chr22_dbSNP.vcf.gz \
	-dn dbSNP \
	-t 16 -sl > ./example_out/human_chr22_case2.log 2>&1

/usr/bin/time --verbose pseudoDB \
	-sp human_chr22 \
	-fa ./example_data/human_chr22.fa \
	-s ./example_data/human_chr22_sample_list.txt \
	-o ./example_out/human_chr22_case3 \
	-db ./example_data/human_chr22_pseudoDB.vcf.gz \
	-dn pseudoDB \
	-t 16 -sl > ./example_out/human_chr22_case3.log 2>&1
