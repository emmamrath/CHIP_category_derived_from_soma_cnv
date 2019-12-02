#!/bin/bash

cohort=$1 # ISKS
infile=$2 # ../ISKS_cannot_model.txt
outdir=$3 # ../results

# Input file contains list of samples whose modelling terminated in error with:
# Negative binomial optimization failed to converge

outfile="${outdir}"/"${cohort}".cannot_model.tsv

:>"${outfile}"

while IFS= read -r sampleid; do

	outline="${sampleid}\tfail\tcannot_model\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t"
	echo -e "${outline}" >> "${outfile}"

done < "$infile"

