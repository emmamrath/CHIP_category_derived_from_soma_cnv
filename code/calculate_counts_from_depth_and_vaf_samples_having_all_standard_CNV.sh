#!/bin/bash

cohort=$1 # ISKS
indir=$2 # ../results
outdir=$3 # ../results

outfile="${outdir}"/"${cohort}".all_standard_CNV.tsv

# chrom	window_id	start_index	end_index	start_pos	end_pos	fit.k1	fit.k2	fit.f	fit.isize	fit.llik
# 1	1:1	1	54658	943687	249136360	1	1	0	0.00122848975433973	-321059.55105003
# 2	2:547	54659	116046	391781	181575281	1	1	0	0.00122848975433973	-360220.252692134
# 3	3:1161	116047	175967	154590	197808975	1	1	0	0.00122848975433973	-352597.668870773
# 4	4:1760	175968	232834	367927	190790246	1	1	0	0.00122848975433973	-335033.386352555

:>"${outfile}"

for infile in "${indir}"/*.nb2_position.tsv; do

	counts=$(grep -v '^chrom' "${infile}" | cut -d$'\t' -f7-8 | grep -P -v '^1\t1$' | wc -l)
	if [ "$counts" -eq "0" ]; then
		echo $infile $counts
		infile_basename=$(basename $infile)
		IFS='.' read -r -a array <<< "$infile_basename"
		sampleid="${array[0]}"
		outline="${sampleid}\tpass\tall_standard_CNV\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t"
		echo -e "${outline}" >> "${outfile}"
	fi
done

