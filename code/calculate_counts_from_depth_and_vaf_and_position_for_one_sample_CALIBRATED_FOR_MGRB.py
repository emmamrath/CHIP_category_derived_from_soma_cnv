#!/usr/bin/python
# python calculate_counts_from_depth_and_vaf_and_position_for_one_sample.py -i indata -i2 indata2 -i3 indata3 -s sampleid -o outdata
# python calculate_counts_from_depth_and_vaf_and_position_for_one_sample.py -i AACMY_LP1000341-NTP_E08_ST-E00152_HYFWFCCXX_6_nb2_counts.tsv -s AACMY -o AACMY_nb2_counts_tallies.tsv


# Example indata:
# "start_chrom"	"end_chrom"	"start_index"	"end_index"	"depth"	"vaf1"	"vaf2"	"measure_of_subclonality_from_model_fitting"
# 1	1	1	54869	0	0.5	0.5	0.0370067675252907
# 2	2	54870	117922	0	0.5	0.5	0.0370067675252907
# 3	3	117923	174747	0	0.5	0.5	0.0370067675252907
# 3	4	174748	225399	0	0.5	0.5	0.0370067675252907
# 4	4	225400	228636	-0.0269448009564896	0.509426106751942	0.490573893248059	0.0370067675252907
# 4	5	228637	278309	0	0.5	0.5	0.0370067675252907
# 5	6	278310	326909	0	0.5	0.5	0.0370067675252907
# 6	7	326910	327099	-0.0544024354577584	0.5	0.5	0.0370067675252907
# 7	7	327100	343699	0	0.5	0.5	0.0370067675252907

# Example indata2:
# chrom	window_id	start_index	end_index	start_pos	end_pos	fit.k1	fit.k2	fit.f	fit.isize	fit.llik
# 1	1:1	1	54869	943687	249126816	1	1	0.0370067675252907	0.001	-305705.953599272
# 2	2:549	54870	117922	107437	242852778	1	1	0.0370067675252907	0.001	-351156.120207228
# 3	3:1180	117923	174747	154590	197808975	1	1	0.0370067675252907	0.001	-316638.472704259
# 4	4:1748	174748	225399	367927	181911457	1	1	0.0370067675252907	0.001	-282278.460617863
# 4	4:2255	225400	228636	181912396	190789536	1	0	0.0370067675252907	0.001	-17869.0188403672
# 5	5:2287	228637	278309	203042	180643732	1	1	0.0370067675252907	0.001	-277223.533536684
# 6	6:2784	278310	326909	230572	170887693	1	1	0.0370067675252907	0.001	-271014.236880734
# 7	7:3270	326910	327099	95621	2706478	0	0	0.0370067675252907	0.001	-1080.67610533181
# 7	7:3272	327100	343699	2713683	54950551	1	1	0.0370067675252907	0.001	-92292.534445401

# Example indata3:
# 4	181912394	190789537	225400	228636	ACSL1
# 4	181912394	190789537	225400	228636	ANKRD37
# 4	181912394	190789537	225400	228636	C4orf47
# 4	181912394	190789537	225400	228636	CASP3
# 4	181912394	190789537	225400	228636	CCDC110
# 4	181912394	190789537	225400	228636	CCDC111
# 4	181912394	190789537	225400	228636	CDKN2AIP
# 4	181912394	190789537	225400	228636	CLDN22
# 4	181912394	190789537	225400	228636	CLDN24
# 4	181912394	190789537	225400	228636	CYP4V2


# outdata contains one line per sample with the counts and tallies found in that sample.


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2018, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import datetime
import math
import random
import commands
import argparse
import re
import numpy
# import subprocess
# from multiprocessing import Pool

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def is_float(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

######################################################
def make_list_unique(this_list):

	new_list = ''

	if (this_list != ''):
		list_fields = this_list.split('|')
		dictionary_of_list = {}
		for this_gene in list_fields:
			dictionary_of_list[this_gene] = this_gene
		for this_gene in dictionary_of_list:
			if (new_list == ''):
				new_list = this_gene
			else:
				new_list = new_list + '|' + this_gene

	return new_list

######################################################
def index_for_chrom(this_chrom):

	this_index = 0
	if (this_chrom == 'X'):
		this_index = 23
	elif (this_chrom == 'Y'):
		this_index = 24
	elif (this_chrom == 'MT'):
		this_index = 25
	elif (this_chrom == 'M'):
		this_index = 25
	elif (is_integer(this_chrom)):
		this_index = int(this_chrom)
	return this_index

######################################################
def does_it_have_diff_chrom_with_same_vaf_values( this_start_chrom, this_end_chrom, this_vaf1, this_vaf2, non_0p5_vaf_start_chrom, non_0p5_vaf_end_chrom, non_0p5_vaf1, non_0p5_vaf2 ):

	has_diff_chrom_with_same_vaf_values = False

	keep_looking = True
	i = 0
	while (keep_looking):
		if ((this_vaf1 == non_0p5_vaf1[i]) and (this_vaf2 == non_0p5_vaf2[i])):
			# these two events have the same vaf, are they on the same chromosome?
			same_chrom = True
			if (this_start_chrom != this_end_chrom):
				same_chrom = False
			elif (non_0p5_vaf_start_chrom[i] != non_0p5_vaf_end_chrom[i]):
				same_chrom = False
			elif ((this_start_chrom != non_0p5_vaf_start_chrom[i]) or (this_end_chrom != non_0p5_vaf_end_chrom[i])):
				same_chrom = False
			if (same_chrom == False):
				has_diff_chrom_with_same_vaf_values = True
				keep_looking = False
		i = i + 1
		if (i >= len(non_0p5_vaf_start_chrom)):
			keep_looking = False

	return has_diff_chrom_with_same_vaf_values

######################################################
def find_closest_telomere_to_this_position( chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides ):

	# We already know that this is a small telomeric event,
	# and the start_chrom is not the same as the end_chrom.
	# If start_chrom is 2 and end_chrom is 3, then we want to know if this is an event for 2q or 3p.
	# If it is a 2q telomere event, then there is only a little bit of chr3 in the event.
	# If it is a 3q telomere event, then there is only a little bit of chr2 in the event.
	# Does this event straddle more of start_chrom q
	# or does it straddle more of end_chrom p.
	# We are only interested in autosomal chromosomes 1 to 22.

	closest_telomere = ''

	telomere_length_in_nucleotides = 300000

	chrom_idx_for_ref_table = 0
	if (is_integer(chrom)):
		chrom_idx_for_ref_table = int(chrom)
	elif (chrom == 'X'):
		chrom_idx_for_ref_table = 23
	elif (chrom == 'Y'):
		chrom_idx_for_ref_table = 24

	start_window_in_nucleotides = telomere_length_in_nucleotides
	end_window_in_nucleotides = chromosome_size_in_nucleotides[chrom_idx_for_ref_table] - telomere_length_in_nucleotides

	if (start_pos <= start_window_in_nucleotides):
		closest_telomere = str(start_chrom) + 'p'
	elif (end_pos >= end_window_in_nucleotides):
		closest_telomere = str(end_chrom) + 'q'

	return closest_telomere

######################################################
def is_this_event_a_focal_telomere( i, chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides, chrom_array ):

	# If this event is within a stone's throw of the ends, and if it is small, then consider it to be a telomere event.

	telomere_position = ''

	close_enough_distance_to_chromosome_end_in_nucleotides = 6000000

	chrom_idx_for_ref_table = index_for_chrom(chrom)

	start_window_in_nucleotides = close_enough_distance_to_chromosome_end_in_nucleotides
	end_window_in_nucleotides = chromosome_size_in_nucleotides[chrom_idx_for_ref_table] - close_enough_distance_to_chromosome_end_in_nucleotides

	if (end_pos <= start_window_in_nucleotides):
		telomere_position = str(start_chrom) + 'p'
	elif (start_pos >= end_window_in_nucleotides):
		telomere_position = str(end_chrom) + 'q'

	return telomere_position

######################################################
def does_this_event_include_telomeres( i, chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides, chrom_array ):

	# This event may or may not be a focal telomere event.
	# Regardless, does it encompass any telomeres?
	# An entire chromosome event will encompass both telomeres.

	telomere_position_p = ''
	telomere_position_q = ''

	if (i == 0):
		telomere_position_p = str(chrom) + 'p'
	elif (chrom != chrom_array[i-1]):
		telomere_position_p = str(chrom) + 'p'
	if (i == (len(chrom_array) - 1)):
		telomere_position_q = str(chrom) + 'q'
	elif (chrom != chrom_array[i+1]):
		telomere_position_q = str(chrom) + 'q'

	return telomere_position_p, telomere_position_q

######################################################
def does_event_span_a_chromosome( i, chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides, chrom_array, start_chrom_array, end_chrom_array ):

	# We know that vaf1 != 0.5 so it is an event.
	# If this event spans from within 10000 bp of one end of the chromosome to within 10000 bp of the other end, 
	# then we consider it to cover the whole chromosome.

	event_spans_this_chromosome = ''

	how_to_determine_entire_chromosome_events = "whether_it_is_first_and_last_event_on_chromosome" # values are "whether_it_is_first_and_last_event_on_chromosome", "by_distance_from_end_of_chromosomes"

	if (how_to_determine_entire_chromosome_events == "by_distance_from_end_of_chromosomes"):

		close_enough_distance_to_chromosome_end_in_nucleotides = 600000

		chrom_idx_for_ref_table = index_for_chrom(chrom)

		start_window_in_nucleotides = close_enough_distance_to_chromosome_end_in_nucleotides
		end_window_in_nucleotides = chromosome_size_in_nucleotides[chrom_idx_for_ref_table] - close_enough_distance_to_chromosome_end_in_nucleotides

		if ((start_pos <= start_window_in_nucleotides) and (end_pos >= end_window_in_nucleotides)):
			event_spans_this_chromosome = chrom

	# If this event is the only event for this chromosome, then it must span the entire chromosome.
	# Execute this logic in both cases of how_to_determine_entire_chromosome_events

	if (i == 0):
		if (chrom_array[i+1] != chrom):
			event_spans_this_chromosome = chrom
	elif (i == (len(start_chrom_array) - 1)):
		if (chrom != chrom_array[i-1]):
			event_spans_this_chromosome = chrom
	elif ((chrom != chrom_array[i-1]) and (chrom != chrom_array[i+1])):
		event_spans_this_chromosome = end_chrom

	return event_spans_this_chromosome

######################################################
def assign_category_of_this_event( i, chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, depth, vaf1, vaf2, fit_f, vaf_dist_from_0p5, this_size_in_bins, fit_k1, fit_k2, fit_f_file2, this_size_in_nucleotides, num_bins_for_this_sample_FLOAT, num_bins_for_this_chromosome_FLOAT, chrom_array, start_chrom_array, end_chrom_array, chromosome_size_in_nucleotides ):

	# no_subclones, focal_subclones, obvious_subclones, telomere_subclones, germline_dup, failed_modelling

	# Here is the algorithm used to assign categories for each event:
	#
	# obvious_subclones : (vaf_dist_from_0p5 >= 0.016) and (this_size_in_bins >= 3000)
	#
	# germline_dup : (depth between 0.2 and 0.6) and (vaf between 0.06 and 0.2) and (size of event is smaller than half of the chromosome)
	#
	# telomere_subclones : ((depth < 0) or (depth > 0)) and (vaf != 0.5) and (event is of a certain size, between 0.005 to 0.3% the size of the genome)
	#
	# focal_subclones : (depth != 0) or (vaf != 0.5)

	# this_event_subcategory will be blank or: negative_depth_and_non_0p5_vaf, positive_depth_and_non_0p5_vaf, zero_depth_and_non_0p5_vaf

	# this_telomere_position will be blank or: 1p, 1q, 2p, 2q...

	# this_chromosome_event_for_chrom and this_chromosome_event_for_event will be filled in only if vaf != 0.5 and the event spans an entire chromosome.
	# this_chromosome_event_for_chrom will be blank or: '1', '2', '3'... 'X', 'Y'
	# this_chromosome_event_for_event will be blank or: "negative_depth_and_non_0p5_vaf", "positive_depth_and_non_0p5_vaf", "zero_depth_and_non_0p5_vaf"

	# Copy number (CN) state is derived directly from fit.k1 and fit.k2:
	# fit.k1 + fit.k2 < 2: CN_loss
	# fit.k1 + fit.k2 > 2: CN_gain
	# fit.k1 + fit.k2 = 2 AND fit.k1*fit.k2 = 0: CN_LoH
	# fit.k1 + fit.k2 = 2 AND fit.k1*fit.k2 = 1: No_CN 

	this_event_category = "no_subclones"
	this_event_CN_category = ''
	this_event_subcategory = ''
	this_telomere_position = ''
	this_telomere_position_p = ''
	this_telomere_position_q = ''
	this_chromosome_event_for_chrom = ''
	this_chromosome_event_for_event = ''
	got_category = False

	event_size_as_a_percentage_of_genome = float(end_index - start_index + 1) / float(num_bins_for_this_sample_FLOAT) * 100
	event_size_as_a_percentage_of_chromosome = float(end_index - start_index + 1) / float(num_bins_for_this_chromosome_FLOAT) * 100

	if (got_category == False):
		if ((vaf_dist_from_0p5 >= 0.016) and (this_size_in_bins >= 3000)):
			this_event_category = "obvious_subclones"
			got_category = True

	if (got_category == False):
		if ((depth >= 0.2) and (depth <= 0.6) and (vaf_dist_from_0p5 >= 0.06) and (vaf_dist_from_0p5 <= 0.2) and (event_size_as_a_percentage_of_chromosome < 50)):
			this_event_category = "germline_dup"
			got_category = True

	how_to_determine_telomere_events = "use_file2" # values are "use_file1", "use_file2"
	if (how_to_determine_telomere_events == "use_file1"):
		if (got_category == False):
			if ((depth != 0) and (vaf_dist_from_0p5 != 0) and (event_size_as_a_percentage_of_genome >= 0.005) and (event_size_as_a_percentage_of_genome <= 0.3)):
				if ( (i == 0) or (i == (len(start_chrom_array) - 1)) ): 
					this_event_category = "telomere_subclones"
					got_category = True
					this_telomere_position = find_closest_telomere_to_this_position( chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides )
				elif ( (start_chrom != end_chrom) and (start_chrom == end_chrom_array[i-1]) ):
					this_event_category = "telomere_subclones"
					got_category = True
					this_telomere_position = find_closest_telomere_to_this_position( chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides )
				elif ( (start_chrom != end_chrom) and (end_chrom == start_chrom_array[i+1]) ):
					this_event_category = "telomere_subclones"
					got_category = True
					this_telomere_position = find_closest_telomere_to_this_position( chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides )
				elif (start_chrom != end_chrom_array[i-1]):
					this_event_category = "telomere_subclones"
					got_category = True
					this_telomere_position = str(start_chrom) + 'p'
				elif (end_chrom != start_chrom_array[i+1]):
					this_event_category = "telomere_subclones"
					got_category = True
					this_telomere_position = str(start_chrom) + 'q'
	else: # how_to_determine_telomere_events == "use_file2"
		if (got_category == False):
			if ((fit_k1 == 1) and (fit_k2 == 1)):
				do_nothing = 1 # this is not an event
			else:
				this_telomere_position = is_this_event_a_focal_telomere( i, chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides, chrom_array )
				if (this_telomere_position != ''):
					this_event_category = "telomere_subclones"
					got_category = True

	if (got_category == False):
		if ((depth != 0) or (vaf1 != 0.5)):
			if (this_size_in_bins > 200):
				this_event_category = "focal_subclones"
				got_category = True

	# determine this_event_subcategory
	if (vaf1 != 0.5):
		if (depth < 0):
			this_event_subcategory = "negative_depth_and_non_0p5_vaf"
		elif (depth > 0):
			this_event_subcategory = "positive_depth_and_non_0p5_vaf"
		else: # (depth == 0)
			this_event_subcategory = "zero_depth_and_non_0p5_vaf"

	# determine this_chromosome_event values
	if (vaf1 != 0.5):
		event_spans_this_chromosome = does_event_span_a_chromosome( i, chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides, chrom_array, start_chrom_array, end_chrom_array )
		if (event_spans_this_chromosome != ''):
			this_chromosome_event_for_chrom = event_spans_this_chromosome
			if (depth < 0):
				this_chromosome_event_for_event = "CN_loss"
			elif (depth > 0):
				this_chromosome_event_for_event = "CN_gain"
			else: # (depth == 0)
				this_chromosome_event_for_event = "CN_LoH"

	# determine this_event_CN_category
	fit_k1_plus_fit_k2 = fit_k1 + fit_k2
	fit_k1_times_fit_k2 = fit_k1 * fit_k2
	if (fit_k1_plus_fit_k2 < 2):
		this_event_CN_category = 'CN_loss' # this corresponds to negative depth
	elif (fit_k1_plus_fit_k2 > 2):
		this_event_CN_category = 'CN_gain' # this corresponds to positive depth
	else: # fit_k1_plus_fit_k2 == 2
		if (fit_k1_times_fit_k2 == 0):
			this_event_CN_category = 'CN_LoH'
		elif (fit_k1_times_fit_k2 == 1):
			this_event_CN_category = 'No_CN'

	# When modelling passes (but not when it fails), this_event_subcategory and this_event_CN_category give the same information.
	# However, this_event_subcategory is more discerning than this_event_CN_category because depth is more discerning than k1 & k2
	# Apart from the difference in sensitivity,
	# CN_gain 	is 	positive_depth_and_non_0p5_vaf
	# CN_loss 	is 	negative_depth_and_non_0p5_vaf
	# CN_LoH 	is 	zero_depth_and_non_0p5_vaf

	# determine whethere events include any telomeres, regardless of event type
	if (vaf1 != 0.5):
		this_telomere_position_p, this_telomere_position_q = does_this_event_include_telomeres( i, chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, chromosome_size_in_nucleotides, chrom_array )


	debug = False
	if (debug):
		if ((vaf1 != 0.5) or (this_event_CN_category != '')):
			print 'subcategory', this_event_subcategory, 'CN_category', this_event_CN_category, 'k1', fit_k1, 'k2', fit_k2, 'k1+k2', fit_k1_plus_fit_k2, 'k1*k2', fit_k1_times_fit_k2, 'fit_f', fit_f, 'fit_f_file2', fit_f_file2

	debug = False
	if (debug):
		if ((depth != 0) and (vaf1 != 0.5)):
			print this_event_category, 'chrom', start_chrom, end_chrom, 'index', start_index, end_index, this_size_in_bins, 'depth', depth, 'vaf', vaf1, vaf2, vaf_dist_from_0p5, 'fit', fit_f, '%chrom', event_size_as_a_percentage_of_chromosome, '%genome', event_size_as_a_percentage_of_genome

	debug = False
	if (debug):
		print start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, this_event_category, this_event_CN_category, this_event_subcategory, this_telomere_position, this_chromosome_event_for_chrom, this_chromosome_event_for_event

	return this_event_category, this_event_CN_category, this_event_subcategory, this_telomere_position, this_telomere_position_p, this_telomere_position_q, this_chromosome_event_for_chrom, this_chromosome_event_for_event

######################################################
def convert_gene_list_to_array_having_standard_index( chrom_array2, start_pos_array2, end_pos_array2, start_index_array2, end_index_array2, chrom_array3, start_pos_array3, end_pos_array3, start_index_array3, end_index_array3, gene_array3 ):

	array_of_genes_hit_by_event = [''] * len(chrom_array2)
	i2 = 0
	this_chrom_array2 = chrom_array2[i2]
	#this_start_pos_array2 = start_pos_array2[i2]
	#this_end_pos_array2 = end_pos_array2[i2]
	this_start_index_array2 = start_index_array2[i2]
	this_end_index_array2 = end_index_array2[i2]

	for i3 in range( 0, len(chrom_array3) ):

		this_chrom_array3 = chrom_array3[i3]
		this_start_pos_array3 = start_pos_array3[i3]
		this_end_pos_array3 = end_pos_array3[i3]
		this_start_index_array3 = start_index_array3[i3]
		this_end_index_array3 = end_index_array3[i3]
		this_gene_array3 = gene_array3[i3]
		
		array2_upto_array3 = False
		while (array2_upto_array3 == False):
			if ((this_chrom_array2 == this_chrom_array3) and (this_start_index_array2 == this_start_index_array3) and (this_end_index_array2 == this_end_index_array3)):
				array2_upto_array3 = True
			else:
				i2 = i2 + 1
				this_chrom_array2 = chrom_array2[i2]
				#this_start_pos_array2 = start_pos_array2[i2]
				#this_end_pos_array2 = end_pos_array2[i2]
				this_start_index_array2 = start_index_array2[i2]
				this_end_index_array2 = end_index_array2[i2]

		if (array_of_genes_hit_by_event[i2] == ''):
			array_of_genes_hit_by_event[i2] = this_gene_array3
		else:
			array_of_genes_hit_by_event[i2] = array_of_genes_hit_by_event[i2] + '|' + this_gene_array3

	return array_of_genes_hit_by_event

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in CNV data derived from modelling from SNP data, and categorise the CNV events - looking for subclonal hematopoeisis')
	parser.add_argument('-i', action="store", dest="indata", required=True, help='Input file of copy number depth and vafs and overall fit derived from nb2 modelling of snps')
	parser.add_argument('-i2', action="store", dest="indata2", required=True, help='Input file 2 of copy number positions and vafs and fits derived from nb2 modelling of snps')
	parser.add_argument('-i3', action="store", dest="indata3", required=True, help='Input file 3 of genes hit in copy number events, one line per gene hit per non-zero CNV event')
	parser.add_argument('-s', action="store", dest="sampleid", required=True, help='Input sample_id')
	parser.add_argument('-o', action="store", dest="outdata", required=True, help='Output file of tallies from copy number depth and vafs derived from nb2 modelling of snps')
	args = parser.parse_args()

	# The model plotting data that is input to this program:
	# "start_chrom"	"end_chrom"	"start_index"	"end_index"	"depth"			"vaf1"			"vaf2"
	# 1		1		1		52895		0 			0.5			0.5
	# 2		2		52896		53299		0.0947089860362449      0.53176943644185	0.46823056355815
	# 2		2		53300		116825		0			0.5			0.5
	# 3		3		116826		117299		0.0947089860362449 	0.53176943644185	0.46823056355815
	# 3		3		117300		173105		0			0.5			0.5
	# 4		4		173106		224799		0			0.5			0.5

	# The model data from which the model plotting data is derived:
	# "chrom" "pos"   "dp"    "ad"    "affinity"      "gc100" "gc200" "gc400" "gc600" "gc800" "pois.lambda"   "prRR"  "prAA"  "het"
	# 1       1031540 31      12      0.95509336123576        0.33    0.415   0.4475  0.436666666667  0.41625 33.8612179648534        1.20606404314962e-14    2.37405127174839e-36    TRUE
	# 1       1033999 29      14      0.938095583505757       0.41    0.425   0.4225  0.456666666667  0.45125 32.8498752609621        1.95008350059784e-18    9.52022158314114e-28    TRUE
	# 1       1123434 30      14      1.06667486722421        0.41    0.375   0.39    0.415   0.41375 37.9156217530781        3.82711658431309e-18    1.11970565041553e-29    TRUE
	# 1       1124663 39      18      1.02375023806503        0.4     0.415   0.37    0.366666666667  0.345   36.6962795733276        2.46260220315049e-22    6.11668474858626e-38    TRUE
	# 1       1130206 40      22      0.958189219225845       0.31    0.385   0.43    0.441666666667  0.4475  33.7503815149709        5.35710573713984e-29    5.18247336690237e-31    TRUE
	# 1       1491251 20      11      0.964912131925299       0.52    0.52    0.53    0.541666666667  0.54625 31.8670433818152        4.44540545905524e-16    3.01502604247324e-17    TRUE
	# 1       1503099 27      12      0.89635832629498        0.48    0.47    0.495   0.48    0.46875 30.4466481668093        1.26409301910224e-15    1.85217227198124e-28    TRUE

	# constants
	defaut_fit_f_from_model = "0.8422122902" # 0.842212290208

	# values to be determined for this sample
	scCNV_QC = 'undetermined' # 'pass' or 'fail'
	num_bins_for_this_sample = 0
	num_events_where_non_0p5_vaf_or_non_zero_depth = 0
	num_events_non_0p5_vaf = 0
	total_size_in_bins_non_0p5_vaf = 0
	total_size_in_nucleotides_non_0p5_vaf = 0
	avg_size_in_bins_non_0p5_vaf = 0
	avg_size_in_nucleotides_non_0p5_vaf = 0
	num_events_same_vaf_diff_chrom = 0
	size_in_bins_same_vaf_diff_chrom = 0
	size_in_nucleotides_same_vaf_diff_chrom = 0
	avg_non_0p5_vaf_dist_from_0p5 = 0
	sum_non_0p5_vaf_dist_from_0p5 = 0
	num_events_non_0p5_vaf = 0
	avg_model_fit_f = 0
	avg_adjusted_model_fit_f = 0
	sum_fit_f = 0
	sum_adjusted_fit_f = 0
	num_fit_f = 0
	num_diff_vaf1_values = 0
	sum_num_events_for_same_vaf1 = 0
	avg_num_events_for_same_vaf1_values = 0
	num_events_non_zero_depth = 0
	total_size_in_bins_non_zero_depth = 0
	total_size_in_nucleotides_non_zero_depth = 0
	avg_size_in_bins_non_zero_depth = 0
	avg_size_in_nucleotides_non_zero_depth = 0
	sum_non_zero_depth = 0
	sum_abs_value_non_zero_depth = 0
	avg_non_zero_depth = 0
	avg_abs_value_non_zero_depth = 0
	num_diff_depth_values = 0
	sum_num_events_for_same_depth = 0
	avg_num_events_for_same_depth_values = 0
	num_non_normal_vaf_and_depth_both = 0
	num_non_normal_vaf_or_depth_but_not_both = 0
	num_telomere_non_normal_vaf_and_depth = 0
        num_non_telomere_non_normal_vaf_and_depth = 0
	non_0p5_vaf_start_chrom = []
	non_0p5_vaf_end_chrom = []
	non_0p5_vaf1 = []
	non_0p5_vaf2 = []
	non_0p5_vaf_size_in_bins = []
	non_0p5_vaf_size_in_nucleotides = []
	list_of_diff_vaf1 = {}
	list_of_diff_depth = {}
	has_default_fit_score = "No"

	# Nucleotide positions for the hg19 genome. Our bins are 4000, so need to divide by 4000.
	# 	chrom	chrom_length	chrom_start	chrom_end
	# 	chr1	249250621	1		249250621
	# 	chr2	243199373	249250622	492449994
	# 	chr3	198022430	492449995	690472424
	# 	chr4	191154276	690472425	881626700
	# 	chr5	180915260	881626701	1062541960
	# 	chr6	171115067	1062541961	1233657027
	# 	chr7	159138663	1233657028	1392795690
	# 	chr8	146364022	1392795691	1539159712
	# 	chr9	141213431	1539159713	1680373143
	# 	chr10	135534747	1680373144	1815907890
	# 	chr11	135006516	1815907891	1950914406
	# 	chr12	133851895	1950914407	2084766301
	# 	chr13	115169878	2084766302	2199936179
	# 	chr14	107349540	2199936180	2307285719
	# 	chr15	102531392	2307285720	2409817111
	# 	chr16	90354753	2409817112	2500171864
	# 	chr17	81195210	2500171865	2581367074
	# 	chr18	78077248	2581367075	2659444322
	# 	chr19	59128983	2659444323	2718573305
	# 	chr20	63025520	2718573306	2781598825
	# 	chr21	48129895	2781598826	2829728720
	# 	chr22	51304566	2829728721	2881033286
	# 	chrX	155270560	2881033287	3036303846
	# 	chrY	59373566	3036303847	3095677412

	chromosome_size_in_nucleotides = [ 0, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566 ]
	# chromosome_bin_start_positions = [ 0, 1, 249250622, 492449995, 690472425, 881626701, 1062541961, 1233657028, 1392795691, 1539159713, 1680373144, 1815907891, 1950914407, 2084766302, 2199936180, 2307285720, 2409817112, 2500171865, 2581367075, 2659444323, 2718573306, 2781598826, 2829728721, 2881033287, 3036303847 ]
	# chromosome_bin_end_positions = [ 0, 249250621, 492449994, 690472424, 881626700, 1062541960, 1233657027, 1392795690, 1539159712, 1680373143, 1815907890, 1950914406, 2084766301, 2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322, 2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412 ] 

	# The two input files have the same number of lines and same line numbers are for the same copy-number event.

	infile = open(args.indata, 'r')
	inlines = infile.readlines() 
	infile.close() 

	infile2 = open(args.indata2, 'r')
	inlines2 = infile2.readlines() 
	infile2.close()

	infile3 = open(args.indata3, 'r')
	inlines3 = infile3.readlines() 
	infile3.close()

	# input values for each event of this sample
	chrom_array = []
	start_chrom_array = []
	end_chrom_array = []
	start_index_array = []
	end_index_array = []
	depth_array = []
	vaf1_array = []
	vaf2_array = []
	fit_f_array = []
	adjusted_fit_f_array = []
	event_category = []
	event_CN_category = []
	event_subcategory = []
	chrom_array2 = []
	window_id_array2 = []
	start_index_array2 = []
	end_index_array2 = []
	start_pos_array2 = []
	end_pos_array2 = []
	fit_k1_array2 = []
	fit_k2_array2 = []
	fit_f_array2 = []
	fit_isize_array2 = []
	fit_llik_array2 = []
	chrom_array3 = []
	start_pos_array3 = []
	end_pos_array3 = []
	start_index_array3 = []
	end_index_array3 = []
	gene_array3 = []

	in_header = True
	for i in range( 0, len(inlines) ):
		inline = inlines[i]
		inline = inline.strip()
		if (inline != ''):
			if (in_header):
				in_header = False
			else:
				infields = inline.split("\t")
				start_chrom = str(infields[0])
				end_chrom = str(infields[1])
				start_index = int(infields[2])
				end_index = int(infields[3])
				depth = float(infields[4])
				vaf1 = float(infields[5])
				vaf2 = float(infields[6])
				fit_f = float(infields[7])
				start_chrom_array.append( start_chrom )
				end_chrom_array.append( end_chrom )
				start_index_array.append( start_index )
				end_index_array.append( end_index )
				depth_array.append( depth )
				vaf1_array.append( vaf1 )
				vaf2_array.append( vaf2 )
				fit_f_array.append( fit_f )
				adjusted_fit_f_array.append( fit_f )
				event_category.append( '' )
				event_CN_category.append( '' )
				event_subcategory.append( '' )
				num_bins_for_this_sample = end_index

	in_header = True
	for i in range( 0, len(inlines2) ):
		inline = inlines2[i]
		inline = inline.strip()
		if (inline != ''):
			if (in_header):
				in_header = False
			else:
				infields = inline.split("\t")
				chrom = str(infields[0])
				window_id = str(infields[1])
				start_index_file2 = int(infields[2])
				end_index_file2 = int(infields[3])
				start_pos = int(infields[4])
				end_pos = int(infields[5])
				fit_k1 = float(infields[6])
				fit_k2 = float(infields[7])
				fit_f_file2 = float(infields[8])
				fit_isize = float(infields[9])
				fit_llik = float(infields[10])
				chrom_array2.append( chrom )
				window_id_array2.append( window_id )
				start_index_array2.append( start_index_file2 )
				end_index_array2.append( end_index_file2 )
				start_pos_array2.append( start_pos )
				end_pos_array2.append( end_pos )
				fit_k1_array2.append( fit_k1 )
				fit_k2_array2.append( fit_k2 )
				fit_f_array2.append( fit_f_file2 )
				fit_isize_array2.append( fit_isize )
				fit_llik_array2.append( fit_llik )

	in_header = False
	for i in range( 0, len(inlines3) ):
		inline = inlines3[i]
		inline = inline.strip()
		if (inline != ''):
			if (in_header):
				in_header = False
			else:
				infields = inline.split("\t")
				chrom_file3 = str(infields[0])
				start_pos_file3 = int(infields[1])
				end_pos_file3 = int(infields[2])
				start_index_file3 = int(infields[3])
				end_index_file3 = int(infields[4])
				gene = str(infields[5])
				# These arrays have one entry per non-zero event per gene hit
				chrom_array3.append( chrom_file3 )
				start_pos_array3.append( start_pos_file3 )
				end_pos_array3.append( end_pos_file3 )
				start_index_array3.append( start_index_file3 )
				end_index_array3.append( end_index_file3 )
				gene_array3.append( gene )

	# This array has one entry per event (whether non-zero or zero event)
	# Each entry contains a list of genes hit by that event. List is delimited by pipe |
	array_of_genes_hit_by_event = convert_gene_list_to_array_having_standard_index( chrom_array2, start_pos_array2, end_pos_array2, start_index_array2, end_index_array2, chrom_array3, start_pos_array3, end_pos_array3, start_index_array3, end_index_array3, gene_array3 )

	#for i in range( 0, len(array_of_genes_hit_by_event) ):
	#	print i, chrom_array2[i], start_index_array2[i], end_index_array2[i], array_of_genes_hit_by_event[i]

	# counts and tallies for this sample
	num_bins_for_this_sample_FLOAT = float(num_bins_for_this_sample)
	num_category_is_no_subclones = 0
	num_category_is_germline_dup = 0
	num_category_is_obvious_subclones = 0
	num_category_is_telomere_subclones = 0
	num_category_is_focal_subclones = 0
	sum_adjusted_fit_f_category_is_no_subclones = float(0)
	sum_adjusted_fit_f_category_is_germline_dup = float(0)
	sum_adjusted_fit_f_category_is_obvious_subclones = float(0)
	sum_adjusted_fit_f_category_is_telomere_subclones = float(0)
	sum_adjusted_fit_f_category_is_focal_subclones = float(0)

	# counts and tallies for this sample
	num_negative_depth_and_non_0p5_vaf_for_obvious_subclones = 0
	total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_obvious_subclones = 0
	num_positive_depth_and_non_0p5_vaf_obvious_subclones = 0
	total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_obvious_subclones = 0
	num_zero_depth_and_non_0p5_vaf_obvious_subclones = 0
	total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_obvious_subclones = 0
	num_negative_depth_and_non_0p5_vaf_for_focal_subclones = 0
	total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_focal_subclones = 0
	num_positive_depth_and_non_0p5_vaf_focal_subclones = 0
	total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_focal_subclones = 0
	num_zero_depth_and_non_0p5_vaf_focal_subclones = 0
	total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_focal_subclones = 0

	num_CN_gain_events = 0
	total_size_in_nucleotides_CN_gain_events = 0
	avg_positive_depth_per_nucleotide_for_CN_gain_events = 0
	sum_positive_depth_for_CN_gain_events = 0

	num_CN_loss_events = 0
	total_size_in_nucleotides_CN_loss_events = 0
	avg_negative_depth_per_nucleotide_for_CN_loss_events = 0
	sum_negative_depth_for_CN_loss_events = 0

	num_CN_LoH_events = 0
	total_size_in_nucleotides_CN_LoH_events = 0

	num_positive_depth_and_non_0p5_vaf_events = 0
	total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_events = 0
	avg_positive_depth_per_nucleotide_for_positive_depth_and_non_0p5_vaf_events = 0
	sum_positive_depth_for_positive_depth_and_non_0p5_vaf_events = 0

	num_negative_depth_and_non_0p5_vaf_events = 0
	total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_events = 0
	avg_negative_depth_per_nucleotide_for_negative_depth_and_non_0p5_vaf_events = 0
	sum_negative_depth_for_negative_depth_and_non_0p5_vaf_events = 0

	num_zero_depth_and_non_0p5_vaf_events = 0
	total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_events = 0

	genes_hit_by_positive_depth_and_non_0p5_vaf_events = ''
	genes_hit_by_negative_depth_and_non_0p5_vaf_events = ''
	genes_hit_by_zero_depth_and_non_0p5_vaf_events = ''
	genes_hit_by_CN_gain_events = ''
	genes_hit_by_CN_loss_events = ''
	genes_hit_by_CN_LoH = ''

	telomere_loss = {}
	telomere_loss['1p'] = 0
	telomere_loss['1q'] = 0
	telomere_loss['2p'] = 0
	telomere_loss['2q'] = 0
	telomere_loss['3p'] = 0
	telomere_loss['3q'] = 0
	telomere_loss['4p'] = 0
	telomere_loss['4q'] = 0
	telomere_loss['5p'] = 0
	telomere_loss['5q'] = 0
	telomere_loss['6p'] = 0
	telomere_loss['6q'] = 0
	telomere_loss['7p'] = 0
	telomere_loss['7q'] = 0
	telomere_loss['8p'] = 0
	telomere_loss['8q'] = 0
	telomere_loss['9p'] = 0
	telomere_loss['9q'] = 0
	telomere_loss['10p'] = 0
	telomere_loss['10q'] = 0
	telomere_loss['11p'] = 0
	telomere_loss['11q'] = 0
	telomere_loss['12p'] = 0
	telomere_loss['12q'] = 0
	telomere_loss['13p'] = 0
	telomere_loss['13q'] = 0
	telomere_loss['14p'] = 0
	telomere_loss['14q'] = 0
	telomere_loss['15p'] = 0
	telomere_loss['15q'] = 0
	telomere_loss['16p'] = 0
	telomere_loss['16q'] = 0
	telomere_loss['17p'] = 0
	telomere_loss['17q'] = 0
	telomere_loss['18p'] = 0
	telomere_loss['18q'] = 0
	telomere_loss['19p'] = 0
	telomere_loss['19q'] = 0
	telomere_loss['20p'] = 0
	telomere_loss['20q'] = 0
	telomere_loss['21p'] = 0
	telomere_loss['21q'] = 0
	telomere_loss['22p'] = 0
	telomere_loss['22q'] = 0

	entire_chromosome_event = {}
	entire_chromosome_event['1'] = 'no'
	entire_chromosome_event['2'] = 'no'
	entire_chromosome_event['3'] = 'no'
	entire_chromosome_event['4'] = 'no'
	entire_chromosome_event['5'] = 'no'
	entire_chromosome_event['6'] = 'no'
	entire_chromosome_event['7'] = 'no'
	entire_chromosome_event['8'] = 'no'
	entire_chromosome_event['9'] = 'no'
	entire_chromosome_event['10'] = 'no'
	entire_chromosome_event['11'] = 'no'
	entire_chromosome_event['12'] = 'no'
	entire_chromosome_event['13'] = 'no'
	entire_chromosome_event['14'] = 'no'
	entire_chromosome_event['15'] = 'no'
	entire_chromosome_event['16'] = 'no'
	entire_chromosome_event['17'] = 'no'
	entire_chromosome_event['18'] = 'no'
	entire_chromosome_event['19'] = 'no'
	entire_chromosome_event['20'] = 'no'
	entire_chromosome_event['21'] = 'no'
	entire_chromosome_event['22'] = 'no'
	entire_chromosome_event['X'] = 'no'
	entire_chromosome_event['Y'] = 'no'

	# bin sizes taken from one of the samples, to be used for all samples
	# 22 chromosomes + 3 more not used
	num_bins_for_each_chromosome = [0, 57110, 57110, 55718, 54195, 50507, 50507, 50507, 42264, 42264, 35315, 35315, 3381, 3381, 22837, 18632, 21026, 15001, 21844, 21844, 15365, 10127, 6088, 0, 0, 0] 
	debug = False
	if (debug):
		num_bins_for_each_chromosome = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] # 22 chromosomes + 3 more not used
		curr_chrom = ''
		curr_start_index = 0
		curr_end_index = 0
		for i in range( 0, len(start_chrom_array) ):
			start_chrom = start_chrom_array[i]
			end_chrom = end_chrom_array[i]
			start_index = start_index_array[i]
			end_index = end_index_array[i]
			if (curr_chrom != start_chrom):
				if (curr_chrom != ''):
					num_bins_for_each_chromosome[ index_for_chrom(curr_chrom) ] = curr_end_index - curr_start_index + 1
				curr_chrom = start_chrom
				curr_start_index = start_index
			curr_end_index = end_index
		if (curr_chrom != ''):
			num_bins_for_each_chromosome[ index_for_chrom(curr_chrom) ] = curr_end_index - curr_start_index + 1

	for i in range( 0, len(start_chrom_array) ):

		# input variables for this event from input file 1
		start_chrom = start_chrom_array[i]
		end_chrom = end_chrom_array[i]
		start_index = start_index_array[i]
		end_index = end_index_array[i]
		depth = depth_array[i]
		vaf1 = vaf1_array[i]
		vaf2 = vaf2_array[i]
		fit_f = fit_f_array[i]
		this_size_in_bins = end_index - start_index + 1
		vaf_dist_from_0p5 = abs(vaf1 - vaf2) / float(2)

		# input variables for this event from input file 2
		chrom = chrom_array2[i]
		window_id = window_id_array2[i]
		start_index_file2 = start_index_array2[i]
		end_index_file2 = end_index_array2[i]
		start_pos = start_pos_array2[i]
		end_pos = end_pos_array2[i]
		fit_k1 = fit_k1_array2[i]
		fit_k2 = fit_k2_array2[i]
		fit_f_file2 = fit_f_array2[i]
		fit_isize = fit_isize_array2[i]
		fit_llik = fit_llik_array2[i]
		this_size_in_nucleotides = end_pos - start_pos + 1

		adjusted_fit_f = fit_f
		fit_f_for_comparison = '{0:.10f}'.format(fit_f)
		if (fit_f_for_comparison == defaut_fit_f_from_model):
			has_default_fit_score = "Yes"
			adjusted_fit_f = 0
		adjusted_fit_f = float(adjusted_fit_f)
		adjusted_fit_f_array[i] = adjusted_fit_f

		num_bins_for_this_chromosome_FLOAT = float( num_bins_for_each_chromosome[ index_for_chrom(start_chrom) ] )

		# start_index = start_index_file2, they are the same value found in 2 different files
		# end_index = end_index_file2, they are the same value found in 2 different files
		# fit_f = fit_f_file2, they are sometimes the same value found in 2 different files. When the modeling fails, they don't have the same values.
		# if (start_index != start_index_file2): 
		#	print 'start_index != start_index_file2', start_index, start_index_file2
		# if (end_index != end_index_file2):
		#	print 'end_index != end_index_file2', end_index, end_index_file2
		# if (fit_f != fit_f_file2):
		#	print 'fit_f != fit_f_file2', fit_f, fit_f_file2 

		this_event_category, this_event_CN_category, this_event_subcategory, this_telomere_position, this_telomere_position_p, this_telomere_position_q, this_chromosome_event_for_chrom, this_chromosome_event_for_event = assign_category_of_this_event( i, chrom, start_chrom, end_chrom, start_pos, end_pos, start_index, end_index, depth, vaf1, vaf2, fit_f, vaf_dist_from_0p5, this_size_in_bins, fit_k1, fit_k2, fit_f_file2, this_size_in_nucleotides, num_bins_for_this_sample_FLOAT, num_bins_for_this_chromosome_FLOAT, chrom_array2, start_chrom_array, end_chrom_array, chromosome_size_in_nucleotides )
		event_category[i] = this_event_category
		event_CN_category[i] = this_event_CN_category
		event_subcategory[i] = this_event_subcategory
		if ((this_telomere_position != '') and (this_event_CN_category == "CN_loss")):
			telomere_loss[this_telomere_position] = 1
		if ((this_telomere_position_p != '') and (this_event_CN_category == "CN_loss")):
			telomere_loss[this_telomere_position_p] = 1
		if ((this_telomere_position_q != '') and (this_event_CN_category == "CN_loss")):
			telomere_loss[this_telomere_position_q] = 1
		if ((this_chromosome_event_for_chrom != '') and (this_chromosome_event_for_event != '')):
			entire_chromosome_event[this_chromosome_event_for_chrom] = this_chromosome_event_for_event

		if (this_event_category == "no_subclones"):
			num_category_is_no_subclones = num_category_is_no_subclones + 1
			sum_adjusted_fit_f_category_is_no_subclones = sum_adjusted_fit_f_category_is_no_subclones + adjusted_fit_f
		elif (this_event_category == "germline_dup"):
			num_category_is_germline_dup = num_category_is_germline_dup + 1
			sum_adjusted_fit_f_category_is_germline_dup = sum_adjusted_fit_f_category_is_germline_dup + adjusted_fit_f
		elif (this_event_category == "obvious_subclones"):
			num_category_is_obvious_subclones = num_category_is_obvious_subclones + 1
			sum_adjusted_fit_f_category_is_obvious_subclones = sum_adjusted_fit_f_category_is_obvious_subclones + adjusted_fit_f
		elif (this_event_category == "telomere_subclones"):
			num_category_is_telomere_subclones = num_category_is_telomere_subclones + 1
			sum_adjusted_fit_f_category_is_telomere_subclones = sum_adjusted_fit_f_category_is_telomere_subclones + adjusted_fit_f
		elif (this_event_category == "focal_subclones"):
			num_category_is_focal_subclones = num_category_is_focal_subclones + 1
			sum_adjusted_fit_f_category_is_focal_subclones = sum_adjusted_fit_f_category_is_focal_subclones + adjusted_fit_f

		if ((vaf1 != 0.5) or (vaf2 != 0.5)):
			num_events_non_0p5_vaf = num_events_non_0p5_vaf + 1
			total_size_in_bins_non_0p5_vaf = total_size_in_bins_non_0p5_vaf + this_size_in_bins
			total_size_in_nucleotides_non_0p5_vaf = total_size_in_nucleotides_non_0p5_vaf + this_size_in_nucleotides
			sum_non_0p5_vaf_dist_from_0p5 = sum_non_0p5_vaf_dist_from_0p5 + (abs(vaf1 - vaf2) / float(2))
			num_event_non_0p5_vaf = num_events_non_0p5_vaf + 1 
			non_0p5_vaf_start_chrom.append(start_chrom)
			non_0p5_vaf_end_chrom.append(end_chrom)
			non_0p5_vaf1.append(vaf1)
			non_0p5_vaf2.append(vaf2)
			non_0p5_vaf_size_in_bins.append( this_size_in_bins )
			non_0p5_vaf_size_in_nucleotides.append( this_size_in_nucleotides )
			vaf1_for_list = '{0:.6f}'.format(vaf1)
			if vaf1_for_list in list_of_diff_vaf1:
				list_of_diff_vaf1[ vaf1_for_list ] = list_of_diff_vaf1[ vaf1_for_list ] + 1
			else:
				list_of_diff_vaf1[ vaf1_for_list ] = 1

		if (depth != 0):
			num_events_non_zero_depth = num_events_non_zero_depth + 1
			total_size_in_bins_non_zero_depth = total_size_in_bins_non_zero_depth + this_size_in_bins
			total_size_in_nucleotides_non_zero_depth = total_size_in_nucleotides_non_zero_depth + this_size_in_nucleotides
			sum_non_zero_depth = sum_non_zero_depth + depth
			sum_abs_value_non_zero_depth = sum_abs_value_non_zero_depth + abs(depth)
			depth_for_list = '{0:.6f}'.format(depth)
			if depth_for_list in list_of_diff_depth:
				list_of_diff_depth[ depth_for_list ] = list_of_diff_depth[ depth_for_list ] + 1
			else:
				list_of_diff_depth[ depth_for_list ] = 1

		if ((vaf1 != 0.5) or (depth != 0)):
			num_events_where_non_0p5_vaf_or_non_zero_depth = num_events_where_non_0p5_vaf_or_non_zero_depth + 1

		if (vaf1 != 0.5):
			if (this_event_category == "obvious_subclones"):
				if (this_event_subcategory == "negative_depth_and_non_0p5_vaf"):
					num_negative_depth_and_non_0p5_vaf_for_obvious_subclones = num_negative_depth_and_non_0p5_vaf_for_obvious_subclones + 1
					total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_obvious_subclones = total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_obvious_subclones + this_size_in_nucleotides
				elif (this_event_subcategory == "positive_depth_and_non_0p5_vaf"):
					num_positive_depth_and_non_0p5_vaf_obvious_subclones = num_positive_depth_and_non_0p5_vaf_obvious_subclones + 1
					total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_obvious_subclones = total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_obvious_subclones + this_size_in_nucleotides
				elif (this_event_subcategory == "zero_depth_and_non_0p5_vaf"):
					num_zero_depth_and_non_0p5_vaf_obvious_subclones = num_zero_depth_and_non_0p5_vaf_obvious_subclones + 1
					total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_obvious_subclones = total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_obvious_subclones + this_size_in_nucleotides
			elif (this_event_category == "focal_subclones"):
				if (this_event_subcategory == "negative_depth_and_non_0p5_vaf"):
					num_negative_depth_and_non_0p5_vaf_for_focal_subclones = num_negative_depth_and_non_0p5_vaf_for_focal_subclones + 1
					total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_focal_subclones = total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_focal_subclones + this_size_in_nucleotides
				elif (this_event_subcategory == "positive_depth_and_non_0p5_vaf"):
					num_positive_depth_and_non_0p5_vaf_focal_subclones = num_positive_depth_and_non_0p5_vaf_focal_subclones + 1
					total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_focal_subclones = total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_focal_subclones + this_size_in_nucleotides
				elif (this_event_subcategory == "zero_depth_and_non_0p5_vaf"):
					num_zero_depth_and_non_0p5_vaf_focal_subclones = num_zero_depth_and_non_0p5_vaf_focal_subclones + 1
					total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_focal_subclones = total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_focal_subclones + this_size_in_nucleotides

		if (this_event_CN_category != ''):
			if (this_event_CN_category == 'CN_gain'):
				num_CN_gain_events = num_CN_gain_events + 1
				total_size_in_nucleotides_CN_gain_events = total_size_in_nucleotides_CN_gain_events + this_size_in_nucleotides
				sum_positive_depth_for_CN_gain_events = sum_positive_depth_for_CN_gain_events + (depth * this_size_in_nucleotides)
				if (genes_hit_by_CN_gain_events == ''):
					genes_hit_by_CN_gain_events = array_of_genes_hit_by_event[i]
				else:
					genes_hit_by_CN_gain_events = genes_hit_by_CN_gain_events + '|' + array_of_genes_hit_by_event[i]
			elif (this_event_CN_category == 'CN_loss'):
				num_CN_loss_events = num_CN_loss_events + 1
				total_size_in_nucleotides_CN_loss_events = total_size_in_nucleotides_CN_loss_events + this_size_in_nucleotides
				sum_negative_depth_for_CN_loss_events = sum_negative_depth_for_CN_loss_events + (depth * this_size_in_nucleotides)
				if (genes_hit_by_CN_loss_events == ''):
					genes_hit_by_CN_loss_events = array_of_genes_hit_by_event[i]
				else:
					genes_hit_by_CN_loss_events = genes_hit_by_negative_depth_and_non_0p5_vaf_events + '|' + array_of_genes_hit_by_event[i]
			elif (this_event_CN_category == 'CN_LoH'):
				num_CN_LoH_events = num_CN_LoH_events + 1
				total_size_in_nucleotides_CN_LoH_events = total_size_in_nucleotides_CN_LoH_events + this_size_in_nucleotides
				if (genes_hit_by_CN_LoH == ''):
					genes_hit_by_CN_LoH = array_of_genes_hit_by_event[i]
				else:
					genes_hit_by_CN_LoH = genes_hit_by_zero_depth_and_non_0p5_vaf_events + '|' + array_of_genes_hit_by_event[i]
			# else: this_event_CN_category == 'No_CN'

		if (this_event_subcategory != ''):
			if (this_event_subcategory == 'positive_depth_and_non_0p5_vaf'):
				num_positive_depth_and_non_0p5_vaf_events = num_positive_depth_and_non_0p5_vaf_events + 1
				total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_events = total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_events + this_size_in_nucleotides
				sum_positive_depth_for_positive_depth_and_non_0p5_vaf_events = sum_positive_depth_for_positive_depth_and_non_0p5_vaf_events + (depth * this_size_in_nucleotides)
				#if (genes_hit_by_positive_depth_and_non_0p5_vaf_events == ''):
				#	genes_hit_by_positive_depth_and_non_0p5_vaf_events = array_of_genes_hit_by_event[i]
				#else:
				#	genes_hit_by_positive_depth_and_non_0p5_vaf_events = genes_hit_by_positive_depth_and_non_0p5_vaf_events + '|' + array_of_genes_hit_by_event[i]
			elif (this_event_subcategory == 'negative_depth_and_non_0p5_vaf'):
				num_negative_depth_and_non_0p5_vaf_events = num_negative_depth_and_non_0p5_vaf_events + 1
				total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_events = total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_events + this_size_in_nucleotides
				sum_negative_depth_for_negative_depth_and_non_0p5_vaf_events = sum_negative_depth_for_negative_depth_and_non_0p5_vaf_events + (depth * this_size_in_nucleotides)
				#if (genes_hit_by_negative_depth_and_non_0p5_vaf_events == ''):
				#	genes_hit_by_negative_depth_and_non_0p5_vaf_events = array_of_genes_hit_by_event[i]
				#else:
				#	genes_hit_by_negative_depth_and_non_0p5_vaf_events = genes_hit_by_negative_depth_and_non_0p5_vaf_events + '|' + array_of_genes_hit_by_event[i]
			elif (this_event_subcategory == 'zero_depth_and_non_0p5_vaf'):
				num_zero_depth_and_non_0p5_vaf_events = num_zero_depth_and_non_0p5_vaf_events + 1
				total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_events = total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_events + this_size_in_nucleotides
				#if (genes_hit_by_zero_depth_and_non_0p5_vaf_events == ''):
				#	genes_hit_by_zero_depth_and_non_0p5_vaf_events = array_of_genes_hit_by_event[i]
				#else:
				#	genes_hit_by_zero_depth_and_non_0p5_vaf_events = genes_hit_by_zero_depth_and_non_0p5_vaf_events + '|' + array_of_genes_hit_by_event[i]

		this_bin_size = end_index - start_index + 1
		sum_fit_f = sum_fit_f + (fit_f * this_bin_size)
		sum_adjusted_fit_f = sum_adjusted_fit_f + (adjusted_fit_f * this_bin_size)
		num_fit_f = num_fit_f + this_bin_size

		if ((vaf1 != 0.5) and (depth != 0)):
			num_non_normal_vaf_and_depth_both = num_non_normal_vaf_and_depth_both + 1
		elif ((vaf1 != 0.5) or (depth != 0)):
			num_non_normal_vaf_or_depth_but_not_both = num_non_normal_vaf_or_depth_but_not_both + 1

		if ((vaf1 != 0.5) and (depth != 0)):
			if (this_event_category == "telomere_subclones"):
				num_telomere_non_normal_vaf_and_depth = num_telomere_non_normal_vaf_and_depth + 1
			else:
				num_non_telomere_non_normal_vaf_and_depth = num_non_telomere_non_normal_vaf_and_depth + 1

		#print i, chrom, start_index_file2, end_index_file2, this_event_subcategory, this_event_CN_category, 'depth', depth

	# Tally the counts
	num_non_telomere_events_where_non_0p5_vaf_or_non_zero_depth = num_events_where_non_0p5_vaf_or_non_zero_depth - num_category_is_telomere_subclones
	if (num_non_telomere_events_where_non_0p5_vaf_or_non_zero_depth < 0):
		num_non_telomere_events_where_non_0p5_vaf_or_non_zero_depth = 0
	avg_adjusted_fit_f_category_is_no_subclones = 0
	if (num_category_is_no_subclones > 0):
		avg_adjusted_fit_f_category_is_no_subclones = sum_adjusted_fit_f_category_is_no_subclones / float(num_category_is_no_subclones)
	avg_adjusted_fit_f_category_is_germline_dup = 0
	if (num_category_is_germline_dup > 0):
		avg_adjusted_fit_f_category_is_germline_dup = sum_adjusted_fit_f_category_is_germline_dup / float(num_category_is_germline_dup)
	avg_adjusted_fit_f_category_is_obvious_subclones = 0
	if (num_category_is_obvious_subclones > 0):
		avg_adjusted_fit_f_category_is_obvious_subclones = sum_adjusted_fit_f_category_is_obvious_subclones / float(num_category_is_obvious_subclones)
	avg_adjusted_fit_f_category_is_telomere_subclones = 0
	if (num_category_is_telomere_subclones > 0):
		avg_adjusted_fit_f_category_is_telomere_subclones = sum_adjusted_fit_f_category_is_telomere_subclones / float(num_category_is_telomere_subclones)
	avg_adjusted_fit_f_category_is_focal_subclones = 0
	if (num_category_is_focal_subclones > 0):
		avg_adjusted_fit_f_category_is_focal_subclones = sum_adjusted_fit_f_category_is_focal_subclones / float(num_category_is_focal_subclones)

	if (total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_events > 0):
		avg_positive_depth_per_nucleotide_for_positive_depth_and_non_0p5_vaf_events = float(sum_positive_depth_for_positive_depth_and_non_0p5_vaf_events) / float(total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_events)
	if (total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_events):
		avg_negative_depth_per_nucleotide_for_negative_depth_and_non_0p5_vaf_events = float(sum_negative_depth_for_negative_depth_and_non_0p5_vaf_events) / float(total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_events)

	if (total_size_in_nucleotides_CN_gain_events > 0):
		avg_positive_depth_per_nucleotide_for_CN_gain_events = float(sum_positive_depth_for_CN_gain_events) / float(total_size_in_nucleotides_CN_gain_events)
	if (total_size_in_nucleotides_CN_loss_events):
		avg_negative_depth_per_nucleotide_for_CN_loss_events = float(sum_negative_depth_for_CN_loss_events) / float(total_size_in_nucleotides_CN_loss_events)

	for i in range( 0, len(non_0p5_vaf_start_chrom) ):
		this_start_chrom = non_0p5_vaf_start_chrom[i]
		this_end_chrom = non_0p5_vaf_end_chrom[i]
		this_vaf1 = non_0p5_vaf1[i]
		this_vaf2 = non_0p5_vaf2[i]
		this_size_in_bins = non_0p5_vaf_size_in_bins[i]
		this_size_in_nucleotides = non_0p5_vaf_size_in_nucleotides[i]

		has_diff_chrom_with_same_vaf_values = does_it_have_diff_chrom_with_same_vaf_values( this_start_chrom, this_end_chrom, this_vaf1, this_vaf2, non_0p5_vaf_start_chrom, non_0p5_vaf_end_chrom, non_0p5_vaf1, non_0p5_vaf2 )
		if (has_diff_chrom_with_same_vaf_values):
			num_events_same_vaf_diff_chrom = num_events_same_vaf_diff_chrom + 1
			size_in_bins_same_vaf_diff_chrom = size_in_bins_same_vaf_diff_chrom + this_size_in_bins
			size_in_nucleotides_same_vaf_diff_chrom = size_in_nucleotides_same_vaf_diff_chrom + this_size_in_nucleotides

	if (num_events_non_0p5_vaf > 0):
		avg_non_0p5_vaf_dist_from_0p5 = float(sum_non_0p5_vaf_dist_from_0p5) / float(num_events_non_0p5_vaf)
		avg_size_in_bins_non_0p5_vaf = int(round( float(total_size_in_bins_non_0p5_vaf) / float(num_events_non_0p5_vaf) ))
		avg_size_in_nucleotides_non_0p5_vaf = int(round( float(total_size_in_nucleotides_non_0p5_vaf) / float(num_events_non_0p5_vaf) ))
	if (num_fit_f > 0):
		avg_model_fit_f = float(sum_fit_f) / float(num_fit_f)
		avg_adjusted_model_fit_f = float(sum_adjusted_fit_f) / float(num_fit_f)

	if (num_events_non_zero_depth > 0):
		avg_non_zero_depth = float(sum_non_zero_depth) / float(num_events_non_zero_depth)
		avg_abs_value_non_zero_depth = float(sum_abs_value_non_zero_depth) / float(num_events_non_zero_depth)
		avg_size_in_bins_non_zero_depth = int(round( float(total_size_in_bins_non_zero_depth) / float(num_events_non_zero_depth) ))
		avg_size_in_nucleotides_non_zero_depth = int(round( float(total_size_in_nucleotides_non_zero_depth) / float(num_events_non_zero_depth) ))

	for this_vaf1 in list_of_diff_vaf1:
		num_diff_vaf1_values = num_diff_vaf1_values + 1
		this_vaf1_num_events = list_of_diff_vaf1[ this_vaf1 ]
		sum_num_events_for_same_vaf1 = sum_num_events_for_same_vaf1 + this_vaf1_num_events
	if (num_diff_vaf1_values > 0):
		avg_num_events_for_same_vaf1_values = float(sum_num_events_for_same_vaf1) / float(num_diff_vaf1_values)

	for this_depth in list_of_diff_depth:
		num_diff_depth_values = num_diff_depth_values + 1
		this_depth_num_events = list_of_diff_depth[ this_depth ]
		sum_num_events_for_same_depth = sum_num_events_for_same_depth + this_depth_num_events
	if (num_diff_depth_values > 0):
		avg_num_events_for_same_depth_values = float(sum_num_events_for_same_depth) / float(num_diff_depth_values)

	num_telomere_non_normal_vaf_and_depth_FRACTION = 0
	if ( float(num_non_normal_vaf_and_depth_both + num_non_normal_vaf_or_depth_but_not_both) > 0 ):
		num_telomere_non_normal_vaf_and_depth_FRACTION = float(num_telomere_non_normal_vaf_and_depth) / float(num_non_normal_vaf_and_depth_both + num_non_normal_vaf_or_depth_but_not_both)

	num_bins_for_this_sample = float(num_bins_for_this_sample)
	std_num_bins_for_a_sample = float(730211)
	num_events_same_vaf_diff_chrom_FRACTION = 0
	if (num_bins_for_this_sample > 0):
		num_events_same_vaf_diff_chrom_FRACTION = float(num_events_same_vaf_diff_chrom) / num_bins_for_this_sample

	# Here is the algorithm used to assign categories for the sample, 
	# derived from data having number of bins of 730211, 705425, 720025, 713588, 709942, 713357, or 704841 (avg 730211) across all the chromosomes.
	#
	# failed_modelling : (total_size_in_bins_non_0p5_vaf >= 700,000 (this is the highway double lines across the vaf chart)) 
	#
	# failed_modelling : (num_events_where_non_0p5_vaf_or_non_zero_depth > 50) and (avg_adjusted_model_fit_f < 0.03)
	#
	# failed_modelling : (num_events_where_non_0p5_vaf_or_non_zero_depth > 20) and (avg_adjusted_model_fit_f >= 0.03)
	#
	# failed_modelling : (num_events_where_non_0p5_vaf_or_non_zero_depth > 20) and (avg_adjusted_model_fit_f > 0.015) and (total_size_in_bins_non_0p5_vaf >= 50,000)
	#
	# obvious_subclones : (there are any obvious_subclones events which are events greater than a minimum size)
	#
	# germline_dup : (there are at least as many germline_dup events as there are telomere events (or all telomere events are germline_dup_events))
	#
	# telomere_subclones : (there are at least as many telomere events as there are non-telomere events)
	#
	# focal_subclones : (any events where vaf != 0.5 or depth != 0) and (avg_adjusted_model_fit_f >= 0.03) and they are not large enough to be classified as obvious_subclones
	#
	# no_subclones: (all events have vaf=0.5 and depth==0)

	this_sample_category = "no_subclones"

	#print 'num_category_is_obvious_subclones', num_category_is_obvious_subclones, '<==='
	#print 'avg_adjusted_fit_f_category_is_obvious_subclones', avg_adjusted_fit_f_category_is_obvious_subclones, '<==='
	#print 'num_category_is_telomere_subclones', num_category_is_telomere_subclones
	#print 'num_non_telomere_non_normal_vaf_and_depth', num_non_telomere_non_normal_vaf_and_depth
	#print 'num_category_is_germline_dup', num_category_is_germline_dup
	#print 'num_non_telomere_events_where_non_0p5_vaf_or_non_zero_depth', num_non_telomere_events_where_non_0p5_vaf_or_non_zero_depth
	#print 'num_events_where_non_0p5_vaf_or_non_zero_depth', num_events_where_non_0p5_vaf_or_non_zero_depth
	#print 'avg_adjusted_model_fit_f', avg_adjusted_model_fit_f
	#print 'num_category_is_focal_subclones', num_category_is_focal_subclones
	#print 'total_size_in_bins_non_0p5_vaf', total_size_in_bins_non_0p5_vaf
	#print 'avg_size_in_bins_non_0p5_vaf', avg_size_in_bins_non_0p5_vaf
	#print 'size_in_bins_same_vaf_diff_chrom', size_in_bins_same_vaf_diff_chrom
	#print 'total_size_in_bins_non_zero_depth', total_size_in_bins_non_zero_depth
	#print 'avg_size_in_bins_non_zero_depth', avg_size_in_bins_non_zero_depth
	#print 'avg_non_0p5_vaf_dist_from_0p5', avg_non_0p5_vaf_dist_from_0p5
	#print 'avg_non_zero_depth', avg_non_zero_depth

	if (total_size_in_bins_non_0p5_vaf >= 700000):
		this_sample_category = "failed_modelling"

	elif ((num_non_telomere_events_where_non_0p5_vaf_or_non_zero_depth > 50) and (avg_adjusted_model_fit_f < 0.03)):
		this_sample_category = "failed_modelling"

	elif ((num_non_telomere_events_where_non_0p5_vaf_or_non_zero_depth > 30) and (avg_adjusted_model_fit_f >= 0.03)):
		this_sample_category = "failed_modelling"

	elif ((num_non_telomere_events_where_non_0p5_vaf_or_non_zero_depth > 20) and (avg_adjusted_model_fit_f > 0.015) and (total_size_in_bins_non_0p5_vaf >= 50000)):
		this_sample_category = "failed_modelling"

	elif (num_category_is_obvious_subclones >= 1):
		this_sample_category = "obvious_subclones"

	elif ( (num_category_is_germline_dup >= 1) and (num_category_is_germline_dup >= num_category_is_telomere_subclones) ):
		this_sample_category = "germline_dup"

	elif ( (num_category_is_telomere_subclones >= 1) and (num_category_is_telomere_subclones >= num_non_telomere_non_normal_vaf_and_depth) ):
		this_sample_category = "telomere_subclones"

	elif ((num_category_is_focal_subclones > 0) and (avg_adjusted_model_fit_f >= 0.03)):
		this_sample_category = "focal_subclones"

	# Assign the algorithm_category
	algorithm_category = this_sample_category
	algorithm_CN_category = ''

	# Assign the scCNV_QC
	scCNV_QC = "pass"
	if (algorithm_category == "failed_modelling"):
		scCNV_QC = "fail"

	# Assign the algorithm_CN_category
	if (scCNV_QC == "pass"):
		if ((algorithm_category == "obvious_subclones") | (algorithm_category == "focal_subclones")):
			if ((total_size_in_nucleotides_CN_gain_events >= total_size_in_nucleotides_CN_loss_events) and (total_size_in_nucleotides_CN_gain_events >= total_size_in_nucleotides_CN_LoH_events)):
				algorithm_CN_category = 'CN_gain'
			elif (total_size_in_nucleotides_CN_loss_events > total_size_in_nucleotides_CN_LoH_events):
				algorithm_CN_category = 'CN_loss'
			elif (total_size_in_nucleotides_CN_LoH_events > 1000):
				algorithm_CN_category = 'CN_LoH'
			else:
				algorithm_CN_category = 'No_CN'

	debug = False
	if (debug):
		print 'algorithm_category', algorithm_category

	# Tidy up the list of genes hit to remove duplicates
	#genes_hit_by_positive_depth_and_non_0p5_vaf_events = make_list_unique( genes_hit_by_positive_depth_and_non_0p5_vaf_events )
	#genes_hit_by_negative_depth_and_non_0p5_vaf_events = make_list_unique( genes_hit_by_negative_depth_and_non_0p5_vaf_events )
	#genes_hit_by_zero_depth_and_non_0p5_vaf_events = make_list_unique( genes_hit_by_zero_depth_and_non_0p5_vaf_events )
	genes_hit_by_CN_gain_events = make_list_unique( genes_hit_by_CN_gain_events )
	genes_hit_by_CN_loss_events = make_list_unique( genes_hit_by_CN_loss_events )
	genes_hit_by_CN_LoH = make_list_unique( genes_hit_by_CN_LoH )

	# Output the counts

	output_file = open(args.outdata, 'w')

	outline_part1 = "sample_id\tscCNV_QC\talgorithm_category\talgorithm_CN_category\tavg_adjusted_model_fit_f\tavg_model_fit_f\tnum_CN_gain_events\ttotal_size_in_nucleotides_CN_gain_events\tavg_positive_depth_per_nucleotide_for_CN_gain_events\tnum_CN_loss_events\ttotal_size_in_nucleotides_CN_loss_events\tavg_negative_depth_per_nucleotide_for_CN_loss_events\tnum_CN_LoH\ttotal_size_in_nucleotides_CN_LoH_events"

	outline_part2 = "num_events_where_non_0p5_vaf_or_non_zero_depth\tnum_events_non_0p5_vaf\ttotal_size_in_nucleotides_non_0p5_vaf\tavg_size_in_nucleotides_non_0p5_vaf\tavg_non_0p5_vaf_dist_from_0p5\tnum_events_same_vaf_diff_chrom\tsize_in_nucleotides_same_vaf_diff_chrom\tnum_diff_vaf1_values\tavg_num_events_for_same_vaf1_values\tnum_events_non_zero_depth\ttotal_size_in_nucleotides_non_zero_depth\tavg_size_in_nucleotides_non_zero_depth\tavg_non_zero_depth\tavg_abs_value_non_zero_depth\tnum_diff_depth_values\tavg_num_events_for_same_depth_values\tnum_non_normal_vaf_and_depth_both\tnum_non_normal_vaf_or_depth_but_not_both\tnum_telomere_non_normal_vaf_and_depth\tnum_non_telomere_non_normal_vaf_and_depth\thas_default_fit_score\tnum_category_is_no_subclones\tavg_adjusted_fit_f_category_is_no_subclones\tnum_category_is_germline_dup\tavg_adjusted_fit_f_category_is_germline_dup\tnum_category_is_obvious_subclones\tavg_adjusted_fit_f_category_is_obvious_subclones\tnum_category_is_telomere_subclones\tavg_adjusted_fit_f_category_is_telomere_subclones\tnum_category_is_focal_subclones\tavg_adjusted_fit_f_category_is_focal_subclones"

	outline_part3 = "num_negative_depth_and_non_0p5_vaf_for_obvious_subclones\ttotal_size_in_nucleotides_negative_depth_and_non_0p5_vaf_obvious_subclones\tnum_positive_depth_and_non_0p5_vaf_obvious_subclones\ttotal_size_in_nucleotides_positive_depth_and_non_0p5_vaf_obvious_subclones\tnum_zero_depth_and_non_0p5_vaf_obvious_subclones\ttotal_size_in_nucleotides_zero_depth_and_non_0p5_vaf_obvious_subclones\tnum_negative_depth_and_non_0p5_vaf_for_focal_subclones\ttotal_size_in_nucleotides_negative_depth_and_non_0p5_vaf_focal_subclones\tnum_positive_depth_and_non_0p5_vaf_focal_subclones\ttotal_size_in_nucleotides_positive_depth_and_non_0p5_vaf_focal_subclones\tnum_zero_depth_and_non_0p5_vaf_focal_subclones\ttotal_size_in_nucleotides_zero_depth_and_non_0p5_vaf_focal_subclones"

	outline_part4 = "telomere_loss_1p\ttelomere_loss_1q\ttelomere_loss_2p\ttelomere_loss_2q\ttelomere_loss_3p\ttelomere_loss_3q\ttelomere_loss_4p\ttelomere_loss_4q\ttelomere_loss_5p\ttelomere_loss_5q\ttelomere_loss_6p\ttelomere_loss_6q\ttelomere_loss_7p\ttelomere_loss_7q\ttelomere_loss_8p\ttelomere_loss_8q\ttelomere_loss_9p\ttelomere_loss_9q\ttelomere_loss_10p\ttelomere_loss_10q\ttelomere_loss_11p\ttelomere_loss_11q\ttelomere_loss_12p\ttelomere_loss_12q\ttelomere_loss_13p\ttelomere_loss_13q\ttelomere_loss_14p\ttelomere_loss_14q\ttelomere_loss_15p\ttelomere_loss_15q\ttelomere_loss_16p\ttelomere_loss_16q\ttelomere_loss_17p\ttelomere_loss_17q\ttelomere_loss_18p\ttelomere_loss_18q\ttelomere_loss_19p\ttelomere_loss_19q\ttelomere_loss_20p\ttelomere_loss_20q\ttelomere_loss_21p\ttelomere_loss_21q\ttelomere_loss_22p\ttelomere_loss_22q"

	outline_part5 = "entire_chromosome_event_1\tentire_chromosome_event_2\tentire_chromosome_event_3\tentire_chromosome_event_4\tentire_chromosome_event_5\tentire_chromosome_event_6\tentire_chromosome_event_7\tentire_chromosome_event_8\tentire_chromosome_event_9\tentire_chromosome_event_10\tentire_chromosome_event_11\tentire_chromosome_event_12\tentire_chromosome_event_13\tentire_chromosome_event_14\tentire_chromosome_event_15\tentire_chromosome_event_16\tentire_chromosome_event_17\tentire_chromosome_event_18\tentire_chromosome_event_19\tentire_chromosome_event_20\tentire_chromosome_event_21\tentire_chromosome_event_22\tentire_chromosome_event_X\tentire_chromosome_event_Y"

	outline_part6 = "genes_hit_by_CN_gain_events\tgenes_hit_by_CN_loss_events\tgenes_hit_by_CN_LoH_events"

	outline = outline_part1 + "\t" + outline_part2 + "\t" + outline_part3 + "\t" + outline_part4 + "\t" + outline_part5 + "\t" + outline_part6 + "\n"

	output_file.write( outline )

	outline_part1 = str(args.sampleid) + "\t" + str(scCNV_QC) + "\t" + str(algorithm_category) + "\t" + str(algorithm_CN_category) + "\t" + str(avg_adjusted_model_fit_f) + "\t" + str(avg_model_fit_f) + "\t" + str(num_CN_gain_events) + "\t" + str(total_size_in_nucleotides_CN_gain_events) + "\t" + str(avg_positive_depth_per_nucleotide_for_CN_gain_events) + "\t" + str(num_CN_loss_events) + "\t" + str(total_size_in_nucleotides_CN_loss_events) + "\t" + str(avg_negative_depth_per_nucleotide_for_CN_loss_events) + "\t" + str(num_CN_LoH_events) + "\t" + str(total_size_in_nucleotides_CN_LoH_events)

	outline_part2 = str(num_events_where_non_0p5_vaf_or_non_zero_depth) + "\t" + str(num_events_non_0p5_vaf) + "\t" + str(total_size_in_nucleotides_non_0p5_vaf) + "\t" + str(avg_size_in_nucleotides_non_0p5_vaf) + "\t" + str(avg_non_0p5_vaf_dist_from_0p5) + "\t" + str(num_events_same_vaf_diff_chrom) + "\t" + str(size_in_nucleotides_same_vaf_diff_chrom) + "\t" + str(num_diff_vaf1_values) + "\t" + str(avg_num_events_for_same_vaf1_values) + "\t" + str(num_events_non_zero_depth) + "\t" + str(total_size_in_nucleotides_non_zero_depth) + "\t" + str(avg_size_in_nucleotides_non_zero_depth) + "\t" + str(avg_non_zero_depth) + "\t" + str(avg_abs_value_non_zero_depth) + "\t" + str(num_diff_depth_values) + "\t" + str(avg_num_events_for_same_depth_values) + "\t" + str(num_non_normal_vaf_and_depth_both) + "\t" + str(num_non_normal_vaf_or_depth_but_not_both) + "\t" + str(num_telomere_non_normal_vaf_and_depth) + "\t" + str(num_non_telomere_non_normal_vaf_and_depth) + "\t" + str(has_default_fit_score) + "\t" + str(num_category_is_no_subclones) + "\t" + str(avg_adjusted_fit_f_category_is_no_subclones) + "\t" + str(num_category_is_germline_dup) + "\t" + str(avg_adjusted_fit_f_category_is_germline_dup) + "\t" + str(num_category_is_obvious_subclones) + "\t" + str(avg_adjusted_fit_f_category_is_obvious_subclones) + "\t" + str(num_category_is_telomere_subclones) + "\t" + str(avg_adjusted_fit_f_category_is_telomere_subclones) + "\t" + str(num_category_is_focal_subclones) + "\t" + str(avg_adjusted_fit_f_category_is_focal_subclones)

	outline_part3 = str(num_negative_depth_and_non_0p5_vaf_for_obvious_subclones) + "\t" + str(total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_obvious_subclones) + "\t" + str(num_positive_depth_and_non_0p5_vaf_obvious_subclones) + "\t" + str(total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_obvious_subclones) + "\t" + str(num_zero_depth_and_non_0p5_vaf_obvious_subclones) + "\t" + str(total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_obvious_subclones) + "\t" + str(num_negative_depth_and_non_0p5_vaf_for_focal_subclones) + "\t" + str(total_size_in_nucleotides_negative_depth_and_non_0p5_vaf_focal_subclones) + "\t" + str(num_positive_depth_and_non_0p5_vaf_focal_subclones) + "\t" + str(total_size_in_nucleotides_positive_depth_and_non_0p5_vaf_focal_subclones) + "\t" + str(num_zero_depth_and_non_0p5_vaf_focal_subclones) + "\t" + str(total_size_in_nucleotides_zero_depth_and_non_0p5_vaf_focal_subclones)

	outline_part4 = str(telomere_loss['1p']) + "\t" + str(telomere_loss['1q']) + "\t" + str(telomere_loss['2p']) + "\t" + str(telomere_loss['2q']) + "\t" + str(telomere_loss['3p']) + "\t" + str(telomere_loss['3q']) + "\t" + str(telomere_loss['4p']) + "\t" + str(telomere_loss['4q']) + "\t" + str(telomere_loss['5p']) + "\t" + str(telomere_loss['5q']) + "\t" + str(telomere_loss['6p']) + "\t" + str(telomere_loss['6q']) + "\t" + str(telomere_loss['7p']) + "\t" + str(telomere_loss['7q']) + "\t" + str(telomere_loss['8p']) + "\t" + str(telomere_loss['8q']) + "\t" + str(telomere_loss['9p']) + "\t" + str(telomere_loss['9q']) + "\t" + str(telomere_loss['10p']) + "\t" + str(telomere_loss['10q']) + "\t" + str(telomere_loss['11p']) + "\t" + str(telomere_loss['11q']) + "\t" + str(telomere_loss['12p']) + "\t" + str(telomere_loss['12q']) + "\t" + str(telomere_loss['13p']) + "\t" + str(telomere_loss['13q']) + "\t" + str(telomere_loss['14p']) + "\t" + str(telomere_loss['14q']) + "\t" + str(telomere_loss['15p']) + "\t" + str(telomere_loss['15q']) + "\t" + str(telomere_loss['16p']) + "\t" + str(telomere_loss['16q']) + "\t" + str(telomere_loss['17p']) + "\t" + str(telomere_loss['17q']) + "\t" + str(telomere_loss['18p']) + "\t" + str(telomere_loss['18q']) + "\t" + str(telomere_loss['19p']) + "\t" + str(telomere_loss['19q']) + "\t" + str(telomere_loss['20p']) + "\t" + str(telomere_loss['20q']) + "\t" + str(telomere_loss['21p']) + "\t" + str(telomere_loss['21q']) + "\t" + str(telomere_loss['22p']) + "\t" + str(telomere_loss['22q'])

	outline_part5 = str(entire_chromosome_event['1']) + "\t" + str(entire_chromosome_event['2']) + "\t" + str(entire_chromosome_event['3']) + "\t" + str(entire_chromosome_event['4']) + "\t" + str(entire_chromosome_event['5']) + "\t" + str(entire_chromosome_event['6']) + "\t" + str(entire_chromosome_event['7']) + "\t" + str(entire_chromosome_event['8']) + "\t" + str(entire_chromosome_event['9']) + "\t" + str(entire_chromosome_event['10']) + "\t" + str(entire_chromosome_event['11']) + "\t" + str(entire_chromosome_event['12']) + "\t" + str(entire_chromosome_event['13']) + "\t" + str(entire_chromosome_event['14']) + "\t" + str(entire_chromosome_event['15']) + "\t" + str(entire_chromosome_event['16']) + "\t" + str(entire_chromosome_event['17']) + "\t" + str(entire_chromosome_event['18']) + "\t" + str(entire_chromosome_event['19']) + "\t" + str(entire_chromosome_event['20']) + "\t" + str(entire_chromosome_event['21']) + "\t" + str(entire_chromosome_event['22']) + "\t" + str(entire_chromosome_event['X']) + "\t" + str(entire_chromosome_event['Y'])

	outline_part6 = str(genes_hit_by_CN_gain_events) + "\t" + str(genes_hit_by_CN_loss_events) + "\t" + str(genes_hit_by_CN_LoH)

	outline = outline_part1 + "\t" + outline_part2 + "\t" + outline_part3 + "\t" + outline_part4 + "\t" + outline_part5 + "\t" + outline_part6 + "\n"

	output_file.write( outline )
	output_file.close()

if __name__=='__main__':
    main()


