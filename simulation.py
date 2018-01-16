#!/usr/bin/env python
import numpy as np, numpy.random
import random

from utils import *



def simulate(patterns, rep_size, is_DMR):
	normal_rep_pat = make_normal_replicate(patterns, rep_size)
	thr = 0.15
	DMR_change_param = 0.9
	replicate_change_param = 0.03


	if is_DMR == 1:
	#DMR region
		DMR_pat = make_DMR_patterns(patterns, thr,  DMR_change_param, replicate_change_param)
		tumor_rep_pat = make_normal_replicate(DMR_pat, rep_size)

	if is_DMR == 0: 
	#non_DMR region
		tumor_rep_pat = make_normal_replicate(patterns, rep_size)
	
	prune_normal_rep_pat = []
	for pats in normal_rep_pat:
		prune_pats = prune_patterns(pats)
		prune_normal_rep_pat.append(prune_pats)

	prune_tumor_rep_pat = []
	for pats in tumor_rep_pat:
		prune_pats = prune_patterns(pats)
		prune_tumor_rep_pat.append(prune_pats)


	return prune_normal_rep_pat, prune_tumor_rep_pat


def prune_patterns(patterns):
	result = []
	for pat in patterns:
		is_new = 1
		for prev_pat in result:
			if pat == prev_pat:
				is_new = 0
				prev_pat.abundance += pat.abundance
				break
		if is_new == 1:
			result.append(pat)
	return result

def make_one_normal_replicate(pat, replicate_array):
	rand = random.random()
	if rand < 0.5:
		replicate_patterns = make_replicate(pat, 0.02,  0.03, 1)
		for replicate_pattern in replicate_patterns:
			replicate_array.append(replicate_pattern)
	else:
		replicate_patterns = make_replicate(pat, 0.02,  0.03, 0.8)
		for replicate_pattern in replicate_patterns:
			replicate_array.append(replicate_pattern)
	return replicate_array



def make_normal_replicate(normal_pat, rep_size):
	normal_patterns = []

	for i in range(rep_size):
		one_replicate = []
		for pat in normal_pat:
			one_replicate = make_one_normal_replicate(pat, one_replicate)

		normal_patterns.append(one_replicate)

	return normal_patterns




def make_replicate(pat, old_rep_param,  new_rep_param, old_ratio):
#old_rep_param = probality of change in meth status for original pattern 
#new_rep_param = probability of change in meth status for replicate patterns
#old_ratio = ratio of abundance of original pattern in making new patterns(replicate)

	replicate_patterns=[]
	rep_param = [old_rep_param, new_rep_param]
	ratio_param = [old_ratio, 1.0 - old_ratio]

	meth = pat.mPat
	#replicate 1(old)
	if str(meth) == "*":
		  return pat

	tokens = str(meth).split(",")

	for rep_rand, ratio in zip(rep_param, ratio_param):
		new_pat=[]
		for token in tokens:
			element = token.split(":")
			methyl = element[1]

			#rand = random.randint(0,101)
			rand = random.random()

			if rand < rep_rand:
				if methyl == 'M':
					new_pat.append(element[0]+':U')
				if methyl == 'U':
					new_pat.append(element[0]+':M')
			else:
				new_pat.append(token)
		if ratio > 0:
			rep_meth = methylPattern(pat.chr, pat.start, pat.end, pat.cid, pat.pid, float(float(pat.abundance)*ratio), ','.join(new_pat), pat.regions, pat.group, pat.patId)
			#rep_meth.setId(0)
			replicate_patterns.append(rep_meth)

	return replicate_patterns


def choose_high_abundant_pattern(normal_pat, thr):
	#high_abundant_patt =[]
	#high_abundant_abd =[]

	abd_sum = 0
	for pat in normal_pat:
		abd_sum += pat.abundance

	abundance_tag = []
	count = 0
	for pat in normal_pat:
		if float(pat.abundance)/abd_sum > thr:
		#	high_abundant_patt.append(pat)
		#	high_abundant_abd.append(abd)
			abundance_tag.append('high')
			count += 1 
		else:
			abundance_tag.append('low')

	if count == 1:
	#case Simple
		#return 0, high_abundant_patt, high_abundant_abd, abundance_tag
		return 'simple', abundance_tag

	if count > 1 and count < 5:
	#case Moderate
		#return 1, high_abundant_patt, high_abundant_abd, abundance_tag
		return 'moderate', abundance_tag
	if count == 0 or count >= 5:
	#case Hard
		#return 2, high_abundant_patt, high_abundant_abd, abundance_tag
		return 'hard', abundance_tag




def make_DMR_patterns(normal_pat, thr,  DMR_change_param, replicate_change_param):
	#choose category
	DMR_pat = []
	DMR_abd = []

	#case, high_abundant_patt, high_abundant_abd, abundance_tag = choose_high_abundant_pattern(normal_pat, abd_pat, thr)
	
	case, abundance_tag = choose_high_abundant_pattern(normal_pat, thr)

	high_abundant_count = 0
	for tag in abundance_tag:
		if tag == 'high':
			high_abundant_count += 1

	max_change_count = high_abundant_count / 2
	change_count = 0
	remaining_change_count = high_abundant_count

	for pat, tag in zip(normal_pat, abundance_tag):
		if tag == 'low':
			DMR_pat = make_one_normal_replicate(pat, DMR_pat)

		else: # tag == 'high'
			if case == 'simple':
			#case Simple
				rand = random.random()
				if rand < float(1.0)/3:
					new_patterns = make_replicate(pat, replicate_change_param, DMR_change_param, 0)
					for pattern in new_patterns:
						DMR_pat.append(pattern)
				if rand >= float(1.0)/3 and rand < float(2.0)/3:	
					new_patterns = make_replicate(pat, replicate_change_param, DMR_change_param, float(1.0)/3)
					for pattern in new_patterns:
						DMR_pat.append(pattern)
				if rand >= float(2.0)/3:
					new_patterns = make_replicate(pat, replicate_change_param, DMR_change_param, float(2.0)/3)
					for pattern in new_patterns:
						DMR_pat.append(pattern)

			if case == 'moderate' or case == 'hard':
				if random.random() < float(max_change_count - change_count)/remaining_change_count:
					change_count += 1
					remaining_change_count -= 1
					rand = random.random()
					if rand < 0.5:
						new_patterns = make_replicate(pat, replicate_change_param, DMR_change_param, 0)	
						for pattern in new_patterns:
							DMR_pat.append(pattern)
					if rand >= 0.5:
						new_patterns = make_replicate(pat, replicate_change_param, DMR_change_param, float(1.0)/3)
						for pattern in new_patterns:
							DMR_pat.append(pattern)
				else:
					remaining_change_count -= 1
					DMR_pat = make_one_normal_replicate(pat, DMR_pat)
	return DMR_pat









