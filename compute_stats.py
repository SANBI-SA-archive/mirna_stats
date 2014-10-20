#!/usr/bin/env python

import argparse
import sys
import re
parser = argparse.ArgumentParser(description='Parser TP/FP stats')
parser.add_argument('--max_energy', type=float, default=-10.0)
parser.add_argument('--min_score', type=float, default=145.0)
#parser.add_argument('--max_score', type=float, default=175.0)
parser.add_argument('known_targets', type=argparse.FileType(), help='Annotated known targets')
parser.add_argument('miranda_output', type=argparse.FileType(), help='Miranda predictions')
parser.add_argument('output_file', type=argparse.FileType('w'), nargs='?', default=sys.stdout, help='Output file')

args = parser.parse_args()

# format of known targets
# hsa-miR-34a-5p corresponds with PLK1 |ENSG00000166851,ENST00000562272,PLK1,HGNC:9077
# hsa-miR-124-3p corresponds with CDK2 |ENSG00000123374,ENST00000266970,CDK2,HGNC:1771

known_targets = dict()
target_re = re.compile(r'^(?P<mirna>\S+) corresponds with (?P<gene_symbol>\S+) \|(?P<ensembl_gene>[^,]+),(?P<ensembl_transcript>[^,]+),\2,(?P<hugo_name>HGNC:\d+)')
for line in args.known_targets:
    match = target_re.match(line) 
    if match:
#        print "got a match:", match.group('mirna')
	mirna_name = match.group('mirna')
	target_hits = known_targets.get(mirna_name, set())
	target_hits.add(match.group('ensembl_gene'))
	known_targets[mirna_name] = target_hits
    else:
        # didn't get a match? let's complain but skip
	sys.stderr.write("Found line that doesn't match: {}".format(line))
#        print "did not match"

# sample miranda output
# mirna_id 	 utr_id 	 score 	 energy 	 mirna_hit_start 	 mirna_hit_end 	 utr_hit_start 	 utr_hit_end 	 aln_length 	 identity 	 similarity 	 aln_mirna 	 aln_map 	 aln_utr 	 seed_mer 
#hsa-miR-16-5p	ENSG00000001617	123.000000	-9.170000	2	21	171	196	23	65.217391	65.217391	GCGGUUAUAAAU----GCACGACGAU	  |||| | |||    | ||||| | 	CTCCAAAACTTAGACCCATGCTGGTC	no_seed_mer	
predicted_targets = dict()
predicted_nontargets = dict()
miranda_iter = iter(args.miranda_output)
for line in miranda_iter:
    if line.startswith('mirna_id'):
        line = miranda_iter.next()
    fields = line.split()
    score = float(fields[2])
    energy = float(fields[3])
    mirna_name = fields[0]
    gene_name = fields[1]
    if score < args.min_score or energy > args.max_energy:
        target_misses = predicted_nontargets.get(mirna_name, set())
	target_misses.add(gene_name)
	predicted_nontargets[mirna_name] = target_misses
    else: 
	# ok, we've passed the threshold for a hit
        target_hits = predicted_targets.get(mirna_name, set())
        target_hits.add(gene_name)
        predicted_targets[mirna_name] = target_hits

#for mirna_name in predicted_targets:
#    for target_hit in predicted_targets[mirna_name]:
#        print "MATCH:", mirna_name, target_hit

TP = FP = FN = TN = 0
for mirna_name in known_targets:
    known_targets_set = known_targets[mirna_name]
    predicted_nontarget_set = predicted_nontargets.get(mirna_name, set())
    predicted_targets_set = predicted_targets.get(mirna_name, set())
    TP_set = known_targets_set.intersection(predicted_targets_set)
    FP_set = predicted_targets_set.difference(known_targets_set)
    FN_set = known_targets_set.difference(predicted_targets_set)
    TN_set = predicted_nontarget_set.difference(predicted_targets_set.union(known_targets_set))
    TP += len(TP_set)
    FP += len(FP_set)
    FN += len(FN_set)
    TN += len(TN_set)
    print "for mirna: {} TP: {} FP: {} TN: {} FN: {}".format(mirna_name, len(TP_set), len(FP_set), len(TN_set), len(FN_set))

sensitivity = float(TP) / (TP + FN)
specificity = float(TN) / (TN + FP)
print "Total TP: {} FP: {} TN: {} FN: {} sensitivity: {} specificity: {}".format(TP, FP, TN, FN, sensitivity, specificity)
