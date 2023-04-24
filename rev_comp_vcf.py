#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:21:30 2018

@author: lukasendler
"""

import sys, re
import os 
#import numpy as np
#from scipy import stats
import argparse
import gzip
import json
import string

parser = argparse.ArgumentParser(description="""takes a VCF file and gives the reverse complement (for SNPgenie) to STDOUT   
""")

parser.add_argument("--in","-i", dest="vcffile", 
                    help="vcf-file, tries stdin if not set use \"STDIN\" or nothing for piping; default \"False\"", default=False)
parser.add_argument("--snpgenie","-s", dest="snpgenie",
                    help="input file is snpgenie result and give chromosome name default \"False\"", default=False)
parser.add_argument("--chroms","-c", dest="chromlen",  type=str,
                    help="dict of chrom names to length; default \'{\"L\":7229,\"S\":3377\}\'", default='{"L":7229,"S":3377}')

args = parser.parse_args()
# get chromlength
chromlen=json.loads(vars(args)['chromlen'])
vcf_file = vars(args)['vcffile']
snpgenie = vars(args)['snpgenie']
revtab=string.maketrans("ACTG","TGAC")
rev_cols= False
# open vcf file
if not vcf_file or vcf_file == "STDIN":
    if not sys.stdin.isatty():
        inf = sys.stdin
    else:
        sys.exit("No vcf file or stdinput given")
else:
    if re.search("\.b?gz",vcf_file):
        inf = gzip.open(vcf_file,'rb')
    else:
        inf = open(vcf_file,"r")

for line in inf:
    line = line.rstrip()
    if (re.match("^\s*\#+",line)): # entry is comment/header
        print line
        continue
    if (re.match("^\s*file",line) and snpgenie != False): # entry is snpgenie header
        entries = line.split("\t")
        if ((entries[-4],entries[-3],entries[-2],entries[-1]) == ("A","C","G","T")):
            rev_cols = True
        print line
        continue
    entries = line.split("\t")
    if (snpgenie == False):
        chrom = entries[0]
        site = 1
    else:
        chrom = snpgenie
        site = 2
        # reverse the acgt counts - last 4 colums normally ACGT checked above, if not do nothing
        if (rev_cols):
            (entries[-4],entries[-3],entries[-2],entries[-1]) = (entries[-1],entries[-2],entries[-3],entries[-4])
    entries[site] = str(chromlen[chrom] - int(entries[site]) + 1)
    if not (snpgenie != False and len(entries[3]) == 3 ):
        # do not do for codons in snpgenie
        entries[3] = entries[3].translate(revtab)[::-1]
    alts = entries[4].split(",")
    entries[4] = ",".join([ x.translate(revtab)[::-1] for x in alts])
    print "\t".join(entries)

    
    



