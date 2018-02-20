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

parser = argparse.ArgumentParser(description="""takes a VCF file and gives the reverse complement (for SNPgenie) to STDOUT   
""")

parser.add_argument("--in","-i", dest="vcffile", 
                    help="vcf-file, tries stdin if not set use \"STDIN\" or nothing for piping; default \"False\"", default=False)
parser.add_argument("--chroms","-c", dest="chromlen",  type=str,
                    help="dict of chrom names to length; default \'{\"L\":7229,\"S\":3377\}\'", default='{"L":7229,"S":3377}')

args = parser.parse_args()
# get chromlength
chromlen=json.loads(vars(args)['chromlen'])
vcf_file = vars(args)['vcffile']
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
    entries = line.split("\t")
    chrom = entries[0]
    if not ("N" in entries[4]):
        print line
        continue



