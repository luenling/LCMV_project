#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 13:47:34 2017

@author: lukasendler
"""

import vcf
import sys, re
import os 
import argparse
from scipy import stats
import select
import gzip
import numpy as np
from collections import defaultdict,OrderedDict,namedtuple

parser = argparse.ArgumentParser(description="""Go through multiple vcf files and try to combine them, change the frequency array to each sample
""")

parser.add_argument("--in","-i", dest="vcffile", help="vcf-file", required=True)
parser.add_argument("--minfreq","-m", dest="minfreq",type=float, help="minimal minor allele frequency to consider (if not reached in at least one sample, allele will be removed)", default=0.01)
parser.add_argument("--sb", dest="minsb",type=float, help="minmal ratio of supporting strands", default=0.1)
parser.add_argument("--fsb", dest="maxfsb",type=float, help="maximal FSB", default=30)
parser.add_argument("--rpp", dest="maxrpp",type=float, help="maximal RPP", default=30)
parser.add_argument("--mc", dest="mc",type=float, help="min. supporting reads", default=3)
parser.add_argument("-v", dest="verb", action="store_true", help="verbose (default: FALSE)", default=False)
args = parser.parse_args()
vcf_file = vars(args)['vcffile']
minfreq=vars(args)['minfreq']
minsb=vars(args)['minsb']
maxfsb=vars(args)['maxfsb']
maxrpp=vars(args)['maxrpp']
mc=vars(args)['mc']
verb=vars(args)['verb']
#maxfreq=1-minfreq
# open vcf file
if vcf_file == "STDIN":
    inf = sys.stdin
    vcf_reader = vcf.Reader(inf)
else:
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    # if AF field not in formats, add it
if not ('AF' in vcf_reader.formats.keys()):
    newForm=vcf_reader.formats[vcf_reader.formats.keys()[0]]
    newForm=newForm._replace(id='AF', num=-1, type='Float', desc='Alternative allele frequency')
    vcf_reader.formats['AF'] = newForm
if not ('FSB' in vcf_reader.infos.keys()):
    newForm=vcf_reader.infos[vcf_reader.infos.keys()[0]]
    newForm=newForm._replace(id='FSB', num=-1, type='Float', desc='PHRED scaled Chi2 test ( or if expected < 50 > 5 with Yates correction or if minimum expected <=5 Fisher test) value for strand bias')
    vcf_reader.infos['FSB'] = newForm
if not ('SB' in vcf_reader.infos.keys()):
    newForm=vcf_reader.infos[vcf_reader.infos.keys()[0]]
    newForm=newForm._replace(id='SB', num=-1, type='Float', desc='ratio of less common strand to more common strand of reads supporting alternative alleles')
    vcf_reader.infos['SB'] = newForm

