def edit_postions(infile,seqs):
    with open(infile, "r") as inf:
        for line in inf:
            if re.match("\#",line) or not line.strip():
                continue
            (chrom,pos,nuc)=line.strip().split()[:3]
            # should be chrom, postition, nucleotide
            pos=int(pos)
            if not chrom in seqs or len(nuc) != 1:
                sys.exit(line.strip()+"contig not in sequences or nucleotide longer than 1")
            seqs[chrom].seq=seqs[chrom].seq[:pos-1] + nuc + seqs[chrom].seq[pos:]
    return seqs

import argparse
import sys
import re
from Bio import SeqIO
#from scipy import stats
#Author: Lukas Endler

parser = argparse.ArgumentParser(description='reads a fasta file with a single contig in stdin and a file with coordinates and nucleotides and replaces the nucs at the given postitions prints to STDOUT')
parser.add_argument("-p","--posfile", dest="posfile", help="file with positions and nucleotides", required=True)
parser.add_argument("-f","--fasta", dest="fasta", help="fastafile", default=False)
args = parser.parse_args()
posfile = vars(args)['posfile']
fasta = vars(args)['fasta']

if fasta:
    inf = open(fasta)
else:
    inf = sys.stdin

seqs=SeqIO.to_dict(SeqIO.parse(inf,"fasta"))
seqs=edit_postions(posfile,seqs)
for i in seqs:
    SeqIO.write(seqs[i],sys.stdout,"fasta")
