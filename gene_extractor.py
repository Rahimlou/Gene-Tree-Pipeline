#! /usr/bin/env python
import argparse
import subprocess
import os
from Bio import SeqIO
def gene_extractor(seqfile, in_path, out_path):
    asmbl = os.listdir(in_path)
    sequences=[i for i in SeqIO.parse(seqfile, 'fasta')]
    if len(sequences) == 1:
        for line in asmbl:
            filename = line
            basename = os.path.splitext(line)[0]
            outfile_name = basename + ".stk"
            in_file = os.path.join(in_path, filename)
            out_file = os.path.join(out_path, outfile_name)
            subprocess.call(["nhmmer", "-A", out_file, seqfile, in_file])
    if len(sequences) > 1:
        for i in sequences:
            sqname=i.name
            in_seq=sqname + ".fasta"
            SeqIO.write(i, in_seq, 'fasta')
            os.mkdir(sqname)
            for line in asmbl:
                filename=line
                in_file=os.path.join(in_path, filename)
                basename = os.path.splitext(line)[0]
                outfile_name = basename + ".stk"
                out_path=os.path.join(sqname, outfile_name)
                subprocess.call(["nhmmer", "-A", out_path, in_seq , in_file])

       
def main():
    parser=argparse.ArgumentParser(description="Making database out of a genome assembly")
    parser.add_argument("-seqfile", help="Query sequence file directory.", type=str, required=True)
    parser.add_argument("-in_path", help="Assembly files directory.", type=str, required=True)
    parser.add_argument("-out_path", help="Output files directory.", type=str, required=True)
    args=parser.parse_args()
        
    gene_extractor(args.seqfile, args.in_path, args.out_path)

if __name__=="__main__":
	main()

