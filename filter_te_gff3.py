#!/usr/bin/env python3

import argparse
import sys


def read_faidx(faidx_path:str)->dict:
    d = {}
    with open(faidx_path,"r") as f:
        for line in f:
            col = line.strip().split("\t")
            d[col[0]] = int(col[1])
    return d

def filter_by_length(gff_path:str,faidx_d:dict,p:float):
    
    with open(gff_path,"r") as f:
        for i,line in enumerate(f):
            if i == 0:
                print(line.strip(), file=sys.stdout)
                continue
            col = line.strip().split("\t")
            te_len = int(col[4]) - int(col[3]) # length of TE copy
            assert te_len >= 0, f"te length is negative, check ${line}"
            if float(te_len/faidx_d[col[0]]) >= p:
                print(f"te_len:{te_len} transcript_len:{faidx_d[col[0]]}", file=sys.stderr)
                print(line.strip(), file=sys.stdout)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-idx","--faidx",type=str, dest="faidx", help="fasta faidx")
    parser.add_argument("-g","--gff",type=str,dest="gff", help="Repeatmasker reformated gff3")
    parser.add_argument("-p","--percent",type=float,dest="percent", help="min. percent of the transcript to keep the TE annotation (0-1)")
    parser.add_argument("-l","--length", action="store_true", help="filter gff3 by transcript length. Work together with -p/--percentage")
    args = parser.parse_args()

    if args.length:
        faidx_d = read_faidx(faidx_path=args.faidx)
        filter_by_length(gff_path=args.gff,faidx_d=faidx_d,p=args.percent)
