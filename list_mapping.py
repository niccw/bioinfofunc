#!/usr/bin/env python3
import argparse

def list_mapping(l,m,s=",",p=1):
    # read map to dict
    d = {}
    with open(m,"r") as f:
        for line in f:
            col = line.strip().split("\t")
            if p == 1:
                d[col[0]] = col[1]
            else:
                d[col[1]] = col[0]
    
    with open(l,"r") as g:
        for line in g:
            mapped = []
            item = line.strip().split(s)
            for k in item:
                match = d[k] if k in d else "NA"
                mapped.append(match)
            print(*mapped,sep=s)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("l", type=str, help="File with each row containing a list to map")
    parser.add_argument("m", type=str, help="2 column tsv for mapping")
    parser.add_argument("-s,--sep", dest= "sep",default=",",type=str, help="List seperator. Default = ','")
    parser.add_argument("-p,--pos", dest= "pos",default=1,type=int, help="Index of column containing the key (i.e. item in the list) (1 base). Default = 1. ")
    args = parser.parse_args()

    list_mapping(l=args.l, m=args.m, s=args.sep, p=args.pos)