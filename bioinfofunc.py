#!/usr/bin/env python3
import sys


def seqLen(fa:str)->dict:
    seql = {}
    with open(fa,"r") as f:
        for line in f:
            if line.startswith(">"):
                seqid = line.rstrip().replace(">","")
                seql[seqid] = 0
            else:
                seql[seqid] = seql[seqid] + len(line.rstrip())
    return seql

if __name__ == "__main__":
    sys.exit("Commonly used function for bioinformatics analysis.")