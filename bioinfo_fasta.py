#!/usr/bin/env python3

import argparse
import sys

class Fasta(object):
    def __init__(self, path: str):
        self.path = path
        seq = {}
        with open(path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    seqid = line.rstrip().replace(">", "")
                    seq[seqid]=""
                else:
                    seq[seqid] = seq[seqid] + line.rstrip()
        self.seqDict = seq

    def labelSeq(self, m)->None:
        for k, v in self.seqDict.items():
            if m.get(k):
                self.seqDict[k + "|" + m[k]] = v
                del self.seqDict[k]

    def seqlen(self)->dict:
        seql = {}
        for k, v in self.seqDict.items():
            seql[k] = len(v)
        return seql

    def parseIdentifier(self)->dict:
        d = {}
        seqName = [*self.seqDict]
        for i in seqName:
            a = i.split("|")
            d[a[0]] = []
            for j in range(1,len(a)):
                d[a[0]].append(a[j])
        return d

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fa", type = str, default="Fasta file.")
    parser.add_argument("--seqlen", action="store_true", help="Calculate seq length.")
    parser.add_argument("--label", action="store_true", help="Add label to ID. Work with label dict")
    parser.add_argument("--labeldict", type=str, help="2 columns tsv containing the label for corresponding seq ID. i.e. <seq1> <label>")

    args = parser.parse_args()

    fa = Fasta(args.fa)
    if args.seqlen:
        fa_len_dict = fa.seqlen()
        for k,v in fa_len_dict.items():
            print(f"{k}\t{v}")

