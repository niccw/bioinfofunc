#!/usr/bin/env python3

import argparse
import sys
from re import finditer

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

    def labelSeq(self, m:dict)->None:
        # m : {'id':'label'}
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

    def getseq(self,query):
        # query is a list, len >= 1
        for k,v in self.seqDict.items():
            if k in query:
                print(f">{k}")
                print(v)

    def gap2bed(self):
        for k,v in self.seqDict.items():
            for m in finditer("N+", str(v)):
                print(*[k,m.span()[0],m.span()[1]], sep="\t")

if __name__ == "__main__":

    help_text= """
    Example:
    - get sequence
    ./bioinfo_fasta.py --getseq sample.fa --id Hm105_s01
    ./bioinfo_fasta.py --getseq sample.fa --id id_list.txt

    - sequence length
    ./bioinfo_fasta.py --seqlen sample.fa

    - Find gap region (Ns) and output as BED
    ./bioinfo_fasta.py --gap2bed sample.fa
    """


    parser = argparse.ArgumentParser(description=help_text,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fa", type = str, default="Fasta file.")
    parser.add_argument("--seqlen", action="store_true", help="Calculate seq length.")
    parser.add_argument("--label", action="store_true", help="Add label to ID. Work with label dict")
    parser.add_argument("--labeldict", type=str, help="2 columns tsv containing the label for corresponding seq ID. i.e. <seq1> <label>")
    parser.add_argument("--getseq", action="store_true", help="Get sequence. Work with --id.")
    parser.add_argument("--id", type=str,dest="id_query", help="ID name/ ID name list (end with .txt)")
    parser.add_argument("--gap2bed", action="store_true", help="Output gap(N) as bed")


    args = parser.parse_args()
    
    if args.seqlen:
        fa = Fasta(args.fa)
        fa_len_dict = fa.seqlen()
        for k,v in fa_len_dict.items():
            print(f"{k}\t{v}")
    
    if args.getseq:
        # check --id is one ID or a ID list
        if args.id_query:
            id_list = []
            if args.id_query.endswith("txt"):
                with open(args.id_query,"r") as f:
                    id_list = f.read().splitlines()
            else:
                id_list.append(args.id_query)
        else:
            sys.exit("Please provide id / id.txt using --id")
        
        fa = Fasta(args.fa)
        fa.getseq(query=id_list)

    if args.gap2bed:
        fa = Fasta(args.fa)
        fa.gap2bed()

    

