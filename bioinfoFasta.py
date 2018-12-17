#!/usr/bin/env python3

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
            for j in range(1,len(a)-1):
                d[a[0]].append(j)
        return d