#!/usr/bin/env python3
import sys
import pandas as pd
import pysam
from array import array

class Fasta(object):
    def __init__(self,path:str):
        self.path = path
        seq={}
        with open(path,"r") as f:
            for line in f:
                if line.startswith(">"):
                    seqid = line.rstrip().replace(">","")
                    seq[seqid]=""
                else:
                    seq[seqid] = seq[seqid] + line.rstrip()
        self.seqDict = seq

    def labelSeq(self,m)->None:
        for k,v in self.seqDict.items():
            if m.get(k):
                self.seqDict[k + "|" + m[k]] = v
                del self.seqDict[k]

    def seqlen(self)->dict:
        seql={}
        for k,v in self.seqDict.items():
            seql[k] = len(v)
        return(seql)

class Gff(object):
    def __init__(self,path):
        self.path = path

class Bed(object):
    def __init__(self,path):
        self.path = path

class Sam(object):
    def __init__(self,path):
        self.path = path
    
    def insertCoverage(self,outpath):
        o = open(outpath+".cov","w+")
        f = pysam.AlignmentFile(self.path)
        # Read contig length
        clen = dict(zip(f.references,f.lengths))

        sl = []  # scaffold list
        clist = []  # coverage list
        prevScaffold = ""
        for read in f.fetch():
            # re-set the covergae dict if this is a new scaffold
            if read.reference_name != prevScaffold:
                print(read.reference_name)
                sl.append(read.reference_name)
                # print previous scaffold if not the head of file
                if prevScaffold != "":
                    idx = [*clen].index(read.reference_name)
                    for i in clist:
                        print(f"{idx}\t{i}",file=o)
                clist = array('L',[0]*clen[read.reference_name])
                prevScaffold = read.reference_name  # update variable
            # skip non unique mapped reads
            if read.has_tag("XA"):
                continue
            # add coverage to position
            epos = read.pos + read.tlen
            for p in range(read.pos,epos):
                clist[p] += 1
    
        # Add empty array for scaffold not present in SAM/BAM
        for k,v in clen.items():
            if k not in sl:
                clist = array('I',[0]*v)
                idx = [*clen].index(k)
                for i in clist:
                    print(f"{idx}\t{i}",file=o)
        f.close()
        o.close()
        omap = open(outpath+".scaffoldmap","w+")
        for i,v in enumerate([*clen]):
            print(f"{v}\t{i}",o=omap)
        print(f"Done. Coverage is written in {outpath}.cov and Scaffold map is written in {outpath}.scaffoldmap")

            

class Eggnog(object):
    def __init__(self,path):
        self.path = path
        row=[]
        with open(path,"r") as f:
            for line in f:
                if not line.startswith("#"):
                    row.append(line.rstrip().split("\t"))
        self.data = pd.DataFrame(row,columns=["query_name","seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score","predicted_gene_name","GO_terms","KEGG_KO","BiGG_Reactions","Annotation_tax_scope","Matching_OGs","best_OG|evalue|score","COG_functional_categories","eggNOG_HMM_model_annotation"])

    def seqOrthologMap(self)->dict:
        d={}
        for _,row in self.data.iterrows():
            if row["eggNOG_HMM_model_annotation"] is not None:
                v = row["seed_eggNOG_ortholog"] + "|" + row["eggNOG_HMM_model_annotation"]
            else:
                v = row["seed_eggNOG_ortholog"]
            d[row["query_name"]] = v
        return d

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
      
if __name__ == "__main__":
    sys.exit("Bioinformatics analysis libray")