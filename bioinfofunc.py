#!/usr/bin/env python3
import sys

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
      
if __name__ == "__main__":
    sys.exit("Bioinformatics analysis libray")