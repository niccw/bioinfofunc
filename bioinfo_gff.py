#!/usr/bin/env python3

import argparse
import pandas as pd

class Gff(object):
    def __init__(self,path):
        self.path = path
        self.gff = Gff.readGff(path)

    
    @staticmethod
    def readGff(gffpath):
        '''
        - gff3 use '#' as header
        - 9 columns:
            1 seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
            2 source - name of the program that generated this feature, or the data source (database or project name)
            3 feature - feature type name, e.g. Gene, Variation, Similarity
            4 start - Start position of the feature, with sequence numbering starting at 1.
            5 end - End position of the feature, with sequence numbering starting at 1.
            6 score - A floating point value.
            7 strand - defined as + (forward) or - (reverse).
            8 frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
            9 attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
        '''
        header = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
        col_list = []
        with open(gffpath,"r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    next(f)
                else:
                    col = line.strip().split("\t")
                    col_list.append(col)
            gff_df = pd.DataFrame(data=col_list,columns = header)
            return gff_df
    

    def reformat_gff(self,m,df,gene_only=True):
        print("##gff-version 3")

        for _,r in self.gff.iterrows():
            if gene_only:
                if r[2] == "mRNA" or r[2] == "CDS":
                    print(*r,sep="\t")
                    continue
            attr_dict = Gff._attribute_to_dict(r[8])
            # if ID is not in the rfdata, just print it out
            if "ID" not in attr_dict:
                print(*r,sep="\t")
                continue
            for j in m:
                if j[0] == "-":
                    # find the corresponding cell in the rfdata df and put it to the attribute d
                    if len(df.loc[df["ID"]==attr_dict["ID"],j[1]]) > 0:
                        attr_dict[j[1]] = df.loc[df["ID"]==attr_dict["ID"],j[1]].to_string(index=False)
                    else:
                        attr_dict[j[1]] = ""
                else:
                    # replace the old attribute value by the new value from rfdata df
                    if len(df.loc[df["ID"]==attr_dict["ID"],j[1]]) > 0:
                        attr_dict[j[0]] = df.loc[df["ID"]==attr_dict["ID"],j[1]].to_string(index=False)
                    else:
                        attr_dict[j[0]] = ""
            # update the attribute list using the dict
            attr = ""
            for k,v in attr_dict.items():
                if attr == "":
                    if len(v) == 0:
                        continue
                    attr = k + "=" + v 
                else:
                    if len(v) == 0:
                        continue
                    attr = attr + ";" + k + "=" + v
            # print the row
            r[8] = attr
            print(*r,sep="\t")
            

    def extract(self,target_list):
        with open(target_list,"r") as f:
            targets = f.read().splitlines()

        print("##gff-version3")
        # extract target from gff
        for _,r in self.gff.iterrows():
            attr_d = Gff._attribute_to_dict(r[8])
            if "ID" in attr_d and attr_d["ID"] in targets:
                print(*r,sep="\t")
            else:
                continue
        
    @staticmethod
    def _attribute_to_dict(a):
            d = {}
            for i in a.split(";"):
                if len(i) > 0:
                    k,v = i.split("=")
                    d[k] = v
            return d

    @staticmethod
    def _merge_row_by_ID(df):
        df = df.groupby("ID").agg(lambda x : ",".join(str(i) for i in x))
        return df

    def repeat_info(self,gene_gff):
        pass


if __name__ == "__main__":
    help_text = """
    Example:
    - Extract gff by ID:
    ./bioinfo_gff.py --extract sample.gff3 --id id.list > result.gff3

    - Reformat attribute using data from tsv based on ID:
    ./bioinfo_gff.py --reformat sample.gff3 --rfmap sample.txt --rfdata sample.tsv  > reformated.gff3

    ##### sample.txt #####
    # To add new attribute, put '-' in the first column (old attribute)
    -   lc2
    -   condition
    -   de
    old_attr   new_attr

    ##### sample.tsv #####
    ID  lc2 condition   de  new_attr
    123   1   0.5       up    new
    ...
    ######################

    """

    parser = argparse.ArgumentParser(description=help_text,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("gff", type = str, help = "GFF3 file.")
    parser.add_argument("--reformat", action="store_true", help = "Reformat the gff attribute. Work with --rfmap and --rfdata.")
    parser.add_argument("--extract", action="store_true", help = "Reformat the gff attribute. Work with --id")
    parser.add_argument("--rfmap", type = str, dest="rfmap",help = "2 columns tsv specific the rule of reformat. <attribite name in gff/- (for new attribute)> <target column name in rfdata>.")
    parser.add_argument("--rfdata", type = str, dest="rfdata",help = "tsv metadata of gff. 'ID' column is used to map to gff 'ID' attribute for row mapping. Make sure ID is unique!!")
    parser.add_argument("--id", type = str, dest= "id",help = "ID list (no header)")
    parser.add_argument("--geneonly", action="store_true", help = "Only process gene, ignore and print out mRNA and CDS and other rows")

    args = parser.parse_args()

    if args.reformat:
        # read rfmap to list
        map_ls = []
        with open(args.rfmap,"r") as f:
            for line in f:
                col = line.strip().split("\t")
                map_ls.append(col)

        # read rfdata to dataframe
        rfdata = pd.read_csv(args.rfdata,sep="\t")
        # make sure ID is unique
        if not len(set(rfdata["ID"])) == len(rfdata["ID"]):
            rfdata = Gff._merge_row_by_ID(rfdata)
            rfdata["ID"] =  rfdata.index
            rfdata.reset_index(drop=True,inplace=True)
        # read gff
        gff = Gff(args.gff)

        gff.reformat_gff(map_ls,rfdata,args.geneonly)

    if args.extract:
        gff = Gff(args.gff)
        gff.extract(args.id)






     

    


