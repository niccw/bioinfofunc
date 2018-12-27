#!/usr/bin/env python3

import pandas as pd

class Bed(object):
    def __init__(self,path,annot=False):
        self.path = path
        self.cluster_dict = None
        self.cluster_param = None

        # read bed to multi-array
        bed_chr = []
        bed_start = []
        bed_end = []
        if annot:
            bed_annot = []
            bed_col = 4
        else:
            bed_col = 3
        
        with open(path,"r") as f:
            for line in f:
                if not line.startswith("#"):
                    row = line.strip().split("\t")[:bed_col]
                    bed_chr.append(row[0])
                    bed_start.append(int(row[1]))
                    bed_end.append(int(row[2]))
                    if annot:
                        bed_annot.append(row[3])
        
        if annot:
            bed_df = pd.DataFrame(data= {
                "chr": bed_chr,
                "start": bed_start,
                "end": bed_end,
                "annot": bed_annot
            } 
            )
        else:
            bed_df = pd.DataFrame(data= {
                "chr": bed_chr,
                "start": bed_start,
                "end": bed_end
            }
            )
        # sort the bed
        self.bed = bed_df.sort_values(by=["chr", "start","end"])
    
    def cluster(self,span=2000):
        # gene cluster for each chr in a dict
        d = {}
        for c in self.bed.chr:
            d[c] = []

        # dict to mark down the features to be compared
        last_row = {
            "chr":"",
            "start":0,
            "end":0,
            "group":None # 0 based
        }

        # check gene one by one
        for _,row in self.bed.iterrows():
            # check if this is new chr
            if row["chr"] != last_row["chr"]:
                last_row["group"] = None
            else:
                # same chr, check the distance with the last gene
                distance = row["start"] - last_row["end"]
                if distance <= span: # group them
                    # if last row is already in a group
                    if last_row["group"] is not None:
                        d[row["chr"]][last_row["group"]].append([row["start"],row["end"]])
                        last_row["group"] = last_row["group"]
                    else:
                        # start a new group
                        d[row["chr"]].append([]) # start a new cluster list in this chr
                        group_idx = len(d[row["chr"]]) - 1 # 0 based
                        # add last row and current row to group
                        d[row["chr"]][group_idx].append([last_row["start"],last_row["end"]])
                        d[row["chr"]][group_idx].append([row["start"],row["end"]])
                        # update group label 
                        last_row["group"] = group_idx
                else: # row not group with last row
                    last_row["group"] = None

            last_row["chr"] = row["chr"]
            last_row["start"] = row["start"]
            last_row["end"] = row["end"]
        self.cluster_dict = d
        self.cluster_param = "span: {}".format(span)
    
    @staticmethod
    def print_cluster(cluster_dict):
        for chrom in [*cluster_dict]:
            #enumerate
            for idx,cluster in enumerate(cluster_dict[chrom]):
                if cluster:
                    features_list = []
                    for feature in cluster:
                        try:
                            features_list.append(":".join([str(x) for x in feature]))
                        except:
                            print(feature)
                    features_str = ",".join(features_list)
                    print(f"{chrom}\t{idx}\t{features_str}\t{len(features_list)}")