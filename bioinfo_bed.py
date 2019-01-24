#!/usr/bin/env python3
import argparse
import pandas as pd

class Bed(object):
    def __init__(self,path,annot=False, header=False):
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
            if header:
                next(f)
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

    @staticmethod
    def read_bed_closest(bed, attr_col:int)->dict:
        c = attr_col - 1 # 1-base to 0-base
        d = {}
        with open(bed,"r") as f:
            for line in f:
                col = line.strip().split("\t")
                k = "|".join([col[0], col[1], col[2]])
                if k in d:
                    d[k] = d[k] + "," + col[c]
                else:
                    d[k] = col[c]
        return d

    @staticmethod
    def label_bed(bed,label_dict,a="annot",header=False):
        with open(bed,"r") as f:
            if header:
                h = f.readline().strip().split("\t")
                h.append(a)
                print(*h,sep="\t")
            for line in f:
                col = line.strip().split("\t")
                k = "|".join([col[0], col[1], col[2]])
                if k in label_dict:
                    col.append(label_dict[k])
                    print(*col, sep="\t")

    def collapse(self):

        last_row = {
            "chr":"",
            "start":"",
            "end":"",
            "annot":""
        }

        for _,row in self.bed.iterrows():
            if row["chr"] == last_row["chr"] and row["start"] == last_row["start"] and row["end"] == last_row["end"]:
                last_row["annot"] =  last_row["annot"] + "," + str(row["annot"])
            else:
                if last_row["chr"] != "":
                    print(f'{last_row["chr"]}\t{last_row["start"]}\t{last_row["end"]}\t{last_row["annot"]}')
                last_row["chr"] = row["chr"]
                last_row["start"] = row["start"]
                last_row["end"] = row["end"]
                last_row["annot"] = row["annot"]
        # print last row
        print(f'{last_row["chr"]}\t{last_row["start"]}\t{last_row["end"]}\t{last_row["annot"]}')

            
if __name__ == "__main__":
    help_text = """
    Example:
    - Annotate bed by column in another bed (dbed)
    ./bioinfo_bed.py --bbannot main.bed --dictbed dict.bed --acol 4
    - Collapse bed (by chr,start,end), desciprtion/annotation are concatanated by ','
    ./bioinfo_bed.py --collapse sample.bed


    """
    parser = argparse.ArgumentParser(description=help_text,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bed", help="bed file")
    parser.add_argument("--header", action="store_true", help="Default = False")
    # Cluster
    parser.add_argument("--cluster", action="store_true", help="Cluster featuresby distance. Work with -d/--distance")
    parser.add_argument("-d,--distance", dest="d", type = int, default = 2000, help="Distance (in bp) to cluster features. Default = 2000")
    # Special annotatation
    parser.add_argument("--bbannot", action="store_true", help="Annotate bed based on another bed")
    parser.add_argument("--dictbed", dest="dbed", type = str, help="BED provide information of annotation. Work with --acol. This bed has no header.")
    parser.add_argument("--acol", dest="acol", type = int, help="Column containing the annotation info (1-base)")
    parser.add_argument("--annot", dest="annot", type = str, help="New column name")
    # Collapse
    parser.add_argument("--collapse", action="store_true", help="Collapse bed, concat. annot by ','")


    args = parser.parse_args()
    
    # Cluster
    if args.cluster:
        bed = Bed(args.bed,header=args.header)
        bed.cluster(span=args.d)
        Bed.print_cluster(bed.cluster_dict)
    
    # Special annotatation
    if args.bbannot:
        d = Bed.read_bed_closest(bed=args.dbed,attr_col=args.acol)
        Bed.label_bed(bed=args.bed,label_dict=d, header=args.header, a=args.annot)

    # Collapse
    if args.collapse:
        bed = Bed(args.bed,annot=True, header=args.header)
        bed.collapse()
