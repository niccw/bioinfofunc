#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

class Rm_out(object):
    def __init__(self,path):
        self.path = path
        self.out = Rm_out.out2df(self.path)
    

    @staticmethod
    def out2df(out):
        header = ["sw_score", "perc_div", "perc_del", "perc_ins", "query_seq", "query_pos_begin", "query_pos_end", "query_left", "strand", "matching_repeat", "repeat_class", "repeat_pos_begin", "repeat_pos_end", "repeat_pos_left", "id", "asterisk"]
        col_list = []
        with open(out,"r") as f:
            # first 3 lines are header, without #
            for i in range(3):
                next(f)
            for line in f:
                 cols = line.rstrip().split()
                 if len(cols) == 15:
                     cols.append("")
                 col_list.append(cols)
            
        bed_df = pd.DataFrame(data=col_list,columns = header)
        return bed_df

    def repeat_class_dict(self):
        d = {}
        for _, r in self.out.iterrows():
            d[r["matching_repeat"]] = r["repeat_class"]
        return d

    def out2gff3(self):
        print("## gff-version 3")
        for _, r in self.out.iterrows():
            if r["strand"] == "C":
                s = "-"
                tstart = r["repeat_pos_end"]
                tend = r["repeat_pos_begin"].replace("(","")
                tend = r["repeat_pos_begin"].replace(")","")
            else:
                s = "+"
                tstart = r["repeat_pos_begin"]
                tend = r["repeat_pos_end"]

            note = "Tstart=" + tstart + ";" + "Tend=" + tend + ";" + "ID=" +  r["matching_repeat"]
            print_col = [r["query_seq"], "RepeatMasker", r["repeat_class"], r["query_pos_begin"], r["query_pos_end"], r["sw_score"], s, ".", note]
            print(*print_col, sep = "\t")
    
    @staticmethod
    def out2_raw_gff3(out_path):
        header = ["sw_score", "perc_div", "perc_del", "perc_ins", "query_seq", "query_pos_begin", "query_pos_end", "query_left", "strand", "matching_repeat", "repeat_class", "repeat_pos_begin", "repeat_pos_end", "repeat_pos_left", "id", "asterisk"]
        with open(out_path,"r") as f:
            print("## gff-version 3")
            for i in range(3):
                    next(f)
            for line in f:
                col = line.rstrip().split()
                if len(col) == 15:
                    col.append("")
                r = dict(zip(header,col))
                if r["strand"] == "C":
                    s = "-"
                    tstart = r["repeat_pos_end"]
                    tend = r["repeat_pos_begin"].replace("(","").replace(")","")
                else:
                    s = "+"
                    tstart = r["repeat_pos_begin"]
                    tend = r["repeat_pos_end"]

                note = "Target=" + r["matching_repeat"] + " " + tstart + " " + tend
                print_col = [r["query_seq"], "RepeatMasker", r["repeat_class"], r["query_pos_begin"], r["query_pos_end"], r["sw_score"], s, ".", note]
                if len(line) > 1:
                    print(*print_col, sep = "\t")

    def out2saf(self):
        print(*["GeneID", "Chr", "Start", "End", "Strand"], sep="\t")
        for _, r in self.out.iterrows():
            # GeneID = class|chr|start|end
            geneid = "|".join([r["matching_repeat"], r["query_seq"], r["query_pos_begin"], r["query_pos_end"]])
            strand = "-" if r["strand"] == "C" else "+"
            print(*[geneid, r["query_seq"],  r["query_pos_begin"], r["query_pos_end"], strand], sep = "\t")

    def out2bed(self): # bed file is 0-base, out and gff3 is 1-base
        """
         * bed file is 0-base, out and gff3 is 1-base
         output: 6 column bed file

        """
        for _, r in self.out.iterrows():
            name = r["matching_repeat"] + ";" + r["repeat_class"]
            print(*[r["query_seq"], r["query_pos_begin"]-1,  r["query_pos_end"], name, "0", r["strand"]], sep = "\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type = str, help="RepeatMasker .out file")
    parser.add_argument("--gff3", action="store_true", help="Convert out file to GFF3 format")
    parser.add_argument("--raw_gff3", action="store_true", help="Convert out file to raw GFF3 format from repeatmasker. The attribute column is Target=<family> start end")
    parser.add_argument("--saf", action="store_true", help="Convert out file to SAF format (for feature count)")
    parser.add_argument("--bed", action="store_true", help="Convert out file to BED format , name = matching_repeat;repeat_class")
    args = parser.parse_args()
    
    
    if args.gff3:
        o = Rm_out(args.out)
        print("Converting .out to gff3",file=sys.stderr)
        o.out2gff3()
    if args.raw_gff3:
        Rm_out.out2_raw_gff3(args.out)
    if args.saf:
        o = Rm_out(args.out)
        print("Converting .out to saf",file=sys.stderr)
        o.out2saf()
    if args.bed:
        o = Rm_out(args.out)
        print("Converting .out to bed",file=sys.stderr)
        o.out2bed()
