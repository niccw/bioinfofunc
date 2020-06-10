#!/usr/bin/env python3
import sys

def delta2bed(delta):
    # .delta is one-based, space-seperated
    # start with >, aligment between 2 pair, mark down the chr(s)
    # if len > 1, start end pos of alignment, read the first 4 digit
    ## sub_start,sub_end,query_start,query_end
    # output: sub_chr sub_start sub_end [query_chr]:[query_start]-[query_end]
    # if len == 1, ignore
    # if meet 0, it represet the end of a indel info, just ignore it

    with open(delta,"r") as f:
        for i in range(2):
            next(f)
        for line in f:
            # header?
            if line.startswith(">"):
                [sub_chr, query_chr, _,_] = line.replace(">","").strip().split(" ")
                continue
            col = line.strip().split(" ")
            if len(col) > 1:
                # aligment position info
                # sub pos is always + strand
                [sub_start, sub_end, query_start, query_end,_,_,_] = col
                print(f"{sub_chr}\t{sub_start}\t{sub_end}\t{query_chr + ':' + query_start + '-' + query_end}")

            else: # indel info
                continue
                
if __name__ == "__main__":
    delta_path = sys.argv[1]
    delta2bed(delta=delta_path)