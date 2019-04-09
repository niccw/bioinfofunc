#!/usr/bin/env python3
"""
Scripts to process BED file.
This script contains functions that do not read the row to memory (while bioinfo_bed.py do)
"""

import argparse

class Bed2(object):
    def __init__(self,path):
        self.path = path

    def filter_bed(self,filter_words_p):
        """
        Check if the name column (i.e. rnd-1_family-1024;DNA/Academ-1) contain the filter word.
        """
        # Mark down the filter_words
        filter_words_list = []
        with open(filter_words_p,"r") as f:
            for line in f:
                filter_words_list.append(line.strip())
        #print(filter_words_list)

        with open(self.path,"r") as f:
            for line in f:
                col = line.strip().split("\t")
                name = col[3]
                rep_class = name.split(";")[1]
                if "/" in rep_class:
                    rep_class_general = rep_class.split("/")
                else:
                    rep_class_general = rep_class
                if rep_class_general in filter_words_list:
                    next
                else:
                    print(*col,sep="\t")
                
    
if __name__ == "__main__":
    help_text = """
    """

    parser = argparse.ArgumentParser(description=help_text,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bed", help="bed file")
    parser.add_argument("--filter", action="store_true", help="Default = False")
    parser.add_argument("--words", dest="fw", type = str, help="txt file contains words to be filter. (compare against name column of BED file)")

    args = parser.parse_args()

    # filter
    if args.filter:
        bed = Bed2(args.bed)
        bed.filter_bed(filter_words_p=args.fw)