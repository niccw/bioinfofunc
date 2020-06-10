#!/usr/bin/env python3
import sys
import argparse

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

class Gff3(object):
    def __init__(self,path):
        self.path = path

    @staticmethod
    def _attribute_to_dict(a):
        d = {}
        for i in a.split(";"):
            if len(i) > 0:
                k,v = i.split("=")
                d[k] = v
        return d

    def rc_gff2bed(self):
        """
        bed file is 0-base; gff3 is 1-base
        """
        with open(self.path,"r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    next(f)
                    continue
                col = line.strip().split("\t")
                chrom = col[0]
                general_class = col[2]
                start = col[3]
                end = col[4]
                strand = col[6]
                # parse attr to dict
                attr_dict = Gff3._attribute_to_dict(col[8])
                name = attr_dict["ID"] + ";" + general_class
                # output it as bed format (six column)
                print(*[chrom,int(start)-1,int(end)-1,name,"0",strand], sep="\t")

if __name__ == "__main__":
    help_text = """
    """
    parser = argparse.ArgumentParser(description=help_text,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("gff", help="gff3 file")
    parser.add_argument("--bed", action="store_true", help="Convert gff3(repeatcraft rmerge) to bed (6 column)")

    args = parser.parse_args()

    # bed
    if args.bed:
        gff = Gff3(args.gff)
        gff.rc_gff2bed()
