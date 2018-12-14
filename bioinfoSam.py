#!/usr/bin/env python3

from array import array
import argparse
import pysam
import sys

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Customised functions for SAM/BAM")
    parser.add_argument("sam", type = str, help="SAM or BAM file")
    parser.add_argument("-i,--insertCoverage",type = str, dest="ic",default=True, help="Calcualte the insert coverage and output .cov and .scaffoldmap")
    parser.add_argument("-o,--out",type = str, dest="o",help="Prefix of output files")

    args = parser.parse_args()

    if args.ic:
        s = Sam(args.sam)
        s.insertCoverage(args.o)
