#!/usr/bin/env python3

from array import array
import argparse
import pysam
import sys

class Sam(object):
    def __init__(self,path):
        self.path = path
    
    def insertCoverage(self,outpath,test=False):
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
                #print(read.reference_name)
                sl.append(read.reference_name)
                # print previous scaffold if not the head of file
                if prevScaffold != "":
                    idx = sorted([*clen]).index(read.reference_name)  # use the sorted index (= chromIdx.tsv)
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

        # print last chrom
        idx = sorted([*clen]).index(read.reference_name)
        for i in clist:
            print(f"{idx}\t{i}",file=o)
            
        # Add empty array for scaffold not present in SAM/BAM
        if not test:
            for k,v in clen.items():
                if k not in sl:
                    clist = array('I',[0]*v)
                    idx = [*clen].index(k)
                    for i in clist:
                        print(f"{idx}\t{i}",file=o)
        f.close()
        o.close()
        omap = open(outpath+".scaffoldmap","w+")
        for i,v in enumerate(sorted([*clen])):
            print(f"{v}\t{i}",file=omap)
        print(f"Done. Coverage is written in {outpath}.cov and Scaffold map is written in {outpath}.scaffoldmap")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Customised functions for SAM/BAM (sort + index)")
    parser.add_argument("sam", type = str, help="SAM or BAM file")
    parser.add_argument("-ic",action="store_true", help=" Calcualte the insert coverage and output .cov and .scaffoldmap")
    parser.add_argument("-o,--out",type = str, dest="o",help="Prefix of output files")
    parser.add_argument("--test",action="store_true", help="Run in test mode.")
    
    args = parser.parse_args()

    if args.ic:
        s = Sam(args.sam)
        if args.test:
            print("Run in test mode.",file=sys.stderr)
        s.insertCoverage(args.o,test=args.test)
     
