#!/usr/bin/env python3

from array import array
import argparse
import pysam
import sys
from pathlib import Path
import bisect


class Sam(object):
    def __init__(self, path):
        self.path = path

    def insertCoverage(self, outpath, test=False):
        o = open(outpath+".cov", "w+")
        f = pysam.AlignmentFile(self.path)
        # Read contig length
        clen = dict(zip(f.references, f.lengths))

        sl = []  # scaffold list
        clist = []  # coverage list
        prevScaffold = ""
        for read in f.fetch():
            # re-set the covergae dict if this is a new scaffold
            if read.reference_name != prevScaffold:
                sl.append(read.reference_name)
                # print previous scaffold if not the head of file
                if prevScaffold != "":
                    # use the sorted index
                    idx = sorted([*clen]).index(prevScaffold)
                    for i in clist:
                        print(f"{idx}\t{i}", file=o)
                # reset array
                clist = array('L', [0]*clen[read.reference_name])
                prevScaffold = read.reference_name  # update variable
            # skip non unique mapped reads
            if read.has_tag("XA"):
                continue
            # add coverage to position
            epos = read.pos + read.tlen
            for p in range(read.pos, epos):
                clist[p] += 1

        # print last chrom
        idx = sorted([*clen]).index(prevScaffold)
        for i in clist:
            print(f"{idx}\t{i}", file=o)

        # Add empty array for scaffold not present in SAM/BAM
        if not test:
            for k, v in clen.items():
                if k not in sl:
                    clist = array('I', [0]*v)
                    idx = sorted([*clen]).index(k)
                    for i in clist:
                        print(f"{idx}\t{i}", file=o)
        f.close()
        o.close()
        omap = open(outpath+".scaffoldmap", "w+")
        for i, v in enumerate(sorted([*clen])):
            print(f"{v}\t{i}", file=omap)
        print(
            f"Done. Coverage is written in {outpath}.cov and Scaffold map is written in {outpath}.scaffoldmap", file=sys.stderr)

    def extract_read(self, qname, flag=None, oname="out.bam"):
        # self only contain the path, not a sam file yet
        # qname is a list
        sam = pysam.AlignmentFile(self.path)
        outfile = pysam.AlignmentFile(oname, "wb", template=sam)

        if flag is None:
            for read in sam.fetch():
                print(read.qname)
                # assume qname is long, use binary search
                if len(qname) > 0: # still have items to look for
                    binary_idx = bisect.bisect(qname,read.qname)-1
                    m = qname[binary_idx] == read.qname
                    if m:
                        outfile.write(read)
                        # assume unique value in qname list
                        qname.pop(binary_idx)
                else:
                    break  # found all qname needed, stop looking at remaining sam/bam
        else:  # extract read with specific flag only
            for read in sam.fetch():
                if len(qname) > 0:
                    binary_idx = bisect.bisect(qname,read.qname)-1
                    m = qname[binary_idx] == read.qname
                    if m and read.flag == flag:
                        outfile.write(read)
                        # assume unique value in qname list
                        qname.pop(binary_idx)
                else:
                    break


if __name__ == "__main__":
    help_text = """
    Customised functions for SAM/BAM (sort + index)self.

    Examples:
    - Calculate insert coverage:
     ./bioinfo_sam.py -ic example.sam -o example
    - Extract reads:
     ./bioinfo_sam.py --extract example.sam -f qname.list --flag 163
    """

    parser = argparse.ArgumentParser(
        description=help_text, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sam", type=str, help="SAM or BAM file (sort + index)")
    parser.add_argument("-ic", action="store_true",
                        help=" Calcualte the insert coverage and output .cov and .scaffoldmap")
    parser.add_argument("-o,--out", type=str, dest="o",
                        help="Prefix of output files")
    parser.add_argument("--test", action="store_true",
                        help="Run in test mode.")
    parser.add_argument("--extract", action="store_true",
                        help="extract read by qname. Work with -f and --flag (optional)")
    parser.add_argument("-f", type=str, dest="f", help="list file")
    parser.add_argument("--flag", type=int, dest="flag",
                        default=None, help="extract read with specific flag")

    args = parser.parse_args()

    if args.ic:
        s = Sam(args.sam)
        if args.test:
            print("Run in test mode.", file=sys.stderr)
        s.insertCoverage(args.o, test=args.test)

    if args.extract:
        with open(args.f, "r") as q:
            qname_list = q.read().splitlines()
        qname_list.sort()
        s = Sam(args.sam)
        # filter bam file name
        outname = Path(args.sam).stem + ".filtered.bam"
        s.extract_read(qname=qname_list, oname=outname, flag=args.flag)
