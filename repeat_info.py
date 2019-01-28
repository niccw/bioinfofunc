#!/usr/bin/env python3
import argparse
import pybedtools
import sys


def repeat_info(repeat_bed,gene_bed,exon_bed):
    r = pybedtools.BedTool(repeat_bed)
    g = pybedtools.BedTool(gene_bed)
    e = pybedtools.BedTool(exon_bed)

    r_in_gene = r.intersect(g,u=True)
    r_not_in_gene = r.intersect(g,v=True)

    r_in_gene_annot = r_in_gene.intersect(g,wo=True) #overlap gene in #8 and overlap bp in #9

    # overlap with exon
    r_in_gene_annot_ei = r_in_gene_annot.intersect(e,wao=True) # #13 > 0 == in exon, else not in exon

    # closest gene of r_not_in_gene
    r_not_in_gene = r_not_in_gene.sort()
    r_not_in_gene_closest = r_not_in_gene.closest(g,d=True) # #8 is the cloest gene #9 is the distance between 

    # filter and save record to list array
    out_bed = []

    # 0   1    2    3     4      5   6   7   
    # chr s    e    re  gene   e/i  neigh  bp
    for i in r_in_gene_annot_ei:
        # i is an pybedtools interval object, can be accessed as a list or dictionary
        ei = "e" if int(i[12]) > 0 else "i"
        out_bed.append([i[0],i[1],i[2],i[3],i[7],ei,"NA","NA"])

    for i in r_not_in_gene_closest:
        out_bed.append([i[0],i[1],i[2],i[3],i[7],"NA",i[7],i[8]])

    # print the result out
    for i in out_bed:
        print(*i, sep="\t")


if __name__ == "__main__":
    help_text = """
    ./repeat_info.py --repeat_bed repeat.bed --gene_bed gene.bed --exon_bed exon.bed
    """
    parser = argparse.ArgumentParser(description=help_text,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r,--repeat_bed", type = str,dest="repeat_bed", help="repeat element in bed format")
    parser.add_argument("-g,--gene_bed", type = str,dest="gene_bed", help="gene.bed")
    parser.add_argument("-e,--exon_bed", type = str,dest="exon_bed", help="exon.bed")
    args = parser.parse_args()

    if len(sys.argv) >=3:
        repeat_info(repeat_bed=args.repeat_bed, gene_bed=args.gene_bed, exon_bed=args.exon_bed)
    else:
        sys.exit(help_text)


