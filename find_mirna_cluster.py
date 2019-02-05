#!/usr/bin/env python3

import argparse
import sys

def find_mirna_cluster(mirna_gff, mrna_gff):

    # mirna bed (sorted)
    # mrna bed (sorted)

    # output:
    # cluster_id    mir-id (sep by ,)   distance_between (sep by ,)

    """
    #1 same orientation
    #2 not seperated by 1) mirna in opposite direction, 2) mRNA
    #3 cluster size >= 2
    """
    last_mirna = {
        "name" : "",
        "chr" : "",
        "end_pos" : 0,
        "orientation" : ""
    }

    def update_dict(name,chr,end_pos:int,orientation):
        last_mirna["Name"] = name
        last_mirna["chr"] = chr
        last_mirna["end_pos"] = end_pos
        last_mirna["orientation"] = orientation

    cluster_mir_id = [] # [[mir-1,mir-2,mir-3],[mir-33,mir-36,mir-38]...]
    cluster_distance = [] # [[150,300],[200,30]]
    cluster_chr = []  #["chr1","chr1","chr2"]

    # mRNA cache
    cache = [] # start,end 

    # open two bed files simutaneouly 
    #chr2L   .       miRNA_primary_transcript        243035  243141  .       -       .       ID=MI0005821;Alias=MI0005821;Name=dme-mir-965
    with open(mirna_gff,"r") as f, open(mrna_gff,"r") as g:
        for line in f:
            # skip header
            if line.startswith("#"):
                continue
            col = line.strip().split("\t")

            # only process full transcript
            if col[2] != "miRNA_primary_transcript":
                continue
            else:
                # parse attribute 
                attr_d = _attribute_to_dict(col[8])
                
                # 1. if not in the same chr, just update and read continue
                if col[0] != last_mirna["chr"]:
                    if last_mirna["chr"] == "":
                        update_dict(attr_d["Name"],col[0],int(col[4]),col[6])
                        continue
                    else:
                        update_dict(attr_d["Name"],col[0],int(col[4]),col[6])
                        continue
                # 2. check if in same orientation as last mirna
                if col[6] != last_mirna["orientation"]:
                    update_dict(attr_d["Name"],col[0],int(col[4]),col[6])
                    continue
                else:
                    mrna_overlap = False 
                    # 2.1 check if there is mrna in between the two mirna
                    span_start = last_mirna["end_pos"]
                    span_end = int(col[3])

                    ####start--read mrna gff ####
                    # The file peacefully closed after the last line since we used with open for file handling
                    # it won't throw StopIteration, so we need to track this manually
                    last_line_flag = True
                    for line in g:
                        last_line_flag = False
                        if line.startswith("#"):
                            continue
                        colg = line.strip().split("\t")

                        # check the cache
                        if len(cache) > 0 :
                            # is the cache at the same chr with current miRNA?
                            if col[0] != cache[2]:
                                # cache not in the same chrom, forget it and keep looping mRNA until they are on the same chr
                                cache = []
                                continue 
                            # is the cahce overlap with the miRNAs inter-distance?
                            if overlap(span_start,cache[0],span_end,cache[1]) > 0:
                                mrna_overlap = True
                                # still cannot forget the cache, since the cache mRNA can be very large
                                # keep the cache and check next time
                                break
                            else:
                                if cache[1] < span_start: # cache_end < span start
                                    # cache is still lag behide the miRNA, keep looking to new line
                                    cache = [] # can forget the cache , keep going and check continue mRNA
                                    continue
                                if cache[0] > span_end: # cache_start > span end
                                    # over-read, back to miRNA
                                    break

                        # process current line in mRNA
                        mrna_start = int(colg[3])
                        mrna_end = int(colg[4])

                        if col[0] != colg[0]:
                            continue

                        if overlap(span_start,mrna_start,span_end,mrna_end) > 0:
                            mrna_overlap = True
                            break
                        if mrna_end < span_start:
                            continue
                        if mrna_start > span_end:
                            # over-read, mark the position to cache and compare cache first in continue loop
                            # ! if over-read and no overlap, then we can keep mrna_overlap = False
                            mrna_chr = colg[0]
                            cache = [mrna_start, mrna_end, mrna_chr]
                            break
                    #### end---read mrna gff ####
                    
                    # check if we have reached the end of mRNA gff
                    if last_line_flag:
                        # is the last line in the same chr as the miRNA?
                        if colg[0] != col[0]:
                            # no mRNA information in the current chr, so they must be a cluster if they fullfill 1 and 2
                            pass
                        else: # check if they overlap
                            if overlap(span_start,mrna_start,span_end,mrna_end) > 0:
                                mrna_overlap = True

                    if mrna_overlap:
                        update_dict(attr_d["Name"],col[0],int(col[4]),col[6])
                        continue
                    else:
                        # can group this miRNA with the last one :)
                        # check if this is the first initiation
                        if len(cluster_mir_id) > 0:  
                            # also check if the chr match ot not !!!!!! ROARRRR
                            if col[0] != cluster_chr[-1]:
                                pass
                            elif last_mirna["Name"] in cluster_mir_id[-1]: # the last miRNA is already in the 'current' cluster
                                cluster_mir_id[-1].append(attr_d["Name"])
                                cluster_distance[-1].append(int(col[3])-last_mirna["end_pos"]) 
                                update_dict(attr_d["Name"],col[0],int(col[4]),col[6])
                                continue
                            else:
                                # this is the first cluster
                                # or this is another miRNA copies on another chromsome ROARRRR i.e. dme-mir-10404
                                pass

                        cluster_mir_id.append([last_mirna["Name"],attr_d["Name"]])
                        cluster_chr.append(col[0])
                        cluster_distance.append([int(col[3])-last_mirna["end_pos"]])
                        update_dict(attr_d["Name"],col[0],int(col[4]),col[6])

    # finish reading files
    for i,_ in enumerate(cluster_mir_id):
        cluster_distance_string = ",".join([str(x) for x in cluster_distance[i]])
        cluster_id_string = ",".join(cluster_mir_id[i])
        print(f"{cluster_chr[i]}\t{cluster_id_string}\t{cluster_distance_string}")

def overlap(s1:int,s2:int,e1:int,e2:int)->int:
    return min(e1,e2)-max(s1,s2)

def _attribute_to_dict(a:str)->dict:
        d = {}
        for i in a.split(";"):
            if len(i) > 0:
                k,v = i.split("=")
                d[k] = v
        return d


if __name__ == "__main__":
    help_text="""
    This script assume only same miRNA only occur once per chr. Work fine in Dme but need to be tested on other fly species.

    Example:
    ./find_mirna_cluster.py --mirna mirna.gff3 --mrna mrna.gff > mirna_cluster.tsv
    """
    parser = argparse.ArgumentParser(description=help_text,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mirna", type = str,dest= "mirna_gff",help="Sorted miRNA gff3. Must contain the 'miRNA_primary_transcript' field")
    parser.add_argument("--mrna", type = str,dest= "mrna_gff",help="Sorted mRNA gff3. ")

    args = parser.parse_args()

    if not args.mirna_gff or not args.mrna_gff:
        parser.print_help()
        sys.exit()
    find_mirna_cluster(mirna_gff=args.mirna_gff, mrna_gff=args.mrna_gff)

    #find_mirna_cluster(mirna_gff=sys.argv[1], mrna_gff=sys.argv[2])