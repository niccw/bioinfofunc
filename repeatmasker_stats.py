#!/usr/bin/env python3

# stats fro class (LINE, SINE) and superfamily (i.e CR1,L2)
# input: reforma gff3, feature($3) is class/superfamily
# count of each class and superfamily
# total bp for each class and superfamily
# average length of each alss and superfamily
# gff3 file

import sys

def repeatmasker_gff_stat(path:str):
    #header = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
    with open(path,"r") as f:
        # initiate dicts for class and superfamily
        # d[class] = [#,bp]
        class_d ={}
        superfamily_d = {}

        for line in f:
            if line.startswith("#"):
                continue
            [_,_,feature,start,end,_,_,_,_] = line.strip().split("\t")
            # first process the id, most of the classes are in class/superfamily format, if no "/" -> treat it as class
            if "/" in feature:
                c = feature.split("/")[0]
                superfamily = feature.split("/")[1]
                f = c + "_" + superfamily
            else:
                c = feature
                f = False
            # calculate the length of the item
            bp = int(end) - int(start)
            
            # fill the class dict 
            if c not in class_d:
                class_d[c] = [1,bp]
            else:
                class_d[c][0] += 1
                class_d[c][1] += bp
            # fill the superfamilt dict
            if f:
                if f not in superfamily_d:
                    superfamily_d[f] = [1,bp]
                else:
                    superfamily_d[f][0] += 1
                    superfamily_d[f][1] += bp
    # print out the result
    # level name # bp avg.length
    print(f"level\tname\tcount\ttotal_bp\tavg.length")
    # print class
    for k,v in class_d.items():
        print(f"class\t{k}\t{v[0]}\t{v[1]}\t{v[1]/v[0]}")
    # print superfamily
    for k,v in superfamily_d.items():
        print(f"superfamily\t{k}\t{v[0]}\t{v[1]}\t{v[1]/v[0]}")

if __name__ == "__main__":
    repeatmasker_gff_stat(sys.argv[1])
