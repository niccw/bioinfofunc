#!/usr/bin/env python3

import argparse
import pandas as pd

class Gff(object):
    def __init__(self,path):
        self.path = path

    
    @staticmethod
    def readGff(g):
        with open(g,"r") as f:
            for line in f:
                pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", type = str, help = "GFF3 file.")
    args = parser.parse_args()

    pass


