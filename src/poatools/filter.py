import re
import os
import sys
from typing import List, Dict, Tuple


dna = list("ATGC-")
def poa_filter(infile, outfile):


    o = open(outfile, "w")
    for line in open(infile):
        tname, tpos, ref, alt, _, cell, zmw, qsub = line.split()
        poa_file = "{}/poa/{}.subset.poa".format(cell, zmw)
        if not os.path.exists(poa_file):
            print("{}\t{}\t{} partial order alignment file does not exist".format(cell, zmw, poa_file))
            continue

        subread_count = 0
        base_count = {base: 0 for base in dna}
        for idx, line in enumerate(open(poa_file).readlines()):
            if idx == 0: # qsub
                arr = line.split("\t")
                jdx = arr[1].index("*")
            elif idx == 1 or idx == 2: # ccs
                continue
            else:
                arr = line.split()
                if len(arr) != 0:
                    subread_count += 1
                    subread = line.split()[1]
                    if not jdx in subread:
                        continue
                    subread_base = subread[jdx]
                    base_count[subread_base] += 1
        subread_ref_count = base_count[ref] 
        subread_alt_count = base_count[alt]
        if subread_count == subread_alt_count:
            o.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:BQ:CELL:ZMW\t{}:{}:{}:{}\n".format(tname, tpos, ref, alt, "./.", "93", cell, zmw))
    o.close()





