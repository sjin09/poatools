#!/usr/bin/env python3
__version__ = "0.0.1"
__author__ = "Sangjin Lee"


# modules
import os
import sys
import logging
import poatools.poa 
import poatools.filter
import poatools.whitelist
from poatools.parse_args import parse_args

def main():
    # parse_args
    parser, options = parse_args(program_version=__version__)
    if options.sub == "poa":
        print("poatools is performing partial order alignment between CCS reads and subreads")
        poatools.poa.get_poa(
            options.cell, 
            options.zmw, 
            options.ccs, # input # bamfile
            options.subreads, # target contigs/scaffolds/chromosomes
            options.qsub, # target contigs/scaffolds/chromosomes fofn
            options.length,
            options.poa, # partial order alignment
            options.pdf, # bool: return pdf of the partial order alignment
        )
    elif options.sub == "filter":
        print("poatools is filtering substitutions based on subread support")
        poatools.filter.poa_filter(
            options.input,
            options.output 
        )
    elif options.sub == "whitelist":
        print("poatools is identifying whitelist ZMWs for partial order alignment")
        poatools.whitelist.get_whitelist(
            options.bam,
            options.vcf, 
            options.min_bq, 
            options.min_mapq, 
            options.min_hq_base_fraction, 
            options.min_sequence_identity, 
            options.whitelist, 
        )

if __name__ == "__main__":
    main()
    sys.exit(0)
