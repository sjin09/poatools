#!/usr/bin/env python3
__version__ = "0.0.1"
__author__ = "Sangjin Lee"


# modules
import os
import sys
import logging
import himut.vcflib 
import himut.mutlib 
import himut.raw_caller
import himut.pon_caller
import himut.tol_caller
import himut.basic_caller
import himut.debug_caller
import himut.somatic_caller
import multiprocessing as mp
from himut.parse_args import parse_args

def main():
    # stdout/stderr
    # sys.stdout = StringIO()
    # sys.stderr = StringIO()
    # sys.stdout.getvalue()
    
    # parse_args
    parser, options = parse_args(program_version=__version__)
    if options.sub == "call":  # call somatic substitutions
        thread_count = mp.cpu_count()
        if options.threads > thread_count:
            options.threads = thread_count
        if options.create_panel_of_normals:
            print("himut is calling substitutions for panel of normal preparation")
            himut.pon_caller.call_pon_substitutions(
                options.reads, # input # bam file
                options.region, # target contigs/scaffolds/chromosomes
                options.region_list, # target contigs/scaffolds/chromosomes fofn
                options.ploidy, # haploid, diploid or polyploid
                options.sample_is_reference, # are the reads and the reference genome from the same sample? # True/False
                30, # int: 0-60
                0.8, # blast sequence identity
                0.3, # proportion of BQ=93 bases
                options.common_snps, # common snps
                options.sample_snps, # germline mutations
                20, # bq ## variant_filter: params
                options.mismatch_window, # mismatch window size
                options.min_mismatch_count, # minimum mismatch count
                options.verbose, # true/false print other metrics
                options.threads, # maxminum number of threads
                __version__,
                options.vcf, # output # himut vcf file
            )
        else:
            if options.tree_of_life_sample:
                print("himut is calling somatic single base substitutions from a tree of life sample")
                himut.tol_caller.call_somatic_substitutions(
                    options.reads, # input # bam file
                    options.region, # target contigs/scaffolds/chromosomes
                    options.region_list, # target contigs/scaffolds/chromosomes fofn
                    options.ploidy, # haploid, diploid or polyploid
                    options.sample_is_reference, # are the reads and the reference genome from the same sample? # True/False
                    options.min_mapq, # int: 0-60
                    options.min_sequence_identity, # blast sequence identity
                    options.min_hq_base_fraction, # proportion of BQ=93 bases
                    options.sample_snps, # germline mutations
                    options.min_bq, # variant_filter: params
                    options.mismatch_window,
                    options.min_mismatch_count,
                    options.min_ref_count, # number of reads supporting the reference base
                    options.min_alt_count, # number of reads supporting the substitution
                    options.min_hap_count, # number of reads supporting h0 and h1 haplotype
                    options.min_hap_proportion, # (h0 reads + h1 reads) / (h0 reads + h1 reaads ... hn reads)
                    options.verbose, # true/false print other metrics
                    options.threads, # maxminum number of threads
                    __version__,
                    options.vcf, # output # vcf file
                )
            else:
                print("himut is calling somatic single base substitutions")
                himut.somatic_caller.call_somatic_substitutions(
                    options.reads, # input # bamfile
                    options.region, # target contigs/scaffolds/chromosomes
                    options.region_list, # target contigs/scaffolds/chromosomes fofn
                    options.ploidy, # haploid, diploid or polyploid
                    options.sample_is_reference, # are the reads and the reference genome from the same sample? # True/False
                    options.min_mapq, # int: 0-60
                    options.min_sequence_identity, # blast sequence identity
                    options.min_hq_base_fraction, # proportion of BQ=93 bases
                    options.common_snps, # common snps
                    options.sample_snps, # germline mutations
                    options.panel_of_normals, # panel of normals
                    options.min_bq, # variant_filter: params
                    options.mismatch_window,
                    options.min_mismatch_count,
                    options.min_ref_count, # number of reads supporting the reference base
                    options.min_alt_count, # number of reads supporting the substitution
                    options.min_hap_count, # number of reads supporting h0 and h1 haplotype
                    options.min_hap_proportion, # (h0 reads + h1 reads) / (h0 reads + h1 reaads ... hn reads)
                    # options.mat, # signatureAnalzyer/sigProfiler pentanucleotide sequence context mutational signature reference
                    options.verbose, # true/false print other metrics
                    options.threads, # maxminum number of threads
                    __version__,
                    options.vcf, # output # himut vcf file
                )
    elif options.sub == "raw":  # TODO: REMOVE
        thread_count = mp.cpu_count()
        if options.threads > thread_count:
            options.threads = thread_count
        print("himut is calling raw substitutions")
        himut.raw_caller.call_raw_substitutions(
            options.reads, # input # bamfile
            options.region, # target contigs/scaffolds/chromosomes
            options.region_list, # target contigs/scaffolds/chromosomes fofn
            options.threads, # maxminum number of threads
            __version__,
            options.vcf, # output # himut vcf file
        )
    elif options.sub == "basic":  # TODO: REMOVE
        thread_count = mp.cpu_count()
        if options.threads > thread_count:
            options.threads = thread_count
        print("himut is calling basic substitutions")
        himut.basic_caller.call_basic_substitutions(
            options.reads, # input # bamfile
            options.region, # target contigs/scaffolds/chromosomes
            options.region_list, # target contigs/scaffolds/chromosomes fofn
            options.ploidy, # haploid, diploid or polyploid
            options.sample_is_reference, # are the reads and the reference genome from the same sample? # True/False
            options.min_mapq, # int: 0-60
            options.min_sequence_identity, # blast sequence identity
            options.min_hq_base_fraction, # proportion of BQ=93 bases
            options.common_snps, # common snps
            options.sample_snps, # germline mutations
            options.panel_of_normals, # panel of normals
            options.min_bq, # variant_filter: params
            options.mismatch_window,
            options.min_mismatch_count,
            options.min_ref_count, # number of reads supporting the reference base
            options.min_alt_count, # number of reads supporting the substitution
            options.threads, # maxminum number of threads
            __version__,
            options.vcf, # output # himut vcf file
        )
    elif options.sub == "debug": # TODO: REMOVE
        print("himut is calling substitutions for each CCS read")
        himut.debug_caller.call_substitutions(
            options.reads, # input # bam file
            options.ploidy, # haploid, diploid or polyploid
            options.sample_is_reference, # are the reads and the reference genome from the same sample? # True/False
            options.sample_snps, # germline mutations
            options.sample_subs, # somatic mutations
            options.min_bq, # min base quality score
            options.mismatch_window, # mismatch window size
            options.min_mismatch_count, # minimum mismatch count
            options.tsv, # output # himut vcf file
        )
    elif options.sub == "merge": # returns
        himut.vcflib.vcf_merge(options.fofn, options.vcf)
    elif options.sub == "sigcounts": # returns sbs96 counts
        himut.mutlib.get_sigcounts(options.input, options.reference, options.output)
    elif options.sub == "normcounts": # returns normalised sbs96 counts
        himut.mutlib.get_normcounts(options.sigcounts, options.ref, options.reads, options.normcounts)
    else:
        logging.warning("The subcommand does not exist!\n")
        parser.print_help()


if __name__ == "__main__":
    main()
    sys.exit(0)
