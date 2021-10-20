# modules
import sys
import argparse

# argparse
def parse_args(program_version, arguments=sys.argv[1:]):
    # command
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="poatools performs partial order alignment between CCS reads and subreads from the same zero mode waveguide",
    )
    # subcommand
    subparsers = parser.add_subparsers(dest="sub")
    # subcommands: whitelist
    parser_whitelist = subparsers.add_parser(
        "whitelist",
        help="return whitelist ZMWs",
    )
    parser_whitelist.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file to read",
    )
    parser_whitelist.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="himut or deepvarinat VCF file to read",
    )
    parser_whitelist.add_argument(
        "--min_bq",
        type=int,
        default=93,
        required=False,
        help="minimum base quality score (default = 93)",
    )
    parser_whitelist.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score (default = 60)",
    )
    parser_whitelist.add_argument(
        "--min_hq_base_fraction",
        type=float,
        default=0.5,
        required=False,
        help="minimum high quality base (BQ=93) proportion (default = 0.5)",
    )
    parser_whitelist.add_argument(
        "--min_sequence_identity",
        type=float,
        default=0.99,
        required=False,
        help="minimum sequence identity threshold (default = 0.99)",
    )
    parser_whitelist.add_argument(
        "-o",
        "--whitelist",
        type=str,
        required=True,
        help="file to write",
    )
    # subcommands: poa
    parser_poa = subparsers.add_parser(
        "poa",
        help="partial order alignment",
    )
    parser_poa.add_argument(
        "--cell",
        type=str,
        required=True,
        help="SMRT cell (e.g m64125_201017_124255)"
    )
    parser_poa.add_argument(
        "--zmw",
        type=str,
        required=True,
        help="Zero-mode waveguide (ZMW) number"
    )
    parser_poa.add_argument(
        "--ccs",
        type=str,
        required=True,
        help="ZMW.ccs.fastq"
    )
    parser_poa.add_argument(
        "--subreads",
        type=str,
        required=True,
        help="ZMW.subreads.fasta"
    )
    parser_poa.add_argument(
        "--qsub",
        type=str,
        required=True,
        help="query substitution (ZMW:POS_ALT/REF) or query substitution list separated by comma"
    )
    parser_poa.add_argument(
        "--length",
        type=int,
        default=100,
        required=True,
        help="multiple sequence alignment length"
    )
    parser_poa.add_argument(
        "--poa",
        type=str,
        required=True,
        help="file to return all the annotated partial order alignment",
    )
    parser_poa.add_argument(
        "--pdf",
        required=False,
        action="store_true",
        help="return partial order alignment as a PDF file (default = False)",
    )
    parser_poa.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )
    # subcommands: filter
    parser_filter = subparsers.add_parser(
        "filter",
        help="filter substitutions based on partial order alignment",
    )
    parser_filter.add_argument(
        "--input",
        type=str,
        required=True,
        help="file to read tsub, SMRTcell and ZMW, qsub"
    )
    parser_filter.add_argument(
        "--output",
        type=str,
        required=True,
        help="file to return filtered substitutions"
    )
    # no arguments
    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
