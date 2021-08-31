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
    # subcommands: poa
    parser_poa = subparsers.add_parser(
        "poa",
        help="partial order alignment.",
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
        help="query substitution (ZMW:POS_ALT/REF)"
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
        "--poa_subset",
        type=str,
        required=True,
        help="file to return subset of the annotated partial order alignment with the qsub",
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
