# modules
import sys
import argparse

# argparse
def parse_args(program_version, arguments=sys.argv[1:]):
    # main_arguments
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="poatools performs partial order alignment between CCS reads and subreads from the same zero mode waveguide",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )
    # subcommands
    subparsers = parser.add_subparsers(dest="sub")

    # subcommands: poa
    parser_poa.add_argument(
        "--ccs",
        type=str,
        required=True,
        help="CCS.fastq"
    )
    parser_poa.add_argument(
        "--subreads",
        type=str,
        required=True,
        help="subreads.fasta"
    )
    parser_poa.add_argument(
        "--sub",
        type=str,
        required=True,
        help="query substitution (ZMW:POS_ALT/REF)"
    )
    parser_poa.add_argument(
        "--out",
        type=str,
        required=True,
        help="file to return partial order alignment",
    )

    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
