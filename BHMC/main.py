#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
from BHMC import BHMC

__author__ = "Fuqiang Gong"
__email__ = "fqgong@stu.xmu.edu.cn"

def main_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(description="""
    BHMC is a simple script to generate the structure of pure metalic cluster.""")

    subparsers = parser.add_subparsers()

    parser_run = subparsers.add_parser(
        "run"
        )
    parser_run.add_argument('PARAM', type=str,
                        help="parameter file, json/yaml format")
    parser_run.set_defaults(func=gen_run)
    return parser

def main():
    info()
    parser = main_parser()
    try:
        import argcomplete
        argcomplete.autocomplete(parser)
    except ImportError:
        # argcomplete not present.
        pass

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == "__main__":
    main()
