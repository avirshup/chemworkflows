#!/usr/bin/env python

import argparse

from . import runapp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('appname')
    parser.add_argument('inputfile')
    parser.add_argument('--outputdir', default=None)

    parser.add_argument('--localdocker', action='store_true')
    parser.add_argument('--here', action='store_true')
    parser.add_argument('--preprocess', action='store_true')
    parser.add_argument('--restart', action='store_true')
    parser.add_argument('--setoutput', nargs='*',
                        help=SETOUTPUTHELP)
    parser.add_argument('--dumptasks', action='store_true')

    args = parser.parse_args()

    runapp.main(args)

SETOUTPUTHELP = ("Set a node's outputs. Used to pass user input from an interactive task."
                 "The outputs must be described in a JSON file, which will be used as "
                 "the node's outputf fields.\n"
                 "   USAGE: --setoutput [taskname1]:[json file1], ...")
