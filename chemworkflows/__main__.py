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

    args = parser.parse_args()

    runapp.runapp(args)
