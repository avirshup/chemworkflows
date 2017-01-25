#!/usr/bin/env python

# NOTE: currently hardcoded to the MMinimize app

import sys
from . import runapp


runapp.runapp(*sys.argv[1:])
