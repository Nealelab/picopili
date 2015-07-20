#! /usr/bin/env python

####################################
# final_file_check.py
# written by Raymond Walters, July 2015
"""
Checks for file existance, and sends success/error to ricopili email address
"""
# wrapper to file_check_email(filename,taskname) from py_helpers.py
#
####################################

import os
import subprocess
import argparse
from py_helpers import *


### parse arguments
parser = argparse.ArgumentParser(prog='final_file_check.py')

parser.add_argument('--filename', 
                    type=str,
                    metavar='FILENAME',
                    help='full path to file to check for existence',
                    required=True)
parser.add_argument('--taskname',
                    type=str,
                    metavar='STRING',
                    help='Task to report succeeded or failed; no spaces',
                    required=True)

args = parser.parse_args()


isSuccess = file_check_email(args.filename, args.taskname)

if isSuccess:
    exit(0)
else:
    exit(1)



