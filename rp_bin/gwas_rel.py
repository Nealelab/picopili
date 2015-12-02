#! /usr/bin/env python

####################################
# gwas_rel.py
# written by Raymond Walters, December 2015
"""
Runs GWAS of plink family data
"""
# Overview:
# 1) split files for parallelization
# 2) gwas chunks
# 3) aggregate and format (add a2, maf, info, etc)
# 4) plots/output
#
####################################


# treat this as driver script
# need to write:
# chunker (with filter on meta-info)
# aggregator (for each chunk: read gwas, read meta info, format/join, append to final outfiles)
# manhattan/qq plots