#! /usr/bin/env python

####################################
# pca_unrel.py
# written by Raymond Walters, July 2015
"""
Runs PCA for GWAS data with related individuals
"""
# Overview:
# 1) Input plink bed/bim/fam
# 2) Run strict QC on the input data
# 2) Define set of unrelated individuals using PRIMUS
# 3) Compute PCA on the unrelated set and projects to remainder
# 4) Plot projected PCs
#
####################################


