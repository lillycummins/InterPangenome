#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 16:57:51 2021

@author: lillycummins

Motive: given an abricate summary file, reformat file with summed partial hits of genes, capped at 100
"""

import argparse
import os
import time
import pandas as pd

parser = argparse.ArgumentParser(description = 'Process partial hits from ABRicate summary.')
parser.add_argument("--input", "-i", type=str, required=True) #abricate summary.tsv file
parser.add_argument("--output", "-o", type=str, required=True) #output file name
args = parser.parse_args()

#check user input
if not os.path.isfile(args.input):
    print('Cannot find input summary file.')
    print('Exiting')

#if output is left empty make file name with time stamp
if args.output == None:
    args.output = os.path.join(os.getcwd(),'summed_capped_'+time.strftime('%Y%m%d_%H%M')) + os.sep
elif not os.path.isdir(args.output):
    os.makedirs(args.output)

with open(args.input) as sumfile:
    df = pd.read_csv(sumfile, sep='\t') #load summary file as a dataframe
    df = df.replace('.',0) #replace entries where no match was found with 0
    headings = df.columns.values.tolist() #take column headers 
    genes = headings[2:] #create list of gene clusters to iterate through
    for cluster in genes: #for each gene cluster in column headings
        cluster_coverage_list = df[cluster].tolist() #take a list of observed coverage of cluster per assembly
        final_coverage_list=[]
        for coverage in cluster_coverage_list: #need to check if there is one entry per coverage reported
            if type(coverage) == str: #this means there are multiple hits reported in one assembly
                new_coverage_list=[]
                split_coverages=coverage.split(";")
                new_coverages = [float(i) for i in split_coverages]
                total_coverage = sum(new_coverages) #add all partial hits together
                final_coverage_list.append(total_coverage)
            else:
                final_coverage_list.append(coverage)     
        for i in range(0,len(final_coverage_list)):
            if (final_coverage_list[i] >100) == True:
                final_coverage_list[i] = 100 #cap %coverage at 100
            else:
                continue
        df[cluster]= final_coverage_list #replace column with list that has summed partial hits
      
df.to_csv(args.output+".tsv", sep="\t", index=False)
