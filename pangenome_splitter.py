#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:20:23 2020

@author: lillycummins

Motive: Split pangenome reference sequences into separate core and accessory lists
using gene presence/absence data.
"""

import os
import time
import csv 
import argparse

#parse command line arguments
parser = argparse.ArgumentParser(description="Separate pangenome reference fasta file into separate core and accessory fastas using gene presence/absence data.")
parser.add_argument("-m","--matrix", help="Input pangenome gene presence/absence matrix", dest="csv_in", type=str, required=True)
parser.add_argument("-f","--fasta", help="Input pangenome fasta file", dest="fasta_in", type=str, required=True)
parser.add_argument("-o", "--output", help="Specify output file name", dest="output", type=str, required=True)
parser.add_argument("-core", "--core_threshold", help="Core genome sample threshold", dest="threshold", type=float, required=True)
args=parser.parse_args()

#check user input
if not os.path.isfile(args.fasta_in):
    print('Cannot find input fasta file.')
    print('Exiting')

if not os.path.isfile(args.csv_in):
    print('Cannot find input csv file.')
    print('Exiting')

#if output is left empty make file name with time stamp
if args.output == None:
    args.output = os.path.join(os.getcwd(),'Pangenome_splitter_'+time.strftime('%Y%m%d_%H%M')) + os.sep
elif not os.path.isdir(args.output):
    os.makedirs(args.output)

#create storage lists for gene names
core_genes = list()
accessory_genes = list()

#assigning genes as core or accessory
print('Classifying genes...')
with open(args.csv_in, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader) #skip header line
    for tabline in reader: #file contains 14 columns before gene presence/absence, presence is marked by genome name: convert to True with bool then convert to 1 with float
        tabline[14:len(tabline)] = [float(bool(x)) for x in tabline[14:len(tabline)]]
        presence_percentage = float(sum(tabline[14:len(tabline)])) / float(len(tabline)-14)
        if presence_percentage < args.threshold: #core threshold 
            accessory_genes.append(tabline[0])
        elif presence_percentage >= args.threshold:
            core_genes.append(tabline[0])

#create storage lists for gene names and sequences
accessory_seqs=list()
core_seqs=list()

#creating fasta style lists        
with open(args.fasta_in, 'r') as pan_fasta:
    fasta_lines=pan_fasta.readlines() 
    is_acc=False
    for line in fasta_lines:
        if line.startswith('>') and line.strip(">\n") in accessory_genes:
            #some way to save all the lines until the next > to accessory_seqs
            accessory_seqs.append(line)
            is_acc=True                    
        elif line.startswith('>') and line.strip(">\n") in core_genes:
            #some way to save all the lines until the next > to core_seqs
            core_seqs.append(line)
            is_acc=False
        else:
            if is_acc:
                accessory_seqs.append(line)
            else:
                core_seqs.append(line)
            
print('Writing sequences to file...')
with open(args.output+'_accessory.fasta', 'w') as f:
    for item in accessory_seqs:
        f.write(item)
with open(args.output+'_core.fasta', 'w') as f:
    for item in core_seqs:
        f.write(item)
        
print('Complete.')