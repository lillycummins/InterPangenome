#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 11:20:37 2021

@author: lillycummins

Motive: count number of genes per COG cat from eggnog annotation files and show 
#COG cat % as heatmap for each ST (optimised for visualisation in spyder v5.0.5)
"""

import os
import csv 
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

parser = argparse.ArgumentParser('Create functional distribution heatmap')
parser.add_argument("--path", "-p", type=str, required=True) #path to directory of eggnog annotation tsvs
parser.add_argument("--output", "-o", type=str, required=True) #path to output matrix
args = parser.parse_args()

possiblecategories=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
                    'P','Q','S','T','U','V','W','Y','Z']
df=pd.DataFrame()

for filename in os.listdir(args.path):
    COG_categories=list()
    ST=filename.split('_')[0]
    with open(args.path+filename) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        next(reader)
        for column in reader: #column 6 for narrowest or 9 for best COG category
            if '#' in column[0]:
                continue
            else:
                category = column[9]
                COG_categories.extend([category])

###COUNT INSTANCES OF EACH COG CATEGORY   
    d=defaultdict(int)
    for x in possiblecategories:
        d[x]=0
            
    for x in COG_categories:
        if x=='-':
            d['?']+=1
        for n in list(x):
            if n in possiblecategories: #check character is category and not a comma or space
                d[n]+=1
            
    ditems=sorted(d.items(), key=lambda x:x[0]) #order categories alphabetically        
    
    x, function = list(zip(*ditems)) #x is d.values, function is d.keys
    
    percent=np.empty(len(ditems))
    for i in range(0,len(ditems)):
        percent[i]=function[i]/len(COG_categories)*100
    percents=percent.tolist()
    row_df = pd.DataFrame([percents], index = [ST], columns=x)
    df=pd.concat([df,row_df])


##REMOVE 'EMPTY' COG CATEGORIES
def find_empty_cats(data_set):
    redundant_cols = list()
    for i in range(0,data_set.shape[1]):
        column = data_set.iloc[:,i] #examine every column
        if all(column.iloc[x] == 0 for x in range(0,data_set.shape[0])) == True: #find empty cat
            redundant_cols.append(column.name)
        else:
            continue
    return redundant_cols

df=df.drop(columns = find_empty_cats(df))

df=df.drop(columns = ['S','?']) #remove groups of unknown function

##PHYLOGROUP ASSIGNMENT 
phylogroup = {'ST10':'A','ST11':'E','ST117':'F','ST12':'B2','ST127':'B2','ST131':'B2',
              'ST14':'B2','ST141':'B2','ST144':'B2','ST167':'A','ST17':'B1','ST21':'B1',
              'ST28':'B2','ST3':'B1','ST372':'B2','ST38':'D','ST410':'A','ST648':'F',
              'ST69':'D','ST73':'B2','ST95':'B2'}

phylo_ser=pd.Series(data=phylogroup,index=['ST10','ST11','ST117','ST12','ST127','ST131',
                                           'ST14','ST141','ST144','ST167','ST17','ST21',
                                           'ST28','ST3','ST372','ST38','ST410','ST648',
                                           'ST69','ST73','ST95'])

cols=['mediumpurple','salmon','pink','orange','lightgreen','c']
lut = dict(zip(phylo_ser.unique(), cols))
row_colors = phylo_ser.map(lut)

##PLOT CLUSTERMAP
plt.figure()
sns.set(font_scale = 1.2)
cmap = sns.color_palette("Blues", as_cmap=True)
cm = sns.clustermap(df, xticklabels=1, tree_kws=dict(linewidths=1.5), cmap=cmap, 
                    row_colors=row_colors, cbar_pos=(.03, .6, .03, .3))
ax=cm.ax_heatmap
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, horizontalalignment='right')
ax.set(ylabel=None)
handles = [Patch(facecolor=lut[name]) for name in sorted(lut)]
plt.legend(handles, sorted(lut), title='Phylogroup', prop={'size': 11},
            bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc = 'upper right')

plt.show()

df.to_csv(args.output+'.csv')