#!/usr/bin/env python3
# coding: utf-8

# In[9]:


#variant_filtering
import pandas as pd
import re
import argparse
import sys
import warnings

# check if you have python 3
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

# initiate the parser
parser = argparse.ArgumentParser(description=" A script to filter vep outputs"
                                )
# add long and short argument
parser.add_argument("--vepfile", help="Input vep output file")
parser.add_argument("--list", help = "Input list of genes sought")


# read arguments from the command line
args = parser.parse_args()
if not args.vepfile:  
    print("A vep file is required \n Use --help for more details")
    sys.exit(1)
if not args.list:
    print("A gene list file is required \n Use --help for more details")
    sys.exit(1)

x = pd.read_csv(args.vepfile, sep = '\t', skiprows = 59)
gene_lists = pd.read_csv(args.list)


f = open(args.vepfile)
name = f.name

outname = re.split("/",name)[-1]


variant_types = x["Consequence"].value_counts()

y = pd.merge(x, gene_lists, on="SYMBOL", how="inner")

filt_types = y["Consequence"].value_counts()


can = y[y["CANONICAL"] == "YES"]

cantypes = can["Consequence"].value_counts()


Missense = can[can["Consequence"] == "missense_variant"]


#output results to file
#add the summaries of variants too
filt_types.to_csv("%s_types_counts_adme_only.txt"%outname, header=True, index=True, sep='\t', mode='a')
cantypes.to_csv("%s_canonical_counts_adme_only.txt"%outname, header=True, index=True, sep='\t', mode='a')
can.to_csv("%s_canonical_adme_only.txt"%outname, header=True, index=None, sep='\t', mode='a')
Missense.to_csv("%s_missense_can_adme_only.txt"%outname, header=True, index=None, sep='\t', mode='a')


#filter some results for those that score highly
filtMissense = Missense[Missense["FATHMM_converted_rankscore"] != "-"]
errorMissense = Missense[Missense["FATHMM_converted_rankscore"] == "-"]
numscores = filtMissense[(pd.to_numeric(filtMissense["MetaLR_rankscore"]) > 0.81) | (pd.to_numeric(filtMissense["FATHMM_converted_rankscore"]) > 0.81) | (pd.to_numeric(filtMissense["MetaSVM_rankscore"]) > 0.81) | (pd.to_numeric(filtMissense["VEST3_rankscore"]) > 0.81)]

numfilters = filtMissense[(pd.to_numeric(filtMissense["MetaLR_rankscore"]) > 0.81) & (pd.to_numeric(filtMissense["FATHMM_converted_rankscore"]) > 0.81) & (pd.to_numeric(filtMissense["MetaSVM_rankscore"]) > 0.81) & (pd.to_numeric(filtMissense["VEST3_rankscore"]) > 0.81)]

#output these to file
numscores.to_csv("%s_missense_can_1of4_only.txt"%outname, header=True, index=None, sep='\t', mode='a')
numfilters.to_csv("%s_missense_can_all4tools_only.txt"%outname, header=True, index=None, sep='\t', mode='a')


# In[ ]:




