#!/usr/bin/env python3
# coding: utf-8

# In[ ]:


#vep_freq_merging
import pandas as pd
import re
import argparse
import sys
import warnings

# check if you have python 3
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

# initiate the parser
parser = argparse.ArgumentParser(description=" A script to filter vep outputs")

# add long and short argument
parser.add_argument("--vepfile", help="Input vep output file")
parser.add_argument("--frqfile", help = "Input plink frq file")
parser.add_argument("--output", help="Path to folder of the output files ")



# read arguments from the command line

args = parser.parse_args()
if not args.vepfile:  
    print("A vep file is required \n Use --help for more details")
    sys.exit(1)
if not args.frqfile:
    print("A plink frq file for the vcf \n Use --help for more details")
    sys.exit(1)
if not args.output:  
    output = './'
else:
    output = args.output
    
#x = pd.read_csv(args.vepfile, sep = '\t', skiprows = 292)

##f = open(args.vepfile)
##name = f.name
##outname = re.split("/",name)[-1]

vepx = pd.read_csv(args.vepfile, sep = '\t', skiprows = 88)
freqx = pd.read_csv(args.frqfile, sep = '\s+')



vepfreq = pd.merge(vepx,freqx, left_on= '#Uploaded_variation', right_on='SNP', how='inner')

#output results of the merge to file
vepfreq.to_csv("%s_vepfreqmerged.txt"%args.output, header=True, index = None, sep='\t', mode='a')

