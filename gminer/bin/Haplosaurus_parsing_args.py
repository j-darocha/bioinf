#!/usr/bin/env python3
# coding: utf-8

# In[ ]:


#haplosaurus parsing and counting
import pandas as pd
import re
import argparse
import sys
import warnings

# check if you have python 3
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

# initiate the parser
parser = argparse.ArgumentParser(description=" A script to parse the haplosaurus output and count the number of haplotypes present"
                                )
# add long and short argument
parser.add_argument("--hapsfile", help="Input haplosaurus output")

# read arguments from the command line
args = parser.parse_args()
if not args.hapsfile:  
    print("A haplosaurus file is required \n Use --help for more details")
    sys.exit(1)


x = pd.read_csv(args.hapsfile, sep = '\t', names = ["Transcript", "Nucleotide_Transcript_Changes", "Annot_info","Protien_transcript_changes", "Sift_Polyphen_Pred", "Gap", "RSIDs", "SampleIDs_genotypes"])
num = x.shape[0]
num
coolist = []
size = list(range(0,num))
for l in size:
    coolist.append("Haplonumber_%s" %l)

f = open(args.hapsfile)
name = f.name

outname = re.split("/",name)[-1]    
    
    
x.insert(0, "Haplo", coolist, True)

x_updated = x.replace(to_replace = ":", value = ",", regex = True)
#x_updated

rowlist = []
for index, rows in x_updated.iterrows():
    my_list =[rows.Haplo , rows.SampleIDs_genotypes] 
    rowlist.append(my_list)
    


#rowlist

#more on the right track
counties = 0 
current_haplonumber = 0
het_counts = []
homo_counts = []
haplonames = []
haplocounts = []
haplohomocounts = []
totcounts = []

for item in rowlist:
    current_haplonumber = item[0].split(",")
    wee = item[1].split(",")
    #print(wee)
    counties = 0
    homocount = 0
    for entries in wee:
        if entries == "1":    
            counties = counties + 1
            #print(current_haplonumber)
            #print(counties)
        elif entries == "2":
            homocount = homocount + 1
    het_counts.append([current_haplonumber,counties])
    haplonames.append(current_haplonumber[0])
    haplocounts.append(counties)
    haplohomocounts.append(homocount)
    totcounts.append(counties + (homocount*2))
        #print(counties)    
#print(counties)

haplodict = {"Haplo": haplonames, "Count_Het": haplocounts, "Count_Homo": haplohomocounts, "Total_Alleles": totcounts}
frame_of_the_counts = pd.DataFrame(haplodict)




merged_dataframes = pd.merge(x_updated, frame_of_the_counts, on = "Haplo")
merged_dataframes.to_csv("%s_haplosaurus_parsed.txt"%outname, header=True, index=None, sep='\t', mode='a')




