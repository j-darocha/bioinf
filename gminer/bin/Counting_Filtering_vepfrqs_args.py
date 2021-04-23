#!/usr/bin/env python3
# coding: utf-8

# In[ ]:


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
parser = argparse.ArgumentParser(description=" A script to filter vep outputs")
# add long and short argument
parser.add_argument("--vepfile", help="Input vep output file")
parser.add_argument("--list", help = "Input list of genes sought")
parser.add_argument("--output", help="Path to folder of the output files ")


# read arguments from the command line
args = parser.parse_args()
if not args.vepfile:  
    print("A vep file is required \n Use --help for more details")
    sys.exit(1)
if not args.list:
    print("A gene list file is required \n Use --help for more details")
    sys.exit(1)
if not args.output:  
    output = './'
else:
    output = args.output

    
    
x = pd.read_csv(args.vepfile, sep = '\t')
gene_lists = pd.read_csv(args.list)

#x = pd.read_csv("/home/jorge/adme_scripts/all_core_mrged_admeonly_proper.txt", sep = '\t')
#gene_lists = pd.read_csv("/home/jorge/adme_scripts/core_extended_list.txt", sep = '\t')



outname = args.output


variant_types = x["Consequence"].value_counts()

#put it here
xnonzero = x[x["MAF"] != 0]
x1p = x[x["MAF"] > 0.01]
x5p = x[x["MAF"] > 0.05]



y = pd.merge(xnonzero, gene_lists, on="SYMBOL", how="inner")
y1p = y[y["MAF"] > 0.01]
y5p = y[y["MAF"] > 0.05]


filt_types = y["Consequence"].value_counts()
y1p_filt_types = y1p["Consequence"].value_counts()
y1p_filt_types = y1p["Consequence"].value_counts()



can = y[y["CANONICAL"] == "YES"]

y1p_can = can[can["MAF"] > 0.01]
y5p_can = can[can["MAF"] > 0.05]

cantypes = can["Consequence"].value_counts()
y1p_cantypes = y1p_can["Consequence"].value_counts()
y5p_cantypes = y5p_can["Consequence"].value_counts()


Missense = can[(can["Consequence"] == "missense_variant") | (can["Consequence"] == "missense_variant,splice_region_variant")]
y1p_Missense = y1p_can[(y1p_can["Consequence"] == "missense_variant") | (y1p_can["Consequence"] == "missense_variant,splice_region_variant") ]
y5p_Missense = y5p_can[(y5p_can["Consequence"] == "missense_variant") | (y5p_can["Consequence"] == "missense_variant,splice_region_variant") ]

#output results to file

#counts of types at levels
cantypes.to_csv("%s_canonical_counts_adme_only.txt"%outname, header=True, index=True, sep='\t', mode='a')
y1p_cantypes.to_csv("%s_canonical_counts_adme_only_1p.txt"%outname, header=True, index=True, sep='\t', mode='a')
y5p_cantypes.to_csv("%s_canonical_counts_adme_only_5p.txt"%outname, header=True, index=True, sep='\t', mode='a')


#canonical variant lists at different levels
can.to_csv("%s_canonical_adme_only.txt"%outname, header=True, index=None, sep='\t', mode='a')
y1p_can.to_csv("%s_canonical_adme_only_1p.txt"%outname, header=True, index=None, sep='\t', mode='a')
y5p_can.to_csv("%s_canonical_adme_only_5p.txt"%outname, header=True, index=None, sep='\t', mode='a')


#canonical missense variant lists at different levels
Missense.to_csv("%s_missense_can_adme_only.txt"%outname, header=True, index=None, sep='\t', mode='a')
y1p_Missense.to_csv("%s_missense_can_adme_only_1p.txt"%outname, header=True, index=None, sep='\t', mode='a')
y5p_Missense.to_csv("%s_missense_can_adme_only_5p.txt"%outname, header=True, index=None, sep='\t', mode='a')


#filter some results for those that score highly
#filtMissense = Missense[Missense["FATHMM_converted_rankscore"] != "-"]
#errorMissense = Missense[Missense["FATHMM_converted_rankscore"] == "-"]
#filtMissense = Missense[(Missense["FATHMM_converted_rankscore"] != "-") | (Missense["MetaLR_rankscore"] != "-") | (Missense["MetaSVM_rankscore"] != "-") | (Missense["VEST3_rankscore"] != "-") ]

errorMissense = Missense[(Missense["FATHMM_converted_rankscore"] == "-") | (Missense["MetaLR_rankscore"] == "-") | (Missense["MetaSVM_rankscore"] == "-") | (Missense["VEST3_rankscore"] == "-") ]

filtMissense1 = Missense[(Missense["FATHMM_converted_rankscore"] != "-")] 
filtMissense2 = filtMissense1[(filtMissense1["MetaLR_rankscore"] != "-")] 
filtMissense3 = filtMissense2[(filtMissense2["MetaSVM_rankscore"] != "-")] 
filtMissense4 = filtMissense3[(filtMissense3["VEST3_rankscore"] != "-")]

numscores = filtMissense4[(pd.to_numeric(filtMissense4["MetaLR_rankscore"]) > 0.81) | (pd.to_numeric(filtMissense4["FATHMM_converted_rankscore"]) > 0.81) | (pd.to_numeric(filtMissense4["MetaSVM_rankscore"]) > 0.81) | (pd.to_numeric(filtMissense4["VEST3_rankscore"]) > 0.81)]
numscores_1p = numscores[numscores["MAF"] > 0.01]
numscores_5p = numscores[numscores["MAF"] > 0.05]


numfilters = filtMissense4[(pd.to_numeric(filtMissense4["MetaLR_rankscore"]) > 0.81) & (pd.to_numeric(filtMissense4["FATHMM_converted_rankscore"]) > 0.81) & (pd.to_numeric(filtMissense4["MetaSVM_rankscore"]) > 0.81) & (pd.to_numeric(filtMissense4["VEST3_rankscore"]) > 0.81)]
numfilters_1p = numfilters[numfilters["MAF"] > 0.01]
numfilters_5p = numfilters[numfilters["MAF"] > 0.05]


#output these to file
numscores.to_csv("%s_missense_can_1of4_only.txt"%outname, header=True, index=None, sep='\t', mode='a')
numscores_1p.to_csv("%s_missense_can_1of4_only_1p.txt"%outname, header=True, index=None, sep='\t', mode='a')
numscores_5p.to_csv("%s_missense_can_1of4_only_5p.txt"%outname, header=True, index=None, sep='\t', mode='a')


numfilters.to_csv("%s_missense_can_all4tools_only.txt"%outname, header=True, index=None, sep='\t', mode='a')
numfilters_1p.to_csv("%s_missense_can_all4tools_only_1p.txt"%outname, header=True, index=None, sep='\t', mode='a')
numfilters_5p.to_csv("%s_missense_can_all4tools_only_5p.txt"%outname, header=True, index=None, sep='\t', mode='a')

#ADME pred model based filters
error_model_Missense = Missense[(Missense["LRT_score"] == "-") | (Missense["MutationAssessor_score_rankscore"] == "-") | (Missense["PROVEAN_converted_rankscore"] == "-") | (Missense["VEST3_rankscore"] == "-") | (Missense["CADD_PHRED"] == "-") ]

model_Missense1 = Missense[(Missense["LRT_score"] != "-")]
model_Missense2 = model_Missense1[(model_Missense1["MutationAssessor_score"] != "-")]
model_Missense3 = model_Missense2[(model_Missense2["PROVEAN_converted_rankscore"] != "-")]
model_Missense4 = model_Missense3[(model_Missense3["VEST3_rankscore"] != "-")]
model_Missense5 = model_Missense4[(model_Missense4["CADD_PHRED"] != "-")]

proper_model_Missense = model_Missense5[(pd.to_numeric(model_Missense5["LRT_score"]) < 0.0025) | (pd.to_numeric(model_Missense5["MutationAssessor_score"]) > 2.05) | (pd.to_numeric(model_Missense5["PROVEAN_converted_rankscore"]) > 0.654) | (pd.to_numeric(model_Missense5["VEST3_rankscore"]) > 0.50489) | (pd.to_numeric(model_Missense5["CADD_PHRED"]) > 19) ]

proper_model_Missense_1p = proper_model_Missense[proper_model_Missense["MAF"] > 0.01]
proper_model_Missense_5p = proper_model_Missense[proper_model_Missense["MAF"] > 0.05]


proper_model_Missense_ALL = model_Missense5[(pd.to_numeric(model_Missense5["LRT_score"]) < 0.0025) & (pd.to_numeric(model_Missense5["MutationAssessor_score"]) > 2.05) & (pd.to_numeric(model_Missense5["PROVEAN_converted_rankscore"]) > 0.654) & (pd.to_numeric(model_Missense5["VEST3_rankscore"]) > 0.50489) & (pd.to_numeric(model_Missense5["CADD_PHRED"]) > 19) ]

proper_model_Missense_ALL_1p = proper_model_Missense_ALL[proper_model_Missense_ALL["MAF"] > 0.01]
proper_model_Missense_ALL_5p = proper_model_Missense_ALL[proper_model_Missense_ALL["MAF"] > 0.05]


#Output Model Values to file
proper_model_Missense.to_csv("%s_ADME_model_Pred_1o5.txt"%outname, header=True, index=None, sep='\t', mode='a')
proper_model_Missense_1p.to_csv("%s_ADME_model_Pred_1o5_1p.txt"%outname, header=True, index=None, sep='\t', mode='a')
proper_model_Missense_5p.to_csv("%s_ADME_model_Pred_1o5_5p.txt"%outname, header=True, index=None, sep='\t', mode='a')

proper_model_Missense_ALL.to_csv("%s_ADME_model_Pred_ALL.txt"%outname, header=True, index=None, sep='\t', mode='a')
proper_model_Missense_ALL_1p.to_csv("%s_ADME_model_Pred_ALL_1p.txt"%outname, header=True, index=None, sep='\t', mode='a')
proper_model_Missense_ALL_5p.to_csv("%s_ADME_model_Pred_ALL_5p.txt"%outname, header=True, index=None, sep='\t', mode='a')

#LoF based filtering

HC_LOF_noncan = y[y["LoF"] == "HC"]
LC_LOF_noncan = y[y["LoF"] == "LC"]

HC_LOF_can = HC_LOF_noncan[HC_LOF_noncan["CANONICAL"] == "YES"]
LC_LOF_can = LC_LOF_noncan[LC_LOF_noncan["CANONICAL"] == "YES"]

HC_LOF_noncan_1p = HC_LOF_noncan[HC_LOF_noncan["MAF"] > 0.01]
HC_LOF_can_1p = HC_LOF_can[HC_LOF_can["MAF"] > 0.01]

#output lof results to file
HC_LOF_noncan.to_csv("%s_HC_LOF_noncan.txt"%outname, header=True, index=None, sep='\t', mode='a')
LC_LOF_noncan.to_csv("%s_LC_LOF_noncan.txt"%outname, header=True, index=None, sep='\t', mode='a')
HC_LOF_can.to_csv("%s_HC_LOF_can.txt"%outname, header=True, index=None, sep='\t', mode='a')
LC_LOF_can.to_csv("%s_LC_LOF_can.txt"%outname, header=True, index=None, sep='\t', mode='a')
HC_LOF_noncan_1p.to_csv("%s_HC_LOF_noncan_1p.txt"%outname, header=True, index=None, sep='\t', mode='a')
HC_LOF_can_1p.to_csv("%s_HC_LOF_can_1p.txt"%outname, header=True, index=None, sep='\t', mode='a')
