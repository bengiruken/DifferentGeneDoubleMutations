
"""
Created on Tue Mar 21 15:50:11 2023

@author: bengi
"""

import pandas as pd

#get cancer  driver mutation list

df_driver = pd.read_csv("Oncogenic_Driver_Mutations_from_CGI.txt",sep="\t")
df_driver.columns
df_driver["GeneMut"] =  df_driver['gene'] +"_" + df_driver['WT_Residue']

Drivers = set(df_driver["GeneMut"].to_list())

##################    
    

##########

df_doubles = pd.read_csv("DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6_FDR_corrected_Q_0_3.txt",sep="\t")
df_doubles.columns
df_doubles = df_doubles[df_doubles['DoubleMut#'] >= 3]
df_doubles.reset_index(inplace=True,drop=True)
###

import numpy as np
df_doubles["Log2_OR4"]=np.log2(df_doubles["OddsR4"])

df_doubles["Tendency"]=["None"]*len(df_doubles)

for ind in df_doubles.index:
    if df_doubles.loc[ind,"OddsR4"]!=0:
        if df_doubles.loc[ind,"Log2_OR4"]>0:
            df_doubles.at[ind,"Tendency"]="Co-occurence"
        else:
            df_doubles.at[ind,"Tendency"]="Mutual Exclusivity"
    if df_doubles.loc[ind,"OddsR4"]==0:
        df_doubles.at[ind,"Tendency"]="NotApplicable"

###########3
#get oncogene and tumor suppressor gene list
OGtsgFile = "OncoKB_CancerGeneList_OG_TSG_Label.txt"

import pandas as pd
df_og_tsg = pd.read_csv(OGtsgFile, sep="\t")
df_og_tsg = df_og_tsg[df_og_tsg['OG_TSG'].isin(["OG","TSG"])]
df_og_tsg.columns
og_tsg_dict = dict(zip(df_og_tsg['Hugo Symbol'], df_og_tsg['OG_TSG'])) 

############
df_doubles["Function1"]=["None"]*len(df_doubles) 
df_doubles["Function2"]=["None"]*len(df_doubles) 

for ind in df_doubles.index:
    g1 = df_doubles.loc[ind,'Mut1'].split("_")[0]
    g2 = df_doubles.loc[ind,'Mut2'].split("_")[0]
    if g1 in og_tsg_dict.keys():
        function1 = og_tsg_dict[g1]
    else:
        function1 = "Not_OG_TSG"
        
    if g2 in og_tsg_dict.keys():
        function2 = og_tsg_dict[g2]
    else:
        function2 = "Not_OG_TSG"
    
    df_doubles.loc[ind,"Function1"] = function1
    df_doubles.loc[ind,"Function2"] = function2
    


##################    
    
df_doubles["MutType1"]=["None"]*len(df_doubles) 
df_doubles["MutType2"]=["None"]*len(df_doubles)  

for ind in df_doubles.index:
    m1 = df_doubles.loc[ind,'Mut1']
    m2 = df_doubles.loc[ind,'Mut2']
    
    if m1 in Drivers:
        tag1 = "Driver"
    else:
        tag1 = "Passenger"
    df_doubles.loc[ind,"MutType1"] = tag1
    
    if m2 in Drivers:
        tag2 = "Driver"
    else:
        tag2 = "Passenger"
        
    df_doubles.loc[ind,"MutType2"] = tag2
    
    
    
#################3
#figure 1: upset plot for og+og, og+tsg, tsg+tsg doublets, stacked bar plot for  driver+passenger info

"""upset plot"""
df_doubles.columns

df_doubles["Function"] = df_doubles['Function1'] +"+"+df_doubles['Function2'] 
set(df_doubles["Function"].to_list())

for ind in df_doubles.index:
    if df_doubles.loc[ind,"Function"]  =="TSG+OG":
        df_doubles.loc[ind,"Function"] = "OG+TSG"
    elif df_doubles.loc[ind,"Function"]  =="Not_OG_TSG+OG":
        df_doubles.loc[ind,"Function"] = "OG+Not_OG_TSG"
    elif df_doubles.loc[ind,"Function"]  =="Not_OG_TSG+TSG":
        df_doubles.loc[ind,"Function"] = "TSG+Not_OG_TSG"
   
set(df_doubles["Function"].to_list())
#OG/TSG kategorisi olanı almadım        
    
        
df_doubles.columns

df_doubles["MutType"] = df_doubles['MutType1'] +"+"+df_doubles['MutType2'] 

for ind in df_doubles.index:
    if df_doubles.loc[ind,"MutType"]  =="Passenger+Driver":
        df_doubles.loc[ind,"MutType"] = "Driver+Passenger"    
df_doubles["Doublet"] = df_doubles["Mut1"]+"+"+df_doubles["Mut2"]
################
df  = df_doubles[["Function","MutType","Doublet"]]

#df.groupby(["Function","MutType"]).size().unstack().plot(kind='bar', stacked=True)


import matplotlib.pyplot as plt
import seaborn as sns 
from matplotlib.font_manager import FontProperties
# Set the desired pastel color palette
sns.set_palette("muted")
# Set the desired figure size
fig, ax = plt.subplots(figsize=(8, 6))
ax = df.groupby(["Function","MutType"]).size().unstack().plot(kind='bar', stacked=True)
# Set the x-axis and y-axis labels
# Set the Arial font
font = FontProperties(family='Arial')
ax.set_xlabel("Doublet constituent gene  function")
ax.set_ylabel("Number of doublets")

ax.figure.savefig("Different_Gene_Doublet_count_MutationType.svg",dpi=1000,bbox_inches="tight")


##

pct = lambda x: 100 * x / x.sum()

out = df.groupby(["Function","MutType"]).count().groupby('Function').apply(pct)
print(out)



fig1, ax1 = plt.subplots(figsize=(8, 6))
ax1 = out.unstack().plot(kind='bar', stacked=True)

ax1.set_xlabel("Doublet constituent gene  function")
ax1.set_ylabel("Fraction of doublets (%)")

ax1.figure.savefig("Different_Gene_Doublet_Percentage_MutationType.svg",dpi=1000,bbox_inches="tight")


##############
###########
"""figure 2"""
#tendency  versus mutType plot 
"""
df_doubles.columns

df_doubles["Function"] = df_doubles['Function1'] +"+"+df_doubles['Function2'] 

for ind in df_doubles.index:
    if df_doubles.loc[ind,"Function"]  =="TSG+OG":
        df_doubles.loc[ind,"Function"] = "OG+TSG"
        
df_doubles.columns
"""

df_doubles["MutType"] = df_doubles['MutType1'] +"+"+df_doubles['MutType2'] 

for ind in df_doubles.index:
    if df_doubles.loc[ind,"MutType"]  =="Passenger+Driver":
        df_doubles.loc[ind,"MutType"] = "Driver+Passenger"    
df_doubles["Doublet"] = df_doubles["Mut1"]+"+"+df_doubles["Mut2"]
################
df  = df_doubles[["Tendency","MutType","Doublet"]]
df.groupby(["Tendency","MutType"]).size().unstack().plot(kind='bar', stacked=True)


pct = lambda x: 100 * x / x.sum()

out = df.groupby(["Tendency","MutType"]).count().groupby('Tendency').apply(pct)
print(out)

###

fig1, ax1 = plt.subplots(figsize=(8, 6))
ax1 = out.unstack().plot(kind='bar', stacked=True)

ax1.set_xlabel("Doublet constituent gene  function")
ax1.set_ylabel("Fraction of doublets (%)")

ax1.figure.savefig("Different_Gene_Doublet_Percentage_Tendency_vs_MutationType.svg",dpi=1000,bbox_inches="tight")

len(df[df["Tendency"]=="Mutual Exclusivity"])

"""signaling pathway info"""


#kegg 2021 downloaded from enrichR
PATH ="DifferentGeneDoublets_new/"

file= "KEGG_2021_Human.txt"

PATH+file
#maybe we can start with signaling pathways only
Pathway_Gene_dict = {}
Genes = set()
with open(PATH+file,"r") as infile:
    for line in infile:
        #if ("Hepatit" not in line) and ("cancer" not in line) and ("disease" not in line) and ("carcinoma" not in line) and ("infection" not in line) and ("virus" not in line):
        if ("signaling" in line) and ("pathway" in line):
            print(line)
            splitted=line.rstrip("\n").split("\t")
            print(splitted)
            Pathway_Gene_dict[splitted[0]] = set(splitted[1:]).difference(set(['']))
            Genes = Genes.union(set(splitted[1:]).difference(set([''])))
                
            
#old set of genes:
["WNT","TP53","TGF-Beta","RTK RAS", "PI3K", "NRFF2","NOTCH","MYC","HIPPO","Cohesin Complex","Cell cycle"]

PathPAth_intersection = {}
for path in Pathway_Gene_dict.keys():
    PathPAth_intersection[path] = {}
    for path2 in Pathway_Gene_dict.keys():
        if path!=path2:
            PathPAth_intersection[path][path2] = len(Pathway_Gene_dict[path].intersection(Pathway_Gene_dict[path2]))
    
################
Gene_PAthway_Dict = {gene:set() for gene in Genes}   

for gene in  Gene_PAthway_Dict.keys():
    for path in Pathway_Gene_dict.keys():
        if gene in Pathway_Gene_dict[path]:
            Gene_PAthway_Dict[gene].add(path)
            
Gene_PAthway_Dict_Count = {gene:len(Gene_PAthway_Dict[gene]) for gene in Genes}   
Gene_PAthway_Dict["PIK3CA"]
###############


df_doubles["SharePathway"] = [None]*len(df_doubles)
#####
Genes = set()
for ind in df_doubles.index:
    mut1 = df_doubles.loc[ind,"Mut1"]
    gene1 = mut1.split("_")[0]
    mut2 = df_doubles.loc[ind,"Mut2"]
    gene2 = mut2.split("_")[0]
    Genes.add(gene1)
    Genes.add(gene2)
    
    
    if (gene1 in Gene_PAthway_Dict.keys()) and (gene2 in Gene_PAthway_Dict.keys()):
        Intersect = Gene_PAthway_Dict[gene1].intersection(Gene_PAthway_Dict[gene2])
        if Intersect!=set():
            tag5 = "ShareSignalingPath"
        else:
            tag5 = "NotShareSignalingPath"
    else:
        tag5 = "AtLeastOneNoPathway"
    
    
    df_doubles.loc[ind,"SharePathway"] = tag5
    
    
    
df  = df_doubles[["SharePathway","Function","Doublet"]]


fig1, ax1 = plt.subplots(figsize=(8, 6))
ax1 = df.groupby(["Function","SharePathway"]).size().unstack().plot(kind='bar', stacked=True)

ax1.set_xlabel("Doublet constituent gene  function")
ax1.set_ylabel("Number of doublets ")

ax1.figure.savefig("Different_Gene_Doublet_Function_vs_Pathway_count.svg",dpi=1000,bbox_inches="tight")



###

pct = lambda x: 100 * x / x.sum()

out = df.groupby(["Function","SharePathway"]).count().groupby('Function').apply(pct)
print(out)



fig1, ax1 = plt.subplots(figsize=(8, 6))
ax1 = out.unstack().plot(kind='bar', stacked=True)

ax1.set_xlabel("Doublet constituent gene  function")
ax1.set_ylabel("Fraction of doublets (%) ")

ax1.figure.savefig("Different_Gene_Doublet_Function_vs_Pathway_Percentage.svg",dpi=1000,bbox_inches="tight")

###

df_doubles.to_csv("DiffGeneAnnotatedDoubles_Genie_vol6_Tcga.txt",sep="\t",index=False)  














###Draw pie chart 
#double mutatn gene distribution according to function
GenesDict = {}
og_tsg_dict

for ind in df_doubles.index:
    g1 = df_doubles.loc[ind,'Mut1'].split("_")[0]
    g2 = df_doubles.loc[ind,'Mut2'].split("_")[0]
    GenesDict[g1] = og_tsg_dict[g1]
    GenesDict[g2] = og_tsg_dict[g2]

GenesDf = pd.DataFrame({"Gene":GenesDict.keys(), "Function":GenesDict.values()})

GenesDf.columns

GenesDf.groupby(["Function"]).count().unstack().plot(kind='pie')

###############
AllTumors = set()
with open("Merged_Genie_Tcga_VAFgreater_0_25_OG_TSG.txt","r") as infile: 
    for line in infile:
        print(line)
        if "Hugo_"  not in line:
            splitted=line.rstrip("\n").split("\t")
            AllTumors.add(splitted[2])
            
########
Gene_Mutation_List={} #each gene keep track of mutations with #tumors>=3 --12,724 key genes with 
#at least one mutation observed on more than 2 tumors
GeneMut_TumorDict={} #{Gene_Mutation:Tumor_set} dictionary, 65,872 mutations in total
with open("Merged_Genie_Tcga_VAFgreater_0_25_OG_TSG.txt","r") as infile: 
    for line in infile:
        if "Unkonwn" not in line:
            splitted=line.rstrip("\n").split("\t")
            gene=splitted[0]
            mutation=splitted[7]
            genemut=gene+"_"+mutation
            patient=splitted[2]
            
            if gene not in Gene_Mutation_List.keys():
                Gene_Mutation_List[gene]=[]
                Gene_Mutation_List[gene].append(mutation)
            elif gene in Gene_Mutation_List.keys():
                Gene_Mutation_List[gene].append(mutation)
                
            if genemut not in GeneMut_TumorDict.keys():
                GeneMut_TumorDict[genemut] = set()
                GeneMut_TumorDict[genemut].add(patient)
            elif genemut in GeneMut_TumorDict.keys():
                GeneMut_TumorDict[genemut].add(patient)

GeneMut_TumorDict
Gene_Mutation_List
##########

GeneMut_TumorDict_GEQ5 = {}  
for genemut in  GeneMut_TumorDict.keys() :
    if len(GeneMut_TumorDict[genemut])  >= 5 :
        GeneMut_TumorDict_GEQ5 [genemut] = GeneMut_TumorDict[genemut]
##########
#her gendeki double mutant tümörler
Gene_Double_Mutants = {k:set() for k  in GenesDict.keys()}

for ind in df_doubles.index:
    g1 = df_doubles.loc[ind,'Mut1'].split("_")[0]
    g2 = df_doubles.loc[ind,'Mut2'].split("_")[0]
    m1 = df_doubles.loc[ind,'Mut1']
    m2 = df_doubles.loc[ind,'Mut2']
    Gene_Double_Mutants[g1] = Gene_Double_Mutants[g1].union(GeneMut_TumorDict_GEQ5[m1].intersection(GeneMut_TumorDict_GEQ5[m2]))
    Gene_Double_Mutants[g2] = Gene_Double_Mutants[g2].union(GeneMut_TumorDict_GEQ5[m1].intersection(GeneMut_TumorDict_GEQ5[m2]))

Gene_Double_Mutants_GEQ = {}
for g in Gene_Double_Mutants.keys():
    if len(Gene_Double_Mutants[g]) > 50:
        Gene_Double_Mutants_GEQ[g] = len(Gene_Double_Mutants[g])
        
################gene double vs single mutant oranına göre filtrelemek için kullan
Gene_All_Mutants = {k:set() for k  in Gene_Double_Mutants.keys()}

for mut in GeneMut_TumorDict_GEQ5.keys():
    gene = mut.split("_")[0]
    if gene in Gene_Double_Mutants.keys():
        Gene_All_Mutants[gene] = Gene_All_Mutants[gene].union(GeneMut_TumorDict_GEQ5[mut])

GeneDoubleAllMutantFraction = {k:0 for k  in Gene_Double_Mutants.keys() if len(Gene_All_Mutants[k]) >100 }   

for g in  GeneDoubleAllMutantFraction.keys():
        
    GeneDoubleAllMutantFraction[g] = (len(Gene_Double_Mutants[g])/ len(Gene_All_Mutants[g]))*100
    
len(Gene_All_Mutants["PAK1"])    
    
##########
PIK3_SpecificGenePartners = {}    

for ind in df_doubles.index:
    g1 = df_doubles.loc[ind,'Mut1'].split("_")[0]
    g2 = df_doubles.loc[ind,'Mut2'].split("_")[0]
    m1 = df_doubles.loc[ind,'Mut1']
    m2 = df_doubles.loc[ind,'Mut2']
    if (g1 == "KRAS"):
        if g2 not in PIK3_SpecificGenePartners.keys():
            PIK3_SpecificGenePartners[g2] = set()
            
            PIK3_SpecificGenePartners[g2] = PIK3_SpecificGenePartners[g2].union(GeneMut_TumorDict_GEQ5[m1].intersection(GeneMut_TumorDict_GEQ5[m2]))
        elif g2 in PIK3_SpecificGenePartners.keys():
           
            PIK3_SpecificGenePartners[g2] = PIK3_SpecificGenePartners[g2].union(GeneMut_TumorDict_GEQ5[m1].intersection(GeneMut_TumorDict_GEQ5[m2]))
               
    elif (g2 == "KRAS"):
        if g1 not in PIK3_SpecificGenePartners.keys():
            PIK3_SpecificGenePartners[g1] = set()
            
            PIK3_SpecificGenePartners[g1] = PIK3_SpecificGenePartners[g1].union(GeneMut_TumorDict_GEQ5[m1].intersection(GeneMut_TumorDict_GEQ5[m2]))
        elif g1 in PIK3_SpecificGenePartners.keys():
           
            PIK3_SpecificGenePartners[g1] = PIK3_SpecificGenePartners[g1].union(GeneMut_TumorDict_GEQ5[m1].intersection(GeneMut_TumorDict_GEQ5[m2]))
                  
    
            
##############

"""oncoprint input of cbioportal"""

#  PIK3CA partners

for tumor in Gene_All_Mutants["KRAS"]:
    gene="KRAS_all_geq5"
    with open("KRAS_PArtners_DoubleMutantGenes_10DoubleMutant_Oncoprint.txt","a") as outfile:
        outfile.write(tumor+"\t"+gene+"\t"+"HOMDEL"+"\t"+"CNA"+"\n")

            
for gene in PIK3_SpecificGenePartners.keys():
        if len(PIK3_SpecificGenePartners[gene]) >= 10:
            for tumor in PIK3_SpecificGenePartners[gene]:
                with open("KRAS_PArtners_DoubleMutantGenes_10DoubleMutant_Oncoprint.txt","a") as outfile:
                    outfile.write(tumor+"\t"+gene+"\t"+"AMP"+"\t"+"CNA"+"\n")
            
            











#oncoprint for gene couples
import pandas as pd
import seaborn as sns


import numpy as np

import pandas as pd
import seaborn as sns

# Create a sample data frame
df = pd.DataFrame({
    'Gene': ['TP53', 'BRCA1', 'BRCA2', 'KRAS', 'EGFR', 'PTEN'],
    'Tumor Sample': ['Sample1', 'Sample1', 'Sample2', 'Sample2', 'Sample3', 'Sample3'],
    'Tissue': ['Lung', 'Breast', 'Lung', 'Pancreas', 'Brain', 'Ovary']
})

# Define a dictionary to map string values to numeric values
tissue_map = {'Lung': 1, 'Breast': 2, 'Prostate': 3, 'Brain': 4, 'Ovary': 5, 'Skin': 6, 'Pancreas': 7}

# Replace the string values in the data frame with their corresponding numeric values
df['Tissue'] = df['Tissue'].replace(tissue_map)

# Create the pivot table using the mode function as the aggfunc
heatmap_data = pd.pivot_table(df, values='Tissue', index='Gene', columns='Tumor Sample', aggfunc=lambda x: x.mode()[0])

# Generate the categorical heatmap
sns.heatmap(heatmap_data, cmap='viridis', linewidths=0.5)

##############alternative 2
import pandas as pd
import seaborn as sns

# Create a sample data frame
df = pd.DataFrame({
    'Gene': ['TP53', 'BRCA1', 'BRCA2', 'KRAS', 'EGFR', 'PTEN'],
    'Tumor Sample': ['Sample1', 'Sample1', 'Sample2', 'Sample2', 'Sample3', 'Sample3'],
    'Tissue': ['Lung', 'Breast', 'Lung', 'Pancreas', 'Brain', 'Ovary']
})

# Define a dictionary to map string values to numeric values
tissue_map = {'Lung': 1, 'Breast': 2, 'Prostate': 3, 'Brain': 4, 'Ovary': 5, 'Skin': 6, 'Pancreas': 7}

# Replace the string values in the data frame with their corresponding numeric values
df['Tissue'] = df['Tissue'].replace(tissue_map)

# Create the pivot table
heatmap_data = pd.pivot_table(df, values='Tissue', index='Gene', columns='Tumor Sample', aggfunc='first')

# Generate the categorical heatmap
cmap = sns.color_palette('pastel', as_cmap=True)
ax = sns.heatmap(heatmap_data, cmap=cmap, linewidths=0.5)

# Add tick labels to the legend
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, verticalalignment='center')

# Define the legend labels using the tissue_map dictionary
legend_labels = [tissue for tissue in tissue_map.keys()]
legend_values = [tissue_map[tissue] for tissue in tissue_map.keys()]

# Add the legend to the heatmap
cbar = ax.collections[0].colorbar
cbar.set_ticks(legend_values)
cbar.set_ticklabels(legend_labels)

############version 3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Create a sample data frame
df = pd.DataFrame({
    'Gene': ['TP53', 'BRCA1', 'BRCA2', 'KRAS', 'EGFR', 'PTEN'],
    'Tumor Sample': ['Sample1', 'Sample1', 'Sample2', 'Sample2', 'Sample3', 'Sample3'],
    'Tissue': ['Lung', 'Lung', 'Breast', 'Breast', 'Lung', 'Lung']
})

# Define a dictionary to map string values to numeric values
tissue_map = {'Lung': 1, 'Breast': 2, 'Prostate': 3, 'Brain': 4, 'Ovary': 5, 'Skin': 6, 'Pancreas': 7}

# Replace the string values in the data frame with their corresponding numeric values
df['Tissue'] = df['Tissue'].replace(tissue_map)

# Create the pivot table
heatmap_data = pd.pivot_table(df, values='Tissue', index='Gene', columns='Tumor Sample', aggfunc='first')

# Generate the categorical heatmap with smaller rectangular cells and frames
cmap = sns.color_palette('pastel', as_cmap=True)
fig, ax = plt.subplots(figsize=(3, 4))
sns.heatmap(heatmap_data, cmap=cmap, linewidths=1, square=False, cbar_kws={"shrink": 0.6}, 
            annot=True, fmt='.0f', annot_kws={'fontsize': 12}, 
            cbar=True, vmin=1, vmax=7, ax=ax)

# Add tick labels to the legend
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, verticalalignment='center')

# Define the legend labels using the tissue_map dictionary
legend_labels = [tissue for tissue in tissue_map.keys()]
legend_values = [tissue_map[tissue] for tissue in tissue_map.keys()]

# Add the legend to the heatmap
cbar = ax.collections[0].colorbar
cbar.set_ticks(legend_values)
cbar.set_ticklabels(legend_labels)

# Add frames to the heatmap
for _, spine in ax.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(1)

# Adjust the layout to prevent the bottom labels from being cut off
fig.tight_layout()


