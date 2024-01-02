
"""
Created on Thu May  5 21:43:49 2023

@author: bengi
"""


#get oncogene and tumor suppressor gene list
OGtsgFile = "OncoKB_CancerGeneList_OG_TSG_Label.txt"

import pandas as pd
df_og_tsg = pd.read_csv(OGtsgFile, sep="\t")
df_og_tsg = df_og_tsg[df_og_tsg['OG_TSG'].isin(["OG","TSG"])]
df_og_tsg.columns
og_tsg_dict = dict(zip(df_og_tsg['Hugo Symbol'], df_og_tsg['OG_TSG']))

og_tsg_list = list(set(df_og_tsg['Hugo Symbol'].to_list()))
##############3
#exclude hypermutated samples
HyperMut  = set()
NotHyperMut = set()
with open("Hypermutated_NotHyperLabellingPatient_Mutation_Counts_for_HyperMutated_parameter=8_IQR.txt","r") as infile:
    for line in infile:
        splitted=line.rstrip("\n").split("\t")
        if splitted[-1] == "hypermutated":
            HyperMut.add(splitted[0])
        elif splitted[-1] == "not_hypermutated":
            NotHyperMut.add(splitted[0])
            
        
3519+59048    






##########
AllTumorData = set()
with open("Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt","r") as infile: 
    for line in infile:
        if "Tumor" not in line:
            splitted = line.rstrip("\n").split("\t")
            AllTumorData.add(splitted[2])

AllTumorCount = len(AllTumorData)

    
#############   


#
#Form the potential doublets from the mutations observed on >=3 tumors
#FolderName: /Users/bengi/Desktop/CommsBioRevision_Final_9Dec/REVISION#2/DifferentGene_TCGA_GENIE_Vol13
Gene_Mutation_List={} #each gene keep track of mutations with #tumors>=3 --12,724 key genes with 
#at least one mutation observed on more than 2 tumors
GeneMut_TumorDict={} #{Gene_Mutation:Tumor_set} dictionary, 65,872 mutations in total
with open("Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt","r") as infile:       #open("Merged_Genie_Tcga_VAFgreater_0_25_OG_TSG.txt","r") as infile: 
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
len(GeneMut_TumorDict["CDKN2A_R80"])  
len(GeneMut_TumorDict["NFE2L2_E79"])          
#just change the filtering 
print("GeneMut_tumorDict  completed")
 
"""Case 1"""
#GEQ5 potential doubles
GeneMut_TumorDict_GEQ3 = {}  
for genemut in  GeneMut_TumorDict.keys() :
    gene = genemut.split("_")[0]
    # if gene in og_tsg_list: #if you want to check only hypermutated samples
    #     print(gene)
    if len(GeneMut_TumorDict[genemut].intersection(NotHyperMut))  >= 3 : #should be in at least 3 nonhper tumors
        GeneMut_TumorDict_GEQ3 [genemut] = GeneMut_TumorDict[genemut].intersection(NotHyperMut )
print("muts in at least 3  nonhyper tumors  completed")
len(GeneMut_TumorDict_GEQ3["CDKN2A_R80"])  
len(GeneMut_TumorDict_GEQ3["NFE2L2_E79"])   

"ABL1_R785"
"BRAF_V600"
GeneMut_TumorDict_GEQ3["ABL1_R785"].intersection(GeneMut_TumorDict_GEQ3["NFE2L2_E79"])
GeneMut_TumorDict["ABL1_R785"].intersection(GeneMut_TumorDict["BRAF_V600"])  

GeneMut_TumorDict_GEQ3["KRAS_G12"].intersection(GeneMut_TumorDict_GEQ3["BRAF_V600"]) #only among hypers
GeneMut_TumorDict_GEQ3["KRAS_G12"].intersection(GeneMut_TumorDict_GEQ3["KEAP1_G417"]) #there are two nonhypers, mention in text

###########
# import pickle
# filehandler = open(b"GeneMutation_tumor_dictionary.p","wb")
# pickle.dump(GeneMut_TumorDict_GEQ3,filehandler)

# file = open("GeneMutation_tumor_dictionary.p",'rb')
# object_file = pickle.load(file)


##########

Mutation_List = sorted(list(GeneMut_TumorDict_GEQ3.keys()))  #
#present in at least 3 tumors
#doubles should be present in at least onne tumor
######        
DoubleCount=0 ###134,346,651

Doublets = []  #8189 genes with >=5 mutations
for mut1 in Mutation_List:
    for mut2 in Mutation_List:
        if (mut1 < mut2) and (mut1.split("_")[0]!=mut2.split("_")[0]):
            if GeneMut_TumorDict_GEQ3[mut1].intersection(GeneMut_TumorDict_GEQ3[mut2]) !=  set():
                DoubleCount+=1
                # with open("DifferentGene_Potential_Doublets_Tcga_GenieVol6_Nonhyper_GEQ3.txt","a") as outfile:
                #     outfile.write(mut1+"\t"+mut2+"\n")
print("potential doubles completed")
print(DoubleCount)

"A1BGJ"<"A2BG"



GeneMut_TumorDict_GEQ3 
        
###no filtering at this point
#Obtain the number of gene mutant tumors (mutations observed on >=3 tumors)
Gene_MutantTumors_Dict={}
for genemut in GeneMut_TumorDict_GEQ3.keys():
        gene=genemut.split("_")[0]
        patient=GeneMut_TumorDict_GEQ3[genemut]
        
        if gene not in Gene_MutantTumors_Dict.keys():
            Gene_MutantTumors_Dict[gene]=set()
            Gene_MutantTumors_Dict[gene]=Gene_MutantTumors_Dict[gene].union(patient)
        elif gene in Gene_MutantTumors_Dict.keys():
            Gene_MutantTumors_Dict[gene]=Gene_MutantTumors_Dict[gene].union(patient)
len(Gene_MutantTumors_Dict["NFE2L2"] )           
### 844 genes geq5
###
#write the statistics for the potential doublets observed on >0 tumors without any further filtering
            #after this step read the file and compute p values for 4 different contingency tables

with open("DifferentGeneDoubletsTumorCounts_Tcga_GenieVol6_Nonhyper_GEQ3.txt","a") as outfile:
    outfile.write("Mut1"+"\t"+"Mut2"+"\t"+"DoubleMut#"+"\t"+"OnlyMut1"+"\t"+"OnlyMut2"+"\t"+"AllTumor#"+"\n")


with open("DifferentGene_Potential_Doublets_Tcga_GenieVol6_Nonhyper_GEQ3.txt","r") as infile:
    for line in infile:
        splitted=line.rstrip("\n").split("\t") 
        doublet = (splitted[0],splitted[1])
        
        mut1=doublet[0]
        mut2=doublet[1]
        doubleMutant=GeneMut_TumorDict_GEQ3[mut1].intersection(GeneMut_TumorDict_GEQ3[mut2])
        only1=GeneMut_TumorDict_GEQ3[mut1].difference(GeneMut_TumorDict_GEQ3[mut2])
        only2=GeneMut_TumorDict_GEQ3[mut2].difference(GeneMut_TumorDict_GEQ3[mut1])
        AllMutant= len(NotHyperMut)   #AllTumorCount
        with open("DifferentGeneDoubletsTumorCounts_Tcga_GenieVol6_Nonhyper_GEQ3.txt","a") as outfile:
            outfile.write(mut1+"\t"+mut2+"\t"+str(len(doubleMutant))+"\t"+\
                         str(len(only1))+"\t"+str(len(only2))+"\t"+str(AllMutant)+"\n")

GeneMut_TumorDict.keys()  

print("contingency table file completed")

#gece buradan başla 10 mart
#####STATISTICAL TEST      
import pandas as pd
from scipy.stats import fisher_exact
import numpy as np

df=pd.read_csv("DifferentGeneDoubletsTumorCounts_Tcga_GenieVol6_Nonhyper_GEQ3.txt",sep="\t")
df.columns
#df["Remaining"] = df['AllTumor#']- df['DoubleMut#']-df['OnlyMut1']-df['OnlyMut2']




len(df.index)
#first write header
with open("DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6.txt","a") as outfile:
    outfile.write("Mut1"+"\t"+"Mut2"+"\t"+"DoubleMut#"+"\t"+"OnlyMut1"+"\t"+"OnlyMut2"+"\t"+"RemainingTumor#"+"\t"+\
                  "AllTumor#VAF>"+"\t"+"p4"+"\t"+"OddsR4"+"\n")
#calculate p-values with fisher exact test for different contingency tables
#for i in df.index:
for i in df.index:
    
    a=df.iloc[i]['DoubleMut#']
    b=df.iloc[i]["OnlyMut1"]
    c=df.iloc[i]["OnlyMut2"]
    
    d4=len(NotHyperMut)-(a+b+c)
    
    table4=np.array([[a,b],[c,d4]]) #values indep., evaluated  among all tumors with vaf>0.25
    
    oddsr4, p4 = fisher_exact(table4, alternative='two-sided')
    try:
        with open("DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6.txt","a") as outfile:
            outfile.write(df.iloc[i]['Mut1']+"\t"+df.iloc[i]['Mut2']+"\t"+str(a)+"\t"+str(b)+"\t"+str(c)+"\t"+str(d4)\
                          +"\t"+str(len(NotHyperMut))+"\t"+str(p4)+"\t"+ str(oddsr4) +"\n")
    except IndexError:
        print(i)
print("Fisher completed")        
   
##############
#multiprocessing
"""check case studies"""



# #FDR correction

# ########
"""FDR correction part""" 

import statsmodels.stats.multitest

import pandas as pd




df4 =  pd.read_csv('DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6.txt',sep="\t")
print(1)



P4=df4["p4"].to_list()
print(2)


q4=statsmodels.stats.multitest.multipletests(P4, alpha=0.5, method='fdr_bh', is_sorted=False, returnsorted=False)
reject4=q4[0]    
qval4=q4[1]
print(3)

df4["Reject4"]=reject4
df4["Qval4"]=qval4
      
Outfile_different_FDR = "DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6.txt"
df4.to_csv(Outfile_different_FDR,sep="\t",index=False)
print("FDR completed")
len(df4)
df4.columns




#5214835 pairs evaluated
df5  = df4[df4["Qval4"]<0.3]
df5.to_csv("DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6_FDR_corrected_Q_0_3.txt",sep="\t",index=False)
print("FDR filtered")
len(df5)  #358008

with open('DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6_FDR_corrected_Q_0_3.txt',"r") as infile:
    for line in infile:
        if  ("KRAS_G12" in line) and (("CDKN2A_" in line) or ("SMAD4_" in line) or ("GNAS_" in line) or ("U2AF1_" in line) or ("KEAP1_" in line) or ("NFE2L2_" in line)):
            with open("kras_g12_tcga_genieVol6.txt","a") as outfile:
                outfile.write(line)
                
with open('DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6_FDR_corrected_Q_0_3.txt',"r") as infile:
    for line in infile:
        if  ("ESR1" in line) and ("PIK3CA" in line) :
            with open("PIK3CA_ESR1_tcga_genieVol6.txt","a") as outfile:
                outfile.write(line)
                                

