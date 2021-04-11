# -*- coding: utf-8 -*-
"""
Created on Tue Apr  18 20:40:15 2020

@author: haris
"""
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import pandas as pd
import numpy as np 


def unique(list1): 
    x = np.array(list1)
    return np.unique(x)


datasetsequence= pd.read_csv('D:\\Fast University\\Semester2\\Bio-informatics\\project\\protein-data-set\\pdb_data_seq.csv',index_col ="structureId")
datasetclass= pd.read_csv('D:\\Fast University\\Semester2\\Bio-informatics\\project\\protein-data-set\\pdb_data_no_dups.csv',index_col ="structureId")

#Merging Data of both files based on structure id
Finaldataset=pd.merge(datasetsequence,datasetclass,on='structureId',how='inner') 
df_filtered_Protein =Finaldataset[Finaldataset['macromoleculeType_x'] == 'Protein']
df_filtered_Protein=df_filtered_Protein.drop_duplicates(subset="sequence")
df_filtered_Protein = df_filtered_Protein.drop(df_filtered_Protein[df_filtered_Protein['sequence'].str.find('X') !=-1].index)
df_filtered_Protein = df_filtered_Protein.drop(df_filtered_Protein[df_filtered_Protein['sequence'].str.find('Z') !=-1].index)
df_filtered_Protein = df_filtered_Protein.drop(df_filtered_Protein[df_filtered_Protein['sequence'].str.find('B') !=-1].index)
df_filtered_Protein = df_filtered_Protein.drop(df_filtered_Protein[df_filtered_Protein['sequence'].str.find('U') !=-1].index)
df_filtered_Protein = df_filtered_Protein.drop(df_filtered_Protein[df_filtered_Protein['sequence'].str.find('O') !=-1].index)

#df_filtered_Protein = df_filtered_Protein.replace(df_filtered_Protein[df_filtered_Protein['sequence'].str.find('HYDROLASE') !=-1].index,'TRANSFERASE')
df_filtered_Protein.to_csv('D:\\Fast University\\Semester2\\Bio-informatics\\project\\protein-data-set\\FinalDS.csv',index = True, header=True)
dd= pd.read_csv('D:\\Fast University\\Semester2\\Bio-informatics\\project\\protein-data-set\\FinalDS.csv')
df=dd[['sequence','densityPercentSol','densityMatthews','resolution','classification']]
#df=dd[['sequence','classification']]
print(len(df))
################preprocessing of classification labels###############################3
a=df['classification']
TRANSFERASE=[]
HYDROLASE=[]
Other=[]
for i in range(len(a)):
    if a[i].upper().find('HYDROLASE') != -1:
        HYDROLASE.append(a[i])
    elif a[i].upper().find('TRANSFERASE') != -1:
        TRANSFERASE.append(a[i])
    else:
        Other.append(a[i])
    
Other=unique(Other)
TRANSFERASE=unique(TRANSFERASE) 
HYDROLASE=unique(HYDROLASE)     
       
for i in range(len(TRANSFERASE)):
    df=df.replace([TRANSFERASE[i]],['TRANSFERASE']) 
for i in range(len(HYDROLASE)):
    df=df.replace([HYDROLASE[i]],['HYDROLASE'])
for i in range(len(Other)):
    df=df.replace([Other[i]],['Others'])
    
df.to_csv('D:\\Fast University\\Semester2\\Bio-informatics\\project\\protein-data-set\\PreprocessedDataset.csv',index = True, header=True)



#################################Feature extraction############################################################################
#df.loc[:,'Acount']=df['sequence'].str.count('A');
#df.loc[:,'Tcount']=df['sequence'].str.count('T');
#df.loc[:,'Ccount']=df['sequence'].str.count('C');
#df.loc[:,'Gcount']=df['sequence'].str.count('G');
#df.loc[:,'Rcount']=df['sequence'].str.count('R');
#df.loc[:,'Ncount']=df['sequence'].str.count('N');
#df.loc[:,'Dcount']=df['sequence'].str.count('D');
#df.loc[:,'Qcount']=df['sequence'].str.count('Q');
#df.loc[:,'Ecount']=df['sequence'].str.count('E');
#df.loc[:,'Hcount']=df['sequence'].str.count('H');
#df.loc[:,'Icount']=df['sequence'].str.count('I');
#df.loc[:,'Lcount']=df['sequence'].str.count('L');
#df.loc[:,'Kcount']=df['sequence'].str.count('K');
#df.loc[:,'Mcount']=df['sequence'].str.count('M');
#df.loc[:,'Fcount']=df['sequence'].str.count('F');
#df.loc[:,'Pcount']=df['sequence'].str.count('P');
#df.loc[:,'Scount']=df['sequence'].str.count('S');
#df.loc[:,'Tcount']=df['sequence'].str.count('T');
#df.loc[:,'Wcount']=df['sequence'].str.count('W');
#df.loc[:,'Vcount']=df['sequence'].str.count('V');

    ###df.loc[i,'instability_index']=analysed_seq.instability_index()
    ##df.loc[i,'flexibility']=analysed_seq.flexibility()

    
    
    
for i, row in df.iterrows():
    analysed_seq=ProteinAnalysis(row['sequence'])
    try:
        df.loc[i,'instability_index']=analysed_seq.instability_index()
        df.loc[i,'Molecularweight']=analysed_seq.molecular_weight()
        df.loc[i,'aromaticity']=analysed_seq.aromaticity()
        df.loc[i,'isoelectric_point']=analysed_seq.isoelectric_point()
        df.loc[i,'Helix']=analysed_seq.secondary_structure_fraction()[0]  
        df.loc[i,'Turn']=analysed_seq.secondary_structure_fraction()[1]
        df.loc[i,'Sheet']=analysed_seq.secondary_structure_fraction()[2]
        df.loc[i,'Amino_Acid_A_percentage']=analysed_seq.get_amino_acids_percent()['A']
        df.loc[i,'Amino_Acid_C_percentage']=analysed_seq.get_amino_acids_percent()['C']
        df.loc[i,'Amino_Acid_D_percentage']=analysed_seq.get_amino_acids_percent()['D']
        df.loc[i,'Amino_Acid_E_percentage']=analysed_seq.get_amino_acids_percent()['E']
        df.loc[i,'Amino_Acid_F_percentage']=analysed_seq.get_amino_acids_percent()['F']
        df.loc[i,'Amino_Acid_G_percentage']=analysed_seq.get_amino_acids_percent()['G']
        df.loc[i,'Amino_Acid_H_percentage']=analysed_seq.get_amino_acids_percent()['H']
        df.loc[i,'Amino_Acid_I_percentage']=analysed_seq.get_amino_acids_percent()['I']
        df.loc[i,'Amino_Acid_K_percentage']=analysed_seq.get_amino_acids_percent()['K']
        df.loc[i,'Amino_Acid_L_percentage']=analysed_seq.get_amino_acids_percent()['L']
        df.loc[i,'Amino_Acid_M_percentage']=analysed_seq.get_amino_acids_percent()['M']
        df.loc[i,'Amino_Acid_N_percentage']=analysed_seq.get_amino_acids_percent()['N']
        df.loc[i,'Amino_Acid_P_percentage']=analysed_seq.get_amino_acids_percent()['P']
        df.loc[i,'Amino_Acid_Q_percentage']=analysed_seq.get_amino_acids_percent()['Q']
        df.loc[i,'Amino_Acid_R_percentage']=analysed_seq.get_amino_acids_percent()['R']
        df.loc[i,'Amino_Acid_S_percentage']=analysed_seq.get_amino_acids_percent()['S']
        df.loc[i,'Amino_Acid_T_percentage']=analysed_seq.get_amino_acids_percent()['T']
        df.loc[i,'Amino_Acid_V_percentage']=analysed_seq.get_amino_acids_percent()['V']
        df.loc[i,'Amino_Acid_W_percentage']=analysed_seq.get_amino_acids_percent()['W']
        df.loc[i,'Amino_Acid_Y_percentage']=analysed_seq.get_amino_acids_percent()['Y']
    except:
        pass
print(df)  

df.to_csv('D:\\Fast University\\Semester2\\Bio-informatics\\project\\protein-data-set\\FeatureExtracted.csv', header=True)



#my_seq = "LEVLLGSGDGSLVFVPSEFSVPSGEKIVFKNNAGFPHNVVFDEDEIPAGVDAVKISMPEEELLNAPGETYVVTLDTKGTYSFYCSPHQGAGMVGKVTVN"
#analysed_seq = ProteinAnalysis(my_seq)
##print(analysed_seq.molecular_weight())
##print(analysed_seq.gravy())
##print(analysed_seq.count_amino_acids())
##print(analysed_seq.get_amino_acids_percent())
##print(analysed_seq.aromaticity())
##print(analysed_seq.instability_index())
##print(analysed_seq.flexibility())
##print(analysed_seq.isoelectric_point())
#print(analysed_seq.secondary_structure_fraction())
#from bs4 import BeautifulSoup
#import urllib3
#http = urllib3.PoolManager()
#response = http.request('GET', 'https://www.ebi.ac.uk/Tools/services/rest/emboss_pepstats/result/emboss_pepstats-E20200407-161911-0574-90494793-p2m/out')
#print(response.status)
##print(response.data)
#soup = BeautifulSoup(response.data, 'html.parser')
#print(soup)
#
#Html_file= open("td","w")
#Html_file.write(str(soup))
#Html_file.close()
#pre = soup.findall('td')
#print(pre)
#html_page = urllib3.urlopen("https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup")
#soup = BeautifulSoup(html_page, 'html.parser')
#print(soup.prettify())
