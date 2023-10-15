#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from Bio import AlignIO


# "We need to return to the very first results of important conserved residues and their identities (when not conserved) and repeat the same analyses for all groups of DCAs and ACAs. A table (in appendix) is needed to report them all and a summary table just to list the residues and their functions and their rate of conservation in our data sets."

# # Metazone ACA

# In[2]:


msafile = "./data/hum_mus_zfish_droso aligned.fasta" # metazoan data
alignment = AlignIO.read(msafile, "fasta")


# In[3]:


# Converting the alignment to a numpy array
align_array = np.array([list(rec) for rec in alignment])
columns_of_interest = [189,272,274,286,299,399] # Defining the columns of interest
id_positions = [] # Initializing empty lists to store the results
letters = []
# Looping over the alignment
for i in range(len(alignment)):
    # Initializing an empty list to store the adjusted positions for each sequence
    lst1 = [alignment[i].id.split("|")[1]]
    # Looping over the columns of interest
    for column in columns_of_interest:
        # Counting the gaps before each column of interest
        count = np.count_nonzero(align_array[i:i+1, 0:column] == "-")
        # Calculating the adjusted position
        position = column - count
        # Appending the position to the list
        lst1.append(position)
    # Appending the list of positions to the final list
    id_positions.append(lst1)
    # Getting the amino acids at the columns of interest
    amino_acids = ([alignment[i].seq[column-1] for column in columns_of_interest])
    # Appending the amino acids to the final list
    letters.append(amino_acids)
# Printing the results as strings
for j in range(len(id_positions)):
    print(f"{id_positions[j][0]}: {', '.join(str(x) for x in id_positions[j][1:])} : {', '.join(letters[j])}")


# In[4]:


counts,values = pd.Series(letters).value_counts().values, pd.Series(letters).value_counts().index
df_results1 = pd.DataFrame(list(zip(values,counts)),columns=["Metazoan ACA","count"])
df_results1


# ## DCAb short 50

# In[5]:


msafile = "./data/DCAb short 50.fasta"
alignment = AlignIO.read(msafile, "fasta")

# Converting the alignment to a numpy array
align_array = np.array([list(rec) for rec in alignment])
columns_of_interest = [109,140,142,146,210,296] # Defining the columns of interest
id_positions = [] # Initializing empty lists to store the results
letters = []
# Looping over the alignment
for i in range(len(alignment)):
    # Initializing an empty list to store the adjusted positions for each sequence
    lst1 = [alignment[i].id.split("|")[1]]
    # Looping over the columns of interest
    for column in columns_of_interest:
        # Counting the gaps before each column of interest
        count = np.count_nonzero(align_array[i:i+1, 0:column] == "-")
        # Calculating the adjusted position
        position = column - count
        # Appending the position to the list
        lst1.append(position)
    # Appending the list of positions to the final list
    id_positions.append(lst1)
    # Getting the amino acids at the columns of interest
    amino_acids = ([alignment[i].seq[column-1] for column in columns_of_interest])
    # Appending the amino acids to the final list
    letters.append(amino_acids)
# Printing the results as strings
for j in range(len(id_positions)):
    print(f"{id_positions[j][0]}: {', '.join(str(x) for x in id_positions[j][1:])} : {', '.join(letters[j])}")


# In[6]:


counts,values = pd.Series(letters).value_counts().values, pd.Series(letters).value_counts().index
df_results2 = pd.DataFrame(list(zip(values,counts)),columns=["Bacterial DCA short","count"])
df_results2


# ## DCAb long 96
# 

# In[7]:


msafile = "./data/DCAb long 96.fasta"
alignment = AlignIO.read(msafile, "fasta")

# Converting the alignment to a numpy array
align_array = np.array([list(rec) for rec in alignment])
columns_of_interest = [109,140,142,146,210,296] # Defining the columns of interest
id_positions = [] # Initializing empty lists to store the results
letters = []
# Looping over the alignment
for i in range(len(alignment)):
    # Initializing an empty list to store the adjusted positions for each sequence
    lst1 = [alignment[i].id.split("|")[1]]
    # Looping over the columns of interest
    for column in columns_of_interest:
        # Counting the gaps before each column of interest
        count = np.count_nonzero(align_array[i:i+1, 0:column] == "-")
        # Calculating the adjusted position
        position = column - count
        # Appending the position to the list
        lst1.append(position)
    # Appending the list of positions to the final list
    id_positions.append(lst1)
    # Getting the amino acids at the columns of interest
    amino_acids = ([alignment[i].seq[column-1] for column in columns_of_interest])
    # Appending the amino acids to the final list
    letters.append(amino_acids)
# Printing the results as strings
for j in range(len(id_positions)):
    print(f"{id_positions[j][0]}: {', '.join(str(x) for x in id_positions[j][1:])} : {', '.join(letters[j])}")


# In[8]:


counts,values = pd.Series(letters).value_counts().values, pd.Series(letters).value_counts().index
df_results3 = pd.DataFrame(list(zip(values,counts)),columns=["Bacterial DCA long","count"])
df_results3


# ## Eukaryotic DCA:

# In[9]:


msafile = "./data/60seq.aln"
alignment = AlignIO.read(msafile, "clustal") 

# Converting the alignment to a numpy array
align_array = np.array([list(rec) for rec in alignment])
columns_of_interest = [408, 450, 452, 456, 566, 687] # Defining the columns of interest
id_positions = [] # Initializing empty lists to store the results
letters = []
# Looping over the alignment
for i in range(len(alignment)):
    # Initializing an empty list to store the adjusted positions for each sequence
    lst1 = [alignment[i].id.split("|")[1]]
    # Looping over the columns of interest
    for column in columns_of_interest:
        # Counting the gaps before each column of interest
        count = np.count_nonzero(align_array[i:i+1, 0:column] == "-")
        # Calculating the adjusted position
        position = column - count
        # Appending the position to the list
        lst1.append(position)
    # Appending the list of positions to the final list
    id_positions.append(lst1)
    # Getting the amino acids at the columns of interest
    amino_acids = ([alignment[i].seq[column-1] for column in columns_of_interest])
    # Appending the amino acids to the final list
    letters.append(amino_acids)
# Printing the results as strings
for j in range(len(id_positions)):
    print(f"{id_positions[j][0]}: {', '.join(str(x) for x in id_positions[j][1:])} : {', '.join(letters[j])}")


# In[10]:


counts,values = pd.Series(letters).value_counts().values, pd.Series(letters).value_counts().index
df_results4 = pd.DataFrame(list(zip(values,counts)),columns=["Eukaryotic DCA","count"])
df_results4


# # Bacterial ACA

# In[11]:


msafile = "./data/UP ACA 994 renamed.fasta" 
alignment = AlignIO.read(msafile, "fasta")

# Converting the alignment to a numpy array
align_array = np.array([list(rec) for rec in alignment])
columns_of_interest = [814,941,943,947,962,1080] # Defining the columns of interest
id_positions = [] # Initializing empty lists to store the results
letters = []
# Looping over the alignment
for i in range(len(alignment)):
    # Initializing an empty list to store the adjusted positions for each sequence
    lst1 = [alignment[i].id]    
    # Looping over the columns of interest
    for column in columns_of_interest:
        # Counting the gaps before each column of interest
        count = np.count_nonzero(align_array[i:i+1, 0:column] == "-")
        # Calculating the adjusted position
        position = column - count
        # Appending the position to the list
        lst1.append(position)
    # Appending the list of positions to the final list
    id_positions.append(lst1)
    # Getting the amino acids at the columns of interest
    amino_acids = ([alignment[i].seq[column-1] for column in columns_of_interest])
    # Appending the amino acids to the final list
    letters.append(amino_acids)
# Printing the results as strings
for j in range(len(id_positions)):
    print(f"{id_positions[j][0]}: {', '.join(str(x) for x in id_positions[j][1:])} : {', '.join(letters[j])}")


# In[12]:


counts,values = pd.Series(letters).value_counts().values, pd.Series(letters).value_counts().index
df_results5 = pd.DataFrame(list(zip(values,counts)),columns=["Bacterial ACA","count"])
df_results5


# In[13]:


# from IPython.display import display_html
# from itertools import chain,cycle

# def display_side_by_side(*args,titles=cycle([''])):
#     html_str=''
#     for df,title in zip(args,chain(titles,cycle(['</br>']))):
#         html_str+='<th style=\"text-align:center\"><td style=\"vertical-align:top\">'
#         html_str+=f'<h2 style=\"text-align: center;\"> {title}</h2>'
#         html_str+=df.to_html().replace('table','table style=\"display:inline\"')
#         html_str+='</td></th>'
#     display_html(html_str,raw=True)
# display_side_by_side(df_results1, df_results2, df_results3, df_results4, df_results5.head(10), titles=['Metazoe', 'Short', 'long', 'Eukaryotic DCA', 'Bacterial ACA'])


# # a summery table

# In[14]:


df_concat = pd.concat([df_results1, df_results2, df_results3, df_results4, df_results5], axis=1)
df_concat.fillna('', inplace=True)
df_concat


# In[15]:


# Create a DataFrame from the two lists
data = {'AC': [id_positions[j][0] for j in range(len(id_positions))],
        'letters': [letters[j] for j in range(len(letters))]}

df = pd.DataFrame(data)

# Define the target letters list
target_letters = ['Q', 'H', 'H', 'E', 'H', 'T']

# Define a custom function to check if a list matches the target
def list_matches_target(lst):
    return lst == target_letters

# Use the apply method to filter rows
df_assumed_active = df[df['letters'].apply(list_matches_target)]
df_assumed_inactive = df[~df['letters'].apply(list_matches_target)] #negation

print(len(df_assumed_active), len(df_assumed_inactive))

