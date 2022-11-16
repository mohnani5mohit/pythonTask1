#import libraries

import matplotlib.pyplot as plt
import pandas as pd


# Flush all the FASTA sequences from file to fasta variable and create a list of each sequence without 1st line in list_file variable
fasta = []
list_file = []
with open(r"C:\Users\mohit\OneDrive\Desktop\BB50242_2022_APC1.2_I_cds.fasta.txt") as file_one:
    for line in file_one:
        line = line.strip()
        if not line:
           continue
        if line.startswith(">"):
            active_sequence_name = line[1:]
            if active_sequence_name not in fasta:
                list_file.append(''.join(fasta))
                fasta = []
            continue
        sequence = line
        fasta.append(sequence)

# Flush the last fasta block to the list_file list
if fasta:
    list_file.append(''.join(fasta))
#Declare all the variables that are required
Number_of_seq=0
count_ATG=0
count_nonATG=0
list_leucine=[]
temp = []
total_leucine=0
e=0
ctt=0
ctc=0
cta=0
tta=0
ttg=0
ctg=0

#Counting of ATG and leucine coding codons
for row in list_file[1:]:
    
    t=row[0:3]
    
    if t=="ATG":
        count_ATG=count_ATG+1
    else:
        count_nonATG=count_nonATG+1
        
    list_leucine=' '.join(row[i:i+3] for i in range(0,len(row),3))
    temp = list_leucine.split()
    for j in temp:
        if j == 'CTT':
            ctt=ctt+1
        elif j == 'CTC':
             ctc=ctc+1
        elif j =='CTA':
             cta=cta+1
        elif j == 'TTA':
             tta=tta+1
        elif j == 'TTG':
             ttg=ttg+1
#Count of CUG that is CTG coding for leucine
        elif j=='CTG':
            ctg=ctg+1
    #temp = []        
    Number_of_seq=Number_of_seq+1

#Total number of Leucine coding codons

total_leucine=ctt+ctc+cta+tta+ttg+ctg

#Opening file again to count how many complement sequences are present

with open(r"C:\Users\mohit\OneDrive\Desktop\BB50242_2022_APC1.2_I_cds.fasta.txt") as f:
    contents = f.read()
    count_comp = contents.count("complement")

#Printing all the counts we got
#Result in word file

print("Total number of sequences are ",Number_of_seq)

print("coding genes",Number_of_seq-count_comp)

print("number of complementary genes",count_comp)

print("total number of leucine codons",total_leucine)
print("total number of CUG coding for leucine",ctg)

#Getting total number of CUG (CTG) codons

leucine_codon=total_leucine-e

print("Number of codons coding for leucine except for CUG",leucine_codon)

print("Number of Sequences starting with ATG",count_ATG)
print("Number of Non-ATG sequences",count_nonATG)


#Creating a DataFrame to plot a graph for different codons to coding for leucine
data1={'Codon':['CTT','CTC','CTA','TTA','TTG','CUG'],'Count':[ctt,ctc,cta,tta,ttg,ctg]}

df = pd.DataFrame(data=data1)


df.plot(kind='bar',x='Codon',y='Count')

#display of graph
#Result in word file
plt.show()
