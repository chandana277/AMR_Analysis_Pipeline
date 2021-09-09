# Imports #

import numpy as np
import pandas as pd
import os
import glob
from Bio import SeqIO
import itertools
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import sys
import warnings



# Collect the path for .rgi files from user and load
path = sys.argv[1]
readfiles=glob.glob(os.path.join(path,"*.txt"))


# Collect the path for .gbk files from user and load
path = sys.argv[2]
gbkfiles=glob.glob(os.path.join(path,"*.gbk"))


################ Start of all defined functions ####################


# Function to extract the required data from gbk files(entire genomes:) using biopython library #  
def extract(infile):
  gene_start=[]
  gene_end=[]
  gene_strand=[]
  gene_name=[]
  loc_tag=[]
  function=[]
  protein_seq=[]
  contig_name=[]
  unique=[]

  for index, record in enumerate(SeqIO.parse(infile, "genbank")):
    #print("index %i, ID = %s, length %i, with %i features"% (index, record.id, len(record.seq), len(record.features)))
    for i in record.features:
        if i.type == "CDS" and "gene" in i.qualifiers:
          locations=i.location
          gene_start.append(locations.start)
          gene_end.append(locations.end)
          gene_strand.append(locations.strand)
          loc_tag.append(i.qualifiers['locus_tag'])
          function.append(i.qualifiers['product'])
          protein_seq.append(str(i.qualifiers['translation']))
          gene_name.append(i.qualifiers['gene'])
          contig_name.append(record.id)
        elif i.type =="CDS":
          locations=i.location
          gene_start.append(locations.start)
          gene_end.append(locations.end)
          gene_strand.append(locations.strand)
          loc_tag.append(i.qualifiers['locus_tag'])
          function.append(i.qualifiers['product'])
          protein_seq.append(str(i.qualifiers['translation']))
          gene_name.append("UID")
          contig_name.append(record.id)
          
  salmonella_gene_frame=pd.DataFrame()
  salmonella_gene_frame['GeneStart']=gene_start
  salmonella_gene_frame['GeneEnd']=gene_end
  salmonella_gene_frame['GeneStrand']=gene_strand
  salmonella_gene_frame['Locus_Tag']=loc_tag
  salmonella_gene_frame['GeneName']=gene_name
  salmonella_gene_frame['Product']=function
  salmonella_gene_frame['ProteinSequence']=protein_seq
  salmonella_gene_frame['contig_name']=contig_name
  
  #print(contig_name)
  for i in contig_name:
    if i not in unique:
      unique.append(i)


  return salmonella_gene_frame,unique


# Function to create groups based on contigs #
def make_groups(frame,column):
  group=frame.groupby(frame[column])
  datasets = {}  
  for groups, data in group:
    datasets[groups] = data
  return datasets

#Creating separate dictionary for all the unique AMR genes of the genome as "key" using the field "Best_Hit_ARO" #
def createeachdict_drug(drugindex):
	temp_dict={}
	for j,k in datadict.items():       
	    temp=k[k['Best_Hit_ARO']==drugindex]
	    if len(temp)>0:
	        temp_dict[j]=temp
	return temp_dict


# A function to assign color to all and any AMR genes in the neighborhood based on CARD color codes for strict, loose and perfect hits #
color_Dict={"Loose":"#DC7633 ","Perfect":"#28B463","Strict":"#F4D03F"}
def checkforRGIinneighborhood(e,genome,rgi_gene):
    for q in range(len(e["GeneEnd"])):
        for p in range(len(datadict[genome]["Stop"])):
            datadict[genome].reset_index(drop=True, inplace=True)
            #print(e["GeneEnd"][q],p)
            if e["GeneEnd"][q] == datadict[genome]["Stop"][p] and e["GeneName"][q]==rgi_gene:
                e["Genecolor"][q]=color_Dict[datadict[genome]["Cut_Off"][p]]
                                    
    return e

# A function to delete keys of those AMR genes who are not present in atleast 25% of the total genomes ##
def delete_keys_lessthan_25percent_instances(dict_element):
    total_genomes= len(gbkfiles)
    minimum_genomes = int((total_genomes /25)*100)
    emptykeyslist=[]
    for i,j in dict_element.items():
        if len(j)<=minimum_genomes:
            emptykeyslist.append(i)
    for i in emptykeyslist:
        del dict_element[i]
    return dict_element

# A function to compute neighbors of each AMR gene based on contig info ##
def find_neighbor(rec,uname,data,number_of_genes,drug,genome,key,instancetype):
      
      neighbor_genes=[]
      rec.reset_index(drop=True, inplace=True)
      contig_flag=0

      for j in range(len(rec['Start'])):
            m=[]
            n=[]
            upwardgenes=[]
            downwardgenes=[]
            recarray=[]
            i=rec.loc[j].Start
            w=rec.loc[j].Stop
            k=rec.loc[j].req_cont
            g=number_of_genes


            if k in uname:
                newlist=data[k]
                newlist.reset_index(drop=True, inplace=True)
                for l in range(len(newlist)):
                    if newlist['GeneStart'][l] > i and newlist['GeneEnd'][l] > w:
                        downwardgenes.append((newlist['GeneStart'][l],newlist['GeneEnd'][l],newlist['GeneStrand'][l],newlist['Locus_Tag'][l],newlist['ProteinSequence'][l],newlist['GeneName'][l],"#A9F1EE","Not_Applicable"))
                    else:
                        upwardgenes.append((newlist['GeneStart'][l],newlist['GeneEnd'][l],newlist['GeneStrand'][l],newlist['Locus_Tag'][l],newlist['ProteinSequence'][l],newlist['GeneName'][l],"#A9F1EE","Not_Applicable"))
                    #print(upwardgenes) 

                newu=pd.DataFrame(upwardgenes,columns=["GeneStart","GeneEnd","Strand","Locus_Tag","ProteinSequence","GeneName","Genecolor","Gene_Cut_Off"])
                m=(newu.iloc[(newu['GeneStart']-i).abs().argsort()[:g+1]]).sort_values(by="GeneStart")
                #print(m)
                newd=pd.DataFrame(downwardgenes,columns=["GeneStart","GeneEnd","Strand","Locus_Tag","ProteinSequence","GeneName","Genecolor","Gene_Cut_Off"])
                n=(newd.iloc[(newd['GeneStart']-i).abs().argsort()[:g]]).sort_values(by="GeneStart")
                #print(n)


                recarray.append((rec['Start'][j],rec['Stop'][j],rec['Orientation'][j],rec['Locus_Tag'][j],rec['Predicted_Protein'][j],rec['Best_Hit_ARO'][j],"#ccccff",rec["Cut_Off"][j]))
                o=pd.DataFrame(recarray,columns=["GeneStart","GeneEnd","Strand","Locus_Tag","ProteinSequence","GeneName","Genecolor","Gene_Cut_Off"])


            else:
                print(k,genome)
                print("contig does not exist----"+"drugclass:"+drug)
                print(rec['Best_Hit_ARO'][j])
      
      
      if len(m)!=0 and len(n)==0:
          m.reset_index(drop=True, inplace=True)
          m=m.drop([len(m)-1])
          e = pd.concat([m,o,n], ignore_index=True)
          e.reset_index(drop=True, inplace=True)
          neighbor_genes.append(e)
              
      if len(m)==0 and len(n)!=0:

          m.reset_index(drop=True, inplace=True)
          n=n.drop([len(n)-1])
          e = pd.concat([m,o,n], ignore_index=True)
          e.reset_index(drop=True, inplace=True)
          neighbor_genes.append(e)
         

      if(len(m)!=0 and len(n)!=0):
        m.reset_index(drop=True, inplace=True)
        m=m.drop([len(m)-1])
        e = pd.concat([m,o,n], ignore_index=True)
        e.reset_index(drop=True, inplace=True)
        neighbor_genes.append(e)
            
      return neighbor_genes,contig_flag


# A function to get separate dicts of locus tags,protein sequences and gene names inorder to compare and write to a fasta file####
def getrequiredgenes(frame,number_of_genes,drug):
    temp_locus_array=[]
    temp_protein_array=[]
    temp_genename=[]
    #print(drug)
    for i in frame:
        temp_locus_array.append(list(i['Locus_Tag']))
        temp_protein_array.append(i['ProteinSequence'])
        temp_genename.append((i['GeneName']))
    
    
    locus_to_protein_dict={}
    for i in frame:
        for j in range(len(i)):
            
            locus_to_protein_dict[i['Locus_Tag'][j].strip()]=i['ProteinSequence'][j]
  
    a=temp_locus_array
    a=list(itertools.chain.from_iterable(a))
   
    b=temp_protein_array
    b=list(itertools.chain.from_iterable(b))
   
    c=temp_genename
    c=list(itertools.chain.from_iterable(c))
    
 
    if number_of_genes==10:
        if len(frame)>1:
            return(a[:10]+a[-10:], b[:10]+b[-10:],c[:10]+c[-10:])## for more than one card genes in a genomes
        else:
            return a,b,c
        #return(a[:5]+a[-5:], b[:5]+b[-5:],c[:5]+c[-5:])## for more than one card genes in a genomes
    if number_of_genes==14:
        
        if len(frame)>1:
            print(drug)
            return(a[:14]+a[-14:], b[:14]+b[-14:],c[:14]+c[-14:])## for more than one card genes in a genomes
        else:
            return a,b,c


# A function to manipulate contig to compare whether the RGI and GBK genes belong to same contig#####

def locus_generator(frame,genome_name):
    tag=[]
    RGI_name_array=[]
        
    frame.reset_index(drop=True, inplace=True)
    for index in range(len(frame)):
        temp_name=frame["Best_Hit_ARO"][index].split(" ")
        if len(temp_name)>1:
            name=(temp_name[0][0]+temp_name[1][0]+"_"+temp_name[2])
            RGI_name_array.append(name)
        else:
            name=(temp_name[0])
            RGI_name_array.append(name)
            
        tag.append(genome_name+"("+name+")"+frame["Cut_Off"][index][0]+"_"+str(frame["Best_Identities"][index]))
        
    return tag,RGI_name_array


########### End of all defined functions ######################


# Read each rgi file(.txt) using pandas and store them in a dataframe#
dataframelist=[]
filenames=[]
datadict={}

for i in sorted(readfiles):   
    filenames.append(os.path.basename(i).split(".")[0])
    dataframelist.append(pd.DataFrame(pd.read_csv(i,sep="\t")))

# Dictionary of dataframes with keys as filenames ###
for i in range(len(filenames)):
    datadict[filenames[i]]=pd.DataFrame(dataframelist[i])  

# Replacing the orientation of rgi dataframe to match the data in .gbk format to help with comparison #
for i,j in datadict.items():
    j["Orientation"]=j["Orientation"].replace("-","-1")
    j["Orientation"]=j["Orientation"].replace("+","+1")

# Adding two new columns 1> Modified Locus tag and 2> contig to extract the right neighbors #
for k,v in datadict.items():
    newcon=[]
    for i in v["Contig"]:
        head,sep,tail=i.partition("_")
        newcon.append(head)
    v['req_cont']=newcon    #adding a new column "req_cont" into dataframe
    v['Locus_Tag'],v["Best_Hit_ARO"]=locus_generator(v,k)


# Read the .gnk file data, store the info in the sorted order in dictionaries #
gbk_names=[]
uniquenames=[]
datasetslist=[]
gbkdict={}
uniquedict={}

for i in sorted(gbkfiles):    
    gbk_names.append(os.path.basename(i).split(".")[0])

for i in gbk_names:
    uniquenames.append(str(i))
    datasetslist.append(str(i))
    
gbkfiles=sorted(gbkfiles) #read the files in the sorted order

for i in range(len(gbk_names)):
    a=gbk_names[i]
    b=uniquenames[i]
    gbkdict[a],uniquedict[b]=extract(gbkfiles[i])

# Convert the Locus tag column from  ['SA200015_0001'] --> SA200015_001 to help while writing to fasta file required to BLAST #
for j,i in gbkdict.items():
    i['Locus_Tag']=i['Locus_Tag'].apply(lambda i:str(i).replace("[","").replace("]","").replace("'",""))
    i['ProteinSequence']=i['ProteinSequence'].str.strip('[]')

# To supress the warning that arise when a single value of dataframe column is manipulated #
def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

pd.set_option('mode.chained_assignment', None)
# The block that causes the chained assignment error #

for i,j in gbkdict.items():
    for value in range(len(j["GeneName"])):
        if isinstance(j["GeneName"].iloc[value], list):
            j["GeneName"].iloc[value] = j["GeneName"].iloc[value][0]

# For each .gbk file, divides the data based on available contigs using groupby and dictionary, used to find the neighbors of same contig# 
datasetdict={}
for i in range(len(gbk_names)):
    for j,k in gbkdict.items():
        datasetdict[j]=make_groups(k,"contig_name") 


# Finds the union of all the AMR genes present in more than 25 percent of the total genomes #
uniquedrugdict={}
unionofdrugclasess=[]

for j,k in datadict.items():
    uniquedrugclasses=[]
    for l in range(len(k)):
        if k['Best_Hit_ARO'][l]  not in uniquedrugclasses:
                uniquedrugclasses.append(k["Best_Hit_ARO"][l])
    uniquedrugdict[j]=uniquedrugclasses
        
for i,j in uniquedrugdict.items():
    for item in j:
        if item not in unionofdrugclasess:
            unionofdrugclasess.append(item)

# Converts the longer AMR gene names into shorter version to make them as dictionary keys and to make it easier to visualize in gene order image#
# Haemophilis influenzea PBP3 --> Hi_PBP3#

listofdrugnames_modified=[]
for k in unionofdrugclasess:
    if len(k.split(" "))>1:
        temp=k.split(" ")
        listofdrugnames_modified.append(temp[0][0]+temp[1][0]+"_"+temp[2])
    else:

        listofdrugnames_modified.append(k.split("; ")[0].split(" ")[0])


# A main dictionary of all the AMR models with all the genomes 
main_dictionary={}
for i in range(len(listofdrugnames_modified)):
    main_dictionary[listofdrugnames_modified[i]]=createeachdict_drug(unionofdrugclasess[i])
#main_dictionary=delete_keys_lessthan_25percent_instances(main_dictionary)
 

# Divide the AMR genes into single instance and multiple instance #            
Dict_multigene_instances={}
Dict_singlegene_instances={}
for i,j in main_dictionary.items():
    temp={}
    flag=0
    for a,b in j.items():
        
        if len(b)>1:
            flag=1
    if flag==1:
        Dict_multigene_instances[i]=j
    else:
        Dict_singlegene_instances[i]=j


# Code to get one to one genome comparisons when there are multiple instances and divide the frames based on locus tags  #
dict_rgi_multiple_occurance={}
g={}
for i,j in Dict_multigene_instances.items():
    temp={}
    l=[]
    deletekeylist=[]
    for a,b in j.items():
        if len(b)>1 or len(b)==1:
            l.append(a)
            b.reset_index(drop=True, inplace=True)            
            x=make_groups(b,"Locus_Tag")
            for e,f in x.items():
                temp[e]=f

    dict_rgi_multiple_occurance[i]=temp
    g[i]=l
 
for i,j in Dict_multigene_instances.items():
    for t in g[i]:
        del j[t]
        
for i,j in dict_rgi_multiple_occurance.items():    
     j.update(Dict_multigene_instances[i])


dict_rgi_single_occurance={}
g={}
for i,j in Dict_singlegene_instances.items():
    temp={}
    l=[]
    deletekeylist=[]
    for a,b in j.items():
        if len(b)==1:
            l.append(a)
            b.reset_index(drop=True, inplace=True)            
            x=make_groups(b,"Locus_Tag")
            for e,f in x.items():
                temp[e]=f

    dict_rgi_single_occurance[i]=temp
    g[i]=l
 
for i,j in Dict_singlegene_instances.items():
    for t in g[i]:
        del j[t]
        
for i,j in dict_rgi_single_occurance.items():    
     j.update(Dict_singlegene_instances[i])




# Create separate single instance multiple instance folders
# Find the neighbors and store them in separate dictionaries

singleinstance_neighboringdict_combined_range10={}
Single_contig_end_flag_dict={}
for i,j in dict_rgi_single_occurance.items():
    Dict_neighboring_genes_range10={}
    temp_contig_flag={}
    for k,l in j.items():
        key=k.split("(")[0]
        temp,contigflag=find_neighbor(j[k],uniquedict[key],datasetdict[key],10,i,k,key,"Single_Instance_Neighborhood")
        if len(temp)>0:
            Dict_neighboring_genes_range10[k]=temp
            temp_contig_flag[k]=contigflag
    singleinstance_neighboringdict_combined_range10[i]=Dict_neighboring_genes_range10
    Single_contig_end_flag_dict[i]=temp_contig_flag


multipleinstance_neighboringdict_combined_range_10={}
multiple_contig_end_flag_dict={}
for i,j in dict_rgi_multiple_occurance.items():
    Dict_neighboring_genes_range10={}
    temp_contig_flag={}
    for k,l in j.items():
        key=k.split("(")[0]
        temp,contigflag=find_neighbor(j[k],uniquedict[key],datasetdict[key],10,i,k,key,"Multiple_Instance_Neighborhood")
        if len(temp)>0:
            Dict_neighboring_genes_range10[key]=temp
            temp_contig_flag[key]=contigflag
    multipleinstance_neighboringdict_combined_range_10[i]=Dict_neighboring_genes_range10
    multiple_contig_end_flag_dict[i]=temp_contig_flag
    



# Code to extract the locus tags, protein sequences and gene names to use to later for comparison #
Drug_singleinstance_range_10_locus_dict10={}
Drug_singleinstance_range_10_protein_dict10={}
Drug_singleinstance_range_10_genename_dict10={}

for i,j in singleinstance_neighboringdict_combined_range10.items():
    locustags_dict_10={}
    protein_dict_10={}
    genename_dict_10={}
    for k,l in j.items():
        locustags_dict_10[k],protein_dict_10[k],genename_dict_10[k]=getrequiredgenes(j[k],10,i)  
    Drug_singleinstance_range_10_locus_dict10[i]=locustags_dict_10
    Drug_singleinstance_range_10_protein_dict10[i]=protein_dict_10
    Drug_singleinstance_range_10_genename_dict10[i]=genename_dict_10


Drug_multiple_instance_range_10_locus_dict10={}
Drug_multiple_instance_range_10_protein_dict10={}
Drug_multiple_instance_range_10_genename_dict10={}

for i,j in multipleinstance_neighboringdict_combined_range_10.items():

    locustags_dict_10={}
    protein_dict_10={}
    genename_dict_10={}
    for k,l in j.items():

        locustags_dict_10[k],protein_dict_10[k],genename_dict_10[k]=getrequiredgenes(j[k],10,i)
    #print(len(locustags_dict_10),len(protein_dict_10),len(genename_dict_10))  
    Drug_multiple_instance_range_10_locus_dict10[i]=locustags_dict_10
    Drug_multiple_instance_range_10_protein_dict10[i]=protein_dict_10
    Drug_multiple_instance_range_10_genename_dict10[i]=genename_dict_10
        



######################## Part 1 ends ##################################

# The BLAST results are generated and stored in .txt files #
# The second part processess the BLAST results to generate various clusters #

import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import average,complete,single,weighted,centroid
import plotly.figure_factory as ff
import scipy.spatial as scs
from scipy.spatial.distance import pdist,squareform
from numpy import savetxt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import plotly.figure_factory as ff
import scipy.cluster.hierarchy as sch
import matplotlib
import seaborn as sns; sns.set()
from scipy.cluster.hierarchy import dendrogram, fcluster, leaves_list,linkage
from scipy.spatial import distance
from scipy.cluster import hierarchy# You can use SciPy one too
import os


# Read the .txt files of BLAST results and store in a dictionary ######
filepath_single=sys.argv[3]
Single_instance_blastdataframelist=[]
Single_instance_blastfilename=[]

readfiles=glob.glob(os.path.join(filepath_single,"*.txt"))
for i in readfiles:    
    Single_instance_blastfilename.append(os.path.basename(i).split(".")[0])
    i=pd.read_csv(i,index_col = False,sep="\t", names=['query_id','sub_id','PI','len','Evalue','bitscore','qseq','sseq'])
    Single_instance_blastdataframelist.append(i)

temp_Single_instance_blastdatadict={}### dictionary of dataframes with keys as filenames ###
for i in range(len(Single_instance_blastfilename)):
    a=Single_instance_blastdataframelist[i]
    b=Single_instance_blastfilename[i]
    temp_Single_instance_blastdatadict[b]=pd.DataFrame(a)


filepath_multiple=sys.argv[4]
Multiple_instance_blastdataframelist=[]
Multiple_instance_blastfilename=[]

readfiles=glob.glob(os.path.join(filepath_multiple,"*.txt"))
for i in readfiles:
    Multiple_instance_blastfilename.append(os.path.basename(i).split(".")[0])
    i=pd.read_csv(i,index_col = False,sep="\t", names=['query_id','sub_id','PI','len','Evalue','bitscore','qseq','sseq'])
    #i=pd.read_csv(i,index_col = False,sep="\t")
    Multiple_instance_blastdataframelist.append(i)

temp_Multiple_instance_blastdatadict={}### dictionary of dataframes with keys as filenames ###
for i in range(len(Multiple_instance_blastfilename)):
    a=Multiple_instance_blastdataframelist[i]
    b=Multiple_instance_blastfilename[i]
    temp_Multiple_instance_blastdatadict[b]=pd.DataFrame(a)


# To sort and filter only those BLAST results with more than 70 percent identity 

Single_instance_blastdatadict={}
for i,j in temp_Single_instance_blastdatadict.items():
    j=j[j["PI"]>70]
    Single_instance_blastdatadict[i]=j

# A function to normalize the bitscore values 
def normalized_bitscore(z): 
  manval=[]
  a=z[z['query_id']==z['sub_id']]
  manval=[]
  a=z[z['query_id']==z['sub_id']]
  a.reset_index(drop=True, inplace=True)
  q={}
  for i in range(len(a)):
    q[a['query_id'][i]]=a['bitscore'][i]

  for i in range(len(z)):
      manval.append(round(float(z['bitscore'][i])/q[z['query_id'][i]],3))

  z['normalized_bitscore']=manval
  return z



Single_instance_drugclass_genomedict={}
for k,l in Drug_singleinstance_range_10_locus_dict10.items():
    u=Single_instance_blastdatadict[k]
    temp_dict={}
    for a,b in l.items():        
        temp=[]
        u.reset_index(drop=True, inplace=True)
        for index in range (len(u)):
            if  u["query_id"][index] in b:
                temp.append((u["query_id"][index],u["sub_id"][index],u["PI"][index],u["bitscore"][index]))
                    
        #frame=pd.DataFrame(temp,columns=["query_id","sub_id","PI","bitscore"])
        temp_dict[a]=pd.DataFrame(temp,columns=["query_id","sub_id","PI","bitscore"])
      
    Single_instance_drugclass_genomedict[k]=temp_dict
    

for i,j in Single_instance_drugclass_genomedict.items():
    for a,b in j.items():
        b=normalized_bitscore(b)

# A function to remove BLAST entry results for same gene comparison by taking the highest bit score 
def remove_duplicates(temp):
    df=pd.DataFrame(temp,columns=["query_id","sub_id","bitscore"])
    d = df.sort_values("bitscore", ascending=False)
    d = d.drop_duplicates(['query_id'])
    d.reset_index(drop=True, inplace=True)
    
    return d


single_similarity_array_dict={}
for k,l in Drug_singleinstance_range_10_locus_dict10.items():
    similarity_array=[]
    for genome1 in Single_instance_drugclass_genomedict[k].keys():  
        u=Single_instance_drugclass_genomedict[k][genome1]
        for genome2,b in l.items():
            temp=[]
            sum1=0
            sum2=0
            u.reset_index(drop=True, inplace=True)
            for index in range (len(u)):
                if u["sub_id"][index] in b:
                    temp.append((u["query_id"][index],u["sub_id"][index],u["normalized_bitscore"][index]))

            tempframe=remove_duplicates(temp)
            sum1=round(tempframe["bitscore"].sum(),3)
            
            if Single_contig_end_flag_dict[k][genome1] ==1:
                length_original_neighborhood=len(Drug_singleinstance_range_10_locus_dict10[k][genome1])
                difference=length_original_neighborhood-len(tempframe)
                sum2= sum1+(21-difference-len(tempframe))
                
            elif Single_contig_end_flag_dict[k][genome2]==1:
                length_original_neighborhood=len(Drug_singleinstance_range_10_locus_dict10[k][genome2])
                difference=length_original_neighborhood-len(tempframe)                                 
                
                sum2= sum1+(21-difference-len(tempframe))
            else:
                sum2=sum1
           
            similarity_array.append((genome1,genome2,round(sum2,3)))
    single_similarity_array_dict[k]=similarity_array
    
# A function to convert the similarity matrix to symmteric matric matrix by taking the average of scores
def converttosymmetrix(j):
    newj=[]
    
    for index1 in range(len(j)):
        
        score1=j[index1][2]
        for index2 in range(len(j)):
            if j[index2][1]==j[index1][0] and j[index2][0]==j[index1][1]:
                avg=round((j[index1][2]+j[index2][2])/2,3)
        newj.append((j[index1][0],j[index1][1],avg))
    return newj


# Use the function to convert to symmetric matrix and store the final matrix in separate dictionaries
Single_instance_matrixframe={}
for i,j in single_similarity_array_dict.items():
    new=converttosymmetrix(j)
    finalresult=pd.DataFrame(new,columns=["Query_id","Sub_id","normalized_bitscore"])
    sim_mat=pd.crosstab(index=finalresult.iloc[:,0], columns=finalresult.iloc[:,1],values=finalresult.iloc[:,2], aggfunc=lambda x: x,colnames=None)
    Y=sim_mat.values
    Single_instance_matrixframe[i]=sim_mat

# Generate a distance matrix and generate UPGMA clusters 
os.mkdir("Outputs/UPGMA_single_clusters")
save_path="Outputs/UPGMA_single_clusters"
dataarray=[]
hierarchy.set_link_color_palette(['r', 'g', 'y', 'm'])
matplotlib.rcParams['lines.linewidth'] = 5
for i,j in Single_instance_matrixframe.items(): 
    if len(j)>1:

        df = pd.DataFrame(j,columns=j.keys())
        r  = df.values
        df = df.transform(lambda x: 1 - x/r.max() )
        df = df.round(decimals=3)
        Y = df.values
        np.fill_diagonal(Y, 0)
        Y = distance.squareform(Y)
        fig = plt.figure(figsize=(15,5))
        Z = hierarchy.linkage(Y, 'average')
        dendrogram(Z, leaf_rotation=90, leaf_font_size=8, labels=j.keys())
        savename=os.path.join(save_path,i+"_.png")
        plt.ylabel("Distance between neighborhoods")
        plt.xlabel("Neighborhoods with their respective Genome IDs")
        plt.savefig(savename,bbox_inches='tight', dpi=100)
        plt.close()


# Execute if only there are multiple instance AMR gene models #

if len(temp_Multiple_instance_blastdatadict.keys())>0:

    Multiple_instance_blastdatadict={}
    for i,j in temp_Multiple_instance_blastdatadict.items():
        j=j[j["PI"]>70]
        j = j.reset_index().drop_duplicates(subset=['query_id','sub_id','bitscore'],keep='first').set_index('index')
        Multiple_instance_blastdatadict[i]=j


    # Divide the mixed BLAST results into their respective AMR gene models and genomes 
    multiple_instance_drugclass_genomedict={}
    for k,l in Drug_multiple_instance_range_10_locus_dict10.items():
        u=Multiple_instance_blastdatadict[k]
        temp_dict={} 
        for a,b in l.items():#          
            temp=[]
            u.reset_index(drop=True, inplace=True)
            for index in range (len(u)):
                if  u["query_id"][index] in b:
                    temp.append((u["query_id"][index],u["sub_id"][index],u["PI"][index],u["bitscore"][index]))
                        
            #frame=pd.DataFrame(temp,columns=["query_id","sub_id","PI","bitscore"])
            temp_dict[a]=pd.DataFrame(temp,columns=["query_id","sub_id","PI","bitscore"])
          
        multiple_instance_drugclass_genomedict[k]=temp_dict

    # using the function and normalizing 
    for i,j in multiple_instance_drugclass_genomedict.items():
        for a,b in j.items():
            b=normalized_bitscore(b)


    # For each gene model compare the genomes and generate a similarity matrix #
    multiple_similarity_array_dict={}
    for k,l in Drug_multiple_instance_range_10_locus_dict10.items():
        similarity_array=[]
        for genome1 in multiple_instance_drugclass_genomedict[k].keys():  
            u=multiple_instance_drugclass_genomedict[k][genome1]
            for genome2,b in l.items():
                temp=[]
                sum1=0
                sum2=0
                u.reset_index(drop=True, inplace=True)
                for index in range (len(u)):
                    if u["sub_id"][index] in b:
                        temp.append((u["query_id"][index],u["sub_id"][index],u["normalized_bitscore"][index]))

                tempframe=remove_duplicates(temp)
                sum1=round(tempframe["bitscore"].sum(),3)
                
                if multiple_contig_end_flag_dict[k][genome1] ==1:
                    length_original_neighborhood=len(Drug_multiple_instance_range_10_locus_dict10[k][genome1])
                    difference=length_original_neighborhood-len(tempframe)
                    sum2= sum1+(21-difference-len(tempframe))
                    
                elif multiple_contig_end_flag_dict[k][genome2]==1:
                    length_original_neighborhood=len(Drug_multiple_instance_range_10_locus_dict10[k][genome2])
                    difference=length_original_neighborhood-len(tempframe)                                 
                    
                    sum2= sum1+(21-difference-len(tempframe))
                else:
                    sum2=sum1
               
                similarity_array.append((genome1,genome2,round(sum1,3)))
        multiple_similarity_array_dict[k]=similarity_array

    multiple_instance_matrixframe={}
    for i,j in multiple_similarity_array_dict.items():
        new=converttosymmetrix(j)
        finalresult=pd.DataFrame(new,columns=["Query_id","Sub_id","normalized_bitscore"])
        sim_mat=pd.crosstab(index=finalresult.iloc[:,0], columns=finalresult.iloc[:,1],values=finalresult.iloc[:,2], aggfunc=lambda x: x,colnames=None)
        Y=sim_mat.values
        multiple_instance_matrixframe[i]=sim_mat

    os.mkdir("Outputs/UPGMA_multiple_clusters")
    save_path="Outputs/UPGMA_multiple_clusters"
    dataarray=[]
    hierarchy.set_link_color_palette(['r', 'g', 'y', 'm'])
    matplotlib.rcParams['lines.linewidth'] = 5
    for i,j in multiple_instance_matrixframe.items(): 
        if len(j)>1:

            df = pd.DataFrame(j,columns=j.keys())
            r  = df.values
            df = df.transform(lambda x: 1 - x/r.max() )
            df = df.round(decimals=3)
            Y = df.values
            np.fill_diagonal(Y, 0)
            Y = distance.squareform(Y)
            fig = plt.figure(figsize=(15,5))
            Z = hierarchy.linkage(Y, 'average')
            dendrogram(Z, leaf_rotation=90, leaf_font_size=8, labels=j.keys())
            savename=os.path.join(save_path,i+"_.png")
            plt.ylabel("Distance between neighborhoods")
            plt.xlabel("Neighborhoods with their respective Genome IDs")
            plt.savefig(savename,bbox_inches='tight', dpi=100)
            plt.close()
