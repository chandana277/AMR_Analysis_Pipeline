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

# A function to visualize the gene order of each neighborhood and save into their respective folders in the form of .jpeg format#
def contigend_visualization(contig_array,end_direction,genome,reverse_term,identity,drugclass,instancetype):
    #if drugclass=="fluoroquinolone_lincosamide":
    RGIgene=genome.split("(")[1].split(")")[0]
    #print(RGIgene)
    font_dict={'size':15,'weight':'bold','family':'Helvetica'}

    save_temp_name="Outputs/All_neighborhoods/"+instancetype+"/"+drugclass

    Features=[]
    temp_array_totrack_length=[]   

    contig_array.reset_index(drop=True, inplace=True)


    if end_direction=="upward":
        b=GraphicFeature(start=contig_array["GeneStart"][0]-500, end=contig_array["GeneStart"][0],color="#0A090A",label="Ends_upward",fontdict=font_dict)
        Features.append(b)  
        temp_array_totrack_length.append((contig_array["GeneStart"][0]-500,contig_array["GeneStart"][0]))


    for i in range(len(contig_array)):
            #print(contig_array["GeneStart"][i])
            if str(contig_array["GeneName"][i])== RGIgene:
                a=GraphicFeature(start=contig_array["GeneStart"][i], end=contig_array["GeneEnd"][i], strand=int(contig_array["Strand"][i]),color=contig_array["Genecolor"][i],label=str(contig_array["GeneName"][i]),fontdict=font_dict,thickness=26,linecolor="#F72808",linewidth=2.5)
            else:
                a=GraphicFeature(start=contig_array["GeneStart"][i], end=contig_array["GeneEnd"][i], strand=int(contig_array["Strand"][i]),color=contig_array["Genecolor"][i],label=str(contig_array["GeneName"][i]),fontdict=font_dict,thickness=26,linewidth=1.7)

            Features.append(a)
            temp_array_totrack_length.append((contig_array["GeneStart"][i],contig_array["GeneEnd"][i]))


    if end_direction=="downward":
        b=GraphicFeature(start=contig_array["GeneEnd"].iloc[-1]+250, end=contig_array["GeneEnd"].iloc[-1]+1500,color="#0A090A",label="Ends_downward",fontdict=font_dict)
        Features.append(b) 
        temp_array_totrack_length.append((contig_array["GeneEnd"].iloc[-1]+250,contig_array["GeneEnd"].iloc[-1]+1500))

    if end_direction=="both":
        #print(contig_array["GeneStart"][0])
        z=GraphicFeature(start=contig_array["GeneStart"][0]-500, end=contig_array["GeneStart"][0],color="#0A090A",label="Ends_upward",fontdict=font_dict)
        Features.append(z)  
        temp_array_totrack_length.insert(0,((contig_array["GeneStart"][0]-500,contig_array["GeneStart"][0])))

        #print(contig_array["GeneEnd"].iloc[-1])
        c=GraphicFeature(start=contig_array["GeneEnd"].iloc[-1]+250, end=contig_array["GeneEnd"].iloc[-1]+1500,color="#0A090A",label="Ends_downward",fontdict=font_dict)
        Features.append(c) 
        temp_array_totrack_length.append((contig_array["GeneEnd"].iloc[-1]+250,contig_array["GeneEnd"].iloc[-1]+1500))


    length=temp_array_totrack_length[-1][1] - temp_array_totrack_length[0][0]

    record = GraphicRecord(first_index=temp_array_totrack_length[0][0],sequence_length=length,features=Features)  

    ax,_=record.plot(figure_width=30,strand_in_label_threshold=7)

    if (reverse_term=="-1"):
        ax.invert_xaxis()

    temp_name=genome+".jpeg" 
    #import matplotlib.pyplot as mpl
    #plt.rcParams['font.size'] = 22
    #csfont = {'fontname':'Comic Sans MS','fontsize':20}
    #ax.set_title(genome,fontname="serif", fontsize=18,color="blue") 

    #ax.figure.set_size_inches(20, 4)
    plt.close('all')
    name = os.path.join(save_temp_name,temp_name)

    ax.figure.savefig(name)
    

# A function to delete keys of those AMR genes who are not present in atleast 25% of the total genomes ##
def delete_keys_lessthan_25percent_instances(dict_element):
    emptykeyslist= []
    total_genomes= len(gbkfiles)
    minimum_genomes = int((total_genomes * 25)/100)
    
    for i,j in dict_element.items():
        if len(j)<14:
            emptykeyslist.append(i)
    for i in emptykeyslist:
        del dict_element[i]
    return dict_element


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
                #print(uname)
                print("contig does not exist----"+"drugclass:"+drug)
                print(rec['Best_Hit_ARO'][j])
      
      
      if len(m)!=0 and len(n)==0:
         # print(drug,len(m),len(m)-1)

          m.reset_index(drop=True, inplace=True)
          #print(len(m),len(m)-1)
          m=m.drop([len(m)-1])
          e = pd.concat([m,o,n], ignore_index=True)
          e.reset_index(drop=True, inplace=True)
          neighbor_genes.append(e)

          e=checkforRGIinneighborhood(e,key,rec["Best_Hit_ARO"][j])

          if len(m)==g:
             e.reset_index(drop=True, inplace=True)
             contig_flag=1
             #neighbor_genes.append(e)
             contigend_visualization(e,"downward",genome,rec["Orientation"][j],rec["Best_Identities"][j],drug,instancetype)
          else:
             contigend_visualization(e,"both",genome,rec["Orientation"][j],rec["Best_Identities"][j],drug,instancetype)


    
      if len(m)==0 and len(n)!=0:

          m.reset_index(drop=True, inplace=True)
          #print(len(m),len(m)-1)
          n=n.drop([len(n)-1])
          e = pd.concat([m,o,n], ignore_index=True)
          e.reset_index(drop=True, inplace=True)
          neighbor_genes.append(e)
          e=checkforRGIinneighborhood(e,key,rec["Best_Hit_ARO"][j])

      if(len(m)!=0 and len(n)!=0):
        m.reset_index(drop=True, inplace=True)
        #print(len(m),len(m)-1)
        m=m.drop([len(m)-1])
        e = pd.concat([m,o,n], ignore_index=True)
        e.reset_index(drop=True, inplace=True)

        e=checkforRGIinneighborhood(e,key,rec["Best_Hit_ARO"][j])
        if(len(n)<g and len(m)==g):
            contig_flag=1
            contigend_visualization(e,"downward",genome,rec["Orientation"][j],rec["Best_Identities"][j],drug,instancetype)
            e.reset_index(drop=True, inplace=True)
            neighbor_genes.append(e)
        elif (len(m)<g and len(n)==g): 
            contig_flag=1
            contigend_visualization(e,"upward",genome,rec["Orientation"][j],rec["Best_Identities"][j],drug,instancetype)
            e.reset_index(drop=True, inplace=True)
            neighbor_genes.append(e)
        elif (len(m)<g and len(n)<g): 
            contig_flag=1
            contigend_visualization(e,"both",genome,rec["Orientation"][j],rec["Best_Identities"][j],drug,instancetype)
            e.reset_index(drop=True, inplace=True)
            neighbor_genes.append(e)

        elif(len(m)==g and len(n)==g):
            e.reset_index(drop=True, inplace=True)
            neighbor_genes.append(e)
            contigend_visualization(e,"noend",genome,rec["Orientation"][j],rec["Best_Identities"][j],drug,instancetype)


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


#   Generates the summary of all the gene models and genomes #

summaryfile=open("Summary of all models.txt","w")
summaryfile.truncate(0)

summaryfile.write("----------------------------------------------------------------------------------------"+"\n")
summaryfile.write("SUMMARY"+"\n")
summaryfile.write("----------------------------------------------------------------------------------------"+"\n")

summaryfile.write("Total Number of Genomes Analysed: "+str(len(filenames))+"\n")

Total_genemodels=len(Dict_multigene_instances.keys())+len(Dict_singlegene_instances.keys())
summaryfile.writelines("Total number of Gene Models: "+str(Total_genemodels)+"\n")
multiple_inst_number=str(len(Dict_multigene_instances.keys()))
single_inst_number=str(len(Dict_singlegene_instances.keys()))
summaryfile.writelines("Number of Multiple Instance Gene Models( More than one gene per gene model per genome): "+multiple_inst_number+"\n")
summaryfile.writelines("Number of Single Instance Gene Models(One gene per gene model per genome): "+single_inst_number+"\n\n")

stat1=[]
for i,j in datadict.items():
  
    strict_count1=len(j[j["Cut_Off"]=="Strict"])
    perfect_count1=len(j[j["Cut_Off"]=="Perfect"])
    
    stat1.append((i,strict_count1,perfect_count1))
temp=pd.DataFrame(stat1,columns=["Genome","Strict_Hits","Perfect_Hits"])
#print(temp)
summaryfile.write(temp.to_string(header = True, index = False))

summaryfile.write("\n\n\n")

summaryfile.write("----------------------------------------------------------------------------------------"+"\n")
summaryfile.write("Multiple Instance Gene Models:"+"\n")
summaryfile.write("----------------------------------------------------------------------------------------"+"\n")


#print((dict_rgi_multiple_occurance.keys()))

for a,b in Dict_multigene_instances.items():
    summaryfile.write("Gene Model: "+a+ "   Total Number of Genomes: "+str(len(b))+"\n")
    stat2=[]
    for i,j in b.items():
        #LooseID,loose_count2=j[j["Cut_Off"]=="Loose"]["Best_Identities"].values,len(j[j["Cut_Off"]=="Loose"])
        StrictID,strict_count2=j[j["Cut_Off"]=="Strict"]["Best_Identities"].values,len(j[j["Cut_Off"]=="Strict"])
        perfect_count2=len(j[j["Cut_Off"]=="Perfect"])
        
        stat2.append((i,strict_count2,perfect_count2,StrictID))
    temp=pd.DataFrame(stat2,columns=["Genome","Strict_Hits","Perfect_Hits","StrictID"])
    summaryfile.write(temp.to_string(header = True, index = False))

    summaryfile.write("\n")

summaryfile.write("----------------------------------------------------------------------------------------"+"\n")
summaryfile.write("Single Instance Gene Models:"+"\n")
summaryfile.write("----------------------------------------------------------------------------------------"+"\n")

#print((dict_rgi_multiple_occurance.keys()))

for a,b in Dict_singlegene_instances.items():
    summaryfile.write("Gene Model: "+a+ "   Total Number of Genomes: "+str(len(b))+"\n")
    stat2=[]
    for i,j in b.items():
        #LooseID,loose_count2=j[j["Cut_Off"]=="Loose"]["Best_Identities"].values,len(j[j["Cut_Off"]=="Loose"])
        StrictID,strict_count2=j[j["Cut_Off"]=="Strict"]["Best_Identities"].values,len(j[j["Cut_Off"]=="Strict"])
        Perfect_count2=len(j[j["Cut_Off"]=="Perfect"])
        
        stat2.append((i,strict_count2,Perfect_count2,StrictID))
    temp=pd.DataFrame(stat2,columns=["Genome","Strict_Hits","Perfect_Hits","StrictID"])
    summaryfile.write(temp.to_string(header = True, index = False))
    
    summaryfile.write("\n")
summaryfile.close()


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
# Visualize each neighborhood and save 

os.mkdir("Outputs")
os.mkdir("Outputs/All_neighborhoods")
os.mkdir("Outputs/All_neighborhoods/Multiple_Instance_Neighborhood")
os.mkdir("Outputs/All_neighborhoods/Single_Instance_Neighborhood")

for eachkey in Dict_multigene_instances.keys():
    n=os.path.join("Outputs/All_neighborhoods/Multiple_Instance_Neighborhood/",eachkey)
    os.mkdir(n)
    
for eachkey in Dict_singlegene_instances.keys():
    n=os.path.join("Outputs/All_neighborhoods/Single_Instance_Neighborhood/",eachkey)
    os.mkdir(n)


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
        
# Create separate folders for single and multiple and generate .fasta files containing all genes for each AMR gene model to BLAST #
os.mkdir("Outputs/Multiple_instance")
save_path = 'Outputs/Multiple_instance/'   

for i,j in multipleinstance_neighboringdict_combined_range_10.items():
    
    name = os.path.join(save_path,i+".fasta") 
    filef=open(name,"w")
    for k in j.keys():
        for a,b in zip(Drug_multiple_instance_range_10_locus_dict10[i][k],Drug_multiple_instance_range_10_protein_dict10[i][k]):
            #print(b)
            filef.writelines(">" + '{}'.format(a))
            filef.write("\n")
            filef.writelines('{}'.format(b))
            filef.write("\n")
    filef.close()


os.mkdir("Outputs/Single_instance")
save_path = 'Outputs/Single_instance/'   

for i,j in singleinstance_neighboringdict_combined_range10.items():
    name=i+".fasta"
    name = os.path.join(save_path,i+".fasta") 
    filef=open(name,"w")
    for k in j.keys():
        for a,b in zip(Drug_singleinstance_range_10_locus_dict10[i][k],Drug_singleinstance_range_10_protein_dict10[i][k]):
            #print(b)
            filef.writelines(">" + '{}'.format(a))
            filef.write("\n")
            filef.writelines('{}'.format(b))
            filef.write("\n")
    filef.close()


######################## Part 1 ends ##################################
