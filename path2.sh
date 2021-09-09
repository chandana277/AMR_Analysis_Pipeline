read -p "Please enter the .fna assembly files path : " Fna_Path_RGI_Annotation
      
      
cd ${Fna_Path_RGI_Annotation}  

mkdir allgbksrequired 
mkdir allfnasrequired
mkdir allrgisrequired


CONDA_BASE=$($CONDA_EXE info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
# Change this to the path of prokka environment
conda activate prokka


for i in *.fasta
do
     name=$(echo $i | cut -d'.' -f1)        
     prokka $i --prefix $name --outdir ${name}_prokka  --locustag $name  


     cp ${name}_prokka/${name}.gbk allgbksrequired
     cp ${name}_prokka/${name}.fna allfnasrequired
 
done


CONDA_BASE=$($CONDA_EXE info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
# Change this to the path of rgi environment
conda activate rgi   


for i in allfnasrequired/*.fna
do

     name=$(basename $i | cut -d '.' -f1)         
     rgi main --clean --input_sequence $i --alignment_tool BLAST  --num_threads 1 --output $name

 
done 

for i in *.txt
do
 
    mv $i allrgisrequired
 
done


cd "$OLDPWD" 
python Neighbor_Generator.py ${Fna_Path_RGI_Annotation}/allrgisrequired ${Fna_Path_RGI_Annotation}/allgbksrequired


cd Output/Multiple_instance
mkdir output_blast_multiple_instance
for i in *.fasta
do    

       name=$(echo $i | cut -d'.' -f1)
       sed -i s/\'//g $i; #########removes single quotes if any from the file
       sed -i 's/[][]//g' $i; ############removes square braces if any from the file

       makeblastdb -dbtype prot -in $i
       blastp -query $i -db $i -outfmt  "6 qseqid sseqid pident length evalue bitscore qseq sseq " -out output_blast_multiple_instance/$name.txt
       #mv $i.txt output_blast/$name.txt
      
done

cd "$OLDPWD"

cd Output/Single_instance
mkdir output_blast_single_instance

for i in *.fasta
do    

       name=$(echo $i | cut -d'.' -f1)
       echo $name
       sed -i s/\'//g $i; #########removes single quotes from the file
       sed -i 's/[][]//g' $i; ############removes square braces from the file

       makeblastdb -dbtype prot -in $i
       blastp -query $i -db $i -outfmt  "6 qseqid sseqid pident length evalue bitscore qseq sseq " -out output_blast_single_instance/$name.txt
       
      
done

cd "$OLDPWD"

python Clusters_Generator.py ${Fna_Path_RGI_Annotation}/allrgisrequired ${Fna_Path_RGI_Annotation}/allgbksrequired Single_instance/output_blast_single_instance Multiple_instance/output_blast_multiple_instance

