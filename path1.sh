read -p "Please enter the RGI file/files path : " RGI_Path
read -p "Please enter the GBK file/files path : " GBK_Path


python Neighbor_Generator.py ${RGI_Path} ${GBK_Path}


cd Output/Multiple_instance
mkdir output_blast_multiple_instance
for i in *.fasta
do    

     name=$(echo $i | cut -d'.' -f1)
     sed -i s/\'//g $i; #########removes single quotes from the file
     sed -i 's/[][]//g' $i; ############removes square braces from the file

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
     sed -i s/\'//g $i; #########removes single quotes if any from the file
     sed -i 's/[][]//g' $i; ############removes square braces if any from the file

     makeblastdb -dbtype prot -in $i
     blastp -query $i -db $i -outfmt  "6 qseqid sseqid pident length evalue bitscore qseq sseq " -out output_blast_single_instance/$name.txt
     
    
done

cd "$OLDPWD"
python Clusters_Generator.py ${RGI_Path} ${GBK_Path} Single_instance/output_blast_single_instance Multiple_instance/output_blast_multiple_instance


