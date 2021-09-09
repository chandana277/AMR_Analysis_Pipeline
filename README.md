# AMR_Analysis_Pipeline
A pipeline written in bash script that combines python and other tools.



Steps to run the pipeline on local genomes:

1>Include all the files in the directory and run " **source setup.sh**" (Please use source for the purpose of virtual env)

2>Follow the insttuctions and choose either of the option (1 or 2)

3>In case of option 1 : Please have the separate folders for rgi and gbk files inside the same directory as code files and provide the path.

  As the result, each separate folders for output are created.
  
  Viualizations for all the gene models can be obtained in the folder "All_neighborhoods"
  
  UPGMA clusters are found in separate folders 
  
4>In case of option 2 : Please have a separate folder for all the .fa files and include the path


*********** Important Information *************

1> Please use source to execute setup.sh

2> To use option 1 of the pipeline,please make sure you have the below softwares installed.

    BLAST 2.9.0
    
3> To use option 2 of the pipeline, please make sure you have the below softwares installed.

    Prokka 1.14.6
    RGI 4.2.2
    BLAST 2.9.0


Base Version Functionalities:

1> Assumes the user is analysing more than 1 genome.

2> Analyses only strict and perfect AMR genes.

3> Computes 10 neighbors upstream and downstream of AMR gene
