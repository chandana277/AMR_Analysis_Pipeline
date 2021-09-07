#!/bin/bash

echo " WELCOME TO CLUSTERING AND ANALAYSIS OF NEIGHBORHOOD GENES"
echo "Executing requirements.txt"

pip install -r requirements.txt

echo "###################################################################################################"

echo "PLEASE CHOOSE FROM THE FOLLOWING OPTIONS FOR FURTHER EXECUTION:"
echo "1)LOCAL .GBK AND .RGI ANNOTATED FILE EXISTS"
echo "2)ONLY LOCAL .FNA FILE EXISTS "

echo "###################################################################################################"

read -p "Please enter the option : " option


case ${option} in 

     1)  sh path1.sh  ;;
	 
	2)  sh path2.sh  
	 	 exit 1 	  ;;
      
  
esac 	
		

