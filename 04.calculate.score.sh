#! /bin/bash

##########part4:calculate score for calculate PRS


#####author:Yang Sui(ysui@broadinsititue.org)
#####Date:2025.03.15

# Print the total number of arguments
echo "Total number of arguments: $#"

# Accessing each argument by its position
echo "Plink Argument: $1"
echo "pfile Argument: $2"
echo "extract Argument: $3"
echo "remove sample id Argument: $4"
echo "score Argument: $5" ###the first columns is ID, the second is effect allele and it has header
echo "score-col-nums Argument: $6" 
echo "out Argument: $7"


#######update the SNP ID
$1 \
--pfile $2 \
--extract $3 \
--remove $4 \
--score $5 1 4 header-read list-variants ignore-dup-ids cols='sid,nallele,dosagesum,scoresums' \
--score-col-nums $6 \
--out $7






