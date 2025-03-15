#! /bin/bash

##########part3:update ID for calculate PRS


#####author:Yang Sui(ysui@broadinsititue.org)
#####Date:2025.03.15

# Print the total number of arguments
echo "Total number of arguments: $#"

# Accessing each argument by its position
echo "Plink Argument: $1"
echo "pfile Argument: $2"
echo "update-name Argument: $3"
echo "out Argument: $4"


#######update the SNP ID
$1 \
--pfile $2 \
--update-name $3 \
--make-pgen \
--out $4






