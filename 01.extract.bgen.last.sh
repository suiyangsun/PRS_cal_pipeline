#! /bin/bash

##########part1:extract the region for calculate PRS
####this is design for MGB, because the bgen file is 'ref-last'

#####author:Yang Sui(ysui@broadinsitute.org)
#####Date:2025.03.15

# Print the total number of arguments
echo "Total number of arguments: $#"

# Accessing each argument by its position
echo "Plink Argument: $1"
echo "bgen Argument: $2"
echo "sample Argument: $3"
echo "extract Argument: $4"
echo "out Argument: $5"


$1 \
--bgen $2 'ref-last' \
--sample $3 \
--extract bed0 $4 \
--make-pgen \
--out $5



