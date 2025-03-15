#! /bin/bash

##########part2:set ID for calculate PRS


#####author:Yang Sui(ysui@broadinsitute.org)
#####Date:2025.03.15

# Print the total number of arguments
echo "Total number of arguments: $#"

# Accessing each argument by its position
echo "Plink Argument: $1"
echo "pfile Argument: $2"
echo "out Argument: $3"


$1 \
--pfile $2 \
--set-all-var-ids 'chr@:#:$r:$a' \
--new-id-max-allele-len 10000 truncate \
--make-pgen \
--out $3






