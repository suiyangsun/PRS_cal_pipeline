#! /bin/bash


source /broad/software/scripts/useuse
reuse -q gzip
reuse -q R-4.1
reuse -q Python-3.9
#$ -cwd

#$ -V

#$ -N GSA_53K_calculate_score

#$ -o /medpop/esp2/yang/project/EC_PRS/script/GSA_53K/log/

#$ -e /medpop/esp2/yang/project/EC_PRS/script/GSA_53K/log/

#$ -pe smp 1 -R y -binding linear:1 -l h_vmem=10g

#$ -l h_rt=100:00:00

########################01.extract the bed########################
#########Extract SNPs from whole variants to save time, especially when need calcuate multiple scores
mkdir -p /medpop/esp2/yang/project/EC_PRS/GSA_53K/all/

cat /medpop/esp2/yang/project/EC_PRS/weight/processed.hg38.all.weight.txt \
| tail -n+2|cut -f 1 |sort |uniq |awk -F ':' '{print $1,$2-1,$2}' \
> /medpop/esp2/yang/project/EC_PRS/GSA_53K/all/hg38.all.bed


im="/medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen"

for i in {1..22}
do
awk -v OFS="\t" -v var="$i" '$1==var' /medpop/esp2/yang/project/EC_PRS/GSA_53K/all/hg38.all.bed \
> /medpop/esp2/yang/project/EC_PRS/GSA_53K/all/chr$i.hg38.all.bed


bgen=$im/GSA_53K.merged.chr${i}.bgen
sample=$im/GSA_53K.merged.chr${i}.sample
extract=/medpop/esp2/yang/project/EC_PRS/GSA_53K/all/chr$i.hg38.all.bed
out=/medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}


######extract the bed file
sh /medpop/esp2/yang/project/PRS_pipeline/01.extract.bgen.last.sh \
$plink2 $bgen sample $extract $out


done


########################02.set the id########################

plink2='/medpop/esp2/yang/software/plink2.231029/plink2'

for i in {1..22}
do

pfile=/medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}
out=/medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.setid


######set the SNP ID
sh /medpop/esp2/yang/project/PRS_pipeline/02.setID.sh \
$plink2 $pfile $out

done



#########prepare the ID
for i in {1..22}
do

cat /medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.setid.pvar \
|tail -n+2 \
|awk -v FS='\t' -v OFS="\t" '{if($4<$5){print $3,$1":"$2":"$4":"$5}else{print $3,$1":"$2":"$5":"$4}}' \
> /medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.SNPID.txt

done


########################03.Update the id########################
for i in {1..22}
do

pfile=/medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.setid
update=/medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.SNPID.txt
out=/medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.upid


sh /medpop/esp2/yang/project/PRS_pipeline/03.updateID.sh \
$plink2 $pfile $update $out


done



########################04.Calculate the score########################
weight="/medpop/esp2/yang/project/EC_PRS/weight"
SNPID="/medpop/esp2/yang/project/EC_PRS/GSA_53K/SNPID"
mkdir -p $SNPID

mkdir -p /medpop/esp2/yang/project/EC_PRS/GSA_53K/sscore_SQ



######There's four scores needed to calculate
awk  -v FS='\t' -v OFS="\t" 'NR == 1 || $7 == "X"' \
/medpop/esp2/yang/project/EC_PRS/weight/processed.hg38.all.weight.txt \
> /medpop/esp2/yang/project/EC_PRS/weight/EC.processed.hg38.all.weight.txt

awk  -v FS='\t' -v OFS="\t" 'NR == 1 || $8 == "X"' \
/medpop/esp2/yang/project/EC_PRS/weight/processed.hg38.all.weight.txt \
> /medpop/esp2/yang/project/EC_PRS/weight/Lipid.processed.hg38.all.weight.txt


awk  -v FS='\t' -v OFS="\t" 'NR == 1 || $9 == "X"' \
/medpop/esp2/yang/project/EC_PRS/weight/processed.hg38.all.weight.txt \
> /medpop/esp2/yang/project/EC_PRS/weight/Non-EC.processed.hg38.all.weight.txt


awk  -v FS='\t' -v OFS="\t" 'NR == 1 || $10 == "X"' \
/medpop/esp2/yang/project/EC_PRS/weight/processed.hg38.all.weight.txt \
> /medpop/esp2/yang/project/EC_PRS/weight/Non-EC_Non-Lipid.processed.hg38.all.weight.txt




for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do

cut -f 1 /medpop/esp2/yang/project/EC_PRS/weight/${pheno}.processed.hg38.all.weight.txt \
|tail -n+2 > /medpop/esp2/yang/project/EC_PRS/GSA_53K/SNPID/${pheno}.SNPID.txt


for i in {1..22}
do

pfile=/medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.upid
extract=/medpop/esp2/yang/project/EC_PRS/GSA_53K/SNPID/${pheno}.SNPID.txt
update=/medpop/esp2/yang/project/EC_PRS/GSA_53K/all/GSA_53K.chr${i}.SNPID.txt
removesample=/medpop/esp2/yang/project/MGB_pheno/GSA_53K_sampleQCfail.txt
score=/medpop/esp2/yang/project/EC_PRS/weight/${pheno}.processed.hg38.all.weight.txt
score_col_nums=6
out=/medpop/esp2/yang/project/EC_PRS/GSA_53K/sscore_SQ/${pheno}_chr${i}


sh /medpop/esp2/yang/project/PRS_pipeline/04.calculate.score.sh \
$plink2 $pfile $extract $removesample $score $score_col_nums $out


done
done


########################05.combine all chr together########################
# ####################all

mkdir -p /medpop/esp2/yang/project/EC_PRS/sscore_clean

weight="/medpop/esp2/yang/project/EC_PRS/weight"


for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do

cat /medpop/esp2/yang/project/EC_PRS/GSA_53K/sscore_SQ/${pheno}_chr*.sscore \
|grep -v "#IID" > /medpop/esp2/yang/project/EC_PRS/GSA_53K/sscore_SQ/chrall.${pheno}.all.sscore

Rscript /medpop/esp2/yang/project/cli_trail/script/combinechr.R \
/medpop/esp2/yang/project/EC_PRS/GSA_53K/sscore_SQ/chrall.${pheno}.all.sscore \
/medpop/esp2/yang/project/EC_PRS/sscore_clean/GSA_53K_sum.chrall.${pheno}.all.sscore

gzip -f /medpop/esp2/yang/project/EC_PRS/sscore_clean/GSA_53K_sum.chrall.${pheno}.all.sscore

########And variants used in this score
cat /medpop/esp2/yang/project/EC_PRS/GSA_53K/sscore_SQ/${pheno}_chr*.sscore.vars \
> /medpop/esp2/yang/project/EC_PRS/sscore_clean/GSA_53K_sum.chrall.${pheno}_all.sscore.vars


rm /medpop/esp2/yang/project/EC_PRS/GSA_53K/sscore_SQ/chrall.${pheno}.all.sscore


done


########################06.count number########################
echo -e "id\torignial\tprocessed\tscore" > /medpop/esp2/yang/project/EC_PRS/sscore_clean/GSA_53K_variants_wc.txt


for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do


orignial=$(cat /medpop/esp2/yang/project/EC_PRS/weight/${pheno}.processed.hg38.all.weight.txt |tail -n+2 |wc -l)


processed=$(cat /medpop/esp2/yang/project/EC_PRS/weight/${pheno}.processed.hg38.all.weight.txt |tail -n+2 |wc -l)


score=$(cat /medpop/esp2/yang/project/EC_PRS/sscore_clean/GSA_53K_sum.chrall.${pheno}_all.sscore.vars |wc -l)

echo -e "$pheno\t$orignial\t$processed\t$score" >> /medpop/esp2/yang/project/EC_PRS/sscore_clean/GSA_53K_variants_wc.txt

done




