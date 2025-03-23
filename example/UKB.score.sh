#! /bin/bash


source /broad/software/scripts/useuse
reuse -q gzip
reuse -q R-4.1
reuse -q Python-3.9
#$ -cwd

#$ -V

#$ -N UKB_calculate_score

#$ -o /medpop/esp2/yang/project/EC_PRS/script/UKB/log/

#$ -e /medpop/esp2/yang/project/EC_PRS/script/UKB/log/

#$ -pe smp 5 -R y -binding linear:5 -l h_vmem=20g

#$ -l h_rt=100:00:00

########################01.extract the bed########################

mkdir -p /medpop/esp2/yang/project/EC_PRS/UKB/all/

cat /medpop/esp2/yang/project/EC_PRS/weight/processed.hg37.all.weight.txt \
| tail -n+2|cut -f 1 |sort |uniq |awk -F ':' '{print $1,$2-1,$2}' \
> /medpop/esp2/yang/project/EC_PRS/UKB/all/hg37.all.bed



im="/broad/ukbb/imputed_v3"

plink2='/medpop/esp2/wallace/tools/plink2/plink2'



for i in {1..22}
do
awk -v OFS="\t" -v var="$i" '$1==var' /medpop/esp2/yang/project/EC_PRS/UKB/all/hg37.all.bed \
> /medpop/esp2/yang/project/EC_PRS/UKB/all/chr$i.hg37.all.bed


bgen=$im/ukb_imp_chr${i}_v3.bgen
sample=/medpop/esp2/projects/UK_Biobank/linkers/app7089/ukb22828_c1_b0_v3_s487207.sample
extract=/medpop/esp2/yang/project/EC_PRS/UKB/all/chr$i.hg37.all.bed
out=/medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}


######extract the bed file
sh /medpop/esp2/yang/project/PRS_pipeline/01.extract.bgen.first.sh \
$plink2 $bgen $sample $extract $out


done


########################02.set the id########################

plink2='/medpop/esp2/yang/software/plink2.231029/plink2'

for i in {1..22}
do

pfile=/medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}
out=/medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}.setid


######set the SNP ID
sh /medpop/esp2/yang/project/PRS_pipeline/02.setID.sh \
$plink2 $pfile $out

done



#########prepare the ID
for i in {1..22}
do

cat /medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}.setid.pvar \
|tail -n+2 \
|awk -v FS='\t' -v OFS="\t" '{if($4<$5){print $3,$1":"$2":"$4":"$5}else{print $3,$1":"$2":"$5":"$4}}' \
> /medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}.SNPID.txt

done


########################03.Update the id########################
for i in {1..22}
do

pfile=/medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}.setid
update=/medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}.SNPID.txt
out=/medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}.upid


sh /medpop/esp2/yang/project/PRS_pipeline/03.updateID.sh \
$plink2 $pfile $update $out


done



########################04.Calculate the score########################
weight="/medpop/esp2/yang/project/EC_PRS/weight"
SNPID="/medpop/esp2/yang/project/EC_PRS/UKB/SNPID"
mkdir -p $SNPID

mkdir -p /medpop/esp2/yang/project/EC_PRS/UKB/sscore_SQ



######There's four scores needed to calculate
awk  -v FS='\t' -v OFS="\t" 'NR == 1 || $7 == "X"' \
/medpop/esp2/yang/project/EC_PRS/weight/processed.hg37.all.weight.txt \
> /medpop/esp2/yang/project/EC_PRS/weight/EC.processed.hg37.all.weight.txt

awk  -v FS='\t' -v OFS="\t" 'NR == 1 || $8 == "X"' \
/medpop/esp2/yang/project/EC_PRS/weight/processed.hg37.all.weight.txt \
> /medpop/esp2/yang/project/EC_PRS/weight/Lipid.processed.hg37.all.weight.txt


awk  -v FS='\t' -v OFS="\t" 'NR == 1 || $9 == "X"' \
/medpop/esp2/yang/project/EC_PRS/weight/processed.hg37.all.weight.txt \
> /medpop/esp2/yang/project/EC_PRS/weight/Non-EC.processed.hg37.all.weight.txt


awk  -v FS='\t' -v OFS="\t" 'NR == 1 || $10 == "X"' \
/medpop/esp2/yang/project/EC_PRS/weight/processed.hg37.all.weight.txt \
> /medpop/esp2/yang/project/EC_PRS/weight/Non-EC_Non-Lipid.processed.hg37.all.weight.txt




for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do

cut -f 1 /medpop/esp2/yang/project/EC_PRS/weight/${pheno}.processed.hg37.all.weight.txt \
|tail -n+2 > /medpop/esp2/yang/project/EC_PRS/UKB/SNPID/${pheno}.SNPID.txt


for i in {1..22}
do

pfile=/medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}.upid
extract=/medpop/esp2/yang/project/EC_PRS/UKB/SNPID/${pheno}.SNPID.txt
update=/medpop/esp2/yang/project/EC_PRS/UKB/all/UKB.chr${i}.SNPID.txt
removesample=/medpop/esp2/yang/project/UKB_pheno/UKB_samplefailQC_FID_IID.txt
score=/medpop/esp2/yang/project/EC_PRS/weight/${pheno}.processed.hg37.all.weight.txt
score_col_nums=6
out=/medpop/esp2/yang/project/EC_PRS/UKB/sscore_SQ/${pheno}_chr${i}


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

cat /medpop/esp2/yang/project/EC_PRS/UKB/sscore_SQ/${pheno}_chr*.sscore \
|grep -v "#IID" > /medpop/esp2/yang/project/EC_PRS/UKB/sscore_SQ/chrall.${pheno}.all.sscore

Rscript /medpop/esp2/yang/project/cli_trail/script/combinechr.R \
/medpop/esp2/yang/project/EC_PRS/UKB/sscore_SQ/chrall.${pheno}.all.sscore \
/medpop/esp2/yang/project/EC_PRS/sscore_clean/UKB_sum.chrall.${pheno}.all.sscore

gzip -f /medpop/esp2/yang/project/EC_PRS/sscore_clean/UKB_sum.chrall.${pheno}.all.sscore

########And variants used in this score
cat /medpop/esp2/yang/project/EC_PRS/UKB/sscore_SQ/${pheno}_chr*.sscore.vars \
> /medpop/esp2/yang/project/EC_PRS/sscore_clean/UKB_sum.chrall.${pheno}_all.sscore.vars


rm /medpop/esp2/yang/project/EC_PRS/UKB/sscore_SQ/chrall.${pheno}.all.sscore


done


########################06.count number########################
echo -e "id\torignial\tprocessed\tscore" > /medpop/esp2/yang/project/EC_PRS/sscore_clean/UKB_variants_wc.txt


for pheno in EC Lipid Non-EC Non-EC_Non-Lipid
do


orignial=$(cat /medpop/esp2/yang/project/EC_PRS/weight/${pheno}.processed.hg37.all.weight.txt |tail -n+2 |wc -l)


processed=$(cat /medpop/esp2/yang/project/EC_PRS/weight/${pheno}.processed.hg37.all.weight.txt |tail -n+2 |wc -l)


score=$(cat /medpop/esp2/yang/project/EC_PRS/sscore_clean/UKB_sum.chrall.${pheno}_all.sscore.vars |wc -l)

echo -e "$pheno\t$orignial\t$processed\t$score" >> /medpop/esp2/yang/project/EC_PRS/sscore_clean/UKB_variants_wc.txt

done




