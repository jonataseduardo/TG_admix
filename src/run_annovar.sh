#!/bin/sh

#PBS -N run_annovar 
#PBS -l nodes=1:ppn=1
#PBS -l mem=24gb
#PBS -l walltime=24:00:00
#PBS -q long
#PBS -j oe
#PBS -o $HOME/TG_admix/data/TG_phase3/log_files/continental_vcf.log

ROOT_P=/raid/genevol/1kg/phase3/data/
VCF_P=phase3_chr/ALL.chr
VCF_S=.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

#ROOT_P=""
#VCF_P=$HOME/TG_admix/src/teste_chr
#VCF_S=.vcf.gz

OUT_P=$HOME/TG_admix/data/TG_phase3/annotations/
POP_P=/raid/genevol/1kg/phase3/data/phase3_pops/pop_noadmix
POP=ALL

mkdir -p $OUT_P

#CHR=$PBS_ARRAYID
CHR=21

NAME_OUT=$OUT_P$POP"_chr"$CHR"_bial_nogen"

#$HOME/bin/bcftools view -G \
#  --max-alleles 2  \
#  --output-file $NAME_OUT.vcf \
#  --output-type v \
#  $ROOT_P$VCF_P$CHR$VCF_S

$HOME/annovar/table_annovar.pl \
  $NAME_OUT.vcf \
  $HOME/annovar/humandb/ \
  -buildver hg19 \
  -out $OUP_P$POP"_chr"$CHR \
  -remove \
  -protocol refGene,cytoBand,avsnp147,dbnsfp30a \
  -dot2underline\
  -operation g,r,f,f \
  -nastring . \
  -vcfinput

#gzip -f $NAME_OUT.vcf
#
#tabix -f $NAME_OUT.vcf.gz
#
#mv $POP"_chr"$CHR.hg19_multianno.txt $OUT_P
#rm $OUP_P$POP"_chr"$CHR*
