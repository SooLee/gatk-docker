#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script runs the VQSR step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_outdir group_prefix ncore mem"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  group_prefix: prefix of the input/output vcf file. This represents a group of samples jointly called. The input file name is assumed to be group_prefix.hc.geno.g.vcf and the output files name will be group_prefix.hc.snp.vqsr.vcf and group_prefix.hc.indel.vqsr.vcf."
  echo -e "  ncore: number of CPUs (recommended: 2)."
  echo -e "  mem: memory (recommended: 15G)."
  echo -e ""
}

# if no arguments supplied, display usage 
if [  $# -le 0 ]
then
  display_usage
  exit 0
fi

# check whether user had supplied -h or --help . If yes display usage 
if [[ $# == 1 && ( $1 == "--help" ||  $1 == "-h" )]]
then
  display_usage
  exit 0
fi

####################

project_indir=$1
project_outdir=$2
group_prefix=$3   # this corresponds to group_prefix
ncore=$4
mem=$5


if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi




### HaplotypeCaller (Simple example)
#### step5. VQSR
# SNP
java -Xms$mem -Xmx$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T SelectVariants -R /resources/human_g1k_v37_decoy.fasta --variant $project_indir/$group_prefix.hc.geno.g.vcf -o $project_outdir/$group_prefix.hc.snp.recal.vcf -selectType SNP -selectType MNP -dt NONE -nt $ncore
java -Xms$mem -Xmx$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T VariantRecalibrator -R /resources/human_g1k_v37_decoy.fasta -input $project_outdir/$group_prefix.hc.snp.recal.vcf --resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 /resources/hapmap_3.3.b37.vcf --resource:omni,VCF,known=false,training=true,truth=true,prior=12.0 /resources/1000G_omni2.5.b37.vcf --resource:1000G,known=false,training=true,truth=false,prior=10.0 /resources/1000G_phase1.snps.high_confidence.b37.vcf --resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 /resources/dbsnp_138.b37.vcf --use_annotation ReadPosRankSum --use_annotation FS --use_annotation MQ --use_annotation DP  --mode SNP  --recal_file $project_outdir/$group_prefix.hc.snp.recal --tranches_file $project_outdir/$group_prefix.hc.snp.recal.tranches --rscript_file $project_outdir/$group_prefix.hc.snp.recal.plots.R -dt NONE -nt $ncore
java -Xms$mem -Xmx$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T ApplyRecalibration -R /resources/human_g1k_v37_decoy.fasta -input $project_outdir/$group_prefix.hc.snp.recal.vcf --ts_filter_level 99.0 --mode SNP --tranches_file $project_outdir/$group_prefix.hc.snp.recal.tranches --recal_file $project_outdir/$group_prefix.hc.snp.recal --out $project_outdir/$group_prefix.hc.snp.vqsr.vcf -dt NONE -nt $ncore

# indel
java -Xms$mem -Xmx$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T SelectVariants -R /resources/human_g1k_v37_decoy.fasta --variant $project_indir/$group_prefix.hc.geno.g.vcf -o $project_outdir/$group_prefix.hc.indel.recal.vcf -selectType INDEL -dt NONE -nt $ncore
java -Xms$mem -Xmx$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T VariantRecalibrator -R /resources/human_g1k_v37_decoy.fasta -input $project_outdir/$group_prefix.hc.indel.recal.vcf --resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 /resources/Mills_and_1000G_gold_standard.indels.b37.vcf --use_annotation QD --use_annotation FS --use_annotation ReadPosRankSum --mode INDEL --recal_file $project_outdir/$group_prefix.hc.indel.recal --tranches_file $project_outdir/$group_prefix.hc.indel.recal.tranches --rscript_file $project_outdir/$group_prefix.hc.indel.recal.plots.R -dt NONE -nt $ncore
java -Xms$mem -Xmx$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T ApplyRecalibration -R /resources/human_g1k_v37_decoy.fasta -input $project_outdir/$group_prefix.hc.indel.recal.vcf --ts_filter_level 99.0 --mode INDEL --tranches_file $project_outdir/$group_prefix.hc.indel.recal.tranches --recal_file $project_outdir/$group_prefix.hc.indel.recal --out $project_outdir/$group_prefix.hc.indel.vqsr.vcf -dt NONE -nt $ncore




