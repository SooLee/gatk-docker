#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script runs the Combine_hc_gvcf step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_outdir group_prefix region_list_str mem"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  group_prefix: prefix of the input/output vcf file. This represents a group of samples jointly called. The input file names are assumed to be group_prefix.hc.geno.region.g.vcf and the output file name will be group_prefix.hc.geno.g.vcf."
  echo -e "  region_list_str: comma-separated regions (eg. 21:1-50000,22:50001-100000,...)."
  echo -e "  mem: memory (recommended: 2G)."
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
group_prefix=$3
region_list_str=$4  ## 2:1-50000,2:50001-100000,...
mem=$5


if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi



region_list=(${region_list_str//,/ }) # convert to array
vcf_list_str=''
for tmpregion in "${region_list[@]}"
do
 vcf_list_str='--variant '$project_indir/$group_prefix.hc.geno.$tmpregion.g.vcf' '$vcf_list_str
done


### HaplotypeCaller (Simple example)
#### step4. combine geno.vcf
java -Xmx$mem -Xms$mem -cp /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants  -R /resources/human_g1k_v37_decoy.fasta $vcf_list_str -out $project_outdir/$group_prefix.hc.geno.g.vcf --assumeSorted -l INFO


