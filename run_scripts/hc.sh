#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script runs the Haplotype Caller step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_outdir prefix_list_str group_prefix region ncore mem"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  prefix_list_str: comma-separated list of prefices (samples). The input file names are assumed to be prefix.hc.raw.g.vcf."
  echo -e "  group_prefix: prefix of the output vcf file. This represents a group of samples to be jointly called. The output file name will be group_prefix.hc.geno.region.g.vcf."
  echo -e "  region: region (eg. 21:1-50000)."
  echo -e "  ncore: number of CPUs (recommended: 2)."
  echo -e "  mem: memory (recommended: 4G)."
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
prefix_list_str=$3  ## comma delimited (each prefix corresponds to a sample)
group_prefix=$4  ## prefix for the group of prefices (samples)
region=$5  ## 2:1-50000
ncore=$6
mem=$7


if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi


prefix_list=(${prefix_list_str//,/ }) # convert to array
vcf_list_str=''
for tmpprefix in "${prefix_list[@]}"
do
 vcf_list_str='--variant '$project_indir/$tmpprefix.hc.raw.g.vcf' '$vcf_list_str
done


### HaplotypeCaller (Simple example)
#### step3. joint genotyping
java -Xmx$mem -Xms$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /resources/human_g1k_v37_decoy.fasta --dbsnp /resources/dbsnp_138.b37.vcf $vcf_list_str -o $project_outdir/$group_prefix.hc.geno.$region.g.vcf -l INFO -stand_call_conf 30.0 -stand_emit_conf 10.0 -dt NONE -rf BadCigar -nt $ncore -L $region

