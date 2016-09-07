#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script runs the BQSR (1) step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_outdir prefix chr_list_str ncore mem"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  prefix: prefix of the input/output bam files. The input file names are assumed to be prefix.indel.chr.bam and the output file name will be prefix.bqsr."
  echo -e "  chr_list_str: comma-separated list of chromosomes (eg. 1,2,3,4,5,6,7,8,...)."
  echo -e "  ncore: number of CPUs (recommended: 2)"
  echo -e "  mem: memory (recommended: 8G)."
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
prefix=$3
chr_list_str=$4  ## 1,2,3,4,...
ncore=$5
mem=$6


if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi



chr_list=(${chr_list_str//,/ }) # convert to array
bam_list_str=''
for tmpchr in "${chr_list[@]}"
do
 bam_list_str='-I '$project_indir/$prefix.indel.$tmpchr.bam' '$bam_list_str
done

### step7 BQSR 
/usr/local/bin/jdk1.7.0_71/bin/java -Xmx$mem -Xms$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T BaseRecalibrator -R /resources/human_g1k_v37_decoy.fasta -knownSites /resources/dbsnp_138.b37.vcf $bam_list_str -o $project_outdir/$prefix.bqsr  -nct $ncore

