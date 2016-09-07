#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script runs the Mark Duplicate step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_outdir split_prefix_list_str prefix mem"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  split_prefix_list_str: comma-separated list of input file prefix (samplename + split_tag) (eg. TEST.MG01HX08_123_HK3GMCCXX_1). Input file names are assumed to be split_prefix.addrg.bam."
  echo -e "  prefix: prefix of the output bam file. The output file name will be prefix.mkdup.bam."
  echo -e "  mem: memory (recommended: 32G)."
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
split_prefix_list_str=$3  ## list of input prefix (sample+split_tag), comma-separated (prefix1,prefix2, ...)
prefix=$4  ## prefix for the output bam file
mem=$5


if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi


split_prefix_list=(${split_prefix_list_str//,/ }) # convert to array
bam_list_str=''
for tmpsplit_prefix in "${split_prefix_list[@]}"
do
  bam_list_str='I='$project_indir/$tmpsplit_prefix.addrg.bam' '$bam_list_str
done


### step5 mark duplicates
/usr/local/bin/jdk1.7.0_71/bin/java -Xmx$mem -Xms$mem -jar /usr/local/bin/picard-1.130/picard.jar MarkDuplicates M=$project_outdir/$prefix.out.dup.matrix COMPRESSION_LEVEL=1 ASSUME_SORTED=True VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false $bam_list_str O=$project_outdir/$prefix.mkdup.bam  TMP_DIR=/tmp
/usr/local/bin/samtools-1.3/bin/samtools index $project_outdir/$prefix.mkdup.bam


