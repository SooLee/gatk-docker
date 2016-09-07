##### These are source files used for a docker container for bwa-gatk-based variant calling pipeline.
##### Authors: Soo Lee (duplexa@gmail.com) & Daniel Kwon (daniel.minseok.kwon@gmail.com)


## Docker image
The docker image is stored as duplexa/gatk_env:v1 on hub.docker.com.


## Prerequisites

You need the resource files in aws S3://maestro-resources/ to be mounted to the docker container as /resources/.


## Steps 

Each of the following shell scripts executes a step of the bwa-gatk-based variant calling pipeline. (In the order displayed)


#### split_fastq.sh
```
This script run the Split_fastq step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/split_fastq.sh fastq_gz project_outdir mate

Arguments:
  fastq_gz: path to the gzipped fastq file (eithe R1 or R2)
  project_outdir: directory in which output files are placed.
  mate: 'R1' or 'R2'. Must match the mate of the fastq file.
```


#### align_sort_addrg.sh 
```
This script run the Alignmet, sort and addRG steps of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/align_sort_addrg.sh project_indir project_outdir split_prefix RGID RGLB RGSM ncore mem

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  split_prefix: input file prefix (samplename + split_tag) (eg. TEST.MG01HX08_123_HK3GMCCXX_1). Input file names are assumed to be split_prefix.R[12].fastq.gz.
  RGID: lane name (eg. L1).
  RGLB: lane ID (eg. MG01HX08_123_HK3GMCCXX_1).
  RGSM: sample ID (eg. S1).
  ncore: number of CPUs (recommended: 4).
  mem: memory (recommended: 2G).
```

#### rmdup.sh
```
This script run the Mark Duplicate step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/rmdup.sh project_indir project_outdir split_prefix_list_str prefix mem

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  split_prefix_list_str: comma-separated list of input file prefix (samplename + split_tag) (eg. TEST.MG01HX08_123_HK3GMCCXX_1). Input file names are assumed to be split_prefix.addrg.bam.
  prefix: prefix of the output bam file. The output file name will be prefix.mkdup.bam.
  mem: memory (recommended: 32G).
```

#### realign.sh
```
This script run the Indel Realignment step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/realign.sh project_indir project_outdir prefix chr ncore mem

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  prefix: prefix of the input/output bam files. The input file name is assumed to be prefix.mkdup.bam and the output file name will be prefix.indel.chr.bam.
  chr: chromosome (eg. 21, or 'decoy', or 'unmapped').
  ncore: number of CPUs (recommended: 2)
  mem: memory (recommended: 8G).
```

#### bqsr1.sh 
```
This script run the BQSR (1) step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/bqsr1.sh project_indir project_outdir prefix chr_list_str ncore mem

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  prefix: prefix of the input/output bam files. The input file names are assumed to be prefix.indel.chr.bam and the output file name will be prefix.bqsr.
  chr_list_str: comma-separated list of chromosomes (eg. 1,2,3,4,5,6,7,8,...).
  ncore: number of CPUs (recommended: 2)
  mem: memory (recommended: 8G).
```

#### bqsr2.sh 
```
This script run the BQSR (2) step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/bqsr2.sh project_indir project_bqsr_indir project_outdir prefix chr ncore mem

Arguments:
  project_indir: directory in which input files are placed.
  project_bqsr_indir: directory in which input bqsr files are placed.
  project_outdir: directory in which output files are placed.
  prefix: prefix of the input/output bam/bqsr files. The input file names are assumed to be prefix.indel.chr.bam and prefix.bqsr and the output file name will be prefix.final.chr.bam.
  chr: chromosome (eg. 21, or 'decoy', or 'unmapped').
  ncore: number of CPUs (recommended: 2)
  mem: memory (recommended: 8G).
```

#### merge_finalbam.sh 
```
This script run the Merge final bam step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/merge_finalbam.sh project_indir project_outdir prefix chr_list_str

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  prefix: prefix of the output bam file. The output file name will be prefix.mkdup.bam.
  chr_list_str: comma-separated list of chromosomes to be merged. (1,2,3,...,MT,decoy,unmapped)
```

#### gvcf.sh
```
This script run the Generate GVCF step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/gvcf.sh project_indir project_outdir prefix chr region ncore mem

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  prefix: prefix of the input/output bam/vcf files. The input file names are assumed to be prefix.final.chr.bam and the output file name will be prefix.region.g.vcf.
  chr: chromosome (eg. 21).
  region: region (eg. 21:1-50000). The region must match the chromosome.
  ncore: number of CPUs (recommended: 2)
  mem: memory (recommended: 4G).
```

#### merge_gvcf.sh 
```
This script run the Merge GVCF step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/merge_gvcf.sh project_indir project_outdir prefix region_list_str mem

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  prefix: prefix of the input/output vcf files. The input file names are assumed to be prefix.region.g.vcf and the output file name will be prefix.hc.raw.g.vcf.
  region_list_str: comma-separated list of regions (eg. 21:1-50000).
  mem: memory (recommended: 2G).
```

#### hc.sh 
```
This script run the Haplotype Caller step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/hc.sh project_indir project_outdir prefix_list_str group_prefix region ncore mem

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  prefix_list_str: comma-separated list of prefices (samples). The input file names are assumed to be prefix.hc.raw.g.vcf.
  group_prefix: prefix of the output vcf file. This represents a group of samples to be jointly called. The output file name will be group_prefix.hc.geno.region.g.vcf.
  region: region (eg. 21:1-50000).
  ncore: number of CPUs (recommended: 2).
  mem: memory (recommended: 4G).
```

#### combine_hc.sh 
```
This script run the Combine_hc_gvcf step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/combine_hc.sh project_indir project_outdir group_prefix region_list_str mem

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  group_prefix: prefix of the input/output vcf file. This represents a group of samples jointly called. The input file names are assumed to be group_prefix.hc.geno.region.g.vcf and the output file name will be group_prefix.hc.geno.g.vcf.
  region_list_str: comma-separated regions (eg. 21:1-50000,22:50001-100000,...).
  mem: memory (recommended: 2G).
```

#### vqsr.sh 
```
This script run the VQSR step of a variant calling pipeline based on bwa-gatk.

Usage:
  /usr/local/bin/run_scripts/vqsr.sh project_indir project_outdir group_prefix ncore mem

Arguments:
  project_indir: directory in which input files are placed.
  project_outdir: directory in which output files are placed.
  group_prefix: prefix of the input/output vcf file. This represents a group of samples jointly called. The input file name is assumed to be group_prefix.hc.geno.g.vcf and the output files name will be group_prefix.hc.snp.vqsr.vcf and group_prefix.hc.indel.vqsr.vcf.
  ncore: number of CPUs (recommended: 2).
  mem: memory (recommended: 15G).
```

