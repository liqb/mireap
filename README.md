mireap
======

Indentify microRNA genes from small RNA sequencing data

Usage

Version

Author


Program: MIREAP (Reap miRNAs from deeply sequenced smRNA library)
Version: 0.2
Contact: Li Qibin <liqb@genomics.org.cn>
  Bioinformatics department, Beijing Genomics Institute


1. Introduction
MIREAP combines small RNA position and depth with a model of
microRNA biogenesis to discover microRNAs from deeply sequenced
small RNA library.

2. Installation
You must have Vienna RNA Package (http://www.tbi.univie.ac.at/RNA)
installed on your computer and make sure that its perl interface
is accessible.

Copy mireap_0.2.tar.gz to a directory (/foo/bar) and unpack
it by command:
  tar -zxvf mireap_0.2.tar.gz

Before running mireap, you need add path /foo/bar/mireap_0.1/lib to
environment variable PERL5LIB:
For csh/tcsh:
  setenv PERL5LIB /foo/bar/mireap_0.2/lib
For sh/ksh/bash:
  export PERL5LIB=/foo/bar/mireap_0.2/lib


3. Usage
mireap.pl -i <smrna.fa> -m <map.txt> -r <reference.fa> -o <outdir>
Options:
-i <file>  Small RNA library, fasta format, forced
-m <file>  Mapping file, tabular format, forced
-r <file>  Reference file, fasta format, forced
-o <dir>   Directory where results produce (current directory)
-t <str>   Sample label (xxx)
-A <int>   Minimal miRNA sequence length (18)
-B <int>   Maximal miRNA sequence length (26)
-a <int>   Minimal miRNA reference sequence length (20)
-b <int>   Maximal miRNA reference sequence length (24)
-u <int>   Maximal copy number of miRNAs on reference (20)
-e <folat> Maximal free energy allowed for a miRNA precursor (-18)
-d <int>   Maximal space between miRNA and miRNA* (35)
-p <int>   Minimal base pairs of miRNA and miRNA*
-v <int>   Maximal bulge of miRNA and miRNA* (4)
-s <int>   Maximal asymmetry of miRNA/miRNA* duplex
-f <int>   Flank sequence length of miRNA precursor (10)
-h         Help

Please convert your small RNA file into fasta format and append
sequencing frequence to sequence Id, just like this entry:
>t0000035 3234
GAATGGATAAGGATTAGCGATGATACA
(t0000035 is read_ID, 3234 is sequencing frequence)

The format of small RNA mapping file should be (delimited by tab or
space):
read_ID,chr_ID,start,end,strand(+/-)

You can make MIREAP run on the test data by execute comand:
perl ../bin/mireap.pl -i rna.fa -m map.txt -r ref.fa


4. Output format
MIREAP produce three files at each run.

*.gff
This file contains miRNA genes discovered by MIREAP, GFF3 format. For
GFF3 format, please refer to http://www.sequenceontology.org/gff3.shtml
Attribute 'Count' denotes the sequenceing frequence.

*.aln
This file contains sequence and structure of the pre-miRNA. Small RNAs
also are aligned to the precursor from which you can get more insights
into the maturation process of miRNAs.

*.log
This log file records parameters, start end time and other informations.

