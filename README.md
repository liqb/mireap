mireap: discover microRNA genes from small RNA sequencing reads
======

## Synopsis

MIREAP combines small RNA position and depth with a model of
microRNA biogenesis to discover microRNAs from deeply sequenced
small RNA libraries.

## Installation
Install Vienna RNA Package (http://www.tbi.univie.ac.at/RNA) on your computer and make sure that its perl interface
(RNA.pm) is accessible.

Copy mireap*.tar.gz to a target directory (/foo/bar) and unpack it:
```
tar -zxvf mireap*.tar.gz
```
Before running mireap, you need add path /foo/bar/mireap*/lib to environment variable PERL5LIB:
```
In csh/tcsh:
setenv PERL5LIB /foo/bar/mireap_0.2/lib
```
```
In sh/ksh/bash:
  export PERL5LIB=/foo/bar/mireap_0.2/lib
```

## Usage
```
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
```
Please convert your small RNA file into fasta format and append
sequencing frequence to sequence Id, just like this entry:
```
>t0000035 3234
GAATGGATAAGGATTAGCGATGATACA
```
t0000035 is tagID, 3234 is read count.

The format of small RNA mapping file should be (delimited by tab or
space):
read_ID, chr_ID, start, end, strand(+/-)

Output format
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


## Tests

You can make a test run using test data by command:
```
perl ../bin/mireap.pl -i rna.fa -m map.txt -r ref.fa
```

## Contributors

QIBIN LI <liqb@genomics.cn>



## License

A short snippet describing the license (MIT, Apache, etc.)


## BUGS

