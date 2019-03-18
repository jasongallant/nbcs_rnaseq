# Neurobiología de la Conducta Social : RNA-seq
Dr. Jason Gallant, Michigan State University, USA

Dr. Vielka Salazar, Cape Breton University, Canada
****

## Getting Set Up:
All of what we need to do requires access to a remote server.  In order to complete the activity, you will need to bring your own laptop with the capability of connecting to via `ssh` to a server, and to run `x-windows` to look at graphical outputs.  Instructions for connecting to the server will be provided by your instructor.  You'll need to bring your laptop to particpate in the course activities.  Please set up your laptop following *in advance* using the instructions appropriate for your operating system below:

Contact Dr. Gallant (jgallant@msu.edu) with difficulties getting computers set up!

### Windows Users:
1. Download and Install:
https://mobaxterm.mobatek.net/download.html
2. Open MobaXTerm and type the following:
```bash
xeyes
```
3. You should see two eyes following your mouse as you move it.

### Mac/Linux Users:
1. Open Terminal.app and type the following:
```bash
xeyes
```
2. You should see two eyes following your mouse as you move it.
3. Mac Users: If you receive an error you need to install X-quartz.  Visit https://www.xquartz.org, download and install.  Repeat 1-2 above.

## Schedule and Course Materials:
*April 1, 2019*
1. [Dr. Gallant's Lecture Notes](introduction.md) - Background on sequencing biochemistry and what RNA-seq is all about
2. [Dr. Salazar's Lecture Notes]() - Nanopore Sequencing

*April 2, 2019*
1. [ONT Manual for Direct RNA Sequencing](direct-rna-sequencing-sqk-rna002-DRS_9080_v2_revB_22Nov2018.pdf)
2. [ONT Protocol for Direct RNA Sequencing](SQK-RNA002_protocol.pdf)

*April 3, 2019*

+ **14:30:** [Computing & Command Line Basics](computing.md) - A Quick Introduction to High Performance Computing and a Tutorial on the Command Line
+ **15:00:** [Practical #1](reads_and_qc.md) - An review of the RNA-seq workflow at "high speed", examining FASTQ Files and evaluating data quality
+ **15:30:**  Break - 15 min
+ **15:45:** [Practical #2](transcriptome_assembly.md) - An overview of de-novo transcriptome assembly using Trinity.  Examination of Trinity transcriptome and evaluation of outputs.
+ **16:15:** [Practical #3](read_alignment.md) - An overview of alignments, aligning Illumina reads and Oxford Nanopore reads
+ **16:45:** Break - 15 min
+ **17:00:** [Practical #4](nanopore_vs_illumina.md) - Comparing Illumina and Oxford Nanopore data.  Coverage, errors and other considerations.
+ **17:30:** [Practical #5](wrap-up.md) - Wrap-up and Review

### Sources and Resources:
+ [ANGUS Titus Brown - MSU](http://ged.msu.edu/angus/index.html)
+ [Eel Pond Tutorial](https://khmer-protocols.readthedocs.org/en/v0.8.4/mrnaseq/index.html)
+ [Trinity RNAseq](http://trinityrnaseq.sourceforge.net)
+ [Broad Institute Genome Bootcamp](http://www.broadinstitute.org/scientific-community/science/platforms/genome-sequencing/broadillumina-genome-analyzer-boot-camp)
+ [Colorado State Computational Biology and Genomics Workshop](https://dbsloan.github.io/TS2018/)
