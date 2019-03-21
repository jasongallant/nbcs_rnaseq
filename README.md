# Neurobiología de la Conducta Social : RNA-seq
Dr. Jason Gallant, Michigan State University, USA

Dr. Vielka Salazar, Cape Breton University, Canada
****

## Getting Set Up:
All of what we need to do requires access to a remote server.  In order to complete the activity, you will need to bring your own laptop with the capability of connecting to via `ssh` to a server, and to run `x-windows` to look at graphical outputs.  Instructions for connecting to the server will be provided by your instructor.  You'll need to bring your laptop to participate in the course activities.  Please set up your laptop following *in advance* using the instructions appropriate for your operating system below:

Contact Dr. Gallant (jgallant@msu.edu) with difficulties getting computers set up!

# Set Up Instructions:
To enable courses to be taught as effectively as possible, you will need to be able to log in to the HPCC resources without difficulty before the class.  If you do not have an HPCC account, please install the appropriate software before the workshop but you will not be able to test that it works before the class. Please try to come 10-15 minutes early to get a temporary account, and test that you can log in.

## File Transfer
If you are taking the Introduction to HPCC/Computational Resources at MSU course you will  need to install a file transfer client. Ex. FileZilla (https://filezilla-project.org/) or  Cyberduck (https://cyberduck.io/?l=en)

## Connecting to HPCC
For Mac users:
Please install the latest version of Xquartz (you can follow the directions here: https://wiki.hpcc.msu.edu/display/hpccdocs/Installing+an+X-server++for+Macs)

## To test that it is working properly:
Open a terminal window (if you do not have the terminal program on your taskbar, you can search for terminal), and
Type:

```bash
ssh -XY [classx]@hpcc.msu.edu #replace classx with the ID of your temporary account [classx]
# type your the password provided for the classx account
```
(type your msu password or the password provided for the temporary account)  This should give you a command prompt.  If you type:

```bash
xeyes
```
You should get a separate window that displays a pair of eyes that follow the mouse around.  If this works, that should be all you need to be ready for the course.  If not, please come several minutes early and we will try to help you get online.


## For WINDOWS users:
We recommend installing Moba Xterm "home edition" (http://mobaxterm.mobatek.net/download-home-edition.html).   Please download the "Home Edition, Installer Version unless you are comfortable installing portable applications yourself.  This program provides both an SSH client (command line) and an X-Windows server system for running Graphical User Interface programs remotely (X11).

Once it is installed, you can open it, and try to connect with the hpcc:

To test that it is working properly:
Open a terminal window (if you do not have the mobaxterm shortcut set up, you can search for mobaxterm), and type:

```bash
ssh -XY yournetid@hpcc.msu.edu # replace yournetid with YOUR ACTUAL NETID or the ID of your temporary account [classx]
# type your the password provided for the classx account
```
This should give you a command prompt.  If you type:

```bash
xeyes
```

You should get a separate window that displays a pair of eyes that follow the mouse around.

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
+ **16:15:** [Practical #3](alignment.md) - An overview of alignments, aligning Illumina reads and Oxford Nanopore reads
+ **16:45:** Break - 15 min
+ **17:00:** [Practical #4](nanopore_vs_illumina.md) - Comparing Illumina and Oxford Nanopore data.  Coverage, errors and other considerations.
+ **17:30:** Wrap-up and Review

### Sources and Resources:
+ [ANGUS Titus Brown - MSU](http://ged.msu.edu/angus/index.html)
+ [Eel Pond Tutorial](https://khmer-protocols.readthedocs.org/en/v0.8.4/mrnaseq/index.html)
+ [Trinity RNAseq](http://trinityrnaseq.sourceforge.net)
+ [Broad Institute Genome Bootcamp](http://www.broadinstitute.org/scientific-community/science/platforms/genome-sequencing/broadillumina-genome-analyzer-boot-camp)
+ [Colorado State Computational Biology and Genomics Workshop](https://dbsloan.github.io/TS2018/)
