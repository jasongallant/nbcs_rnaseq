# Practical on RNA-seq: From Test Tube to FASTQ File

## Review: How Illumina Data is Generated (5-10 min)

Though we already talked about it in lecture, this is a quick review of how Illumina data is created for sequencing:

1.	We start with the RNA of the tissue/organisim/developmental stage of interest
2.	Various quality checks of the RNA are performed  
3.	The RNA is fragmented  
4.	The RNA is converted to cDNA    
5.	The cDNA is ligated with adapter + a unique "barcode" (if pooling samples).  We later can use this barcode to unmix pooled samples  
6.	The cDNA is "size selected"-- for RNAseq this is ~ 250bp so that the 150bp reads are 'overlapping' (see below)  
7.	Various quality checks of the library (the name for the result of the above process) are run  
8.	Cluster generation occurs (cDNAs are hybridized to membrane, and bridge amplification leads to millions of clusters of clonally related DNAs)  
9.	Sequencing begins on first "strand" direction, 150bp (or however many cycles)  
10.	Sequencing begins on second "strand" direction, 150bp (or however many cycles)  
11.	Importantly the data is a series of "pictures" of the flow cell - must be converted into A's T's G's and C's by basecalling software (performed by sequencing center).  
12.	If pools are used, the data must be "demultiplexed", where data is then separated by barcodes  
13.	The "output" of this process are compressed, "FASTQ" files.


## FASTQ Files & Quality Control (30 min)

### Step 1: Login to the HPCC

To get the graphics displays to visualize the output, we need to connect to the HPCC in X-windows mode.

```bash
ssh -XY [classx]@hpcc.msu.edu
#enter password
ssh dev-intel18 # connect to an interactive 'development' node
cd module_1
```

### Step 2: Examining FASTQC Files
Let's take a look at one of these files.  The files are large, so we are only going to peek using a command line trick on our newly sequenced data that looks only at the first few lines of the data, using a new command syntax, combining `gzip` the pipe `|` and `head`.  Because there are so many reads, and because we often pay for storage space, compressing data is preferrable to uncompressed data.  Hence files are 'gzipped'.  In order to look inside of these files, we need to unzip them using `gzip -cd`.  The pipe character `|` is a way of taking the output of one command and feeding it to another command, in this case `head` which gives us the first few lines of a file:


```bash
gzip -cd 74_brain_S26_L002_R1_001.fastq.gz | head
```

This command returns something that looks like the following:

```
@K00392:151:H2CTWBBXY:2:1101:30949:1209 1:N:0:NGCCTCAT+NCTCTACT
NGCCTCGGCTGCTCAAAGCTGGGAACTCATCCAGGAACATCTGAAGGAAGTCCTGAGAAAGAATGAGTTGCCCATGGAATCTCATCACGAAGACTCTGAGGAAGCCCAAGAACTCGTTCTCGCCCACCTCACTGGCCAGGAACCTCAGTA
+
#AAF-FFJJ<JJ-<FJJFJ<FFFJFJ-<F7AF<F7JJJJFJA7JJAFJ<FAJJ<AJJJJJAJJJFJJ-AJJFJJ7FFJF-A-<A-<FJAAJFJJ7F77FFJJJFJJJJJJFFJ-AF-7FAJ-JJFAFJJJJJF7FFJJ-AFJJ<<FF7-7
@K00392:151:H2CTWBBXY:2:1101:22455:1226 1:N:0:NGCCTCAT+NCTCTACT
NTGATGTGCACGTTTCCCGTGTTTGGGCCTGAGCACCCCTCTTTGTGTTTAAAGTGCTTGAATGCCAGATTAAGGACTTCATGGACTAATGCCTCTGGAACAGGATGAAGAGGAATCTGTTTTAAAACTTCCACTGAAACTAAACAAAAG
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJFJJJJJJJFJJJJJJJJJJJFJJJJJJJJJJJJJJJJJFJJJJJJJJFJJJJJFJJJFJJFJJJJJJJJJJJJJ
@K00392:151:H2CTWBBXY:2:1101:22780:1226 1:N:0:NGCCTCAT+NCTCTACT
NGCGAGAGCCAGCTCCAGGTCCTGCTGCTCCTGGAGGCGGAGCCTACGAAGACGCTCCGCCTCCTCCCTCTTACTCTCACGCTTGGCCCATTCCAGCTGCTCCGTTTCACTTCCGAAGCCTCTAAGGTCGGTAAACGCGGCGGGAGCCAT
```

Things to Note about the FASTQ Format:
+ The first line always starts with @.  This contains the "metadata" describing the sequence (what cluster it comes from, what the barcode was, etc.)
+ The second line is the sequence of the read.
+ It is always followed by a + sign on the third line.
+ Finally, the fourth line describes the quality scores of each base.  The quality scores are on a scale from 0-93.  They are encoded by a single ASCII character that corresponds to a number from 0-93 to save space.

Indeed, there are always four lines for every read in a FASTQ file.  Thus, if we count the lines and divide by 4, we can quickly determine how many reads are in our file, which is useful for determining how big our job will be down the line.  

```bash
echo $(zcat 74_brain_S26_L002_R1_001.fastq.gz | wc -l)/4 | bc
```
+ Question: How many reads are in this file?  Can you explain how this command works?

Here's another example of how we can use pipes effectively.  In this case, we want all of the lines of the file, which is supplied by `zcat (cat for gzipped files)`, then we pipe this to `wc -l/4`.  wc counts normally counts words in a file or string, but with the `-l` it counts lines.  We divide this by 4 to get the number of reads.  `bc` is a program that forces the shell to do arithmetic rather than spit out the text of the calculation...

#### Quality Checks

To run Fastqc type:

```bash
module load FastQC/0.11.5-Java-1.8.0_162
fastqc 74_brain_S26_L002_R1_001.fastq.gz
```

The software will update you on its progress as it makes its way through the file.
Beware, this might take a while on very large files.  When FastQC is done with a file, it creates a directory (by default in the same directory as the input file, you can change this).  The directory created contains many output files:

```bash
ls 74_brain_pre_fastqc/
```
produces:

```
74_brain_S26_L002_R1_001_fastqc.html  74_brain_S26_L002_R1_001_fastqc.zip
```

The fastqc_report.html contains a graphical summary which you can open in a web browser:

[74_brain_S26_L002_R1_001_fastqc.html](https://efishgenomics.integrativebiology.msu.edu/2019workshop/74_brain_pre_fastqc/74_brain_S26_L002_R1_001_fastqc.html)


This produces a lot of information!  Let's walk through some of the highlights.  Examine the report and try to answer the following questions:

1. How many sequences are there?  Did we get it right with our quick calulation before?
2. How long were our sequences reads?
3. Which bases had the highest quality scores?
4. Which bases had the lowest quality scores?
5.  Noice the links on the left-- there are red 'Xs' , orange '!s' and green 'checks'.  What do these mean?
5. Is there anything of potential concern?


We have a couple of things to worry about:
1. Adapter content (mild at worst)
2. Per base sequence content and Kmer content (not really a concern, artifact of library construction, see link below)


Not sure why something looks wrong?  Read the manual!
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/

#### Review
+ FASTQC that there are some issues that need to be resolved: Quality scores are not universally great, and we also see some adapter contamination.  We can fix this by filtering our reads in the next section
+ Rather than repeat this time-consuming process, I've already done it for the other file.  Know that you *have* to run FASTQC seperately for both the file of forward and reverse reads

Looks like we have some work to do on these particular files!  Let's do some quality filtering!  First, a digression:


### What are adapters, anyway?


![Cluster Generation Overview](http://tucf-genomics.tufts.edu/images/faq02_pic01.jpg?1378237298)

Adapters are sequences that are "attached" to the cDNA sequences (representing the original mRNA sequences) by a process called <i>ligation</i>, which enable the molecule to hybridize to the flowcell as shown in the picture.  Typically, adapters also contain in-line indicies (barcodes), which enable pooled samples to later be computationally separated.  Indicies have already been automatically removed by the demultiplexing algorithm when your data was pre-processed after sequencing.  The adapters still have to be removed, and their sequence depend partly on the sequencing technology used.  We are typically performing our sequencing with TrueSeq chemistry on the HiSeq Illumina Platform. This uses adapters that typically look like this:

	>PrefixPE/1
	TACACTCTTTCCCTACACGACGCTCTTCCGATCT
	>PrefixPE/2
	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

Separating these out of your sequences could take a while!  Fortunately, there are several applications to choose from to help us with this:

To do this, we'll use the software [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). Note that this software is  baked into Trinity now, for convenience.  We're going to run this on our reads and inspect the results so we can peek under the hood.

### Run Trimmomatic to Remove Adapters

```bash
trimmomatic Command
````
#### Run Quality Checks on Quality Trimmed Data

To run Fastqc type:

```bash
module load FastQC/0.11.5-Java-1.8.0_162
fastqc 74_brain_S26_L002_R1_001.fastq.gz.PwU.qtrim.fq -o ./74_brain_post_fastqc/
```

[74_brain_S26_L002_R1_001_fastqc.html](https://efishgenomics.integrativebiology.msu.edu/2019workshop/74_brain_post_fastqc/74_brain_S26_L002_R1_001.fastq.gz.PwU.qtrim_fastqc.html)

How did we do?  does adapter contamination look better?


### Module 1 Review:
You just learned how to perform quality control!  There is a lot of art and discussion about the best procedures for quality control for RNA-seq data, and it largely depends on the application.  See MacManes (2014) for a nice discussion of this issue: https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full
