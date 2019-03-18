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

	ssh -X USERNAME@hpcc.msu.edu
	**enter password**

### Step 2: Examining FASTQC Files
Let's take a look at one of these files.  The files are HUGE, so we are only going to peek using a command line trick on our newly sequenced data that looks only at the first few lines of the data, using a new command syntax, combining `gzip` the pipe `|` and `head`.  Because there are so many reads, and because we often pay for storage space, compressing data is preferrable to uncompressed data.  Hence files are 'gzipped'.  In order to look inside of these files, we need to unzip them using `gzip -cd`.  The pipe character `|` is a way of taking the output of one command and feeding it to another command, in this case `head` which gives us the first few lines of a file:

	gzip -cd 468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz | head

This command returns something that looks like the following:

	@HWI-ST1348:49:D1C9EACXX:3:1101:1246:1978 1:N:0:ATCACG
	NCAACAACATCCTAAGTAATTCCCATCAGATACCGATGTCCTGCACGGGTGCGAATCCCCCCATCAGTCCCCCCACTGGGACCCTGCTGGACAGGAAGGC
	+
	#4BDFFFFHHHHHJJJHIJJJJJJJJJJJJIJJJJJJIJJJJJJJJJJJ@HIIJJHHHHFFDDDDDDEDDDDDDDDDDDDDDDDDDDDDDDBDDDDDDDD
	@HWI-ST1348:49:D1C9EACXX:3:1101:1850:1980 1:N:0:ATCACG
	NGGGGTTTGGTATTGGGAGATTGCTGGGGGTTTTATGTTGATGATTGTGGTGATGAAGTTGATAGAGCCAAGGATGGAGGAAACACCGGCTAAGTGGAGG
	+
	#4=DDADDHHDFHIJJJGIGIJJJJJJJJJDHIJGIIIJIEIJIJIJIJJ@DEHHFHHHHFDFFFFFEECEDDBDDDDDDBBDDDDDDDDDDDDCDDDDD
	@HWI-ST1348:49:D1C9EACXX:3:1101:1940:1994 1:N:0:ATCACG
	NTTCTGGAGTTGGTCCTTGAGAATTTCGTGTATCCCTGGTACAGGTATTTGTCGCAGCTTTCCTGTGATTTGTGTGCCTCTCAATACACACAGATGTTAT

Things to Note about the FASTQ Format:
+ The first line always starts with @.  This contains the "metadata" describing the sequence (what cluster it comes from, what the barcode was, etc.)
+ The second line is the sequence of the read.
+ It is always followed by a + sign on the third line.
+ Finally, the fourth line describes the quality scores of each base.  The quality scores are on a scale from 0-93.  They are encoded by a single ASCII character that corresponds to a number from 0-93 to save space.

Indeed, there are always four lines for every read in a FASTQ file.  Thus, if we count the lines and divide by 4, we can quickly determine how many reads are in our file, which is useful for determining how big our job will be down the line.  

	cat 468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz | echo $((`wc -l`/4))

Here's another example of how we can use pipes effectively.  In this case, we want all of the lines of the file, which is supplied by `cat`, then we pipe this to `wc -l/4`.  wc counts normally counts words in a file or string, but with the `-l` it counts lines.  We divide this by 4 to get the number of reads.  `echo` is the LINUX print command.


#### Quality Checks

To get the graphics displays to visualize the output, we need to connect to the HPCC in X-windows mode.

	ssh -X hpcc
	**enter password**

To run Fastqc type:

	fastqc yourfile.fastq.gz

The software will update you on its progress as it makes its way through the file.
Beware, this might take a while on very large files.  When FastQC is done with a file, it creates a directory (by default in the same directory as the input file, you can change this).  The directory created contains many output files:

	fastqc_data.txt
	fastqc_report.html
	Icons
	Images
	summary.txt

The fastqc_report.html contains a graphical summary which you can open in a web browser.

	firefox fastqc_report.html

Looks like we have some work to do on these particular files!  This is a good way of scanning to see if there are massive problems with your library.  Typically, I don't get too worried about this, because we still have a lot of cleaning to do!

**Question:**
+ How do we figure out how to tell fastqc where to put the files if we had many reports to run?
+ If it isn't obvious, how would you figure it out?

#### Review
FASTQC that there are some issues that need to be resolved: Quality scores are not universally great, and we also see some adapter contamination.

### What are adapters, anyway?

![Cluster Generation Overview](http://tucf-genomics.tufts.edu/images/faq02_pic01.jpg?1378237298)

Adapters are sequences that are "attached" to the cDNA sequences (representing the original mRNA sequences) by a process called <i>ligation</i>, which enable the molecule to hybridize to the flowcell as shown in the picture.  Typically, adapters also contain in-line indicies (barcodes), which enable pooled samples to later be computationally separated.  Indicies have already been automatically removed by the demultiplexing algorithm when your data was pre-processed after sequencing.  The adapters still have to be removed, and their sequence depend partly on the sequencing technology used.  We are typically performing our sequencing with TrueSeq chemistry on the HiSeq Illumina Platform. This uses adapters that typically look like this:

	>PrefixPE/1
	TACACTCTTTCCCTACACGACGCTCTTCCGATCT
	>PrefixPE/2
	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

Separating these out of your sequences could take a while!  Fortunately, there are several applications to choose from to help us with this:

A few that we've used with success in the lab are:
1. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) A Java Based Trimmer
2. [Sickle](https://github.com/najoshi/sickle)
3. [Scythe](https://github.com/vsbuffalo/scythe)
4. [FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)

For this tutorial, we are going to use a Trimmomatic to clean up our reads.

### Run Trimmomatic to Remove Adapters
