# Practical on RNA-seq: From Test Tube to FASTQ File

## How Illumina Data is Generated
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
12.	If pools are used, the data must be "demultiplexed", where data is then seperated by barcodes  
13.	The "output" of this process are compressed, "FASTQ" files.

## FASTQ Files
In the interests of time and practical constraints, we have had steps 1-13 performed for us (though we saw last time how they were performed, when we prepared our Oxford Nanopore run)

Let's take a look at one of these files.  The files are HUGE, so we are only going to peek using a commandline trick on our newly sequenced data that looks only at the first few lines of the data:

	gzip -cd 468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz | head

This command returns the following information:

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

Indeed, there are always four lines for every read in a FASTQ file.  Thus, if we count the lines and divide by 4, we can quickly determine how many reads are in our file, which is useful for determining how big our job will be down the line:

	cat 468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz | echo $((`wc -l`/4))


## Quality Checks

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

Or you can simply look at the results (without graphics) by peeking at the summary.txt

	cat summary.txt

Which will look like this:

	PASS	Basic Statistics	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	PASS	Per base sequence quality	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	PASS	Per sequence quality scores	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	FAIL	Per base sequence content	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	FAIL	Per base GC content	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	PASS	Per sequence GC content	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	WARN	Per base N content	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	PASS	Sequence Length Distribution	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	FAIL	Sequence Duplication Levels	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	WARN	Overrepresented sequences	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz
	WARN	Kmer Content	468_3676_4097_N_MNIG_6656_ATCACG_R1.fastq.gz

Looks like we have some work to do on these particular files!  This is a good way of scanning to see if there are massive problems with your library.  Typically, I don't get too worried about this, because we still have a lot of cleaning to do!

<b> Teachable Moment: </b> How do we figure out how to tell fastqc where to put the files?  Most (well written) software has a nifty manual built right in!  Typically you can type 'man' in front of most commands to see it, while others require you to type an -h flag after the command or -help or --help.  In the case of fastqc:

	fastqc --help

Gives us what we are looking for:

	FastQC - A high throughput sequence QC analysis tool

	SYNOPSIS

		fastqc seqfile1 seqfile2 .. seqfileN

	    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
	           [-c contaminant file] seqfile1 .. seqfileN

	 ...

So, we now can see that if we type

	fastqc yourfile.fastq.gz -o /your/directory/for/output/

FASTQC will put the files where we like them, rather than where it feels like.

### Review
Yesterday, we discovered with FASTQC that there are some issues that need to be resolved: sometimes bases on the ends of the reads are of lower quality.  We've also seen some adapter contamination as well.

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

For this tutorial, we are going to use a combination of Trimmomatic and FASTX toolkit, which are pretty simple to operate.

### Run Trimmomatic to Remove Adapters

### Run FASTX Tools to Quality Filter Reads

### Remove Orphaned Reads from Filtering Process

Now, we should be fairly well set up to run Trimmomatic on our data files.  Here's where our earlier organization scheme becomes relevant.  Trimmomatic needs to know which files go together in order to operate-- you can run Trimmomatic manually on each pair, but that becomes fairly tedious to do.  If you are very careful about your filenaming scheme, you should be able to write a shell script that will pair the files together.  My preferred method is the textfile that I created earlier.

### Review
We've now learned how sequencing data is generated, and what potential problems result from the sequencing process.  In the last tutorial, we learned about adapter trimming and quality filtering.  Now you should have a nice, clean dataset to work with!  Today, we are going to proceed with performing the actual assembly of the transcriptome data.

### Transcriptome Assembly

What exactly is a <i> transcriptome </i> anyway?  Gene expression is a very dynamic process-- a transcriptome can be thought of as a
"snapshot" of this process, or as a collection of all of the genes that are expressed in a particular tissue at a particular time.  An example might be brain during early development, or a skeletal muscle after exercise.  The sky is the limit!  Our sequenced mRNA is a sample of the expressed genes.

Remember that technologies do not presently exist to sequence entire mRNA molecules, so we have to break them into small (presently ~300bp) sequences and sequence them in vast quantities.  This certainly allows us to generate huge amounts of data!  We ultimately hope to quantify the gene expression in or original sample.

This is pretty easy to do if you have a genome for your organisim of interest-- you simply align the cleaned reads to the genome and voila,  the number of reads aligning to each gene sequence correspond to its expression level.  More on this later.

But what if you don't know what the genes look like?  Unfortunately, this is the case for many organisims that haven't had a genome sequenced yet.  We are left with the enormous problem of putting everything back together again (creating a "reference" transcriptome), before we can perform our quantification. This requires us to perform what is called a <i>de-novo</i> transcriptome assembly.

A  computationally efficient method of putting the transcriptome back together is also a method that appeals to common sense-- this is called the "DeBruin Graph" approach.  This method has several steps, but we will walk through the basic process:

![Assembly Overview](http://www.nature.com/nrg/journal/v12/n10/images/nrg3068-f3.jpg)

1. Each read is broken down into overlapping strings of length "k"  (FYI, a k-mer is the name for a string of length k.)
2. An ordered list of kmers is generated (hash table)-- each k-mer overlaps by (k-1 bases)
3. Each <b>unqiue</b> kmer is assigned a node on a "graph"-- overlapping k-mers are adjacent to one-another.
4. Graphs are pruned; unbranched sets of nodes are compressed into larger nodes
5. Paths are traversed through graph to generate sequences (partically informed by known insert sizes)

Check out this video for a more detailed explanation: http://www.broadinstitute.org/partnerships/education/broade/trinity-screencast

### Memory Issues

This process, while logically straightforward is difficult to implement on the large scales of our transcriptome data.  Effectively, the computer has to keep track of every unique k-mer in the dataset in RAM.  This leads to fantastic amounts of RAM needed (hundreds of GB to terabytes) in order to perform the assembly.  It is relatively easy to gather a dataset that we don't have the computer power to assemble!

Fortunately, many smart people have been working to improve algorithms to make this process more efficient-- one of these is called "digital normalization", which acts to reduce the number of k-mers which need to be tracked by the assembler, reducing the memory footprint and increasing the speed with which the transcriptome can be assembled.    There are two sources of "superflous" k-mers:

1) *Redundant coverage of highly abundant genes* As the RNA-seq datasets represent a population of mRNAs, highly abundnant genes are represented with more reads than lowly abudant genes.  But, these additional reads provide no more additonal <i>information</i> about the sequence, only the amount of that sequence present in the data.

2) *Sequencing Errors*  Next generation sequences have about 1-2% erorr rates (1 bp in every 100bp is incorrect)-- the sequence reported is not the same as the actual sequence.  Because deBruin graphs keep track of every unique k-mer, many of the nodes in the de-bruin graph result from sequencing errors, and are therefore are not informative.  Because sequencing errors are random,

The way to handle this is algorithmically simple-- the reads in a dataset are iterated through, and k-mers in each read are counted.  If the k-mers have been seen enough times previously (the number of times is a threshold set by user), the read is removed from the dataset.  If the k-mers have not been seen before, the read is retained.

This has two practical outcomes:  the most obvious is that "redundant" data from high-coverage genes is removed from the dataset.  The second outcome, is less intiutitve, but important on memory usage.   Because of the 1-2% error rate, we  expect each read has an average of 1 sequencing error, which will introduce k number of erroneous k-mers.   Therefore, for every read that is removed, we remove k erronenous k-mers.  This apparently has minimal impact on the quality of the assembly (and in some cases, actually improves it!).  The interested can check out this paper by [Titus Brown](http://arxiv.org/pdf/1203.4802v2.pdf) to get more details.

### Trinity

Trinity is one of several software packages designed for de-novo RNA sequence assembly.  It is by far the easiest to use, but with that ease of use comes its own limitations-- for instance, you can't pick which size k-mers to use.  The practical consequences of this seem to vary between bioinformaticians-- I like to think of this difference of opinon as the folks that like to use Linux vs. MacOSX.  Sure there are more settings, but what helps you get your job done quicker?

Surely, the simplest method is a great place to start-- and if you get to the point where you need that fine degree of control, you'll have learned a lot!

Some caveats about Trinity:

1.  It is a total memory pig.  You need a lot of memory and a lot of time on the cluster to make this work.  Odds are your first submissions will fail.  Don't get frustrated, the pig can be tamed!
2. Trinity can be restarted!  Yay!  If you get to a point in the process where Trinity suddenly dies, don't delete your work!  You can restart Trinity and it starts more or less where it left off.

General guidelines for Trinity are about 1GB per 1 Million reads (Trinity Website). This is a very loose rule of thumb.  Here are some statistics/parameter combinations on some of my previous assemblies:

    Campylomormyrus (not normalized)
    EO => 43097741 reads
	SM => 34963283 reads

	Trinity Parameters:   --seqType fq --single ../campy_EO.fastq ../campy_SM.fastq --CPU 4 --JM 12G
	Runtime
	=======
	  Start:       Wed Feb 26 23:20:52 EST 2014
	  End:         Thu Feb 27 08:27:09 EST 2014
	  Trinity.pl   32777 seconds (9.5 hours!)
	    Inchworm   3559 seconds
	    Chrysalis  20553 seconds
	    Butterfly  8244 seconds
	    Rest       421 seconds


### How to Run Trinity

Trinity has a lot of options, let's walk through the essential ones:

	module load Trinity
	module load bowtie
	Trinity.pl --seqType fq --left left.fq --right right.fq --seqType fq --CPU 4 --JM 12G --output trinity_dir > trinity.log


--SeqType specifies the type of input file (in our case Fastq).
--JM specifies the amount of memory to use for Jellyfish (This will be proportional to the amount of reads, and has to be determined somewhat empirically.  Good to set high to ensure everything stays running!)  Also note that this should be slightly less than the total memory requested for the job, as Trinity runs programs in the background to monitor the status of things.
--left and --right specify the read1 files and read2 files.  I typically stick the single ended reads in with the left reads.  Note that you can combine multiple files (for instance if you want to do SM and EO together...) by listing multiple files.
--output specifies the name of the directory
--CPU number of CPUs that Butterfly (and other paraallel sub programs) can use.  Pick this carefully, too many and the job won't run!

All the output is printed to the screen, you can direct this to a file so you can monitor the progress during the job run

Optional arguements:
--normalize_reads  will implement normalization on the dataset automagically without user intervention.  We should benchmark assemblies with and without this option to compare the performance.


### Sources and Resources:
[ANGUS Titus Brown - MSU](http://ged.msu.edu/angus/index.html)
[Eel Pond Tutorial](https://khmer-protocols.readthedocs.org/en/v0.8.4/mrnaseq/index.html)
[Trinity RNAseq](http://trinityrnaseq.sourceforge.net)
[Broad Institute Genome Bootcamp](http://www.broadinstitute.org/scientific-community/science/platforms/genome-sequencing/broadillumina-genome-analyzer-boot-camp)
