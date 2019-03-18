### Transcriptome Assembly

What exactly is a <i> transcriptome </i> anyway?  Gene expression is a very dynamic process-- a transcriptome can be thought of as a "snapshot" of this process, or as a collection of all of the genes that are expressed in a particular tissue at a particular time.  An example might be brain during early development, or a skeletal muscle after exercise.  Our sequenced mRNA is a sample of the expressed genes.

Remember that technologies do not presently exist to sequence entire mRNA molecules, so we have to break them into small (presently ~300bp) sequences and sequence them in vast quantities.  This certainly allows us to generate huge amounts of data!  We ultimately hope to quantify the gene expression in or original sample.

This is pretty easy to do if you have a genome for your organism of interest-- you simply align the cleaned reads to the genome and voila,  the number of reads aligning to each gene sequence correspond to its expression level.  More on this later.

But what if you don't know what the genes look like?  Unfortunately, this is the case for many organisms that haven't had a genome sequenced yet.  We are left with the enormous problem of putting everything back together again (creating a "reference" transcriptome), before we can perform our quantification. This requires us to perform what is called a <i>de-novo</i> transcriptome assembly.

A  computationally efficient method of putting the transcriptome back together is also a method that appeals to common sense-- this is called the "DeBruin Graph" approach.  This method has several steps, but we will walk through the basic process:

![Assembly Overview](http://www.nature.com/nrg/journal/v12/n10/images/nrg3068-f3.jpg)

1. Each read is broken down into overlapping strings of length "k"  (FYI, a k-mer is the name for a string of length k.)
2. An ordered list of kmers is generated (hash table)-- each k-mer overlaps by (k-1 bases)
3. Each <b>unique</b> kmer is assigned a node on a "graph"-- overlapping k-mers are adjacent to one-another.
4. Graphs are pruned; unbranched sets of nodes are compressed into larger nodes
5. Paths are traversed through graph to generate sequences (partially informed by known insert sizes)

Check out this video for a more detailed explanation: http://www.broadinstitute.org/partnerships/education/broade/trinity-screencast

### Memory Issues

This process, while logically straightforward is difficult to implement on the large scales of our transcriptome data.  Effectively, the computer has to keep track of every unique k-mer in the dataset in RAM.  This leads to fantastic amounts of RAM needed (hundreds of GB to terabytes) in order to perform the assembly.  It is relatively easy to gather a dataset that we don't have the computer power to assemble!

Fortunately, many smart people have been working to improve algorithms to make this process more efficient-- one of these is called "digital normalization", which acts to reduce the number of k-mers which need to be tracked by the assembler, reducing the memory footprint and increasing the speed with which the transcriptome can be assembled.    There are two sources of "superflous" k-mers:

1) *Redundant coverage of highly abundant genes* As the RNA-seq datasets represent a population of mRNAs, highly abundant genes are represented with more reads than lowly abundant genes.  But, these additional reads provide no more additional <i>information</i> about the sequence, only the amount of that sequence present in the data.

2) *Sequencing Errors*  Next generation sequences have about 1-2% error rates (1-2 bp in every 100bp is incorrect)-- the sequence reported is not the same as the actual sequence.  Because deBruin graphs keep track of every unique k-mer, many of the nodes in the de-bruin graph result from sequencing errors, and are therefore are not informative.  Because sequencing errors are random, the way to handle this is algorithmically simple-- the reads in a dataset are iterated through, and k-mers in each read are counted.  If the k-mers have been seen enough times previously (the number of times is a threshold set by user), the read is removed from the dataset.  If the k-mers have not been seen before, the read is retained.

This has two practical outcomes:  the most obvious is that "redundant" data from high-coverage genes is removed from the dataset.  The second outcome, is less intuitive, but important on memory usage.   Because of the 1-2% error rate, we  expect each read has an average of 1 sequencing error, which will introduce k number of erroneous k-mers.   Therefore, for every read that is removed, we remove k erroneous k-mers.  This apparently has minimal impact on the quality of the assembly (and in some cases, actually improves it!).  The interested can check out this paper by [Titus Brown](http://arxiv.org/pdf/1203.4802v2.pdf) to get more details.

### Trinity

Trinity is one of several software packages designed for de-novo RNA sequence assembly.  It is by far the easiest to use, but with that ease of use comes its own limitations-- for instance, you can't pick which size k-mers to use.  The practical consequences of this seem to vary between bioinformaticians-- I like to think of this difference of opinion as the folks that like to use Linux vs. MacOSX.  Sure there are more settings, but what helps you get your job done quicker?

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
--CPU number of CPUs that Butterfly (and other parallel sub programs) can use.  Pick this carefully, too many and the job won't run!

All the output is printed to the screen, you can direct this to a file so you can monitor the progress during the job run

Optional arguments:
--normalize_reads  will implement normalization on the dataset automagically without user intervention.  We should benchmark assemblies with and without this option to compare the performance.
