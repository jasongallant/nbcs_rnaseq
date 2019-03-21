# Transcriptome Assembly

## Background Information (5-10 min)
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

## Memory Issues

This process, while logically straightforward is difficult to implement on the large scales of our transcriptome data.  Effectively, the computer has to keep track of every unique k-mer in the dataset in RAM.  This leads to fantastic amounts of RAM needed (hundreds of GB to terabytes) in order to perform the assembly.  It is relatively easy to gather a dataset that we don't have the computer power to assemble!

Fortunately, many smart people have been working to improve algorithms to make this process more efficient-- one of these is called "digital normalization", which acts to reduce the number of k-mers which need to be tracked by the assembler, reducing the memory footprint and increasing the speed with which the transcriptome can be assembled.    There are two sources of "superflous" k-mers:

1) *Redundant coverage of highly abundant genes* As the RNA-seq datasets represent a population of mRNAs, highly abundant genes are represented with more reads than lowly abundant genes.  But, these additional reads provide no more additional <i>information</i> about the sequence, only the amount of that sequence present in the data.

2) *Sequencing Errors*  Next generation sequences have about 1-2% error rates (1-2 bp in every 100bp is incorrect)-- the sequence reported is not the same as the actual sequence.  Because deBruin graphs keep track of every unique k-mer, many of the nodes in the de-bruin graph result from sequencing errors, and are therefore are not informative.  Because sequencing errors are random, the way to handle this is algorithmically simple-- the reads in a dataset are iterated through, and k-mers in each read are counted.  If the k-mers have been seen enough times previously (the number of times is a threshold set by user), the read is removed from the dataset.  If the k-mers have not been seen before, the read is retained.

This has two practical outcomes:  the most obvious is that "redundant" data from high-coverage genes is removed from the dataset.  The second outcome, is less intuitive, but important on memory usage.   Because of the 1-2% error rate, we  expect each read has an average of 1 sequencing error, which will introduce k number of erroneous k-mers.   Therefore, for every read that is removed, we remove k erroneous k-mers.  This apparently has minimal impact on the quality of the assembly (and in some cases, actually improves it!).  The interested can check out this paper by [Titus Brown](http://arxiv.org/pdf/1203.4802v2.pdf) to get more details.

## Running Trinity (~ 25 min)
Trinity is one of several software packages designed for de-novo RNA sequence assembly.  It is easy to use, but with that ease of use comes its own limitations-- for instance, you can't pick which size k-mers to use.  The practical consequences of are beyond this introductory course-- but see https://oyster-river-protocol.readthedocs.io/en/latest/ for a really nice, thorough exploration of this and many other issues regarding transcriptome assembly.

For now, the simplest method is a great place to start-- and if you get to the point where you need that fine degree of control, you'll have learned a lot!  We won't have time to construct a useful assembly, as they take many hours.  Instead, we'll examine the steps from a pre-calculated assembly based on the data that you performed QC on.

Some caveats about Trinity:

1.  It is a total memory pig.  You need a lot of memory and a lot of time on the cluster to make this work.  Odds are your first submissions will fail.  Don't get frustrated, the pig can be tamed!
2. Trinity can be restarted!  Yay!  If you get to a point in the process where Trinity suddenly dies, don't delete your work!  You can restart Trinity and it starts more or less where it left off.


### How to Run Trinity
Trinity has a lot of options, here's the "submission script":

```bash
#trinity_run.sb
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=08:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=30           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=160G                  # memory required per node - amount of memory (in bytes)
#SBATCH --job-name trinity_run     # you can give your job a name for easier identification (same as -J)

########## Command Lines to Run ##########
cd ${SLURM_SUBMIT_DIR}
module load Trinity/2.8.4
Trinity.pl --seqType fq --max_memory 150G --left 74_brain_S26_L002_R1_001.fastq.gz --right 74_brain_S26_L002_R2_001.fastq.gz --CPU 25 --trimmomatic
```

+ `--SeqType` specifies the type of input file (in our case Fastq).
+ `--max_memory` specifies the amount of RAM to use, selected proportional to the amount of reads, and has to be determined somewhat empirically.  Good to set high to ensure everything stays running!)  Also note that this should be slightly less than the total memory requested for the job, as Trinity runs programs in the background to monitor the status of things.
+ `--left` and `--right` specify the R1 files and R2 files from the sequencer.  Note the `--trimmomatic` flag which automagically applies MacManes' reccomended filters by default!
+ `--output` specifies the name of the directory
+ `--CPU` specfies the number of CPUs that can be utilized to run calculations in parallel can use.  Pick this carefully, too many and the job won't run!

Here are the runtime statistics for our job:

```
Runtime
=======
Number of reads: 23437550
Output data: Trinity.fasta 219 MByte
Trinity   24136 seconds (~7h)
  Inchworm (phase 1 - read clustering)  3161 seconds (< 1h)
  Chrysalis (phase 1 - read clustering) 19142 seconds (5.3 h)
  Rest (phase 2 - parallel assembly)       1833 seconds (30 min)
```

Notice that the submission script for `Trinity` is written in a very strange seeming language that we haven't seen before: this is called [SLURM](https://slurm.schedmd.com), and is reserved for computationally intensive tasks.  Imagine you all ran this job at once, but there was only 50 CPUs up for grabs.  How would this work?  Now imagine hundreds or even thousands of people.  In order to ensure resources are used fairly, a rule system is enforced by the HPCC manager.  This sets up a queue to run jobs.  It's far from "first come first served", but integrates the length of the job you need to run, how often you use the cluster, how many resources you need, and prioritizes your job for the next available slot using algorithms.  During "high use" or when you need resource intensive jobs, your job may wait for days to run.  Most often, it will run fairly soon (minutes to hours after submission).   For a quick "peek" at the queue, type:

```bash
squeue
```

Impressive, eh?

### Examining Trinity output: (TODO)

## Module 2 Review: (TODO)
