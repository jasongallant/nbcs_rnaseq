# An Introduction to RNA-sequencing

### A Quick Primer on Sequencing Biochemistry
*****

#### 1.  How did sequencing work before 'next generation' sequencing?
To know what we're doing, we need to know where we've been.  The methods we will employ are the "next generation" of DNA/RNA sequencing technologies (often called Next Generation Sequencing or NGS).  To understand why this is such a big deal, you should understand how sequencing used to occur.  This method, called Sanger sequencing was slow, unreliable, and did not operate "in parallel", limiting the amount of data that could be output (typically 800-1000bp per reaction).  Move through the following tutorial to get a sense of how this works:

[Click here to see a tutorial about Sanger Sequencing](http://www.wiley.com/college/pratt/0471393878/instructor/animations/dna_sequencing/index.html)

#### 2. How does next-generation sequencing work?

#### Illumina Sequencing:

There are many "types" of next-generation sequencing.  One of the most widely utilized approaches is the Illumina sequencing platform, which uses "sequencing by synthesis" chemical reaction.  The benefit over "traditional" DNA sequencing is that many millions of sequencing reactions can occur in parallel.  The typical output of an Illumina Sequencing Machine Run (as of 2018) is 100 billion base pairs (assuming you are reading the now aging HiSeq4000).  This is enough to sequence the whole human genome ~ 33 times!  Watch the video below to see how this works.

[![Illumina Sequencing](http://img.youtube.com/vi/HMyCqWhwB8E/hqdefault.jpg)](https://www.youtube.com/watch?v=HMyCqWhwB8E "Illumina Sequencing")

The trade-off for this vast amount of data is that reads of DNA are typically about 300bp.  Meaning that the RNA or DNA that you put in comes out as a "soup" of small reads that need to be assembled.

#### Nanopore Sequencing:
New developments have brought new technologies-- so called "third generation" technologies, such as Oxford Nanopore, that target and sequence entire molecules.  This technique uses synthetic protein 'nanopores' set in an electrically resistant polymer membrane. An ionic current is passed through the nanopore by setting a voltage across this membrane. If an analyte passes through the pore or near its aperture, this event creates a characteristic disruption in current (as shown in the diagram below). Measurement of that current makes it possible to identify the molecule in question, and provides sequences.
![Nanopore](https://nanoporetech.com/sites/default/files/s3/inline-images/sequencing-animated_0.gif)

These "third generation" techniques theoretically allow for the analysis of entire molecules of DNA or RNA in one go-- a big improvement to say the least.  Compared to Illumina, which is fairly mature, "third generation" approaches are are more cutting edge, which means rapid improvement.  Throughput (as of 2018) is lower than Illumina sequencing (as of 2018, it is 10–30 GB per flowcell or about 3-10 human genomes), though this number varies widely from experimenter to experimenter.  Error rates in sequencing are also quite high: the latest chemistry (as of 2018) suggests about 2-13% errors.  This compares to about 0.25% in Illumina reads (as of 2018).  Watch the video below to see how this all works:
[![Nanopore Sequencing](https://i.vimeocdn.com/video/734700212_640.jpg)](https://vimeo.com/297106166 "Nanopore Sequencing")


### What is RNAseq?
****
RNA-seq is a technique that considers the abundance of RNAs in a tissue sample.  To understand the meaning of what RNA-seq is really doing, you have to know a little bit about the significance of RNA in the biology of cells.  Biologists have a model of how cells utilize their genomic information, which is termed "The Central Dogma", and is summarized nicely in this figure:

![The Central Dogma](http://compbio.pbworks.com/f/central_dogma.jpg)

The instructions to build an organism are encoded in the <i>genome</i> of the organism.  In <i>eukaryotic</i> organisms, this information is organized into chromosomes and stored in the <i>nucleus</i>.  In tissues, the entire genome is not <i>expressed</i>, rather just a subset of these genes.  This subset of genes is referred to as the <i> transcriptome</i>.  RNAseq is a methodological approach under the "umbrella" of what is called "transcriptome sequencing".  The transcriptome is the collection of expressed genes in a particular tissue, or at a particular time point, etc, and should be thought of as a "snapshot" of a dynamic process.  The expressed genes, which are <i>encoded</i> in the DNA are copied  into RNA by  enzymes called <i>RNA polymerases</i>.  The RNA is therefore <i>transcribed</i> from DNA, and is called a transcript.  This is only one type of RNA, the messenger RNA (mRNA).  The RNA transcript (or mRNA) is then routed to the <i>ribosome</i> which reads the instructions encoded by this transcript to build a protein out of amino acids.  Proteins are the functional units of cells that make them unique from other cell types.

#### Why is RNAseq useful?
This process is central to understanding many things about a particular tissue.  We can determine what differentiates two tissues from each other by understanding differences in their transcriptome.  For instance, let's say we were interested in the genes that made muscle cells different than skin cells.  Knowing the "contents" of the transcriptome would tell us a number of things:

+ The relative abundance of genes that are expressed in the two tissues.  In this example, genes related to contraction would be much more abundant in muscle cells.  RNAseq is effectively a "sampling" technique.  Consider the photo below.  If there were only three people in this photo, you would have very little ability to assess how common tall, short and medium were among all human beings.  An easy way to address this is by increasing our sample size.  The more people we include, the closer the "sample" approaches the real population (in statistics this is called the central limit theorem).  The same premise works with RNA-seq.  In this case, because of the fantastic depth to which we are sequencing the transcriptome, we start to observe the same molecules multiple times-- the number of times a particular sequence is observed is directly proportional to its abundance.  Thus, the technique can be used to quantify RNA as well as determine its sequence.

![Height Distribution](http://bioluliaes.files.wordpress.com/2013/01/487668.jpg)

+ The structure of genes expressed in the two tissues.  The protein-encoding portions of genes in the genome, called <i>exons</i>, are typically are surrounded by molecular "packing peanuts" called introns.  Adjacent exons are not always assembled sequentially and can be "skipped", leading to multiple versions of the same gene called <i> alternative splice variants </i>.  Alternative splice variants can perform different tasks with the same genetic material.  Perhaps there are differences in the "splicing" between the two tissues of interest.

+ Finally, many organisms do not have genomes sequenced.  We can utilize the "transcriptome" of organisms to determine the sequence of genes that have never been sequenced.  This is called <i> de-novo</i> assembly of transcriptomes.  This could be very valuable for understanding the evolution of particular genes, or evaluating the function of these genes.  Compared to sequencing the whole genome of an organism, sequencing expressed genes can be done for fractional costs and provide a lot of useful information.

#### Practical limitations, and why we need bioinformaticians!

I'm a biologist! -- *why do I need to sit in front of a computer all day?*

+ Unfortunately, NGS sequencing technologies cannot (yet) read "the whole molecule" in a completely error free way. Thus, for all forms of sequencing, some degree of quality control, data manipulation, and other tasks are necessary for inferring the sequences of samples.  

+ All NGS sequencing technologies generate *a lot* of data.  Manipulating large amounts of data efficiently is central to success in using this technique.
