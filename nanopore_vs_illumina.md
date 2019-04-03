# Comparing Oxford Nanopore and Illumina RNA-seq Data

New developments have brought new technologies-- so called "third generation" technologies, such as Oxford Nanopore, that target and sequence entire molecules.  This technique uses synthetic protein 'nanopores' set in an electrically resistant polymer membrane. An ionic current is passed through the nanopore by setting a voltage across this membrane. If an analyte passes through the pore or near its aperture, this event creates a characteristic disruption in current (as shown in the diagram below). Measurement of that current makes it possible to identify the molecule in question, and provides sequences.
![Nanopore](https://nanoporetech.com/sites/default/files/s3/inline-images/sequencing-animated_0.gif)

These "third generation" techniques theoretically allow for the analysis of entire molecules of DNA or RNA in one go-- a big improvement to say the least.  Compared to Illumina, which is fairly mature, "third generation" approaches are are more cutting edge, which means rapid improvement.  Throughput (as of 2018) is lower than Illumina sequencing (as of 2018, it is 10–30 GB per flowcell or about 3-10 human genomes), though this number varies widely from experimenter to experimenter.  Error rates in sequencing are also quite high: the latest chemistry (as of 2018) suggests about 2-13% errors.  This compares to about 0.25% in Illumina reads (as of 2018).  Watch the video below to see how this all works:
[![Nanopore Sequencing](https://i.vimeocdn.com/video/734700212_640.jpg)](https://vimeo.com/297106166 "Nanopore Sequencing")

### What to Expect When You're Expecting (Nanopore Data)
Theoretically, a MinION run should put out the same type of data as the de-novo Trinity Transcriptome assembly-- but with redundancy, as it will sequence transcripts proportionally to their abundance.  For de-novo applications, this is really exciting.  However, a major concern might be error rates, which for a 2,000bp transcript may be as many as 200bp.  

### QC

We'll be using a program called [MinIONQC](https://github.com/roblanf/minion_qc) to perform QC.  For reasons that aren't entirely clear to your instructor, this step does not run on the HPCC, but runs great on his laptop!  I've gone ahead and pre-computed the results for you so that we can explore them together, simply.  It's fairly straightforward to run, here's the code below (for your reference)

```bash
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load R/3.5.1-X11-20180604
Rscript MinIONQC.R -i example_input_minion/ -o ~/Desktop/example_pore_qc
```
Let's look at the results together.  Click [here](minionqc/README.md)

### Notes about Error Correction:

One way of dealing with this, is an error correction step:  this could be "Self Correction" or "Hybrid Correction".  In the hybrid correction model, a high depth sequencing platform with lower error rates (e.g. Illumina) is used to "correct" the reads and reduce the overall error rate.  This can get pretty expensive-- for your money and time, you may be better just sequencing with Illumina.

Another way of doing this would be "self correction".  The logic behind this is fairly straightforward-- since errors are random, and you have potentially multiple transcrpipts per gene, the "consensus" of the reads could be gathered, and errors can be removed.  This is an issue still under active development.  For the purposes of this exercise, we'll simply ignore this problem, but know it is something that should be considered.

### Alignments
First, we've performed an alignment of the raw Nanopore data to the B. gauderio genome (as we did previously with the Illumina data).  This should look similar to you, except we are using the program `minimap2` instead of `bowtie2` to perform the alignments, as it copes with the higher error rates of Oxford Nanopore data better (we are told!)

This alignment takes a while, so we won't be running it due to time constraints, but let's have a look at the code that generates the alignment.

```bash
~/minimap2/minimap2 -L -ax splice -k14 -uf ../../bgaud_genome_data/bgaud_genome.fa FAK56435.fastq > bgaud_ont.aln.sam
#note the -L flag fixes cigar strings for conversion to bam
module load samtools
samtools view -O BAM -o bgaud_ont.bam bgaud_ont.aln.sam
samtools sort bgaud_ont.bam -o bgaud_ont.sorted.bam
samtools index bgaud_ont.sorted.bam
```

### Using IGV to Examine Alignments
Let's again look at the alignments that we created in a slightly more graphical way than we have been.

1. Download [igv_session_ont.xml](igv_session_ont.xml) to your local computer by right clicking the link and selecting "Download Selected File As..."
2. Open IGV and load `igv_session_ont.xml`

**Questions**
1. How does the alignment compare between the two platforms for a particular gene of interest?


### RSEM Analysis on Two Alignments

### Plot Outputs

```bash
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load R/3.5.1-X11-20180604
R
```
### Counting Reads with HTSeq
Again, we'll run `HTSeq`, to count our nanopore reads per transcript.  In the interests of time, we won't run this script as it will take a while...

```bash
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=04:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=2           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=12G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name htseq_ont_brain74      # you can give your job a name for easier identification (same as -J)

########## Command Lines to Run ##########
cd ${SLURM_SUBMIT_DIR}
source htseq/bin/activate
htseq-count -s no -r pos -t gene -i ID bgaud_ont.aln.sam ../../bgaud_genome_data/bgaud_genome.genesonly.gff > brain74_ont.counts
```

### The Grand Finale: Comparing ONT and Illumina Empircially

Now that we've obtained Read Counts for the two datasets-- how do they match up?  We'll use R on the HPCC to perform this calculation:

```bash
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load R/3.5.1-X11-20180604
R
```

Run the following code:

```R
library(data.table)
dat<-fread('bgaud_counts_and_gene_info.illumina.and.ont.txt')

plot(log10(dat$Counts),log10(dat$ONT_Counts))

reads.lm = lm(log10(Counts) ~ log10(ONT_Counts), data=dat)
summary(reads.lm)$r.squared

```

**Questions**
Discussion:
1. How strong is the correlation between the two datasets?
2. Are there points of agreement between the two datasets?
3. Are there points of disagreement between the two datasets?
4. How could we empirically determine which was giving us "correct" data?
5. Knowing these limitations, how could you improve a future experiment?

**Sources:**
[Vincent Lacroix](http://www.genoscope.cns.fr/externe/rna_workshop/slides/Vincent_Lacroix.pdf)
