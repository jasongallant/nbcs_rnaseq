# Read Alignment

## Background Information (5-10 min)

![Alignment Overview](http://www.nature.com/nrg/journal/v12/n10/images/nrg3068-f3.jpg)

1.

Check out this video for a more detailed explanation: http://www.broadinstitute.org/partnerships/education/broade/trinity-screencast


## Performing Alignments with Bowtie2 (~ 25 min)

```bash
## submit_align.sb

#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=04:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=12G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name align_brain_74      # you can give your job a name for easier identification (same as -J)

########## Command Lines to Run ##########

module load bowtie2
module load samtools

cd ${SLURM_SUBMIT_DIR}

bowtie2 -p 16  -x bgaud -1 74_brain_S26_L002_R1_001.fastq.gz.PwU.qtrim.fq -2 74_brain_S26_L002_R2_001.fastq.gz.PwU.qtrim.fq | samtools view -bS - > 74_brain.bam
```

After the alignment is complete, we need to do a few little things to fix it up for the next steps:

```bash
module load samtools
samtools sort 74_brain.bam > 74_brain.sorted.bam
samtools index 74_brain.sorted.bam

```

Next, we run an algorithm called HTSeq, which counts the reads algorithmically that overlap genes in the genome.  This is a very basic method of quantifying expression, but is central to the premise of RNAseq.  Dr. Hoke will be elaborating on this point more, but we'll need this data for our comparision of Illumina and Oxford Nanopore reads in the final section.

```bash
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=04:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=2           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=12G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name htseq_brain74      # you can give your job a name for easier identification (same as -J)

########## Command Lines to Run ##########


cd ${SLURM_SUBMIT_DIR}

source htseq/bin/activate

htseq-count -s no -r pos -t gene -i ID -f bam 74_brain.sorted.bam bgaud_genome.genesonly.gff > brain74.counts

```

Now that we've got the count data, we need to have a look at it.  We're going to fire up R and have a quick look:

```bash
module load R #insert appropriate Version
R

```R
require(data.table)
dat<-fread('./brain74.counts')
hist(log10(dat$V2), 100)
```

The data will look like this: ![histogram](images/hist.jpg)
