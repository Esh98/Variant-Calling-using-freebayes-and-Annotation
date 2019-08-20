# Variant Calling using freebayes and Annotation    
  
This repository is a usable, publicly available tutorial for introduction to basics of variant calling. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bio Informatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format) and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.  


## Variant Calling Workflow

The variant calling workflow begins with quality control and alignment, similar to the other NGS applications. Alignment is followed by alignment clean-up to prepare data for variant calling. Then, variant calling is performed, followed by filtering and annotation of the variant calls.

![var_calling_workflow](/img/variant_calling_workflow.png)

## Set-up

Before we start with variant calling, we need to set-up our directory structure, and make sure the tools are readily available. 

Login to Xanadu and start an interactive session with one processors: 
> NOTE 
> If you are using your own user account in the xanadu cluster plese use:
> -p general --qos=general  
> Otherwise for workshop perposes please use the following:

```
$ srun --pty -p mcbstudent --qos=mcbstudent --mem 8G  bash
```

If the folder structure has not been created for you by the admin, please create the following directory structure for the variant calling project in your home directory:

```bash
Variant-Calling-using-freebayes-and-Annotation/
    ├── 01_raw_data
    ├── 02_reference_data
    ├── 03_fastqc/
    ├── 04_trimmed_reads/
    ├── 05_align/
    ├── 06_variants/
    ├── 07_annotation/
```

With the `-p` option of the `mkdir` command, we create the above structure very quickly:

```bash
$ mkdir -p Variant-Calling-using-freebayes-and-Annotation
$ cd Variant-Calling-using-freebayes-and-Annotation

$ mkdir -p 01_raw_data 02_reference_data 03_fastqc 04_trimmed_reads 05_align 06_variants 07_annotation
```

Now that we have the directory structure created, let's copy over the data to perform our quality control and alignment, including our fastq files and reference data files:

```bash
$ cp /UCHC/PublicShare/VariantWorkshop/data/*fq ./01_raw_data/

$ cp /UCHC/PublicShare/VariantWorkshop/reference/chr20.fa ./02_reference_data/
```

Once you have copied the above data the above two folders will look like:  
```
01_raw_data/
├── NA12878_20_paired_1.fq
└── NA12878_20_paired_2.fq

02_reference_data/
└── chr20.fa
```   


## Dataset

To explore the variant calling workflow, we will be using a subset of a human WGS dataset attained from the [Genome in a Bottle Consortium (GIAB)](http://jimb.stanford.edu/giab). 

![](/img/genome_in_a_bottle.jpeg)

GIAB was initiated in 2011 by the National Institute of Standards and Technology "to develop the technical infrastructure (reference standards, reference methods, and reference data) to enable translation of whole human genome sequencing to clinical practice" [[1](http://jimb.stanford.edu/giab/)].

The human WGS dataset we will be using in class was completed by GIAB and is "essentially the **first complete human genome to have been extensively sequenced and re-sequenced by multiple techniques**, with the results weighted and analyzed to eliminate as much variation and error as possible" [[2](http://www.nist.gov/mml/bbd/dna-022514.cfm)]. To minimize bias from any specific DNA sequencing method, the dataset was sequenced separately by 14 different sequencing experiments and 5 different platforms [[3](http://www.nature.com/nbt/journal/v32/n3/full/nbt.2835.html)]. 

**The dataset acts as a 'truth set' for variation in the human genome to be used as a genotype reference set to compare variant calls against.** Additionally, the DNA is available for validating new sequencing technologies / analysis methods, and ~8300 vials of DNA from a homogenized large batch of the sample cells is [available](https://www-s.nist.gov/srmors/view_detail.cfm?srm=8398) for distribution to other labs [[2](http://www.nist.gov/mml/bbd/ppgenomeinabottle2.cfm)]. Detailed information on the data and methods have been published, and the project information, data and analyses are available on Github (https://github.com/genome-in-a-bottle) [[1](http://jimb.stanford.edu/giab/), [5](http://www.nature.com/articles/sdata201625)].

The source DNA, known as NA12878, was taken from a single person: the daughter in a father-mother-child 'trio' (she is also mother to 11 children of her own) [[4](http://www.nature.com/nmeth/journal/v12/n10/fig_tab/nmeth.3505_SF4.html)]. Father-mother-child 'trios' are often sequenced to utilize genetic links between family members. 

<img src="/img/na12878_tree.jpg" width="700">

While the sample NA12878 was sequenced at a depth of 300x, we will only be using a **subset of the dataset aligning to chromosome 20**. The sequencing files we will be using for NA12878 sample will have a total of **~4 million paired-end reads**. 

## QC and Alignment

In our workflow, we will be going over the quality control steps in order to remove any adapter or vector contamination. We will be using  *FastQC* to ensure there are no obvious problems with our samples and no adapter or vector contamination. 

Choice of alignment tool is often determined by the type of NGS application being conducted. For variant calling we will use [BWA (Burrows-Wheeler Aligner)](http://bio-bwa.sourceforge.net) for alignment. 

BWA is generally slower than Bowtie2 with similar sensitivity and both tools can perform gapped alignment for the identification of indels and can effectively map paired-end reads. However, BWA is a bit more accurate and provides information on which alignments are trustworthy. Small numbers of bad alignments can result in many false variant calls, so accuracy is paramount, and is the basis for choosing BWA.

### Quality evaluation using FASTQC
In order to evalaute the general quality of reads in the file we will be using FASTQC package.  The command can take multiple files as input and outputs a quality report in html format. Once the file is generated you have to transfer the file to your local computer to open it and examine the results carefully.

In order to evalaute the general quality of reads in the file we will be using FASTQC package.  The command can take multiple files as input and outputs a quality report in html format. Once the file is generated you have to transfer the file to your local computer to open it and examine the results carefully.


<pre style="color: silver; background: black;">
module load fastqc 
fastqc -o fastqc -o . ../01_raw_data/NA12878_20_paired_1.fq ../01_raw_data/NA12878_20_paired_2.fq  
</strong>
</pre>

The full slurm scrip is called [fastqc.sh](/03_fastqc/fastqc.sh) which can be found in **03_fastqc/** folder.   

To view the HTML file generated by the FASTQC program, you may have to copy the files to your local computer using the `transfer.cam.uchc.edu` node which is specifically there to transfer files to and out of the Xanadu cluster. To copy the files using a terminal window use the following command, using a new terminal window which is opened in your local computer.   
```bash
scp USER_NAME@transfer.cam.uchc.edu:/--FULL_PATH_TO_YOUR_FILES--/NA12878_20_paired_1_fastqc.html .
```   

As an example:
```bash
scp USER_NAME@transfer.cam.uchc.edu:/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Calling-using-freebayes-and-Annotation/03_fastqc/NA12878_20_paired_1_fastqc.html  
```   
 


### Quality control using sickle  
Now change the working directory to **04_trimmed_reads**:  

Sickle performs quality control on illumina paired-end and single-end short read data using a sliding window. As the window slides along the fastq file, the average score of all the reads contained in the window is calculated. Should the average window score fall beneath a set threshold, <a href="https://github.com/najoshi/sickle/blob/master/README.md">sickle</a> determines the reads responsible and removes them from the run. After visiting the SRA pages for our data, we see that our data are single end reads. Let's find out what sickle can do with these:

<pre style="color: silver; background: black;">-bash-4.2$ module load sickle

-bash-4.2$ sickle

<strong>Usage</strong>: sickle <command> [options]

<strong>Command</strong>:
pe	paired-end sequence trimming
se	single-end sequence trimming

--help, display this help and exit
--version, output version  Information and exit</pre>

We have single-end sequences. 

<pre style="color: silver; background: black;">-bash-4.2$ sickle se

Usage: sickle se [options] -f <fastq sequence file> -t <quality type> -o <trimmed fastq file>

Options:
-f, --fastq-file, Input fastq file (required)
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)
-o, --output-file, Output trimmed fastq file (required)
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.
-x, --no-fiveprime, Don't do five prime trimming.
-n, --trunc-n, Truncate sequences at position of first N.
-g, --gzip-output, Output gzipped files.
--quiet, Don't print out any trimming information
--help, display this help and exit
--version, output version information and exit</pre>




The quality may be any score from 0 to 40. The default of 20 is much too low for a robust analysis. We want to select only reads with a quality of 35 or better. Additionally, the desired length of each read is 75bp. Again, we see that a default of 20 is much too low for analysis confidence. We want to select only reads whose lengths exceed 45bp. 

Once we have performed data trimming we will recheck the quality of data using fastqc.
 

Let's put all of this together for our sickle script using our downloaded fastq files:

```
nano sickle.sh
```  

<pre style="color: silver; background: black;">
module load sickle

module load sickle
sickle pe -t sanger \
        -f ../01_raw_data/NA12878_20_paired_1.fq -r ../01_raw_data/NA12878_20_paired_2.fq \
        -o trimmed_NA12878_20_paired_1.fq -p trimmed_NA12878_20_paired_2.fq \
        -l 45 \
        -q 25 \
        -s singles_NA12878_20_paired_2.fq  

module load fastqc
fastqc -o ../03_fastqc/ trimmed_NA12878_20_paired_1.fq trimmed_NA12878_20_paired_2.fq 

</pre>
<br>
<pre style="color: silver; background: black;">-bash-4.2$ sbatch sickle.sh </pre>  
  

The full slurm script is called [sickle.sh](/04_trimmed_reads/sickle.sh)  is located in **04_trimmed_reads/** folder.   




### BWA modes

Depending on read length, BWA has different modes optimized for different sequence lengths:

- **BWA-backtrack:** designed for Illumina sequence reads up to 100bp (3-step)

- **BWA-SW:** designed for longer sequences ranging from 70bp to 1Mbp, long-read support and split alignment

- **BWA-MEM:** shares similar features to BWA-SW, but BWA-MEM is the latest, and is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

### Aligning reads with BWA-MEM

Change directories into the `02_reference_data` directory:

```bash
$ cd ../02_reference_data
```

#### Creating BWA-MEM index

Similar to the other alignment tools we have used, the first step in the BWA alignment is to create an index for the reference genome. Similar to Bowtie2, BWA indexes the genome with an FM Index based on the Burrows-Wheeler Transform to keep memory requirements low for the alignment process. 

The basic options for indexing the genome using BWA are:

* `-p`: prefix for all index files

```bash
module load bwa/0.7.17
bwa index -p chr20 chr20.fa
```  

At the same time, we will create a fasta file index for the reference genome using samtools:  
```bash
module load samtools/1.7
samtools faidx chr20.fa
```

The full slurm script is called [index.sh](/02_reference_data/index.sh) which can be found in the **02_reference_data** folder.  


#### Aligning reads with BWA-MEM

Now that we have our indexes created, we can get started with read alignment. Change directories to the **05_align** folder:

```bash
$ cd ../05_align/
```

We will perform alignment on our paired-end reads for sample `na12878`. Details on BWA and its functionality can be found in the [user manual](http://bio-bwa.sourceforge.net/bwa.shtml); we encourage you to peruse through to get familiar with all available options.

The basic options for aligning reads to the genome using BWA-MEM are:

* `-t`: number of threads / cores
* `-M`: mark shorter split hits as secondary
	> This is optional for Picard compatibility as MarkDuplicates can directly process BWA's alignment, whether or not the alignment marks secondary hits. However, if we want MergeBamAlignment to reassign proper pair alignments, to generate data comparable to that produced by the Broad Genomics Platform, then we must mark secondary alignments [[GATK discussion forum](https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently#step3)].

Additionally we will specify:

* the path to genome indexes including prefix
* FASTQ files for paired-end reads
* `.err`: save standard error to file
* `-o`: save alignment output to a SAM file

 >**NOTE:** BWA will soft-clip poor quality sequences from the ends of the reads by default, so we do not need to specify a parameter to perform soft clipping.

```bash
module load bwa/0.7.17
bwa mem -M -t 2 \
        ../02_reference_data/chr20 \
        ../04_trimmed_reads/trimmed_NA12878_20_paired_1.fq ../04_trimmed_reads/trimmed_NA12878_20_paired_2.fq \
        -o na12878.sam
```

The complete slurm script is called [bwa_align](/05_align/bwa_align.sh) which can be found in the **05_align/** folder.  


### Alignment clean-up

The last stage of the alignment phase is marking duplicates, and it is usually only required for variant calling. We need to find reads that are likely artifacts from the PCR amplification as they can bias variant calls.

![align_cleanup](/img/workflow_cleanup.png)

If duplicates aren't marked, then the PCR-based errors will be picked up again and again as false positive variant calls. Duplicates are easy to detect: since they have the same mapping information and CIGAR string:  

![dedup1](/img/dedup_begin.png)

Marking duplicates with tools such as *Picard* or *samblaster* will result in the variant caller ignoring these PCR-based errors, and instead seeing:

![dedup1](/img/dedup_end.png)

The variant caller will be more likely to discard the error, instead of calling it as a variant.

We will be using the [Picard](http://broadinstitute.github.io/picard/) suite of tools from the Broad Institute to sort the alignment SAM file and mark duplicates. The documentation for the tools and their usage and options is available in the [user manual](http://broadinstitute.github.io/picard/command-line-overview.html#Tools).

Using the Picard suite on O2 is a little different from tools we have used this far, let's see what information module spider shows us: 

In addition to usual information it gives some information about how to use it. 

```
To use, type
      java -jar $PICARD [options]
```

Java tools usually have a `.jar` executable file and it needs to be run using `java -jar` as well as the full path to the executable file. *You can check what is stored in the `$PICARD` environment variable, before and after you load the module.*

```bash
$ module load picard/2.9.2
```

Let's check what options or specific tools are available to us with *Picard*:

```bash
$ java -jar $PICARD
```

#### Sorting SAM by coordinates

The *Picard* tool, `SortSam`, sorts an input SAM or BAM file by coordinate, queryname, etc. Input and output formats (SAM or BAM) are determined by the file extension.

The description of base options for the `SortSam` tool:

* `INPUT`:	The BAM or SAM file to sort. Required.
* `OUTPUT`:	The sorted BAM or SAM output file. Required.
* `SORT_ORDER`:	Sort order of output file Required. Possible values: {unsorted, queryname, coordinate, duplicate}
* `VALIDATION_STRINGENCY`: Validation stringency for all SAM files read by this program. Possible values: {STRICT, LENIENT, SILENT}
	
> **NOTE:** BWA can produce SAM records that are marked as unmapped but have non-zero MAPQ and/or non-"\*" CIGAR. Typically this is because BWA found an alignment for the read that hangs off the end of the reference sequence. Picard considers such input to be invalid. In general, this error can be suppressed in Picard programs by passing VALIDATION_STRINGENCY=LENIENT or VALIDATION_STRINGENCY=SILENT [[3](https://sourceforge.net/p/picard/wiki/Main_Page/)]. 

working directory will be **05_align/** 

```bash
module load picard/2.9.2
java -Xmx8G -jar $PICARD SortSam \
        INPUT=na12878.sam \
        OUTPUT=na12878_sort.sam \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=SILENT

```
The full slrum script is called [picard_sort.sh](/05_align/picard_sort.sh) and it can be found in the **05_align** directory.  


#### Marking duplicates
The *Picard* tool, `MarkDuplicates`, can locate and tag duplicate reads (both PCR and optical/sequencing-driven) in a BAM or SAM file, where duplicate reads are defined as originating from the same original fragment of DNA. Explanation of the process of determining duplicate reads is provided in the [user manual](http://broadinstitute.github.io/picard/command-line-overview.html#Tools).

The basic options for marking duplicates are:

* `INPUT`:	The sorted BAM or SAM file to sort. Required.
* `OUTPUT`:	The BAM or SAM output file. Required.
* `METRICS_FILE`: File to write duplication metrics to Required.
* `ASSUME_SORTED`: If true, assume that the input file is coordinate sorted even if the header says otherwise. Default value: false. Possible values: {true, false}
* `VALIDATION_STRINGENCY`: Validation stringency for all SAM files read by this program. Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}

```bash
module load picard/2.9.2
java -Xmx8G -jar $PICARD MarkDuplicates \
        INPUT=na12878_sort.sam \
        OUTPUT=na12878_sort_marked.bam \
        METRICS_FILE=metrics.txt \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT 

```

> We use `java -Xmx8G` in the command above to make sure that Java stays within the memory limits we have asked SLURM for. If you are marking duplicates in a large file, it is not unheard of to set up your script or interactive session with over 40G of memory.

The slurm script is called [picard_markduplicates.sh](/05_align/picard_markduplicates.sh) which can also be found in **05_align/** folder.  
  

#### Creating index for BAM file

Now that we have a sorted BAM file that has duplicates marked, let's index it for visualization with IGV. As we have done in previous sessions, we will use *Samtools* to create the index. We will first need to the load the module:

```bash
$ module load samtools/1.9

$ samtools index na12878_sorted_marked.bam na12878_sort_marked.bami
```

The full script is called [samtools_index.sh](/05_align/samtools_index.sh) which can also be found in **05_align/** folder.  



## Variant Calling

We have the aligned and cleaned up the data, and have a BAM file ready for calling variants. 

![](/img/variant_calling_workflow_2.png)
<img src="../img/variant_calling_workflow_2.png" width="450">

Some of the more popular tools for calling variants include [SAMtools mpileup](http://samtools.sourceforge.net/mpileup.shtml), [the GATK suite](https://www.broadinstitute.org/gatk/about/) and [FreeBayes](https://github.com/ekg/freebayes#freebayes-a-haplotype-based-variant-detector) ([Garrison and Marth, 2012](http://arxiv.org/abs/1207.3907)). While it can be useful to work through the [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/) we will be using FreeBayes in this module as it is just as sensitive and precise, but has no license restrictions. After calling variants, we will filter out low quality variants using *[vcftools](https://vcftools.github.io/index.html)*, a toolkit designed to work with Variant Call Format or VCF files.

## Freebayes

*FreeBayes* is a **haplotype-based** variant detector and is a great tool for calling variants from a population. 

> "FreeBayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment."

> "FreeBayes is haplotype-based, in the sense that it calls variants based on the literal sequences of reads aligned to a particular target, not their precise alignment. This model is a straightforward generalization of previous ones (e.g. PolyBayes, samtools, GATK) which detect or report variants based on alignments. This method avoids one of the core problems with alignment-based variant detection--- that identical sequences may have multiple possible alignments:"

<img src="/img/freebayes_2.png" width="600">
---
<img src="/img/freebayes_1.png" width="200">

> "FreeBayes uses short-read alignments (BAM files with Phred+33 encoded quality scores, now standard) for any number of individuals from a population and a reference genome (in FASTA format) to determine the most-likely combination of genotypes for the population at each position in the reference. It reports positions which it finds putatively polymorphic in variant call file (VCF) format. It can also use an input set of variants (VCF) as a source of prior information, and a copy number variant map (BED) to define non-uniform ploidy variation across the samples under analysis."

### Running FreeBayes

Change the working directory to the **06_variants/** folder, if this is not created make sure to create it using the `mkdir` command, and then change the directory to that folder.


```bash
cd ../06_variants/

module load freebayes/1.1.0

freebayes -h
```

```bash
module load freebayes/1.1.0

freebayes -f ../02_reference_data/chr20.fa \
        ../05_align/na12878_sort_marked.bam > na12878.vcf 
```

The full slurm script is called [freebayes_variants.sh](/06_variants/freebayes_variants.sh) which can be found in **06_variants/** folder.  


### Variant Call Format (VCF)

VCF is a text format. It usually has several header lines before the actual data; the header lines start with `##`. There is usually only 1 VCF file generated for all the samples in an experiment. Variants are represented in the rows, and each sample has a column with the status of a given variant:

```
##format=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
#CHROM  POS     ID        REF   ALT    QUAL  FILTER  INFO                                 FORMAT       NA00001         NA00002         
20      14370   rs6054257 G     A      29    0       NS=55;DP=255;AF=0.768;DB;H2          GT:GQ:DP:HQ  0|0:48:1:51,51  1|0:48:8:51,51  
20      13330   .         T     A      3     q10     NS=55;DP=202;AF=0.024                GT:GQ:DP:HQ  0|0:49:3:58,50  0|1:3:5:65,3    
20      1110696 rs6040355 A     G,T    67    0       NS=55;DP=276;AF=0.421,0.579;AA=T;DB  GT:GQ:DP:HQ  1|2:21:6:23,27  2|1:2:0:18,2    
20      10237   .         T     .      47    0       NS=57;DP=257;AA=T                    GT:GQ:DP:HQ  0|0:54:7:56,60  0|0:48:4:51,51  
20      123456  microsat1 G     D4,IGA 50    0       NS=55;DP=250;AA=G                    GT:GQ:DP     0/1:35:4        0/2:17:2        
```

Often the header lines will have some explanation about the various columns in the VCF, including the confusing looking INFO column. Here's an explanation of the INFO column for the first entry in the example above (the example below is representing the same variant as above, "rs6054257", but the VCF was excerpted from a much larger experiment):

<img src="/img/vcf_3.png" width="600">

Below is another example with slightly different fields in the INFO column:

<img src="/img/vcf_2.png" width="600">

Now, let's take a look at the one we just generated:

```bash
$ less na12878.vcf
```

How does this compare to the 2 examples we have seen so far? How does the ID column compare?

## Filtering VCFs

It is very important to filter out low-quality variants before moving to the assessment and annotation steps. Low quality variants usually represent sequencing errors (low quality bases). Freebayes variant quality determination is done as described here: [https://github.com/ekg/freebayes#observation-qualities](https://github.com/ekg/freebayes#observation-qualities).

Today we are going to use `vcftools` to remove entries that have calls with a quality score of lower than 20.

```bash
$ module load gcc/6.2.0 vcftools/0.1.15
```

The manual for `vcftools` is [available here](https://vcftools.github.io/man_latest.html), let's take a quick look at it.

So the most basic options you need to specify are input `--vcf <name>` and output `--out <name-filtered>`. There are many different criteria that can be used for filtering the input vcf, below are a few examples.

> Include/exclude specific sites by chromosome:

	--chr 20 
	--not-chr 20
	
> No two sites are within specified distance to one another:

	--thin 5
	
> Specify minimum depth for each site:

	--minDP 10
	
> Filter by variant type:

	--keep-only-indels 
	--remove-indels 
	
> Include SNPs with specific ID (i.e. dbSNP, this is information we will be adding in the annotation section):

	--snps <string>

We are going to stick with using only the quality score for today's class:
	
```bash
$ vcftools --vcf na12878.vcf --minQ 20 --recode --recode-INFO-all --out na12878_q20  
``` 
Full slurm script is called [vcftools_filter.sh](/06_variants/vcftools_filter.sh) which can be found in **06_variants/** folder.  


> "`--recode` : These options are used to generate a new file in either VCF or BCF from the input VCF or BCF file after applying the filtering options specified by the user. The output file has the suffix ".recode.vcf" or ".recode.bcf". By default, the INFO fields are removed from the output file, as the INFO values may be invalidated by the recoding (e.g. the total depth may need to be recalculated if individuals are removed). This behavior may be overriden by the following options. By default, BCF files are written out as BGZF compressed files."
> 
> "`--recode-INFO-all` : These options can be used with the above recode options to define an INFO key name to keep in the output file. This option can be used multiple times to keep more of the INFO fields. The second option is used to keep all INFO values in the original file."
> 
> Information about `recode` adapted from the [VCFtools manual](https://vcftools.github.io/man_latest.html).

Now we are *(almost)* ready to annotate this VCF with known information from dbSNP, and add functional annotation information to enable variant prioritization.

## Annotating variants

Variant annotation is a crucial step in linking sequence variants with changes in phenotype. Annotation results can have a strong influence on the ultimate conclusions of disease studies. Incorrect or incomplete annotations can cause researchers both to overlook potentially disease-relevant DNA variants and to dilute interesting variants in a pool of false positives. 

<img src="/img/variant_calling_workflow_3.png" width="450">

At this stage, we have a large tab-delimited file containing loci at which a variation was found in the sample DNA sequence relative to the reference. We have filtered out these variations (also referred to as 'variant calls') to keep only those we are highly confident in, and now need to find out more. We can do this by **comparing our variants against known variants, and also use genome annotations to help predict information about our variants.** 

<img src="/img/prioritize.png" width="200">


### Setting up

For this section we are going to need to copy over some reference data required for annotation. Start an interactive session and move into **Variant-Calling-using-freebayes-and-Annotation/** directory. Then copy over the required data.

```
$ srun --pty -p mcbstudent --qos=mcbstudent --mem=20G bash 

$ cp /UCHC/PublicShare/VariantWorkshop/reference/reference_chr20.vcf.gz* ./02_reference_data/

```

Let's also create a new directory for the results of our annotation steps called **07_annotation/** if not created for you:

```
$ mkdir 07_annotation
```

So now the folder structure will look like:  
```
Variant-Calling-using-freebayes-and-Annotation/
├── 01_raw_data
├── 02_reference_data
├── 03_fastqc
├── 04_trimmed_reads
├── 05_align
├── 06_variants
└── 07_annotation
```   



## Annotation with known variants 

Variant annotation is the process of assigning information to DNA variants. There are many different types of information that can be associated with variants, and a first commonly used resource is using databases which contain variants that have previously been described. One popular example is [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/),a free, public archive for genetic variation within and across different species. It is hosted by NCBI in collaboration with NHGRI and although the name implies SNPs; it actually includes range of molecular variation.

<img src="/img/dbsnp.png" width="500">

To add dbSNP information you need to download the organism specific data using their FTP download. **We have already done this for you** and was the zip file that you copied over into your `reference_data` folder. 

> For the full dbSNP dataset that corresponds with our data you can download it via the FTP site: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/ . You may notice is that there are alot of files, and the README does not provide enough information. To find out more on how the datasets differ you can access the [NCBI human variation docs](http://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/).  

To annotate our data with dbSNP information we wil be using [`bcftools`](https://samtools.github.io/bcftools/), a command-line utility for variant calling and manipulating VCF files and its binary counterpart BCF files. It is a part of the `samtools` project, a tool that we are by now pretty familiar with. 

The `bcftools annotate` command allows the user to **add or remove annotations**. 

```bash
$ module load bcftools/1.6

$ bcftools annotate --help
```

The annotation we wish to add and the file we are annotating must be a Bgzip-compressed and tabix-indexed file (usually VCF or BED format). Tabix indexes a tab-delimited genome position file and creates an index file (.tbi), which facilitates quick retrieval of data lines overlapping regions. *NOTE: this has already been done for our dbSNP file*

```bash
$ bgzip ../06_variants/na12878_q20.recode.vcf 
$ tabix ../06_variants/na12878_q20.recode.vcf.gz
```

> Both `bgzip` and `tabix` are not available on Xanadu as standalone modules, but its provided through `samtools` installation.

When running `bcftools annotate`, we also need to specify the column(s) to carry over from the annotation file, which in our case is ID.

```bash
bcftools annotate -c ID \
        -a ../02_reference_data/reference_chr20.vcf.gz ../06_variants/na12878_q20.recode.vcf.gz \
        > na12878_q20_annot.vcf
```

Take a quick peek at the new VCF file that was generated using `less`. You should now see in the ID column `rs` ids which correspond to identifiers from dbSNP. For the variants that are not known, you will find the `.` in place of an ID indicating novelty of that sequence change.


In the inital vcf file of `na12878_q20.recode.vcf` which is before annotation:  
```
#CHROM  POS     ID      REF     ALT     QUAL
chr20   61795   .       G       T       95.0616
```  

The vcf file `na12878_q20_annot.vcf`, which is after annotation:  
```
#CHROM  POS     ID      	REF     ALT     QUAL
chr20   61795   rs4814683       G       T       95.0616
```  
 
   
## Functional annotation with SnpEff

One fundamental level of variant annotation involves categorising each variant based on its relationship to coding sequences in the genome and how it may change the coding sequence and affect the gene product. To do this we will be using a tool called [SnpEff](http://snpeff.sourceforge.net/), a **variant effect predictor program**. 

Our understanding of the protein-coding sequences in the genome is summarised in the set of transcripts we believe to exist. Thus, **variant annotation depends on the set of transcripts used as the basis for annotation**. The widely used annotation databases and browsers – ENSEMBL, RefSeq, UCSC – contain sets of transcripts that can be used for variant annotation, as well as a wealth of information of many other kinds as well, such as ENCODE data about the function of non-coding regions of the genome. SnpEff will take information from the provided annotation database and populate our VCF file by adding it into the `INFO` field name `ANN`. Data fields are encoded separated by pipe sign "\|"; and the order of fields is written in the VCF header.

<img src="/img/snpeff.png" width="700">

Some **common annotations** are listed below, but [the manual](http://snpeff.sourceforge.net/SnpEff_manual.html#input)
provides a more comprehensive list.

* Putative_impact/impact: A simple estimation of putative impact / deleteriousness : (HIGH, MODERATE, LOW, MODIFIER)
* Gene Name: Common gene name (HGNC). Optional: use closest gene when the variant is “intergenic”.
* Feature type: Which type of feature (e.g. transcript, motif, miRNA, etc.). It is preferred to use Sequence Ontology (SO) terms, but ‘custom’ (user defined) are allowed. 
* Feature ID: Depending on the annotation sources, this may be: Transcript ID (preferably using version number), Motif ID, miRNA, ChipSeq peak, Histone mark, etc. Note: Some features may not have ID (e.g. histone marks from custom Chip-Seq experiments may not have a unique ID).
* Biotype: The bare minimum is at least a description on whether the transcript is (“Coding”, “Noncoding”). Whenever possible, use ENSEMBL biotypes.

Take a look at the options available. We will be using `snpEff` to annotate our variants.

```bash

$ module load snpEff/4.3q

$ java -jar $SNPEFF -h
```

> `snpEff` is also a java based tool, similar to `picard` and we have to use a similar syntax to run it.

To run the snpEff command we will need to specify two things:

1. The appropriate genome
2. The VCF file we want to annotate
	
An additional parameter to add to our command is `Xmx8G`, a Java parameter to define available memory. Since we are in an interactive session with 8GB, if we had requested more before starting the session we could increase the number here.

The final command will look like this:

```bash
module load snpEff/4.3q

java -Xmx8G -jar $SNPEFF eff \
        -dataDir /isg/shared/databases/SnpEff/v4_3/data/ \
        hg19 na12878_q20_annot.vcf > na12878_q20_annot_snpEff.vcf
```	

By default snpEff downloads and install databases automatically (since version 4.0) for the organism that is specified. To see what databases are available for human you can use the `databases` command:

```bash
$ java -jar $SNPEFF/snpEff.jar databases | grep Homo_sapiens
```

In Xanadu cluster we do have databases created and it can be used by users with out having to download on there own. So if you have any requried databases that you would like to use please send us a request and we will download it for you.  

The complete slurm script is called [snpEff.sh](/07_annotation/snpEff.sh) which can be found in **07_annotation/** directory.  
  


### SnpEff Output

SnpEff produces three output files:

1. an annotated VCF file 
2. an HTML file containing summary statistics about the variants and their annotations
3. a text file summarizing the number of variant types per gene

Let's take a look at the **text file**:

```bash
$ head snpEff_genes.txt
```
Each row corresponds to a gene, and each column coresponds to a different variant type. This gives you a resource for quickly interrogating genes of interest and see what types of variants they contain, if any.

To look at the **HTML file**, we will need to use `scp` or FileZilla to bring it over to our local machine. 

A sample report is provided "snpEff_summary.html" if for some reason you fail to produce it.

Let's scroll through the report. The first part of the report is a summary, which outlines what was run and what was found.

<img src="/img/snpeff_summary_f.png">

As we scroll through the report, we can obtain more details on the categories of variants in our file. 

There is a section **summarizing variant by type**:

<img src="/img/snpeff_bytype_f.png">

These different types are defined as follows:

|Type  |  What it means  |  Example|
| ------------- |:-------------:| -----:|
|SNP  |  Single-Nucleotide Polymorphism  |  Reference = 'A', Sample = 'C'|
|Ins  |  Insertion  |  Reference = 'A', Sample = 'AGT'|
|Del  |  Deletion  |  Reference = 'AC', Sample = 'C'|
|MNP  |  Multiple-nucleotide polymorphism  |  Reference = 'ATA', Sample = 'GTC'|
|MIXED  |  Multiple-nucleotide and an InDel  |  Reference = 'ATA', Sample = 'GTCAGT'|


Additionally, variants are **categorized by their 'impact'**: {High, Moderate, Low, Modifier}. These impact levels are [pre-defined categories](http://snpeff.sourceforge.net/SnpEff_manual.html#input) based on the 'Effect' of the variant, to help users find more significant variants. 

<img src="/img/snpeff_byimpact_f.png">

***

**Exercise**

Use the HTML report to answer the following questions:

1. The majority of variants idenified are classified as SNPs. How many insertions and deletions were found?
2. How many of our variants are novel (not in dbSNP)?
3. How many variants were found in exonic regions?
4. The Ts/Tv ratio (the transition/transversion rate) tends to be around 2-2.1 for the human genome, although it changes between different genomic regions. What is the ratio reported for our sample? 



