#  NGSTools : Next-Generation sequencing toolkits

Author: Czh3 <zhangchao3@hotmail.com>

  High-throughput sequencing technology is repaidly becoming the standrd method for genomics, transcriptomic and epigenetics. The down stream data analysis is sophisticated because of the unprecedented data throughput. This NGSTools helps you to do this easily by writing a few line of scripts.

###For illumina sequencing platform:
* Hiseq 2000
* Hiseq 2500
* HiseqX-Ten/Five
* Hiseq 3000
* Hiseq 4000
* Miseq

###Based on:
* python::Bio
* fastqc
* fastx
* cutadapt
* SAMtools
* picard
* GATK
* HTseq
* cufflinks
* DEseq2


This package help you to build your own pipeline easily.

Here are several examples of NGS pipelines.

mRNA sequencing
```bash
python2.7 ~/bin/NGSTools/RNA_pipeline.py --sampleList sample.list\
	-d raw  -o pipe_out -c ~/.mouse.cfg\
	-a 1,2
```

```bash
python RNA_pipeline.py -h
usage: RNA_pipeline.py [-h] -s SAMPLELIST [-o OUTDIR] [-d {raw,clean}]
                       [-a ANALYSIS] [-c CONFIG] [--debug DEBUG]

A pipeline of RNA_seq data analysis. <zhangchao3@hotmail.com>

optional arguments:
  -h, --help            show this help message and exit
  -s SAMPLELIST, --sampleList SAMPLELIST
                        sample list for RNA samples information.
                         A file each line contains:
                        sampleName	sampleCondition	fastq1Path	fastq2Path
  -o OUTDIR, --outDir OUTDIR
                        The pipeline output dir
  -d {raw,clean}, --dataType {raw,clean}
                        fastq data type:
                        raw data or clean data.
                         if (clean data): not run cutadapter
  -a ANALYSIS, --analysis ANALYSIS
                        analysis of the pipeline to do.
                        Here is some software to choose to analy
                        [1:QC, quality control
                         2:Mapping, align the reads to reference genome
                         3:Cufflinks, assemble with cufflinkes 
						 4:DESeq2, call DEGs(different expression genes) using DESeq2 package
                         5:DEXSeq, call DEUs(different exon usages) using DEXSeq package
                         6:GATK, call SNP on mRNA using GATK]
  -c CONFIG, --config CONFIG
                        the config file of NGSTools package.
  --debug DEBUG         debug mode
```

BS sequence pipeline
```bash
python ~/bin/NGSTools/BS_seq.py -s sample.list 	\
	-o result 	\
	-c /home/zhangc/.mouse.cfg 	\
	-r -a 1,2,7
```

```bash
python ~/bin/NGSTools/BS_seq.py -h
usage: BS_seq.py [-h] -s SAMPLELIST [-o OUTDIR] [-c CONFIG] [-r] [-a ANALYSIS]
                 [--debug DEBUG]

A pipeline of Bisulfite-Seq(WGBS and RRBS) data analysis. <zhangchao3@hotmail.com>

optional arguments:
  -h, --help            show this help message and exit
  -s SAMPLELIST, --sampleList SAMPLELIST
                        sample list for RNA samples information.
                         A file each line contains:
                        sampleID	sampleName	fastq1Path	fastq2Path
  -o OUTDIR, --outDir OUTDIR
                        The pipeline output dir
  -c CONFIG, --config CONFIG
                        the config file of NGSTools package.
  -r, --rrbs            for RRBS library
  -a ANALYSIS, --analysis ANALYSIS
                        analysis of the pipeline to do.
                        Here is some software to choose to analy
                        [1:QC, quality control
                         2:bismark, align the reads to reference genome
                         3:bs_seeker2, align the reads to reference genome
                         4:picard_rmdup, remove PCR duplicates using picard
                         5:BS_SEEK2_methylation_extractor, call methylation levels using the tools in bismark
                         6:swDMR, call methylation levels using swDMR]
                         7:bismark_methylation_extractor
  --debug DEBUG         debug mode
```

MeDIP pipeline
```bash
python ~/bin/NGSTools/MeDIP.py -s sample.list	\
	-o result	\
	-c ~/.mouse.cfg
```

```bash
python ~/bin/NGSTools/MeDIP.py -h
usage: MeDIP.py [-h] -s SAMPLELIST [-o OUTDIR] [-c CONFIG] [-a ANALYSIS]
                [--debug DEBUG]

A pipeline for MeDIP data analysis. <zhangchao3@hotmail.com>

optional arguments:
  -h, --help            show this help message and exit
  -s SAMPLELIST, --sampleList SAMPLELIST
                        sample list for RNA samples information.
                         A file each line contains:
                        sampleID	sampleName	fastq1Path	fastq2Path
  -o OUTDIR, --outDir OUTDIR
                        The pipeline output dir
  -c CONFIG, --config CONFIG
                        the config file of NGSTools package.
  -a ANALYSIS, --analysis ANALYSIS
                        analysis of the pipeline to do.
                        Here is some software to choose to analy
                        [1:QC, quality control
                         2:BWA, align the reads to reference genome
                         3:Filter Bam, get properly matched reads 4:picard_rmdup, remove PCR duplicates using picard
  --debug DEBUG         debug mode
```
