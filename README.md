# Bacterial Intergenic Region Analysis Pipeline (BIRAP)

##Description

**B**acterial **I**ntergenic **R**egion **A**nalysis **P**ipeline (**BIRAP**) is an open source, easy to use Perl pipeline that can be used to re-annotate bacterial genomes using experimental data. The tool integrates expression profile derived from RNA-seq and/or proteogenomics, compares it with existing in silico annotation and helps validate annotation, identify novel protein coding regions, putative non-coding RNA as well as help correct er-rors in the existing annotation. The pipeline requires "pileup" output from SAMtools for RNA-seq data and peptide/ePST locus file (GFF) from the Proteogenomic Mapping tool for prote-ogenomics data along with the genome and existing in silico annotation. When information regarding the locus of promot-ers and terminators in the genome (in .coords format) is provided, the pipeline will also help identify non-coding RNA.


##Usage

**Bacterial Intergenic Region Analysis Pipeline Help**

###Basic usage:-
> $ perl birap.pl \-\-analysis [1 or 2 or 3] \-\-fasta [genome.fasta] \-\-pileup [file1 file2 ... or *.pileup] \-\-pgm [proteogenomic\.gff] \-\-anno [annotation.gff OR annotation.ptt] \-\-output [output filename base]  \-\-pro [promoters\.coords] \-\-term [terminators.coords] \-\-intlen [value (default:70)]

*********************************************************
###Description of parameters:

**Mandatory arguments**: \-\-analysis, \-\-fasta, \-\-pileup and/or \-\-pgm based on "analysis", \-\-anno, \-\-output

**\-\-analysis**    (*Choose the type of analysis you wish to perform from among the following options:*)  
  **1** for analysis of intergenic regions using **RNA-seq** data  
  **2** for analysis of intergenic regions using **Proteomics** data  
  **3** for analysis of intergenic regions using **integrated omics** approach    
  
**\-\-fasta** genome file  
**\-\-pileup**        one or more bowtie alignment output files in pileup format (obtained using samtools 'mpileup')  
        (please provide space-separated absolute path for each file with respect to current working directory.   'pathtodir/*.pileup' can also be used)
**\-\-pgm**   GFF file from Proteogenomic mapping tool or similar tools providing the locus of expressed peptides in GFF format  
**\-\-anno**  annotation file in .gff or .ptt format  
**\-\-output**        base filename to be used for all output files  

**Optional arguments:**  

**\-\-pro**   file containing loci of promoters in the genome in .coords format (see below)\*\*  
**\-\-term**  file containing loci of terminators in the genome in .coords format (see below)\*\*  
**\-\-intlen**  minimum length of an intergenic region to be considered as a region of interest for downstream RNA-seq data analysis (default cutoff is **70bp**)  

\*\* *.coords FORMAT:*  
\<START\> \<STOP\> \<STRAND\> \<DESCRIPTION\>    
*Sample entry for promoter locus in .coords file:*    
14024   14052   +       sigmaA\_15bp  

*Sample entry for terminator locus in .coords file:*    
3046632 3046596 -       TERM\_001  

If files containing locus of promoters and terminators in the genome are provided, the program will identify putative non-coding RNA based on the existing annotation and the expression profile generated from RNA-seq  

*********************************************************


##Input

RNA-seq input files: 'mpileup' output from SAM tools generated after aligning the reads to the genome
Protein expression: '.gff' output from Proteogenomic Mapping Tool containing locus of peptides and ePSTs. PGM tool manual [here](http://www.agbase.msstate.edu/tools/pgm/).   

###Generating SAMtools output:
Start with Fastq files, map them to the genome using Bowtie aligner as described in the Bowtie tutorial [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). Steps:  

1. Index the bacterial genome using 'bowtie-build'  
2. Align reads to the genome using Bowtie2 and generate the alignment output in '.sam' format  
3. Convert the .sam to .bam file using SAMtools  
4. Sort the '.bam' file using SAMtools.  
*Do not forget to run the following command before going to step 5*   
"samtools faidx genome.fasta" where "genome.fasta" is your bacterial genome.
5. Generate pileup file for reads per base expression as described   [here](http://samtools.sourceforge.net/samtools.shtml) using the following command  
"samtools mpileup -f genome.fasta sorted1.bam,sorted2.bam,. output.pileup"   
OR  
"samtools mpileup -f genome.fasta sorted1.bam output1.pileup"   
"samtools mpileup -f genome.fasta sorted2.bam output2.pileup" and so on..  
  
*NOTE: Do not forget to save the result in Step 5 as "something.pileup".*   
You can generate a single pileup file containing all your replicates when you provide a comma separated list of all 'sorted.bam' files or generate an individual pileup output file for each of your samples. If you are low are computing power, we recommend generating separate pileup files for each of your samples. Just put them all in one folder before you run BIRAP. Pileup files are usually large and take up a lot of memory.   
Once you have the .gff and/or .pileup files, run BIRAP.pl as detailed above in the BIRAP "Usage" section.  

##Output

All result files are either '.tab' or '.fasta' files. Following is the description of results files depending on the choice of analysis. 

###Analysis choice 1: RNA-seq analysis

1. Tab delimited file containing a list of all annotated regions and whether they were expressed in the data set (***GenEx.tab***).   
\<gene\_start\> \<gene\_stop\> \<gene\_strand\> \<gene\_description\> \<expression\>  
2. Tab delimited file containing list of all expressed intergenic regions (***EIR.tab***) in the following format:   \<EIR\_Description\> \<EIR\_start\> \<EIR\_stop\> \<EIR\_association\>  
The "EIR\_association" column details the relationship of the expressed region with respect to a neighboring annotated region.  
3. Fasta file containing the sequences of all expressed intergenic regions (***EIR.fasta***).  

If loci of promoters and terminators in the genome were provided, two additional result files are generated if any putative non-coding RNA are identified. They are:

4. Tab delimited file containing loci of all putative non-coding RNA (***putative\_sRNA.tab***) in the following format  
\<EIR\_Description\> \<EIR\_start\> \<EIR\_stop\> \<signal\_start\> \<signal\_stop\> \<signal\_description\> \<signal\_start\> \<signal\_stop\> \<signal\_description\>.   
Here "signal" associated with the EIR could be a promoter or a terminator or both.  
5. Fasta file containing sequences of putative non-coding RNA (***putative\_sRNA.fasta***).   

###Analysis choice 2: Proteomic analysis  

1. Tab delimited file containing a list of all annotated regions and the number of peptides that mapped to each (***ProtEx.tab***).  
\<protein\_start\> \<protein\_stop\> \<strand\> \<description\> \<peptide\_count\> \<expression\>  
2. Tab delimited file containing loci of all expressed intergenic peptides (***EIP.tab***) and their associate to annotated regions   
\<EIP\_start\> \<EIP\_stop\> \<EIP\_strand\> \<EIP\_description\> \<EIP\_association\>  
3. Fasta file containing the sequences of ePSTs of all intergenic peptides that code for putative novel protein coding regions (***Intergenic\_ePST.fasta***).  

###Analysis choice 3: Integrated analysis  

Results files mentioned above generated from RNA-seq and Proteogenomic analysis as well as the following:  

1. Tab delimited file containing list of all annotated regions and whether they were expressed in the RNA-seq and proteomics data sets (***IntegratedEx.tab***).  
2. Tab delimited file containing list of all expressed peptides and ePSTs (from PGM tool) and whether they were expressed in the RNA-seq data (***PGMEx.tab***).   

##Possible downstream analysis of results from BIRAP:  

All expressed intergenic regions/peptides/ePSTs (tab delimited) files can be converted into '.gff' format and be loaded into annotation tools such as Artemis or IGV for updating annotation.   
Fasta sequences of expressed intergenic regions can be searched in BlastX to identify novel protein coding regions and/or update existing annotation.   
Fasta sequences of putative sRNA can be searched against Rfam or sRNAdb to annotate them.   
If expressed regions are not conserved in any of the above mentioned databases, it is up to the end user to interpret and annotate their findings.  