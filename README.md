# Identification of Genetic Mutations Causing Antibiotic Resistance in E. coli
## Step 1. Manual inspection
After collecting all neccessary items of not resistant to antibiotics E. coli strain (*_genomic.fna.gz., *_genomic.gff.gz  from [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/)) and resistant one (amp_res_1.fastq.gz, amp_res_2.fastq.gz) from [here](https://doi.org/10.6084/m9.figshare.10006541.v3). Let's pick inside and exemine if there is what we expacted.


`zless amp_res_1.fastq.gz`

```console
@SRR1363257.37 GWZHISEQ01:153:C1W31ACXX:5:1101:14027:2198 length=101
GGTTGCAGATTCGCAGTGTCGCTGTTCCAGCGCATCACATCTTTGATGTTCACGCCGTGGCGTTTAGCAATGCTTGAAAGCGAATCGCCTTTGCCCACACG
+
@?:=:;DBFADH;CAECEE@@E:FFHGAE4?C?DE<BFGEC>?>FHE4BFFIIFHIBABEECA83;>>@>@CCCDC9@@CC08<@?@BB@9:CC#######
```

1. Let's find out before trimming how many reads are in each fastq file:

`wc -l amp_res_1.fastq.g`
```console
160713 amp_res_1.fastq.gz
```
`wc -l amp_res_2.fastq.g`

```console
158169 amp_res_2.fastq.gz
```

2. Let's get detailed statistics of out fastq files using [seqkit tool](https://bioinf.shenwei.me/seqkit/) (seqkit v2.9.0), previously installed. 


`seqkit stats amp_res_1.fastq.gz`

```console
file                format  type  num_seqs     sum_len  min_len  avg_len  max_len
amp_res_1.fastq.gz  FASTQ   DNA    455,876  46,043,476      101      101      101
```
`seqkit stats amp_res_2.fastq.gz`
```console 
file                format  type  num_seqs     sum_len  min_len  avg_len  max_len
amp_res_2.fastq.gz  FASTQ   DNA    455,876  46,043,476      101      101      101
```
## Step 2. Coming closer - working with quality of raw reads
We are going to check quality of our raw reads, inspecting charachteristics and comparing them before and after trimming process.
1. For this purpose we will use common tool [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (FastQC v0.12.1), that was previously installed in its own environment, so now we just need to activate it:


`fastqc raw_data/*fastq.gz -o QC`

![alt text](<Screenshot 2025-03-18 132229.png>)

2. Removing low quality base calls can improve downstream steps of the analysis. There are long list of different tolls, we ve choose [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (v0.39).

`trimmomatic PE -phred33 raw_data/*fastq.gz \
-baseout trimmed \
LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20`

```console
Input Read Pairs: 455876 Both Surviving: **446259** (97.89%) Forward Only Surviving: 9216 (2.02%) Reverse Only Surviving: 273 (0.06%) Dropped: 128 (0.03%)
```

(trimmomatic) ktsel@LPT-pc1:~/HWGenomics/HW2/trimmed$ wc -l trimmed_1P 
1785036 trimmed_1P
1785036/4 = **446 259**

 `fastqc trimmed/*P -o QC_trimmed`

3. **Modify Quality Score to 30**


 `trimmomatic PE -phred33 raw_data/*fastq.gz \
-baseout trimmed30 \
LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 MINLEN:20`

Input Read Pairs: 455876 Both Surviving: **376340** (82.55%) Forward Only Surviving: 33836 (7.42%) Reverse Only Surviving: 25307 (5.55%) Dropped: 20393 (4.47%)

(trimmomatic) ktsel@LPT-pc1:~/HWGenomics/HW2/QC_T_30$ wc -l trimmed30_1P
1505360 trimmed30_1P

1505360/4 = **376 340**

`fastqc QC_T_30/*P -o QC_T_30/`

![alt text](<Screenshot 2025-03-19 152315.png>)
- Initial Raw Reads (amp_res_1_fastqc)

Per Base Sequence Quality: Failed, indicating low-quality scores at some positions. Generally, quality scores decrease towards the end of reads​.

- Trimming with Quality Threshold (LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20) (trimmed_1P_fastqc)

Per Base Sequence Quality: Passed. Your trimming parameters significantly improved the read quality by removing low-quality bases at both ends and filtering out reads with poor sliding window quality​.

- Trimming with Stricter Quality Threshold (LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 MINLEN:20) (trimmed30_1P_fastqc)

Per Base Sequence Quality: Passed as well, indicating a further improvement. Reads show consistently higher quality across all bases compared to the previous trimming stage, reflecting stricter trimming parameters​.


    The original reads required trimming due to quality concerns.
    Both trimming methods improved quality, but the stricter method (LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30) provides the best overall quality.
    We can proceed with the reads trimmed using the stricter settings (trimmed30_1P) for any sensitive downstream analyses to ensure accuracy, but consider that it might significantly reduce the read length or total number of retained reads.

## Step 3. Mapping

[bwa](https://github.com/lh3/bwa) (v 0.7.18-r1243-dirty)

1. Indexing

`bwa index GCF_000005845.2_ASM584v2_genomic.fna.gz`

as a result it gave us several files:

```console
GCF_000005845.2_ASM584v2_genomic.fna.gz.ann  GCF_000005845.2_ASM584v2_genomic.fna.gz.pac
GCF_000005845.2_ASM584v2_genomic.fna.gz.amb  GCF_000005845.2_ASM584v2_genomic.fna.gz.bwt  GCF_000005845.2_ASM584v2_genomic.fna.gz.sa 
```
2. Alignment 

`bwa mem raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz trimmed/trimmed_1P trimmed/trimmed_2P > alignment.sam`

 3. Compress SAM file to BAM (binary) samtools (v 1.21)

`samtools view -S -b alignment.sam > alignment.bam`

here -S for SAM input, -b for BAM output

To  get some basic statistics:

`samtools flagstat alignment.bam`


```console
892776 + 0 in total (QC-passed reads + QC-failed reads)
892518 + 0 primary
0 + 0 secondary
258 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
891649 + 0 mapped (99.87% : N/A)
891391 + 0 primary mapped (99.87% : N/A)
892518 + 0 paired in sequencing
446259 + 0 read1
446259 + 0 read2
888554 + 0 properly paired (99.56% : N/A)
890412 + 0 with itself and mate mapped
979 + 0 singletons (0.11% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
Summary of Key Metrics:

High Mapping Rate: The mapping rate of 99.87%. The high percentage of properly paired reads (99.56%). The low number of singletons (0.11%) indicates that most reads have their mate mapped, which is expected for a well-prepared paired-end library.No Duplicates. No Cross-Chromosome Mappings.

`samtools sort alignment.bam -o alignment_sorted.bam`

`samtools index alignment_sorted.bam`

## Step 4. Variant calling

The goal now is to go through our data, and for each position in the reference genome, see how many reads have a mutation at the same position

`samtools mpileup -f raw_data/GCF_000005845.2_ASM584v2_genomic.fna/GCF_000005845.2_ASM584v2_genomic.fna alignment_sorted.bam >  my.mpileup`

To call actual variants, we will be using a program called [VarScan](https://dkoboldt.github.io/varscan/) (v2.4.6)

`varscan mpileup2snp my.mpileup --min-var-freq 0.5 --variants --output-vcf 1 > VarScan_results.vcf`

```console
Only SNPs will be reported
Warning: No p-value threshold provided, so p-values will not be calculated
Min coverage:   8
Min reads2:     2
Min var freq:   0.5
Min avg qual:   15
P-value thresh: 0.01
Reading input from my.mpileup
4641343 bases in pileup file
9 variant positions (6 SNP, 3 indel)
0 were failed by the strand-filter
6 variant positions reported (6 SNP, 0 indel)
```

---

## Step 4. IGV - Variant effect prediction 

1. st SNP Position: 4,390,754
```
ID: .
Chr: NC_000913.3
Position: 4,390,754
Reference: G*
Alternate: T
Qual: .
Type: SNP
Is Filtered Out: No
Alleles:
Alternate Alleles: T
Variant Attributes
NC: 0
ADP: 15
WT: 0
HET: 0
HOM: 1
```
Ген: rsgA (protein_coding);
Продукт гена: ribosome small subunit-dependent GTPase A;
G → T;
Mutation type: A → A synonymous

![alt text](<Screenshot 2025-03-18 223701.png>)

2. nd SNP Position: 3,535,147
```
ID: .
Chr: NC_000913.3
Position: 3,535,147
Reference: A*
Alternate: C
Qual: .
Type: SNP
Is Filtered Out: No
Alleles:
Alternate Alleles: C
Variant Attributes
NC: 0
ADP: 17
WT: 0
HET: 0
HOM: 1
```
Ген: envZ (protein_coding);
Продукт гена: sensor histidine kinase EnvZ;
A → C;
Mutation type: V → Y missence 
![alt text](<Screenshot 2025-03-18 224542.png>)

3. d SNP Position: 1,905,761
```
ID: .
Chr: NC_000913.3
Position: 1,905,761
Reference: G*
Alternate: A
Qual: .
Type: SNP
Is Filtered Out: No
Alleles:
Alternate Alleles: A
Variant Attributes
NC: 0
ADP: 13
WT: 0
HET: 0
HOM: 1
```
Ген: mntP (protein_coding);
Продукт гена: Mn(2+) exporter;
G → A;
Mutation type: G → R missence


![alt text](<Screenshot 2025-03-18 225156.png>)

4. th SNP Position: 852,762

```
ID: .
Chr: NC_000913.3
Position: 852,762
Reference: A*
Alternate: G
Qual: .
Type: SNP
Is Filtered Out: No
Alleles:
Alternate Alleles: G
Variant Attributes
NC: 0
ADP: 14
WT: 0
HET: 0
HOM: 1   
```
Ген: rybA (ncRNA);
Продукт гена: small RNA RybA;
A → G;
Mutation type: Non-coding protein, different classification

![alt text](<Screenshot 2025-03-18 225839.png>)


5. th SNP Position: 482,698
   
```
ID: .
Chr: NC_000913.3
Position: 482,698
Reference: T*
Alternate: A
Qual: .
Type: SNP
Is Filtered Out: No
Alleles:
Alternate Alleles: A
Variant Attributes
NC: 0
ADP: 16
WT: 0
HET: 0
HOM: 1
```
Ген: acrB(protein coding);
Продукт гена: multidrug efflux pump RND permease AcrB;
T → A;
Mutation type: аминокислота меняется с L (лейцин) на Q (глутамин) - missense 

![alt text](<Screenshot 2025-03-18 230222-1.png>) 

6. th SNP Position: 93,043

```
ID: .
Chr: NC_000913.3
Position: 93,043
Reference: C*
Alternate: G
Qual: .
Type: SNP
Is Filtered Out: No
Alleles:
Alternate Alleles: G
Variant Attributes
NC: 0
ADP: 17
WT: 0
HET: 0
HOM: 1
```
Ген: ftsI (protein coding);
Продукт гена: peptidoglycan DD-transpeptidase FtsI;
C → G;
Mutation type: аланин - synonymous

![alt text](<Screenshot 2025-03-18 230834.png>) 


## Step 5.  Automatic SNP annotation

[SnpEff](https://pcingola.github.io/SnpEff/) (v 5.2) 

- Create snpEff.config File
   
`echo "k12.genome : ecoli_K12" > snpEff.config`
- Create Folder for the Database 
  
`mkdir -p data/k12`
- Prepare the .gbk File
```
gunzip GCF_000005845.2_ASM584v2_genomic.gbff.gz
cp GCF_000005845.2_ASM584v2_genomic.gbff data/k12/genes.gbk
```
- Build the Database
  
  `snpEff build -genbank -v k12`
- Annotation
  
  `snpEff ann k12 VarScan_results.vcf > VarScan_results_annotated.vcf`

- Parsing

Since SnpEff was installed via conda, there is no package to parse it, but it should be if installation was via wget. So installing SnpSift (v 5.2) into the same environment. Be carefull vith calling it, cuz naming sucks...Runnig this command we will get extracted annotation (ANN) in TSV file:

`SnpSift extractFields VarScan_results_annotated.vcf ANN > extracted_ann.tsv`

![alt text](<Screenshot 2025-03-19 135218.png>)

We want to split each annotation into separate lines and add headers, so we can use the Python script (parse_ann.py), to run it:

`python parse_ann.py snpeff/extracted_ann.tsv parsed_ann.tsv`

now we have human-readable tsv:
![alt text](<Screenshot 2025-03-19 135611.png>)
