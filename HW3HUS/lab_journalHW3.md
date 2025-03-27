## Step 1. Quality  assesment of raw reads
We are going to check quality of our raw reads, inspecting charachteristics and comparing them before and after trimming process.
1. For this purpose we will use common tool [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (FastQC v0.12.1), that was previously installed in its own environment, so now we just need to activate it:
   
`fastqc input/*fastq.gz -o QC`   

**name | sequence number**   
SRR292678sub_S1_L001_R1_001 |  5499346  
SRR292678sub_S1_L001_R2_001 |  5499346  
SRR292770_S1_L001_R1_001 | 5102041  
SRR292770_S1_L001_R2_001 | 5102041  
SRR292862_S2_L001_R1_001 | 5102041  
SRR292862_S2_L001_R2_001 | 5102041  

## Step 2.1 K-mer profile and genome size estimation 

[jellyfish](https://github.com/gmarcais/Jellyfish/blob/master/doc/Readme.md) (v2.3.0), [tutor](https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/) 

`jellyfish count -m 31 -C -s 1000 -o output_k31.jf SRR292678sub_S1_L001_R1_001.fastq SRR292678sub_S1_L001_R2_001.fastq` 

The “jellyfish count” command takes the following options:  
- -m or “mer” specifies the length   
- -C tells it to ignore directionality (it treats each read the same as its reverse complement).   
- -s is an initial estimate for the size of the hash table jellyfish uses, set > genome size   
- -o specifies the name of the output file. choose a name with the k-mer length in it.

Making a histogram file:

`jellyfish histo -o kmer_31.histo output_k31.jf`

Checking what's inside of our file (head command):

```
1 13894353
2 822668
3 154910
4 53078
5 26464
6 15473
7 9911
8 7249
9 6073
10 4854
```
 On the left is a list of the bins (the number of times a k-mer occurs or its ‘depth’), and on the right is the count for the number of k-mers in the data that fit into that category.

 We can plot it in any tool, here is example using R:

![alt text](Rplot01.png)

`Peak:
  Frequency: 125
  Count: 46189`
### Calculate coverage depth (N) 
From FastQC report:

- Total Bases (T) = 494.9 Mbp (for R1) + 494.9 Mbp (for R2) = 989.8 Mbp  
- Average Read Length (L) = 90 bp (both)   
- K-mer size (K) = 31
- Peak K-mer Depth (M) = Obtained from kmer_31.histo = 124 
```
M=124
L=90
K=31
N=$(echo "($M * $L) / ($L - $K + 1)" | bc -l)
echo "Coverage depth (N): $N"
> Coverage depth (N): 186
```
### Calculate genome size
```
T=989800000  # 989.8 Mbp in bases
genome_size=$(echo "$T / $N" | bc -l)
echo "Estimated genome size: $genome_size bp"
> Estimated genome size: 5321505.37634408602150537634 bp = 5.32Mbp
```

## Step 3. Assembling *E. coli* X genome from paired reads 

Checking the version of software
`spades.py --version`
[SPAdes](https://ablab.github.io/spades/) (v4.0.0)

`spades.py -1 input/SRR292678sub_S1_L001_R1_001.fastq -2 input/SRR292678sub_S1_L001_R2_001.fastq -o spades`

**output** :

 scaffolds.fasta - resulting scaffolds (recommended for use as resulting sequences)    
contigs.fasta - resulting contigs    
assembly_graph.fastg - assembly graph    
contigs.paths - contigs paths in the assembly graph    
scaffolds.paths - scaffolds paths in the assembly graph   
before_rr.fasta - contigs before repeat resolution   
corrected/ - files from read error correction   
configs/ - configuration files for read error correction   
corrected.yaml - internal configuration file   
Output files with corrected reads   
params.txt - information about SPAdes parameters in this run   
spades.log - SPAdes log   
dataset.info - internal configuration file  
input_dataset.yaml - internal YAML data set file   
K<##>/ - directory containing intermediate files from the run with K=<##>.    
These files should not be used as assembly results; we will use resulting **contigs/scaffolds** in files mentioned above

## Step 4. Assembly quality assessment

[QUAST](https://quast.sourceforge.net/) (v5.2.0)

`quast.py spades/contigs.fasta spades/scaffolds.fasta -o quast`

output:
```
Assembly	contigs	scaffolds
# contigs (>= 0 bp)	517	499
# contigs (>= 1000 bp)	151	147
# contigs (>= 5000 bp)	83	82
# contigs (>= 10000 bp)	68	66
# contigs (>= 25000 bp)	50	50
# contigs (>= 50000 bp)	32	33
Total length (>= 0 bp)	5315160	5316783
Total length (>= 1000 bp)	5204929	5205502
Total length (>= 5000 bp)	5038908	5045822
Total length (>= 10000 bp)	4935231	4935654
Total length (>= 25000 bp)	4639340	4685368
Total length (>= 50000 bp)	3960781	4051686
# contigs	206	214
Largest contig	300784	300784
Total length	5244705	5253645
GC (%)	50.54	50.51
N50	105346	105346
N90	21026	21421
auN	128088.5	128949.0
L50	15	15
L90	54	52
# N's per 100 kbp	0.00	32.23
```
## Step 2.2
Repeat the k-mer profile plotting step with corrected reads from spades ouput. Before running the command below, don't forget to `qzip -d [FILE]`.

`jellyfish count -m 31 -C -s 1000 -o output_cor_k31.jf spades/corrected/SRR292678sub_S1_L001_R1_00100.0_0.cor.fastq spades/corrected/SRR292678sub_S1_L001_R2_00100.0_0.cor.fastq`

`jellyfish histo -o kmer_cor_31.histo output_cor_k31.jf`

![alt text](Rplot.png)

`Peak:
  Frequency: 125
  Count: 46189`

## Step 5.Three libraries quality asessment
Initially, we got three libraries - pair-end, and 2 mate-pair. This condition could improve result of assemling. We will figure this out in this step: took precomputed result od assembling with 3 libraries, the code for that could be:

```
  spades.py \
  --pe1-1 input/SRR292678sub_S1_L001_R1_001.fastq \
  --pe1-2 input/SRR292678sub_S1_L001_R2_001.fastq \
  --mp1-1 input/SRR292770_S1_L001_R1_001.fastq.gz \
  --mp1-2 input/SRR292770_S1_L001_R2_001.fastq.gz \
  --mp2-1 input/SRR292862_S2_L001_R1_001.fastq.gz \
  --mp2-2 input/SRR292862_S2_L001_R2_001.fastq.gz \
  -o three_libs_spades_out \
  -t 4 \
  -m 2 \
  --isolate
```

`quast.py input/three_libs_spades_out/three_libs_spades_out/contigs.fasta input/three_libs_spades_out/three_libs_spades_out/scaffolds.fasta -o quast_3lib`

output for 3 libraries:

```
Assembly	contigs	scaffolds
# contigs (>= 0 bp)	369	327
# contigs (>= 1000 bp)	79	54
# contigs (>= 5000 bp)	33	16
# contigs (>= 10000 bp)	30	13
# contigs (>= 25000 bp)	26	10
# contigs (>= 50000 bp)	22	10
Total length (>= 0 bp)	5403327	5437160
Total length (>= 1000 bp)	5331230	5365719
Total length (>= 5000 bp)	5202939	5258076
Total length (>= 10000 bp)	5183802	5238939
Total length (>= 25000 bp)	5133691	5200270
Total length (>= 50000 bp)	4975501	5200270
# contigs	105	90
Largest contig	698474	2815616
Total length	5350156	5391554
GC (%)	50.59	50.57
N50	335515	2815616
N90	79998	180369
auN	319603.4	1633387.0
L50	6	1
L90	20	7
# N's per 100 kbp	0.00	627.52
```

Single-library assembly: N50 = 105 kbp, 206 contigs.

Three-library assembly: N50 = 335 kbp (contigs) / 2.8 Mbp (scaffolds), 105 contigs.

Improvement: 3.2× higher N50, 51% fewer contigs.

## Step 6. Annotation
`prokka --outdir prokka input/three_libs_spades_out/three_libs_spades_out/scaffolds.fasta`

prokka statistics:

```
organism: Genus species strain 
contigs: 205
bases: 5424432
CDS: 5057
gene: 5160
rRNA: 23
repeat_region: 1
tRNA: 79
tmRNA: 1

```

## Step 7. Finding the closest relative of E. coli X

First, we need to locate 16S rRNA in the assembled E. coli X genome. [Barrnap](https://github.com/tseemann/barrnap) (v 0.9) predicts the location of ribosomal RNA genes in genomes.

`barrnap PROKKA_03252025.fna -outseq barrnap_output.gff`

The lenght of 16S rRNA region ~ 1538 with coverege ~ 6

***BLAST Escherichia coli 55989, complete sequence
NCBI Reference Sequence: NC_011748.1***  

Downloading refernce sequence of 55989 as fasta file for futher analysis

## Step 8. What is the genetic cause of HUS?

"To understand the genetic cause of HUS, we will perform a genome-wide comparison with the reference genome and will analyze the regions where these strains differ from each other. If we find a region where E. coli X encodes a new virulence factor or a new gene responsible for antibiotic resistance, it may shed light on the genetic cause of HUS.

Mauve, which visualizes an alignment as a series of conserved segments called Locally Collinear Blocks (LCBs), which are similar to synteny blocks."

[Mauve](https://darlinglab.org/mauve/user-guide/installing.html) is bloody ancient tool, so for Windows installation we previosly install [Java](https://www.java.com/en/download/) (v8)

'Then: File” → “Align with progressiveMauve...”. Press “Add sequences” and select the reference genome, then the annotated E. coli X genome (“scaffolds.gbk” or “scaffolds.gbf”, depending on version), and start the alignment.
After alignment we need to find our gene of interest (shiga toxin-related genes). go to the Sequence Navigator in Mauve (the icon of the binoculars in the upper line). Select “Product” in the left window and enter the name of the desired gene (or its function) in the right window.'

![alt text](<Screenshot 2025-03-26 121640.png>)

stxB Shiga toxin subunit B, and longer one - stxA Shiga toxin subunit A.

## Step 9. Tracing the source of toxin genes in E. coli X
Here what we discovered close to our gene of interest: from the left side of it - unknown protein encoded within prophage CP-933V, on the right side of it - tRNA-Arg-TCT, tRNA-Arg-TCG, tRNA-Met-CAT, Phage DNA adenine methylase(EC 2.1.1.72)
Shiga-like toxin genes in E. coli X were almost certainly acquired through horizontal gene transfer (HGT) from a bacteriophage (prophage), specifically a CP-933V-like phage.

## Step 10. Antibiotic resistance detection
To search for genes responsible for antibiotic resistance, we will use [ResFinder](http://genepi.food.dtu.dk/resfinder)

reference:   
**Antimicrobial	Class	WGS-predicted phenotype	Genetic background**    
tetracycline	tetracycline	Resistant	tet(B);;2;;AF326777   
doxycycline	tetracycline	Resistant	tet(B);;2;;AF326777   
minocycline	tetracycline	Resistant	tet(B);;2;;AF326777   

 E. coli X:   
streptomycin	aminoglycoside	Resistant	aph(6)-Id;;1;;M28829   
amoxicillin	beta-lactam	Resistant	blaCTX-M-15;;1;;AY044436  
ampicillin	beta-lactam	Resistant	blaCTX-M-15;;1;;AY044436   
cefepime	beta-lactam	Resistant	blaCTX-M-15;;1;;AY044436   
cefotaxime	beta-lactam	Resistant	blaCTX-M-15;;1;;AY044436   
ceftazidime	beta-lactam	Resistant	blaCTX-M-15;;1;;AY044436   
piperacillin	beta-lactam	Resistant	blaCTX-M-15;;1;;AY044436  
aztreonam	beta-lactam	Resistant	blaCTX-M-15;;1;;AY044436  
ticarcillin	beta-lactam	Resistant	blaCTX-M-15;;1;;AY044436  
ceftriaxone	beta-lactam	Resistant	blaCTX-M-15;;1;;AY044436   
cephalothin	beta-lactam	Resistant	blaTEM-1B;;1;;AY458016   
sulfamethoxazole	folate pathway antagonist	Resistant	sul2;;3;;HQ840942  
trimethoprim	folate pathway antagonist	Resistant	dfrA7;;1;;AB161450  

tetracycline	tetracycline	Resistant	tet(A);;4;;AJ517790  
doxycycline	tetracycline	Resistant	tet(A);;4;;AJ517790  

Going bac to Mauve to find  antibiotics resistance region. Do the search on 'bla' genes and side part of them.


*E. coli* X likely acquired β-lactamase-mediated antibiotic resistance through plasmid. Left side: Fumarate reductase subunits (frdD/frdC). These are chromosomal genes in *E. coli*, but their presence near resistance genes suggests integration of foreign DNA. Right side: Small Multidrug Resistance (SMR) transporter (SugE) + Lipocalin (Bic). SugE is a plasmid-borne efflux pump for quaternary ammonium compounds, often linked to antibiotic resistance cassettes. Bic (Bleomycin-induced protein) is associated with stress response and found in genomic islands.