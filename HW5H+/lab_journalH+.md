# H+

## Step 1. Data download and preparing

We will work with raw 23andMe data. “Fix your teacher”: my own data (yes, I decided to release it for class purposes) - 23andMe and Genotek. Our raw data looks like this for now:
```
# rsid	chromosome	position	genotype
rs12564807	1	734462	AA
rs3131972	1	752721	AG
rs148828841	1	760998	CC
rs12124819	1	776546	AA
```
"Each line corresponds to a single SNP.  For each SNP, we provide its identifier  (an rsid or an internal id), its location on the reference human genome, and the  genotype call oriented with respect to the plus strand on the human reference sequence."

Firstly we need to convert raw data into vcf format, removing all SNPs corresponding to deletions and insertions, using [plink](https://www.cog-genomics.org/plink/) 

`plink --23file SNP_raw_v4_Full_20170514175358.txt --recode vcf --out snps_clean --output-chr MT --snps-only just-acgt`

Check mitochondrial (MT) and Y-chromosome SNPs:

```
grep -w "Y" snps_clean.vcf | head   # Paternal haplogroup (Y-DNA)
> 595401

grep -w "MT" snps_clean.vcf | head  # Maternal haplogroup (mtDNA)
MT      3       i4001200        T       .       .       .       PR      GT      0/0
MT      7       i4001110        A       .       .       .       PR      GT      0/0
MT      9       i4001358        G       .       .       .       PR      GT      0/0
MT      26      i4000553        C       .       .       .       PR      GT      0/0
MT      40      i4001079        T       .       .       .       PR      GT      0/0
MT      41      i4001190        C       .       .       .       PR      GT      0/0
MT      43      i4000964        C       .       .       .       PR      GT      0/0
MT      46      i4001177        T       .       .       .       PR      GT      0/0
MT      49      i4000987        A       .       .       .       PR      GT      0/0
```
The presence of chromosome Y sequences (contig ID=Y) in the VCF file explicitly confirms the biological sex as male.

## Step 2.  Determine Ancestry and Haplogroups
Determining maternal (mtDNA) Haplogroup using [James Lick mtDNA Haplogroup Analysis](https://dna.jameslick.com/mthap/).
![alt text](<mtDNA Haplogroup.png>)
Summary of an output:

- Most likely haplogroup: H(T152C)
- Key markers present in your mtDNA: 152C, 263G, 750G, 1438G, 4769G, 8860G
- Missing marker (untested): 15326
- Haplogroup H is the most common mtDNA haplogroup in Europe. It is prevalent in Western Europe and associated historically with populations migrating into Europe during and after the last Ice Age.

Determining paternal (Y chromosome) Haplogroup: Use Y-chromosome analysis tools [MorleyDNA](https://ytree.morleydna.com/extractFromAutosomal)
![alt text](MorleyDNA.png)
Summary of an output:

- Most Likely Y-Haplogroup: R1a1a (R1a-M17, R1a-M198)
- R1a1a is common in Eastern European and South Asian populations, historically associated with Indo-European migrations.
- The presence of Y-specific SNP markers (M17, M198, and M417) directly confirms the biological sex as male.
  
## Step 3. Phenotype Annotation (Sex, Eye Color, Traits)
we will do eye-color analysis based on the provided [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694299/). The relevant SNPs are clearly stated:

    rs12913832 (HERC2)

    rs16891982 (SLC45A2)

    rs12203592 (IRF4)

    rs6119471 (ASIP)

    rs12896399 (SLC24A4) (the new SNP added in 8-plex system)

 Prediction method:
 1. rs12913832 (HERC2):

    - A/A or G/A → NOT blue (brown or green)
    - G/G → NOT brown (green or blue)
2. Brown eyes: (A/A or G/A at rs12913832) AND (G/G at rs6119471 OR C/C at rs16891982)
    - Green eyes: (G/G at rs12913832 and C/C at rs16891982) OR (G/A at rs12913832 and T/T at rs12203592)
    - Blue eyes: G/G at rs12913832 and T/T at rs12203592
    - Brown eyes if A/A at rs12913832 AND G/G at rs12896399
    - Blue eyes if G/G at rs12913832 AND T/T at rs12896399


`grep -w -E "rs12913832|rs16891982|rs12203592|rs6119471|rs12896399" snps_clean.vcf`

```
5       33951693        rs16891982      C       G       .       .       PR      GT      0/1
6       396321  rs12203592      C       T       .       .       PR      GT      0/1
14      92773663        rs12896399      G       .       .       .       PR      GT      0/0
15      28365618        rs12913832      A       G       .       .       PR      GT      0/1
```
Summary genotype:
- rs12913832 (HERC2): A/G (0/1)
- rs16891982 (SLC45A2): C/G (0/1)
- rs12203592 (IRF4): C/T (0/1)
- rs12896399 (SLC24A4): G/G (0/0)
- rs6119471 (ASIP): Not provided

Lets use simple python [script](<eye_predict.py>) to analyse our SNP presented and applying rule.

> Predicted eye color: Not blue (Brown or Green)
> 
Well, based on data it's complicated to be sure on 100% exact eye color. It's **not blue** bu hard to distinguish between brown or green eyes because it's heterozygous at key SNP positions, and rs6119471 is missing. With help of ChatGPT [(Sora)](https://sora.chatgpt.com/library) we can imagine the eye color:

![alt text](eye_color.png)

## Step 4. Clinical Annotation and SNP Exploration

java -jar ./snpEff.jar GRCh37.75 snps_clean.vcf > snps_snpeff.vcf
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
sudo apt-get install tabix
tabix -p vcf clinvar.vcf.gz
java -jar SnpSift.jar annotate clinvar.vcf.gz snps_clean.vcf > snps_clean_clinvar.vcf


Our annotated file [(snps_clean_clinvar.vcf)]() contains clinical annotations (CLNDN) describing diseases or clinical traits. We will save it in separate file for further analysis, having this structure: Chromosome, Position, SNP ID (rsID), and clinical description (CLNDN).
`grep -v "^#" snps_clean_clinvar.vcf | grep "CLNDN" | cut -f1,2,3,8 > clinically_relevant_snps.txt`

There are a lot of "not provided" entries, so we will filter them and keep only informative one:
`grep "CLNDN" snps_clean_clinvar.vcf | grep -v "not_provided" > clinically_relevant_snps_filtered.txt`
GWAS Catalog helps associate SNPs with traits identified through Genome-Wide Association Studies:
1. Download GWAS Catalog data:   
   `wget https://www.ebi.ac.uk/gwas/api/search/downloads/full -O gwas_catalog.tsv`
2. Annotate your VCF file using SnpSift:   
   `java -jar SnpSift.jar gwasCat -db gwas_catalog.tsv snps_clean.vcf > snps_clean_gwascat.vcf`
3. Extract interesting GWAS associations:    
   `grep -v "^#" snps_clean_gwascat.vcf | grep "GWASCAT_TRAIT" | cut -f1,2,3,8 > gwas_trait_associations.txt`
