# Differential RNA expression analysis - (Ba)king (Br)ead

This  project focuses on analyzing differential RNA expression in yeast cells during the transition from aerobic respiration to anaerobic fermentation. Yeasts switch metabolism depending on oxygen availability, producing ATP aerobically or ethanol and CO₂ anaerobically. These metabolic changes have practical applications, notably in biotechnology and baking.

## Step 1. Get RNA-seq data
Our meta-data consists info from yeast before and during fermentation (0 and 30 mins) in 2 replicates:

SRR941816: fermentation 0 minutes replicate 1 (413 Mb)   
SRR941817: fermentation 0 minutes replicate 2 (455 Mb)   
SRR941818: fermentation 30 minutes replicate 1 (79.3 Mb)  
SRR941819: fermentation 30 minutes replicate 2 (282 Mb)

We also get reference genome as *Saccharomyces cerevisiae* (strain S288c) and annotation(GFF)
## Step 2. Alignment and quantification
[HISAT2](https://daehwankimlab.github.io/hisat2/)(v 2.2.1) + featureCounts pipeline
1. Build HISAT2 index    
   `hisat2-build GCF_000146045.2_R64_genomic.fna genome_index`
2. Align reads (single-end mode). Repeat for each FASTQ file separately.   
   `hisat2 -p 4 -x genome_index -U SRR941816.fastq.gz | samtools sort > sample16.bam`
3. Convert GFF to GTF. For this purpose use gffread (v0.12.7)
   `gffread GCF_000146045.2_R64_genomic.gff -T -o genome.gtf`
4. Quantify reads with featureCounts (subread v2.0.6)
   `featureCounts -g gene_id -a genome.gtf -o counts.txt *.bam`

   Here we have an error, because some files doesn't have "gene_id" in line, let's figure out which are they:
   ```
   grep -v gene_id genome.gtf     
   NC_001224.1     RefSeq  exon    58009   60724   .       +       .       transcript_id "gene-Q0158"; gene_name "21S_RRNA";     
   NC_001224.1     RefSeq  exon    61868   62447   .       +       .       transcript_id "gene-Q0158"; gene_name "21S_RRNA";

   ```
   Delete lines without gene_id and create a new GTF file:   
   `grep 'gene_id' genome.gtf > genome_fixed.gtf`

   and then again try with new file:     
   `featureCounts -g gene_id -a genome_fixed.gtf -o counts.txt *.bam`
   here is the summary output:
    ```
   Status	sample16.bam	sample17.bam	sample18.bam	sample19.bam
    Assigned	7305701	8016370	1404578	4983577
    Unassigned_Unmapped	512971	505208	65074	229681
    Unassigned_Read_Type	0	0	0	0
    Unassigned_Singleton	0	0	0	0
    Unassigned_MappingQuality	0	0	0	0
    Unassigned_Chimera	0	0	0	0
    Unassigned_FragmentLength	0	0	0	0
    Unassigned_Duplicate	0	0	0	0
    Unassigned_MultiMapping	1330267	1682084	312462	1202454
    Unassigned_Secondary	0	0	0	0
    Unassigned_NonSplit	0	0	0	0
    Unassigned_NoFeatures	587358	591324	99101	368759
    Unassigned_Overlapping_Length	0	0	0	0
    Unassigned_Ambiguity	37538	37704	4327	15800
    ```
5. Simplify counts output   
   `cat counts.txt | cut -f 1,7-10 > simple_counts.txt`
    head:
    ```
    # Program:featureCounts v2.0.6; Command:"featureCounts" "-g" "gene_id" "-a" "genome_fixed.gtf" "-o" "counts.txt" "sample16.bam" "sample17.bam" "sample18.bam" "sample19.bam" 
    Geneid  sample16.bam    sample17.bam    sample18.bam    sample19.bam
    gene-YAL068C    14      16      2       6
    gene-YAL067W-A  0       0       0       0
    gene-YAL067C    116     69      5       11
    gene-YAL065C    11      3       2       3
    gene-YAL064W-B  1       3       0       0
    gene-YAL064C-A  0       1       1       0
    gene-YAL064W    0       0       0       0
    gene-YAL063C-A  0       0       0       0
    ```
## Step 3. Differential Expression Analysis (DESeq2)
Now we move to Rstudio and use provided script [(deseq2.r)](https://doi.org/10.6084/m9.figshare.14239304.v1)

"This script generates following files:
- [result.txt]() will contain calculated metrics for our genes
- [norm-matrix-deseq2.txt]() will contain normalised counts that we will use in visualisation"
  
utilising very last file we could draw [heatmap](output.pdf) via script from the same source.

The heatmap clearly shows significant differences in gene expression between yeast samples before fermentation (NoFermentation_1 and 2) and during fermentation (Fermentation_1 and 2). The distinct clustering indicates strong differential gene expression associated with fermentation conditions. Many genes are upregulated (shown in red) during fermentation, whereas a substantial set of genes are downregulated (shown in green) compared to the non-fermented conditions​.

We have at least 50 significantly differentially expressed genes shown in your heatmap results

## Step 4. Functional Enrichment (GO terms)
Select top 50 genes with lowest adjusted p-value   
`head -n 50 result.txt | cut -f 1 | cut -d "-" -f 2 > genes.txt`

GO analysis

Use [Yeast GO Slim mapper](http://www.yeastgenome.org/cgi-bin/GO/goSlimMapper.pl) (version 2025-03-16) => Upload genes.txt, choose "Yeast GO-Slim: Process", select all terms, and click "Search". Resultunt file is [here](GO_out.txt)

The GO enrichment analysis highlighted several significantly enriched GO terms, notably involving ribosomal biogenesis, rRNA processing, carbohydrate metabolism, and transcription​.

Hypotheses Based on GO Terms:
Upregulated process example:

    Carbohydrate metabolic process (GO:0005975):

        During fermentation, yeast cells shift their metabolism from oxidative phosphorylation to glycolysis, increasing expression of carbohydrate metabolism-related genes (e.g., YBR105C, YER062C, YKR097W, YOL136C). This is logical as fermentation demands rapid sugar breakdown and ethanol production under anaerobic conditions.

Downregulated process example:

    Ribosomal large subunit biogenesis (GO:0042273):

        Ribosome biogenesis (genes like YCR072C, YDL063C, YHR066W, etc.) typically requires significant energy. Under fermentation conditions, yeast may reduce ribosome synthesis as an energy-saving measure, shifting metabolic priorities from growth and protein synthesis towards quick energy production from sugars.

These hypotheses align biologically with yeast transitioning from growth-oriented aerobic respiration to rapid, energy-efficient anaerobic fermentation, reflecting metabolic adaptation during bread dough fermentation.