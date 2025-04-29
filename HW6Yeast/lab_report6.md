# RNA-seq differential expression analysis in yeast during dough fermentation

## Abstract

This study investigates how gene expression in yeast cells changes during the transition from aerobic respiration to anaerobic fermentation, which occurs naturally during bread making. RNA sequencing was performed on yeast samples before and during fermentation to identify differentially expressed genes. Results revealed substantial changes in expression, with many genes upregulated or downregulated in response to fermentation conditions. Functional analysis highlighted significant shifts in metabolic and biosynthetic pathways, illustrating the genetic mechanisms that enable yeast to adapt to anaerobic conditions during bread dough fermentation.

## Introduction
Yeast (*Saccharomyces cerevisiae*) is a widely utilized eukaryotic model organism, particularly valuable in biotechnology and baking. Yeasts can switch their metabolism between aerobic respiration and anaerobic fermentation depending on oxygen availability. Under anaerobic conditions, yeasts ferment sugars to ethanol and carbon dioxide, the latter essential for dough leavening in baking [1]. Understanding the gene expression changes associated with this metabolic shift can help optimize fermentation processes and inform broader biological insights into metabolic regulation.

Differential gene expression analysis via RNA-seq is essential for identifying molecular mechanisms behind physiological changes. By quantifying RNA levels under different conditions, it provides insights into gene regulatory processes and cellular responses to environmental changes [2].

## Methods
RNA-seq data were obtained from yeast samples at two conditions: before fermentation (0 minutes) and during fermentation (30 minutes), each with two replicates (SRR941816, SRR941817, SRR941818, SRR941819) [3]. The reference genome (*Saccharomyces cerevisiae* strain S288c assembly R64) and its annotation file were acquired from NCBI [4].

Reads were aligned using HISAT2 (v2.2.1) [5] in single-end mode, and alignment was sorted with SAMtools. The genome annotation file (GFF) was converted to GTF format using gffread (v0.12.7) [6]. Gene counts were quantified with featureCounts (subread v2.0.6) [7]. Differential expression analysis was performed using DESeq2 [8] in R, and significantly expressed genes were visualized through a heatmap, utilizing provided scripts [9].

Functional enrichment analysis using Gene Ontology (GO) terms was conducted via Yeast GO Slim Mapper [10] to categorize the biological processes significantly impacted during fermentation.

For more detailed instructions please refer to the Lab Journal (Supplementary material 1).

## Results

Alignment and quantification of sequencing reads provided robust data to identify numerous genes significantly altered in expression between pre-fermentation and fermentation conditions. Differential expression analysis revealed substantial changes across a broad set of genes, clearly visualized in a heatmap demonstrating distinct clustering of fermentation versus pre-fermentation samples (Figure 1). Fermentation samples (Fermentation_1 and Fermentation_2) form a distinct cluster on the left side, predominantly colored red, indicating a general pattern of gene upregulation. Pre-fermentation samples (NoFermentation_1 and NoFermentation_2) cluster together on the right side, predominantly colored green, indicating downregulation of many of the same genes.From these significantly differentially expressed genes, the top 50 with the most significant changes were selected for functional enrichment analysis using Gene Ontology (GO) terms.
![alt text](<Heatmap_expression.png>)

GO analysis (result could be found [here](GO_out.txt)) indicated significant enrichment in processes such as ribosomal biogenesis, rRNA processing, carbohydrate metabolism, and transcription. Notably, carbohydrate metabolism (GO:0005975) showed prominent upregulation, whereas ribosomal large subunit biogenesis (GO:0042273) displayed significant downregulation.

## Discussion
Differential analysis revealed extensive changes in gene expression associated with the yeast transition from aerobic respiration to anaerobic fermentation, with at least 50 significantly differentially expressed genes identified. GO analysis highlighted shifts in several biological processes, notably carbohydrate metabolism and ribosomal biogenesis.

The upregulated carbohydrate metabolism genes (YBR105C, YER062C, YKR097W, YOL136C) likely reflect increased glycolytic activity needed for rapid sugar breakdown during fermentation, producing ATP in the absence of oxygen and generating CO₂ essential for bread dough rising.
Conversely, downregulated ribosomal biogenesis genes (YCR072C, YDL063C, YHR066W, among others) likely indicate a metabolic shift where energy-intensive processes like ribosome synthesis are reduced, reallocating resources toward anaerobic energy production under fermentation conditions.

These results clearly illustrate metabolic adaptation mechanisms in yeast, emphasizing how cells dynamically respond to changes in oxygen availability, optimizing survival and growth during bread dough fermentation.

## References
1. Maicas, S. (2020). The role of yeasts in fermentation processes. Microorganisms, 8(8), 1142. https://doi.org/10.3390/microorganisms8081142
2. Conesa, A., Madrigal, P., Tarazona, S., Gomez-Cabrero, D., Cervera, A., McPherson, A., ... & Mortazavi, A. (2016). A survey of best practices for RNA-seq data analysis. Genome Biology, 17(1), 13. https://doi.org/10.1186/s13059-016-0881-8
3. National Center for Biotechnology Information (NCBI). (n.d.). BioProject PRJNA212389. Retrieved April 29, 2025, from https://www.ncbi.nlm.nih.gov/bioproject/PRJNA212389
4. National Center for Biotechnology Information (NCBI). (n.d.). BioProject PRJNA43747. Retrieved April 29, 2025, from https://www.ncbi.nlm.nih.gov/bioproject/PRJNA43747
5. Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (n.d.). HISAT2. Retrieved April 29, 2025, from https://daehwankimlab.github.io/hisat2/
6. Pertea, G. (n.d.). Gffread. Retrieved April 29, 2025, from https://github.com/gpertea/gffread
7. Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: An efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923–930. https://doi.org/10.1093/bioinformatics/btt656
8. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
9. Raiko, M. (2021). Scripts for RNA-seq project. Figshare. https://doi.org/10.6084/m9.figshare.14239304.v1
10. Saccharomyces Genome Database. (n.d.). Yeast GO-Slim Mapper. Retrieved April 29, 2025, from https://www.yeastgenome.org/goSlimMapper
    
## Supplementary material 1
Complete step-by-step protocol and scripts available in [Lab Journal](lab_journal6.md).
