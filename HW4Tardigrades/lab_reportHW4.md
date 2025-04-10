### Identification of DNA-associated Proteins in *Ramazzottius varieornatus*

**Abstract:**

Tardigrades (*Ramazzottius varieornatus*) are microscopic extremophiles noted for extraordinary stress tolerance, particularly resistance to DNA-damaging conditions like radiation. We analyzed the R. varieornatus genome to identify proteins potentially involved in DNA repair or protection. By matching peptide sequences identified from DNA-associated fractions to predicted proteins, we pinpointed two candidate proteins: g4106.t1, associated with the endoplasmic reticulum, and g12510.t1, associated with membranes. These findings highlight potential novel mechanisms of DNA protection in tardigrades.

---

**Introduction**    
Tardigrades, also known as water bears, are renowned for their resilience against extreme environmental conditions, including intense radiation, dehydration, and high pressures. They are the first known animals capable of surviving exposure to outer space. Understanding the genomic basis of their stress tolerance, particularly mechanisms for efficient DNA repair, holds promise for advances in biotechnology and radiation protection.

 
Recently, Boothby et al. [1] proposed extensive horizontal gene transfer (HGT) as a key mechanism underlying tardigrades’ remarkable stress resistance. However, Koutsovoulos et al. [2] subsequently refuted this claim, demonstrating minimal actual HGT and highlighting significant genomic contamination issues. Clarifying the exact genetic basis for tardigrade stress resistance thus remains critical, motivating us to investigate the genome of *Ramazzottius varieornatus* through rigorous bioinformatics approaches.

---

**Methods**    
The assembled genome of *R. varieornatus* was downloaded from NCBI (Accession: GCA_001949185.1_Rvar_4.0) [3] and obtained protein predictions from precomputed AUGUSTUS outputs. Predicted proteins were extracted using the provided script, you can find it in Lab Journal (Supplementary material 1). 

For identifying DNA-associated proteins, peptide sequences from mass spectrometry analysis were matched against predicted proteins using MMseqs2 with default parameters, except for a minimum sequence identity threshold of 80% to ensure high-confidence matches.

Localization predictions and functional annotation employed WoLF PSORT and TargetP-2.0 webservers with default parameters for animal-type organisms. Domain identification utilized HMMER (v3.4) against the latest Pfam-A database. The hmmscan analysis was run locally, employing default significance thresholds. Comprehensive BLAST searches was applied with default parameters (using "All non-redundant GenBank CDS translations+PDB+SwissProt+PIR+PRF excluding environmental samples from WGS projects" database).

For more detailed instructions and program references please refer to the Lab Journal (Supplementary material 1).

---

**Results**    
Out of **16,435 proteins** predicted in *R. varieornatus*, MMseqs2 alignment identified two candidate DNA-associated proteins (g4106.t1 and g12510.t1) matched by peptide data with 100% sequence identity and coverage.

Localization and Domain Analysis: 

**Protein g4106.t1** was predicted by WoLF PSORT primarily to localize to the **endoplasmic reticulum (ER)** (score: 14.5). Detailed analysis revealed a cleavable N-terminal signal peptide (positions 1–24), a transmembrane domain, a clear bipartite nuclear localization signal (NLS), and a peroxisomal targeting motif (SKL-like motif: "GRA"). However, TargetP analysis did not indicate strong mitochondrial or secretory targeting signals. Importantly, no significant functional domains were identified in Pfam analysis for this protein, suggesting a novel or highly specialized functional role.    
**Protein g12510.t1** was localized primarily to the **plasma membrane** (WoLF PSORT score: 29), supported by the identification of a significant **MARVEL membrane-associating domain** (MARVEL domain is known as a MAL and related proteins for vesicle trafficking and membrane link) [4] in Pfam analysis (**E-value: 3.3×10⁻⁷**). Moreover, WoLF PSORT highlighted its structure as an integral membrane protein with seven predicted transmembrane segments. This protein also contains a leucine zipper motif, suggesting potential involvement in protein-protein interactions or DNA-associated transcriptional regulation. Like g4106.t1, g12510.t1 lacked strong mitochondrial or secretory signals in TargetP analysis.

BLAST Database Search:   

BLAST searches yielded no significant matches for either protein, underscoring their novelty and possible tardigrade-specific functions. However, targeted BLAST analysis against available tardigrade genomic sequences (NCBI) confirmed strong tardigrade-specific matches:    
-	g4106.t1: Perfect match (E-value: 0) to the hypothetical protein RvY_04937 from *Ramazzottius varieornatus*, indicating species-specific evolution.     
-	g12510.t1: Highly significant matches (E-value ≤ 2e-149) to hypothetical proteins (RvY_14932 isoforms) in *Ramazzottius varieornatus*, strongly suggesting a conserved tardigrade-specific role.

 There is also weak hit in the g4106.t1 with Sarcoplasmic reticulum histidine-rich calcium-binding protein (47% identity, E-value: 0.001).
 
 BLASTP Results for g4106.t1 protein presented on Figure 1
![alt text](<BLASTP Results for g4106.t1.png>) Figure 1

 BLASTP Results for g12510.t1  protein presented on Figure 2
![alt text](<BLASTP Results for g12510.t1.png.png>) Figure 2

Table 1 summarizes the integrated analysis of identified proteins:


| Protein ID | Localization (WoLF PSORT)                     | Localization (TargetP) | Pfam Domain | E-value (Pfam) | BLAST Top Hit                                      | BLAST E-value |
|------------|-----------------------------------------------|------------------------|-------------|----------------|----------------------------------------------------|---------------|
| g4106.t1   | Endoplasmic reticulum<br>Nuclear (NLS)<br>Peroxisomal (SKL-like) | Other                  | None        | N/A            | Hypothetical protein RvY_04937 (<i>R. varieornatus</i>) | 0             |
| g12510.t1  | Plasma membrane<br>(Integral, MARVEL domain)  | Other                  | MARVEL      | 3.3×10⁻⁷       | Hypothetical protein RvY_14932 isoforms (<i>R. varieornatus</i>) | ≤2e-149       |

These two proteins, given their novelty, specificity, and strong evidence of membrane and nuclear localization, represent robust candidates for further experimental validation as tardigrade-specific proteins involved in DNA-protective or stress response mechanisms.

---

**Discussion**    
Our analysis identified two novel candidate proteins (g4106.t1 and g12510.t1) potentially involved in tardigrade-specific stress response pathways, with indirect links to DNA protection mechanisms. Both proteins lack homologs in standard protein databases (UniProtKB/Swiss-Prot), suggesting they represent tardigrade-specific evolutionary innovations, consistent with the minimal horizontal gene transfer observed by Koutsovoulos et al.[2].  

Protein g4106.t1, localized predominantly to the endoplasmic reticulum (ER), possesses a nuclear localization signal (NLS) and peroxisomal targeting motifs. Its ER-nuclear localization and peptide-based DNA association suggest a stress-induced nuclear translocation, potentially affecting chromatin stability or nuclear stress signaling. 

Protein g12510.t1 is membrane-associated with a clear MARVEL domain and leucine zipper motif, indicative of roles in membrane integrity, signaling, or transcriptional regulation rather than direct DNA repair. This aligns with previous findings by Hashimoto et al. [5], where identified proteins directly bound DNA, suggesting multiple protective mechanisms might operate in parallel. 

Future experimental validation is recommended:   
-	For g4106.t1: Confirm nuclear localization and chromatin association under stress conditions via microscopy and ChIP assays.
-	For g12510.t1: Validate membrane localization and explore leucine zipper-mediated interactions and signaling roles using biochemical interaction assays and functional knockout studies. 


---

**References**

1. Boothby, T. C., Tenlen, J. R., Smith, F. W., Wang, J. R., Patanella, K. A., Nishimura, E. O., Tintori, S. C., Li, Q., Jones, C. D., Yandell, M., Messina, D. N., Glasscock, J., & Goldstein, B. (2015). Evidence for extensive horizontal gene transfer from the draft genome of a tardigrade. Proceedings of the National Academy of Sciences, 112(52), 15976–15981. https://doi.org/10.1073/pnas.1510461112

2. Koutsovoulos, G., Kumar, S., Laetsch, D. R., Stevens, L., Daub, J., Conlon, C., Maroon, H., Thomas, F., Aboobaker, A. A., & Blaxter, M. (2016). No evidence for extensive horizontal gene transfer in the genome of the tardigrade Hypsibius dujardini. Proceedings of the National Academy of Sciences, 113(18), 5053–5058. https://doi.org/10.1073/pnas.1600338113

3. Schoch, C. L., Ciufo, S., Domrachev, M., Hotton, C. L., Kannan, S., Khovanskaya, R., Leipe, D., McVeigh, R., O’Neill, K., Robbertse, B., Sharma, S., Soussov, V., Sullivan, J. P., Sun, L., Turner, S., & Karsch-Mizrachi, I. (2020). NCBI Taxonomy: A comprehensive update on curation, resources and tools. Database, 2020, baaa062. https://doi.org/10.1093/database/baaa062

4. Sanchez-Pulido, L., Martin-Belmonte, F., Valencia, A., & Alonso, M. A. (2002). MARVEL: A conserved domain involved in membrane apposition events. Trends in Biochemical Sciences, 27(12), 599–601. https://doi.org/10.1016/S0968-0004(02)02229-6

5. Hashimoto, T., Horikawa, D. D., Saito, Y., Kuwahara, H., Kozuka-Hata, H., Shin-I, T., Minakuchi, Y., Ohishi, K., Motoyama, A., Aizu, T., Enomoto, A., Kondo, K., Tanaka, S., Hara, Y., Koshikawa, S., Sagara, H., Miura, T., Yokobori, S. I., Miyagawa, K., … Kunieda, T. (2016). Extremotolerant tardigrade genome and improved radiotolerance of human cultured cells by tardigrade-unique protein. Nature Communications, 7, 12808. https://doi.org/10.1038/ncomms12808
---

**Supplementary material 1**



Complete step-by-step protocol and scripts available in [Lab Journal](lab_journalHW4.md).



