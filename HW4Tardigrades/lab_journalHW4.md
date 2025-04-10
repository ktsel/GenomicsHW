# Tardigrade Genome Analysis Project 
This journal perform step-by-step analyzis of the tardigrade genome to identify potential DNA-repair related proteins.
You can find laboratory report [here](lab_reportHW4.md).

## Step 1. Obtain the genome data

```
# Download the assembled genome
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/949/185/GCA_001949185.1_Rvar_4.0/GCA_001949185.1_Rvar_4.0_genomic.fna.gz

# Uncompress the file
gunzip GCA_001949185.1_Rvar_4.0_genomic.fna.gz
```
`grep -c ">" GCA_001949185.1_Rvar_4.0_genomic.fna`  
200
## Step 2: Obtain protein predictions

On this step we will use precomputed data from AUGUSTUS results (protein FASTA and GFF)
```
# Protein FASTA
wget --no-check-certificate 'https://drive.google.com/uc?export=download&id=1hCEywBlqNzTrIpQsZTVuZk1S9qKzqQAq' -O augustus.whole.faa

# GFF file
wget --no-check-certificate 'https://drive.google.com/uc?export=download&id=12ShwrgLkvJIYQV2p1UlXklmxSOOxyxj4' -O augustus.whole.gff
```

"First, we need to extract protein sequences (fasta) from the prediction output (usually GFF/GTF). If you are using AUGUSTUS results, you can do it with the [getAnnoFasta.pl]() script"   
Make the script executable   
`chmod +x getAnnoFasta.pl`  
run   
`./getAnnoFasta.pl augustus.whole.gff`  

Count the number of predicted proteins   
`grep -c ">" augustus.whole.aa`

The amount of obtained proteins: 16435   
You can have a look on the output file [augustus.whole.aa]()

## Step 3. Identify DNA-associated proteins
First, download the peptide data peptides.fa. Physical Localization using [MMseqs2](https://github.com/soedinglab/MMseqs2) (high-performance sequence search):  

Create separate folder mmseqs2, to perform analysis there
```
~/HWGenomics/HW4/input$ mkdir mmseqs2
~/HWGenomics/HW4/input$ cp augustus.whole.aa mmseqs2/
```
1. prepare db with the predicted proteins   
```
mmseqs createdb augustus.whole.aa tardigrade_db      
mmseqs createindex tardigrade_db tmp_index
```
2. Prepare query database with peptide sequences (['peptides.fa']())
   
`mmseqs createdb peptides.fa peptides_query`   

3. Perform alignment search  
   
`mmseqs search peptides_query tardigrade_db mmseqs_results tmp_search --search-type 3 --min-seq-id 0.8`

*Parameters*:    
--search-type 3 for protein search (blastp equivalent).   
--min-seq-id 0.8 Matches required at least 80% sequence identity for inclusion.

4. Convert results to readable format 
   
`mmseqs convertalis peptides_query tardigrade_db mmseqs_results peptides_alignment.tsv --format-output "query,target,evalue,pident,qcov"`

5. Extract identified proteins
   
```
cut -f2 peptides_alignment.tsv | sort | uniq > protein_ids.txt
seqtk subseq augustus.whole.aa protein_ids.txt > identified_proteins.fa
``` 
The peptide sequences obtained from tandem mass spectrometry were aligned against predicted proteins from the Tardigrade genome using MMseqs2. The alignment identified 2 unique proteins matched by peptides (with pident, qcov = 100.000, 1.000). These identified proteins (g4106.t1, g12510.t1) are considered potential DNA-associated candidates for further localization and functional annotation.

## Step 4. Localization prediction
Upload [identified_proteins.fa]() to the servers below:

here is short report from [WoLF PSORT](https://wolfpsort.hgc.jp/) Prediction:

g4106.t1 details E.R.: 14.5, E.R._golg: 9.5, extr: 7, golg: 3.5, lyso: 3, pero: 2, plas: 1, mito: 1
g12510.t1 details plas: 29, cyto: 3

**g4106.t1**   
Primary: E.R. (Endoplasmic Reticulum) with score 14.5.   
Secondary: E.R._golg (9.5), extr (7), etc.    
Conclusion: The protein (g4106.t1) is likely localized to the endoplasmic reticulum.   

**g12510.t1**   
Primary: Plasma membrane (plas: 29)   
Weak signal: Cytoplasm (3)   
Conclusion: Likely a membrane protein (unrelated to DNA repair).  

*Interpretation*
```
Candidate Protein: g4106.t1 (642 aa)

Localization and structural highlights:

    Signal Peptide: Present (cleavable, positions 1–24), suggesting protein secretion or membrane localization.

    Transmembrane domain: Single predicted TMS (positions 62–78).

    Bipartite nuclear localization signal (NLS) identified at residue 477 ("KRQTGPLSAAASRRDAR").

    Peroxisomal targeting signal (SKL): GRA motif at C-terminus.

Functional implications:

    The presence of a cleavable signal peptide and single transmembrane domain suggests initial membrane targeting or secretion. The presence of a nuclear localization signal (NLS) indicates potential nuclear functions.

    While no known DNA-binding domains or motifs were directly found, the NLS strongly suggests a role within the nucleus, potentially related to nuclear processes, such as DNA repair or chromatin remodeling.

    The presence of a peroxisomal targeting signal (GRA) at the C-terminus could suggest a secondary role in reactive oxygen species (ROS) management, indirectly protecting DNA by mitigating oxidative damage.

Hypothesized DNA-related role:

    Possibly a stress-responsive nuclear protein, initially targeted to the membrane or secretory pathway, then cleaved and translocated to the nucleus under specific stress conditions (e.g., dehydration, radiation). It may aid in DNA repair regulation, chromatin remodeling, or ROS protection.

Candidate Protein: g12510.t1 (513 aa)

Localization and structural highlights:

    No signal peptide, indicating it's not secreted and primarily intracellular.

    Multiple (7) transmembrane segments suggest membrane embedding.

    Leucine zipper (PS00029) identified at positions 191 and 350. Leucine zippers typically mediate protein-protein or protein-DNA interactions.

    Ribosomal protein L30 signature (PS00634) identified, which typically associates with ribosomal functions.

Functional implications:

    The leucine zipper motifs strongly indicate possible roles in protein dimerization or DNA binding. These motifs are frequently observed in DNA-binding proteins such as transcription factors and DNA repair enzymes.

    The ribosomal signature (L30-like) could suggest a role linking ribosome biogenesis or translation regulation under stress conditions, potentially affecting the synthesis of proteins required for DNA repair or protection pathways.

    The integral membrane topology with multiple domains suggests embedding in membranes, possibly functioning as a membrane-associated sensor or signaling component in response to DNA damage or stress signals.

Hypothesized DNA-related role:

    Potentially acts as a membrane-associated stress sensor or signaling protein, responding to DNA-damaging conditions (e.g., radiation, dehydration). Its leucine zipper could mediate downstream signaling or transcriptional regulation of DNA repair genes.
```

Prediction from [TargetP-2.0](https://services.healthtech.dtu.dk/service.php?TargetP-2.0)

**g4106.t1**    
Other: 0.7297     
Signal peptide: 0.2669   
Mitochondrial transfer peptide: 0.0034  

**g12510.t1**   
Other: 0.9997     
Signal peptide: 0.0001   
Mitochondrial transfer peptide: 0.0002

## Step 5. BLAST search

**BLASTP search against the “UniProtKB/Swiss-Prot” database and excluding *Ramazzottius varieornatus* gave no hit**
We did simpliest search with adding no additional parameters (using "All non-redundant GenBank CDS translations+PDB+SwissProt+PIR+PRF excluding environmental samples from WGS projects" db), in that case we succeeded we smth close to our organism tho..

 BLASTP Results for g4106.t1  protein presented on Figure 1
![alt text](<BLASTP Results for g4106.t1.png>) Figure 1

 BLASTP Results for g12510.t1   protein presented on Figure 2
![alt text](<BLASTP Results for g12510.t1.png.png>) Figure 2

## Step 6. Pfam prediction
 Pfam domain prediction locally using [HMMER](https://www.ebi.ac.uk/Tools/hmmer/) (v 3.4) 

`conda install -c bioconda hmmer`

Download latest Pfam database (large file ~1GB)    
`wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`
`gzip -d Pfam-A.hmm.gz`

Prepare Pfam database for HMMER   
`hmmpress Pfam-A.hmm`

run HUMMer   
`hmmscan --cpu 4 --domtblout pfam_domains.tsv Pfam-A.hmm identified_proteins.fa`   

*Parameters*:    
--cpu 4: uses 4 cores (adjust according to your CPU).   
--domtblout: detailed domain table format.   

The output file [pfam_domains.tsv]() lists proteins and identified domains. Check brifely with `head`
```
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
MARVEL               PF01284.28   135 g12510.t1            -            513   3.3e-07   30.8   6.0   1   4      0.93   2.2e+04   -4.2   0.1    12    23    80    91    76    93 0.64 Membrane-associating domain
MARVEL               PF01284.28   135 g12510.t1            -            513   3.3e-07   30.8   6.0   2   4   1.4e-11   3.3e-07   30.8   6.0     4   135   147   295   144   295 0.88 Membrane-associating domain
MARVEL               PF01284.28   135 g12510.t1            -            513   3.3e-07   30.8   6.0   3   4      0.17     4e+03   -1.8  14.6    17   101   323   466   321   500 0.60 Membrane-associating domain
MARVEL               PF01284.28   135 g12510.t1            -            513   3.3e-07   30.8   6.0   4   4     0.072   1.7e+03   -0.7   0.1    80    99   487   506   481   510 0.63 Membrane-associating domain
```

The provided Pfam results indicate that the protein g12510.t1 contains the MARVEL (PF01284) domain with significant hits (E-value: 3.3e-07). The [MARVEL domain is known as a MAL and related proteins for vesicle trafficking and membrane link](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR008253/), suggesting that this protein is likely membrane-associated.
Summary:

    Protein ID: g12510.t1

    Significant Pfam domain: MARVEL (PF01284)

    Domain function: Membrane-associating domain

    E-value: 3.3e-07 (highly significant)


## Step 7: Integrated Analysis 

### Protein: g12510.t1

- Subcellular Localization:
    - WoLF PSORT: Plasma membrane (confidence score: 29/100). No N-terminal signal peptide, seven transmembrane domains (integral membrane protein). Leucine zipper motif (PROSITE) identified, suggesting involvement in protein-protein/DNA interactions. Ribosomal protein L30 motif identified, possibly linking translation regulation under stress.
    - TargetP: No strong targeting signal ("Other" with low probabilities).
- BLAST Homology:
    - Best hit: Hypothetical protein RvY_14932-2 from Ramazzottius varieornatus (E-value: 2 × 10<sup>-149</sup>, 99.5% identity).
- Pfam Domains:
    - MARVEL (PF01284): Membrane-associating domain (E-value: 3.3 × 10<sup>-7</sup>).
- Conclusion:
    - Likely a membrane-associated protein, suggesting this protein is embedded in membranes and involved in membrane-associated signaling pathways. Leucine zipper indicates potential DNA-binding or transcriptional regulation capabilities, indirectly influencing DNA repair and stress response.

### Protein: g4106.t1
- Subcellular Localization:
    - WoLF PSORT: Endoplasmic reticulum (score: 14.5/100). Cleavable N-terminal signal peptide indicating initial secretion or membrane targeting. Single predicted transmembrane domain. Bipartite Nuclear Localization Signal (NLS) identified strongly suggests nuclear functions. Peroxisomal targeting (SKL) motif ("GRA"), hinting at peroxisome-related stress response pathways.
    - TargetP: No strong targeting signal ("Other" with low probabilities).
- BLAST Homology:
    - Best hit: Hypothetical protein RvY_04937 from R. varieornatus (E-value: 0.0, 100% identity).
    - Weak hit: Sarcoplasmic reticulum calcium-binding protein (E-value: 0.001, 47% identity).
- Pfam Domains:
    - None detected
- Conclusion:
    - Likely has dual localization, first targeted to membranes, then cleaved, and transported into the nucleus under stress conditions. Nuclear localization suggests potential roles in chromatin remodeling or regulation of DNA repair pathways.The presence of peroxisomal targeting signals hints at a possible indirect DNA protection mechanism through reactive oxygen species (ROS) mitigation.