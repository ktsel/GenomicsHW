# Microbiome in pediatric inflammatory bowel disease 

Скачивание данных из [SRA](<https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR33112824&display=metadata>) базы данных  и распаковка в форвард и реверс ридов. Тул - sratools

```
prefetch SRR33112824
fastq-dump --split-files SRR33112824/
```
less SRR33112824_1.fastq
```
@SRR33112824.1 NB501763:150:HC2M5BGX5:1:11101:20027:1040 length=151
ACTGTNTACAGCAAACTAACGGACAAATNACTCTACTGCTTTANTCCGTTTTAGCTAAAACATGGTTAAAATTGTACNTCTATGTTCAATTTANCCNCAAGCCCTNTAGNCACAATGTCGTAAAGCGTGGAAAGANTAAGGTTGNTCCCTT
+SRR33112824.1 NB501763:150:HC2M5BGX5:1:11101:20027:1040 length=151
AAAAA#EEEEAEEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEE#EEEEEEEEEE6EEEEEEEEEAAA<EEEEAEE/E#AEA/<<EAAEEE<EE#<E#EE/A6E/E#AE/#6EEEEE</EEEEE/E<AE<EEAEEE#/EA/AEEE#EE</A/
@SRR33112824.2 NB501763:150:HC2M5BGX5:1:11101:7297:1040 length=151
GCCATNAAGACAATGCCAGTAAAGATGGNACCTTTGTAATTACNGGTTCTGTCGGCAATAATACCACCTATCAATGCNAGAATGTAAATAGAANAANAAAATGTTNAGTNGATCAAACCTGCCTCTTTTCCATCCNGTCCGAACNTGGCTT
```
# Оценка качества 

![alt text](<Screenshot 2025-06-02 144856.png>)

Прямое и обратное примерно одинаковы по качеству, обрежем используя trimmomatic.

# trimming
```
trimmomatic PE -threads 5 SRR33112824_1.fastq SRR33112824_2.fastq R1_trimmed.fastq R1_solo.fastq R2_trimmed.fastq R2_solo.fastq HEADCROP:14 SLIDI
NGWINDOW:10:20 LEADING:10 TRAILING:10
TrimmomaticPE: Started with arguments:
 -threads 5 SRR33112824_1.fastq SRR33112824_2.fastq R1_trimmed.fastq R1_solo.fastq R2_trimmed.fastq R2_solo.fastq HEADCROP:14 SLIDINGWINDOW:10:20 LEADING:10 TRAILING:10
Quality encoding detected as phred33
Input Read Pairs: 3931119 Both Surviving: 3914906 (99.59%) Forward Only Surviving: 15321 (0.39%) Reverse Only Surviving: 830 (0.02%) Dropped: 62 (0.00%)
TrimmomaticPE: Completed successfully
```
А также проверим качество после тримминга:

![alt text](<Screenshot 2025-06-02 145137.png>)

# Dehosting with BWA + Samtools

Для дехостинга нужно скачать хоста (в нашем случе - человека), проиндексировать, 
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz    

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 
```


`bwa mem -t 10 ref\ index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna R1_trimmed.fastq R2_trimmed.fastq | tee >(samtools flagstat -> flagstat.txt) | samtools view -b -f 4 - | samtools fastq -1 nonhost_R1.fastq -2 nonfost_R2.fastq -`

Команда не сохраняет промежуточные SAM/BAM, сразу дает очищенные FASTQ и сохраняет статистику выравнивания
```
вот результат файла flagstat.txt
7829821 + 0 in total (QC-passed reads + QC-failed reads)
7829812 + 0 primary
0 + 0 secondary
9 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
3295765 + 0 mapped (42.09% : N/A)
3295756 + 0 primary mapped (42.09% : N/A)
7829812 + 0 paired in sequencing
3914906 + 0 read1
3914906 + 0 read2
3290700 + 0 properly paired (42.03% : N/A)
3294212 + 0 with itself and mate mapped
1544 + 0 singletons (0.02% : N/A)
954 + 0 with mate mapped to a different chr
568 + 0 with mate mapped to a different chr (mapQ>=5)
```
# Проанализировали таксономический состав исходных ридов

`kraken2 --confidence 0.5 --paired nonfost_R2.fastq nonhost_R1.fastq --threads 60 --report report1gr_st_r.txt --use-names --db /mnt/disk1/db/k2_db`

Просмотр [аутпута](<report1gr_st_r.txt>) с помощью Pavian Tool - [здесь](<sankey-report1gr_st_r.txt.html>) 

![Screenshot 2025-06-02 154618](https://github.com/user-attachments/assets/fe6007d5-1cba-43b7-b235-be9ab3d837c8)


# Сборка с помощью megahit

`megahit -1 trimm_1P.fastq.gz -2 trimm_2P.fastq.gz -t 60 -o megahit`

# Проверка качества сборки

`quast.py final.contigs.fa -o quast`

![alt text](<Screenshot 2025-06-02 144027.png>)
```
N50 = 1,673 bp - это хороший результат для кишечного метагенома
Largest contig = 33.9 kb - отличный показатель, вероятно представляет почти полные бактериальные гены/опероны
GC = 46.19% - типично для кишечной микробиоты (обычно 40-50%)

```
# Предсказание генов из контигов
`prodigal -i contigs.fasta -o genes.gff -a proteins.faa -d genes.fna -p meta`

Товарищ выдает координаты генов, но без функциональной аннотации, поэтому идем в reCOGnizer

`recognizer -f proteins.faa -o annotation`

![alt text](<Screenshot 2025-06-02 111021.png>)

Вот он нам выдал кучу всего хорошего. 

* reCOGnizer_results.tsv - главный файл с функциональной аннотацией
* reCOGnizer_report.xlsx - сводный отчет в Excel

Функциональные категории:

* COG_quantification.html/tsv - количественный анализ COG функций
* KOG_quantification.html/tsv - количественный анализ KOG функций
* COG_report.tsv / KOG_report.tsv - детальные отчеты

Специализированные анализы:

* NCBI_Curated_report.tsv - курированные NCBI аннотации
* PFam_report.tsv - домены белков Pfam
* PRK_report.tsv - Protein clusters
* SMART_report.tsv - структурные домены
* TIGR_report.tsv - TIGR families

Базово мы можем посмотреть красивую картинку в Krona

![alt text](<Screenshot 2025-06-02 150613.png>)

А далее работать с [файлом](<reCOGnizer_results.xlsx>)
