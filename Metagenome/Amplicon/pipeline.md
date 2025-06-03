
# Atlantic salmon [microbiome](<https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1102691/>)
1 ILLUMINA (Illumina MiSeq) run: 253,513 spots, 115.5M bases
Дизайн: Multiple samples were pooled together in equal proportions based on their molecular weight and DNA concentrations. Pooled samples were purified using calibrated Ampure XP beads, with pooled samples and purified PCR product being used to prepare an Illumina DNA library.
## Скачивание данных с SRA

```
prefetch SRR28755808
fastq-dump SRR28755808 --split-files
gzip *.fastq
```
 ## Работа в R
Полнвй файл работы можно найти [здесь](<pipeline.R>).

На первоначальных этапах работы в R происходит подготовка рабочей директории, устпановка необходимых библиотек. Далее переходим к оценке качества и вот здесь притормозим:

 ![alt text](<Screenshot 2025-06-02 224303.png>)
 ![alt text](<Screenshot 2025-06-02 224326.png>)

Forward reads (_1): Хорошее качество до ~200-220 позиции; качество остается выше Q30 до позиции ~200; длина ридов ~250 bp

Reverse reads (_2): Качество падает раньше, после ~150-180 позиции; сильное падение качества в конце; длина ридов ~250 bp

Идем обрезать, но у нас вылазит ошибка о разности количество ридов в прямом и обратном. Фиксим с помощью BBTools:

```
repair.sh in1=SRR28755808/SRR28755808_1.fastq.gz \
          in2=SRR28755808/SRR28755808_2.fastq.gz \
          out1=SRR28755808/SRR28755808_1_repaired.fastq.gz \
          out2=SRR28755808/SRR28755808_2_repaired.fastq.gz \
          outs=SRR28755808/SRR28755808_singletons.fastq.gz

# Для всех образцов автоматически:
for dir in SRR*; do
    echo "Repairing $dir..."
    repair.sh in1="$dir/${dir}_1.fastq.gz" \
              in2="$dir/${dir}_2.fastq.gz" \
              out1="$dir/${dir}_1_repaired.fastq.gz" \
              out2="$dir/${dir}_2_repaired.fastq.gz" \
              outs="$dir/${dir}_singletons.fastq.gz"
done

```
![alt text](<Screenshot 2025-06-03 212951.png>)
С обновленными файлами продолжаем работать в R.


![alt text](<Screenshot 2025-06-03 180329.png>)
![alt text](<Screenshot 2025-06-03 180348.png>)

## Таксономическая классификация

=== АНАЛИЗ КЛЮЧЕВЫХ ТАКСОНОВ ЛОСОСЯ ===    
❌ Carnobacterium: не найдено   
❌ Aeromonas: не найдено   
❌ Mycoplasma: не найдено   
✅ Clostridiales: найдено 1 ASV   
❌ Lactobacillus: не найдено   
❌ Weissella: не найдено   
❌ Flavobacterium: не найдено   
✅ Psychrobacter: найдено 1 ASV   
✅ Pseudoalteromonas: найдено 1 ASV   
✅ Vibrio: найдено 49 ASV   
✅ Proteobacteria: найдено 335 ASV  
✅ Firmicutes: найдено 27 ASV   
❌ Bacteroidetes: не найдено   

Неклассифицированные последовательности по уровням:
Kingdom  Phylum   Class   Order  Family   Genus Species 
      0       9      20      65     155     362     656


>50% ASV без видовой классификации
