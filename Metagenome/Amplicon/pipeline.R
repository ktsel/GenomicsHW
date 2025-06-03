# Загрузка библиотек

# Установка пакетов 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("dada2", "phyloseq", "DECIPHER", "phangorn"))
install.packages(c("ggplot2", "vegan", "tidyverse", "plotly", "heatmaply"))

library(dada2)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(vegan)
library(tidyverse)

setwd('/Users/79509/Desktop/metagenome/amplicon')

'''
# Скачать весь проект одной командой
prefetch PRJNA1102691

# Конвертировать все SRA файлы в FASTQ
fasterq-dump SRR28755808 --split-files
gzip *.fastq
'''


dir.create("raw_data")
dir.create("filtered")
dir.create("results")
dir.create("plots")


# Путь к файлам (после скачивания)
path <- "raw_data"
list.files(path)

fnFs <- list.files("raw_data", pattern = "_1_repaired.fastq.gz", 
                   full.names = TRUE, recursive = TRUE)
fnRs <- list.files("raw_data", pattern = "_2_repaired.fastq.gz", 
                   full.names = TRUE, recursive = TRUE)


# Остальной код без изменений
fnFs <- sort(fnFs)
fnRs <- sort(fnRs)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Проверка качества
plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])

'''
Forward reads (_1):

Хорошее качество до ~200-220 позиции
Качество остается выше Q30 до позиции ~200
Длина ридов ~250 bp

Reverse reads (_2):

Качество падает раньше, после ~150-180 позиции
Сильное падение качества в конце (типично для Illumina)
Длина ридов ~250 bp
'''

# Фильтрация и обрезка
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Параметры фильтрации (настроить под ваши данные)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(200, 150),  # Forward: 200 bp, Reverse: 150 bp
                     maxN = 0,                # Никаких неопределенных нуклеотидов
                     maxEE = c(2, 3),         # Forward: max 2 ошибки, Reverse: max 3 ошибки
                     truncQ = 2,              # Обрезать хвосты с качеством < Q2
                     rm.phix = TRUE,          # Удалить контрольные последовательности PhiX
                           # Не обрезать начало (праймеры уже удалены)
                     compress = TRUE)

# Проверить результаты фильтрации
head(out)
'''
reads.in reads.out
SRR28755808_1_repaired.fastq.gz   253513    241266
SRR28755809_1_repaired.fastq.gz   171964    163025
SRR28755810_1_repaired.fastq.gz   213095    204609
'''
summary(out)

# Статистика сохранения ридов
retained_percent <- round(out[,2]/out[,1] * 100, 1)
cat("Процент сохраненных ридов:\n")
print(retained_percent)
cat("Средний процент:", mean(retained_percent), "%\n")
cat("Образцов с > 70% сохранением:", sum(retained_percent > 70), "\n")

'''
Процент сохраненных ридов:
SRR28755808_1_repaired.fastq.gz 
                           95.2 
SRR28755809_1_repaired.fastq.gz 
                           94.8 
SRR28755810_1_repaired.fastq.gz 
                           96.0
                           
Средний процент: 95.33333 %
Образцов с > 70% сохранением: 3
'''

# Проверить качество после фильтрации (первые 3 образца)
plotQualityProfile(filtFs[1:3])
plotQualityProfile(filtRs[1:3])

##### oбучение модели ошибок и деноизинг#####
# Обучение модели ошибок
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

'''
121780000 total bases in 608900 reads from 3 samples will be used for learning the error rates
91335000 total bases in 608900 reads from 3 samples will be used for learning the error rates.
'''

# Визуализация ошибок
plotErrors(errF, nominalQ = TRUE)

# Деноизинг (создание ASVs)
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

'''
Forward
Sample 1 - 241266 reads in 43105 unique sequences.
Sample 2 - 163025 reads in 16978 unique sequences.
Sample 3 - 204609 reads in 36099 unique sequences

Reverse
Sample 1 - 241266 reads in 35410 unique sequences.
Sample 2 - 163025 reads in 13726 unique sequences.
Sample 3 - 204609 reads in 30191 unique sequences
'''
# Слияние paired-end ридов
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

'''
232327 paired-reads (in 474 unique pairings) successfully merged out of 239786 (in 1111 pairings) input.
160389 paired-reads (in 228 unique pairings) successfully merged out of 162049 (in 653 pairings) input.
199542 paired-reads (in 468 unique pairings) successfully merged out of 203403 (in 998 pairings) input.

'''

# Создание таблицы ASV (Amplicon Sequence Variants)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Распределение длин последовательностей
table(nchar(getSequences(seqtab)))

'''
Всего ASV: 889 уникальных последовательностей
Доминирующая длина: 253 bp (616 ASV = 69% от всех)
Основной диапазон: 252-274 bp (большинство последовательностей)

Узкое распределение - большинство последовательностей сконцентрировано в диапазоне 252-274 bp
Ожидаемая длина - после обрезки (200 bp forward + 150 bp reverse) и перекрытия получается ~250-270 bp
Четкий пик на 253 bp указывает на качественное объединение ридов

Короткие (200-250 bp): возможно, плохо объединенные риды
Длинные (>280 bp): могут быть артефактами или химерами

'''

# Анализ распределения длин последовательностей
length_dist <- table(nchar(getSequences(seqtab)))
print("Распределение длин ASV:")
print(length_dist)

# Визуализация распределения длин
library(ggplot2)
length_df <- data.frame(
  Length = as.numeric(names(length_dist)),
  Count = as.numeric(length_dist)
)

ggplot(length_df, aes(x = Length, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Распределение длин ASV последовательностей",
       x = "Длина последовательности (bp)",
       y = "Количество ASV") +
  geom_vline(xintercept = c(250, 275), color = "red", linetype = "dashed",
             alpha = 0.7) +
  annotate("text", x = 262, y = max(length_df$Count) * 0.8, 
           label = "Основной диапазон\n250-275 bp", color = "red")

ggsave("plots/asv_length_distribution.png", width = 10, height = 6)



#### Удаление химер ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                    multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Доля сохраненных ридов = 0.9913433

# Функция для отслеживания ридов на каждом этапе
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- sample.names
print(track)

# Сохранение таблицы ASV
write.csv(seqtab.nochim, "results/asv_table.csv")

##### ТАКСОНОМИЧЕСКАЯ КЛАССИФИКАЦИЯ #####

# Скачивание SILVA v138.1 (наиболее подходящая для рыб)
if (!file.exists("silva_nr99_v138.1_train_set.fa.gz")) {
  cat("Скачивание SILVA v138.1 training set...\n")
  download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
                "silva_nr99_v138.1_train_set.fa.gz")
}

if (!file.exists("silva_species_assignment_v138.1.fa.gz")) {
  cat("Скачивание SILVA species assignment...\n")
  download.file("https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz",
                "silva_species_assignment_v138.1.fa.gz")
}


cat("=== КЛАССИФИКАЦИЯ С SILVA ===\n")
taxa_silva <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", 
                             multithread = TRUE, tryRC = TRUE)

# Добавление видовых названий
taxa_silva <- addSpecies(taxa_silva, "silva_species_assignment_v138.1.fa.gz")


# Ключевые таксоны в микробиоме Atlantic salmon (по литературе):
salmon_key_taxa <- c(
  # Кишечник:
  "Carnobacterium",    # Доминирующий род в кишечнике
  "Aeromonas",         # Часто встречается
  "Mycoplasma",        # Важный род
  "Clostridiales",     # Порядок (анаэробы)
  "Lactobacillus",     # Пробиотические бактерии
  "Weissella",         # Молочнокислые бактерии
  
  # Кожа/слизь:
  "Flavobacterium",    # Доминирует на коже
  "Psychrobacter",     # Холодолюбивые бактерии
  "Pseudoalteromonas", # Морские бактерии
  "Vibrio",            # Потенциальные патогены
  
  # Общие:
  "Proteobacteria",    # Доминирующий тип
  "Firmicutes",        # Второй по важности тип
  "Bacteroidetes"      # Часто встречается
)

# Функция для анализа специфичных таксонов
analyze_salmon_taxa <- function(tax_table, key_taxa = salmon_key_taxa) {
  
  cat("=== АНАЛИЗ КЛЮЧЕВЫХ ТАКСОНОВ ЛОСОСЯ ===\n")
  
  # Поиск ключевых таксонов в результатах
  found_taxa <- list()
  
  for (taxon in key_taxa) {
    # Поиск во всех таксономических уровнях
    matches <- apply(tax_table, 1, function(x) any(grepl(taxon, x, ignore.case = TRUE)))
    if (sum(matches) > 0) {
      found_taxa[[taxon]] <- sum(matches)
      cat(sprintf("✅ %s: найдено %d ASV\n", taxon, sum(matches)))
    } else {
      cat(sprintf("❌ %s: не найдено\n", taxon))
    }
  }
  
  return(found_taxa)
}

# Анализ результатов SILVA
salmon_taxa_silva <- analyze_salmon_taxa(taxa_silva)


'''
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

'''

cat("\n=== ПОТЕНЦИАЛЬНЫЕ ПРОБЛЕМЫ ===\n")

# 1. Неклассифицированные последовательности
unclassified_counts <- apply(taxa_silva, 2, function(x) sum(is.na(x)))
cat("Неклассифицированные последовательности по уровням:\n")
print(unclassified_counts)

'''
Kingdom  Phylum   Class   Order  Family   Genus Species 
      0       9      20      65     155     362     656
'''
# 2. Процент неизвестных таксонов
unknown_percent <- round(unclassified_counts / nrow(taxa_silva) * 100, 1)
cat("\nПроцент неизвестных:\n")
print(unknown_percent)

'''
Kingdom  Phylum   Class   Order  Family   Genus Species 
    0.0     1.3     2.9     9.4    22.4    52.2    94.7 
    
'''

if (unknown_percent["Species"] > 50) {
  cat("⚠️  ВНИМАНИЕ: >50% ASV без видовой классификации\n")
  cat("   Это нормально для микробиома рыб \n")
}

#### файлы для PICRUSt2: ####

library(Biostrings)

# Экспорт ASV последовательностей
sequences <- getSequences(seqtab.nochim)
names(sequences) <- paste0("ASV_", 1:length(sequences))

# Сохранение FASTA файла
writeXStringSet(DNAStringSet(sequences), "asv_sequences.fasta")

# Создание таблицы ASV для PICRUSt2
asv_table_picrust <- seqtab.nochim
colnames(asv_table_picrust) <- names(sequences)

# Сохранение в TSV формате (требуется для PICRUSt2)
write.table(asv_table_picrust, "asv_table.tsv", 
            sep = "\t", quote = FALSE, col.names = NA)

