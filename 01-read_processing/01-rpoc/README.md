# Raw data processing

## Setup

Download raw data from ////

```bash
cd /home/allie/rpoc_comparison/01-read_processing/01-rpoc # on pickles
mkdir raw && cd raw
wget ////
cd ..
```

### 1. Install R packages (v4.1.0)

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("magrittr")
install.packages("stringr")
BiocManager::install("dada2")
install.packages("data.table")
install.packages("broom")
install.packages("qualpalr")
install.packages("seqinr")
```

NOTE: Problems with the stringi package? Try first installing it with conda (with the domhain environment activated)

```bash
conda install -c r r-stringi
```

Then update the conda environment

```bash
conda --update all
```

Load up R and try to load the stringi package

```R
library(stringi)
```

### 2. Load required libraries

```R
library(dada2)
library(stringr)
library(data.table)
library(broom)
library(qualpalr)
library(ShortRead)
library(Biostrings)
library(seqinr)
```

### 3. File path setup

```R
rawpath <- "raw"
wdpath <- "/home/allie/rpoc_comparison/01-read_processing/01-rpoc/" # change to where git repository was cloned
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
```

```text
 [1] "2L13-PD1"       "2L13-PE1"       "2L13-PF1"       "2L17-PE1"      
 [5] "2L17-PF1"       "2L2-PD1"        "2L2-PE1"        "2L23-PD1"      
 [9] "2L23-PE1"       "2L3-PF1"        "2L4-PF1"        "2L46-PE1"      
[13] "2L46-PF1"       "2L47-PE1"       "2L47-PF1"       "2L49-PF1"      
[17] "2L52-PF1"       "2L59-PF1"       "2L60-PE1"       "2L60-PF1"      
[21] "2L61-PD1"       "2L61-PE1"       "2L61-PF1"       "2L62-PF1"      
[25] "2L63-PE1"       "2L63-PF1"       "2L64-PD1"       "2L64-PE1"      
[29] "2L67-PF1"       "2L69-PF1"       "2L70-PF1"       "2L72-PE1"      
[33] "2L72-PF1"       "2L75-PD1"       "2L75PE1"        "4L59-PF1"      
[37] "Ecoli"          "Mock20Strain1"  "Mock20Strain2"  "Mock20Strain3" 
[41] "MockCommunity"  "MockCommunity2" "MockCommunity3" "PCRBlank16"    
[45] "PCRBlank9"      "Saureus"      
```

### 4. Plot quality scores

```R
system("mkdir img")
pdf(paste(wdpath, "img/", "forward_quality_plot.pdf", sep=""))
plotQualityProfile(fnFs[10:20])
dev.off()
pdf(paste(wdpath, "img/", "reverse_quality_plot.pdf", sep=""))
plotQualityProfile(fnRs[10:20])
dev.off()
```

### 5. Preliminary filter (removes sequences with N's)

```R
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)
```

### 6. Primer removal

```R
cutadapt <- as.character(system("which cutadapt", intern=T))
system("cutadapt --version")
path.cut <- file.path(rawpath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc("MAYGARAARMGNATGYTNCARGA")
REV.RC <- dada2:::rc("GMCATYTGRTCNCCRTCRAA")
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", "MAYGARAARMGNATGYTNCARGA", "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", "GMCATYTGRTCNCCRTCRAA", "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,"-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))
```

### 7. Filter and trim reads

```R
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimRight=25, maxN=c(0,0), maxEE=c(4,6), rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
retained
```

```text
                                       reads.in reads.out percentage_retained
2L13-PD1_S1_L001_R1_001.fastq.gz          53293     43742            82.07832
2L13-PE1_S7_L001_R1_001.fastq.gz          51802     41454            80.02394
2L13-PF1_S13_L001_R1_001.fastq.gz         50704     44991            88.73264
2L17-PE1_S19_L001_R1_001.fastq.gz         56847     50544            88.91234
2L17-PF1_S25_L001_R1_001.fastq.gz         56009     49541            88.45186
2L2-PD1_S4_L001_R1_001.fastq.gz           50309     40624            80.74897
2L2-PE1_S10_L001_R1_001.fastq.gz          62524     54860            87.74231
2L23-PD1_S5_L001_R1_001.fastq.gz          40489     33113            81.78271
2L23-PE1_S11_L001_R1_001.fastq.gz         60126     48884            81.30260
2L3-PF1_S34_L001_R1_001.fastq.gz          38961     33907            87.02805
2L4-PF1_S35_L001_R1_001.fastq.gz          36064     31520            87.40018
2L46-PE1_S20_L001_R1_001.fastq.gz         61738     51051            82.68975
2L46-PF1_S26_L001_R1_001.fastq.gz         68764     61084            88.83137
2L47-PE1_S21_L001_R1_001.fastq.gz         70837     62151            87.73805
2L47-PF1_S27_L001_R1_001.fastq.gz         55775     48134            86.30031
2L49-PF1_S36_L001_R1_001.fastq.gz         40668     34890            85.79227
2L52-PF1_S31_L001_R1_001.fastq.gz         59780     50357            84.23720
2L59-PF1_S32_L001_R1_001.fastq.gz         40894     36049            88.15230
2L60-PE1_S22_L001_R1_001.fastq.gz         54777     47508            86.72983
2L60-PF1_S28_L001_R1_001.fastq.gz         57276     48543            84.75278
2L61-PD1_S2_L001_R1_001.fastq.gz          50025     40404            80.76762
2L61-PE1_S8_L001_R1_001.fastq.gz          64138     49533            77.22879
2L61-PF1_S14_L001_R1_001.fastq.gz         62070     51385            82.78556
2L62-PF1_S15_L001_R1_001.fastq.gz         71605     63677            88.92815
2L63-PE1_S23_L001_R1_001.fastq.gz         74208     62830            84.66742
2L63-PF1_S29_L001_R1_001.fastq.gz         68862     57124            82.95431
2L64-PD1_S6_L001_R1_001.fastq.gz          39499     32855            83.17932
2L64-PE1_S12_L001_R1_001.fastq.gz         61432     49547            80.65341
2L67-PF1_S16_L001_R1_001.fastq.gz         65863     57544            87.36924
2L69-PF1_S17_L001_R1_001.fastq.gz         67115     58595            87.30537
2L70-PF1_S18_L001_R1_001.fastq.gz         65607     57322            87.37177
2L72-PE1_S24_L001_R1_001.fastq.gz         71292     62928            88.26797
2L72-PF1_S30_L001_R1_001.fastq.gz         52618     46624            88.60846
2L75-PD1_S3_L001_R1_001.fastq.gz          38170     30681            80.37988
2L75PE1_S9_L001_R1_001.fastq.gz           60637     51506            84.94154
4L59-PF1_S33_L001_R1_001.fastq.gz         42618     38122            89.45047
Ecoli_S39_L001_R1_001.fastq.gz            26162     22615            86.44217
Mock20Strain1_S3_L001_R1_001.fastq.gz      4573      3880            84.84583
Mock20Strain2_S4_L001_R1_001.fastq.gz      4579      3608            78.79450
Mock20Strain3_S5_L001_R1_001.fastq.gz      4350      3709            85.26437
MockCommunity_S37_L001_R1_001.fastq.gz    34030     32012            94.06994
MockCommunity2_S1_L001_R1_001.fastq.gz     4140      3577            86.40097
MockCommunity3_S2_L001_R1_001.fastq.gz     3695      3257            88.14614
PCRBlank16_S32_L001_R1_001.fastq.gz         152        60            39.47368
PCRBlank9_S40_L001_R1_001.fastq.gz          219        90            41.09589
Saureus_S38_L001_R1_001.fastq.gz          35241     34221            97.10564
```

### 8. Learn and plot error rates

```R
set.seed(12349)
errF <- learnErrors(filtFs, multithread=T, random=T)
errR <- learnErrors(filtRs, multithread=T, random=T)
png(paste(wdpath, "img/", "error_plot.png", sep=""))
plotErrors(errF, nominalQ=TRUE) 
dev.off()
```

### 9. Dereplication

```R
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# reassign sample names
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

### 10. Sample inference

```R
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

### 11. Filter out samples with fewer than 2500 reads

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 2500
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)]
paste(samples_to_remove)
```

```text
[1] "PCRBlank16" "PCRBlank9" 
```

### 12. Merge paired end reads

```R
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=T, minOverlap=10, trimOverhang=TRUE)
```

### 13. Construct sequence table

```R
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

```text
[1]   44 5149
```

### 14. Sequence length distribution plot

```R
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab))))
png(paste(wdpath, "img/", "length_hist.png", sep=""))
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()
```

### 15. Filter sequences by length

```R
minlen <- 400
maxlen <- 600 
seqlens <- nchar(getSequences(seqtab))
seqtab <- seqtab[,seqlens >= minlen & seqlens <= maxlen]
dim(seqtab)
```

```text
[1]   44 3539
```

### 15. Remove chimeras

```R
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```

```text
[1]   44 2427
[1] 0.9746875
```

### 16. Processing summary

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "non_chimeric")
rownames(track) <- sample.names[samples_to_keep]
track
```

### 17. Save output

```R
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set.fa")
system("awk '/^>/{print \">ASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set.fa")
```

### 18. Clean up ASV names

```R
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.merged.txt", row.names=FALSE, quote=F, sep="\t")
```

### 19. Assign taxonomy

Run kraken2

```bash
kraken2 --db ~/refdb/kraken_rpoc/ --threads 8 --use-names --output rep_set.kraken.out rep_set.fa --unclassified-out rep_set.unclassified.kraken.out --confidence 0.01
```

Unassigned reads may be primarily off target sequences, sequences with no close reference in the database. To get taxonomy for these, use kraken and the nt database (not used in diversity analyses, only for trees)

```bash
kraken2 --db ~/refdb/kraken_nt/ --threads 8 --use-names --output unassigned.kraken.out rep_set.unclassified.kraken.out --confidence 0.01
```

Get taxonomic lineage information from tax ids

```bash
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
# cd refdb && mkdir ncbi_taxonomy
# cd ncbi_taxonomy
# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
# tar xzf new_taxdump.tar.gz
# cd ../..
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > query
```

<!-- Add unidentified to taxonomy list (only do this once)

```bash
sed -i '1 i \0       |       unidentified    |               |' ~/refdb/ncbi_taxonomy/fullnamelineage.dmp
```

We also have to change taxid 81850 to 33958 (was merged on March 11, 2021). Open query file in nano and do CTL+W and replace string (only found once in file). NOTE! Depending on the version of the ncbi taxonomy file you pull and kraken version, you may have some issues with taxon names not being in the file. -->

Use queries to pull full lineage information (takes a little bit of time)

```bash
cat query | while read line; do grep -w -m 1 ^$line ~/refdb/ncbi_taxonomy/fullnamelineage.dmp | awk -F"|" '{print $3, $2}' | sed 's/\t//g' | sed 's/  / /g' | sed 's/cellular organisms; //' | sed 's/; /;/g' | sed 's/ /_/g'; done > lineages
```

Merge with ASV name to get taxonomy table

```bash
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids lineages > taxonomy.txt
```
Remove unwanted ASVs from taxonomy, sequence table

```bash
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python ~/rpoc_comparison/00-scripts/fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```

Build a tree from the filtered representative sequences

```bash
mafft rep_set.filt.fa > rep_set.align.fa
fasttree -nt rep_set.align.fa > rep_set.align.tre
```
