## Arielome Proccesing Pipeline Documentation

Here, I will describe the process carried out by the pipeline. I will choose a toy example of input files: `f1.fq f2.fq`, that required to be processed.

#### 1. BarcodeSpliting 

The first step is to split the reads by barcode, each specie contains a unique barcode. The barcodes can be found at *data/barcodes_spitting.txt*. 

Barcode splitter is not run for `r2.fq`, We use `bioawk` to retrive mate 2 based on `specie-r1.fq`  .

```
cat r1.fq | fastx_barcode_splitter.pl --bcfile data/barcodes_spitting.txt --prefix out/ --suffix _01.fq --bol
bioawk -cfastx '{print $name}' test/pombe_01.fq | seqtk subseq r2.fq - > test/pombe_02.f
bioawk -cfastx '{print $name}' test/planaria_01.fq | seqtk subseq r2.fq - > test/planaria_02.fq
bioawk -cfastx '{print $name}' test/fish_01.fq | seqtk subseq r2.fq - > test/fish_02.fq
bioawk -cfastx '{print $name}' test/fly_01.fq | seqtk subseq r2.fq - > test/fly_02.fq
bioawk -cfastx '{print $name}' test/minigene_01.fq | seqtk subseq r2.fq - > test/minigene_02.fq
```

#### 2. RemoveAdapters

Each specie has a unique adapter, (adapter can be found at the helper script *helper.py*). Before mapping adapters are removed from reads.

From now, I will show the commands that are run for only one speci, but keep in mind that this is run for the four species and in particular cases for the minigenes.


```
cutadapt -g TCTAACGGCGAAATGGC test/fish_01.fq -o test/fish_01_trimmed.fq
cutadapt -g TTAGTCACCTA test/fish_02.fq -o test/fish_02_trimmed.fq
```

#### 3. Mapping

We use bowtie to map the reads to their correspinding transcriptome depending on species (transcriptome files: *data/transcriptome_references/*)

```
bowtie -n 1 --seedlen 10 -I 200 -X 800 --threads 10 data/transcriptome_references/fish \
-1 test/fish_01_trimmed.fq -2 test/fish_02_trimmed.fq -S test/fish.sam
```

#### 4. FilterAndSort

We need to retrive the exact mapping sequence, from the reference transcriptome. Before that reads are sorted by read name.

```
samtools view -f 0x02 -Sb test/fish.sam | samtools sort -n -o test/fish_filter.bam
```

#### 5. ExtractSeqs

```
bamToBed -bedpe -mate1 -bed12 -i test/fish_filter.bam | cut -f 1,2,6,7 | awk '$2<$3 {print ;}' > test/fish.bed
bedtools getfasta -fi data/transcriptome_references/fish.fa -bed test/fish.bed -name > test/fish.fasta
```

