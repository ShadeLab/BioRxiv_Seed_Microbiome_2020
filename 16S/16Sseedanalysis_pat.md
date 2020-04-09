# Raw Sequences Analysis of 16S rRNA gene from Pat data

## Analysis of 16S Miseq Data
Here we used the USEARCH pipeline (v10.0.240_i86linux64) for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). Additional analyses were conducted using QIIME

# Part I: Clustering

usearch v10.0.240_i86linux64, 16.3Gb RAM, 4 cores
(C) Copyright 2013-17 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

## 1) Merge Paired End Reads
```
# decompress the reads
gunzip *.gz

# make directory called "mergedfastq"
mkdir mergedfastq

# merge paired end reads
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_mergepairs seed/raw_reads/*R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout newseedtest/mergedfastq/merged.fq -tabbedout newseedtest/mergedfastq/merged.report.txt -alnout newseedtest/mergedfastq/merged_aln.txt

# -fastq_maxdiffs 10: Allow 10 max differences in overlap region

### Output ###

 1109691  Pairs (1.1M)
    962109  Merged (962.1k, 86.70%)
    477092  Alignments with zero diffs (42.99%)
    143642  Too many diffs (> 10) (12.94%)
      3940  No alignment found (0.36%)
         0  Alignment too short (< 16) (0.00%)
      1832  Staggered pairs (0.17%) merged & trimmed
    246.81  Mean alignment length
    252.80  Mean merged length
      0.28  Mean fwd expected errors
      1.15  Mean rev expected errors
      0.09  Mean merged expected errors
```
## 2) Check Sequence Quality of Merged Seqs
```
mkdir fastq_info
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 mergedfastq/merged.fq -output fastq_info/eestats.txt

### output ###
962109 reads, max len 439, avg 252.8

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50     957985( 99.6%)     961922(100.0%)     962068(100.0%)
   100     948055( 98.5%)     960567( 99.8%)     962023(100.0%)
   150     941546( 97.9%)     957702( 99.5%)     960720( 99.9%)
   200     936264( 97.3%)     956549( 99.4%)     960632( 99.8%)
   250     920142( 95.6%)     951342( 98.9%)     960063( 99.8%)
   300         87(  0.0%)         98(  0.0%)        111(  0.0%)
   350         24(  0.0%)         30(  0.0%)         38(  0.0%)
   400          5(  0.0%)          7(  0.0%)         10(  0.0%)
```
## 3) Filter and Truncate the Merged Seqs
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter mergedfastq/merged.fq -fastq_maxee 1 -fastq_trunclen 250 -fastqout filtered_merged.fq

### output ###
100.0% Filtering, 98.9% passed
    962109  Reads (962.1k)                  
      1326  Discarded reads length < 250
      9441  Discarded reads with expected errs > 1.00
    951342  Filtered reads (951.3k, 98.9%)
```
## 4) Dereplicate Sequences
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_uniques filtered_merged.fq -fastqout uniques_filtered_merged.fastq -sizeout

### output ###
951342 seqs, 73768 uniques, 47786 singletons (64.8%)
Min size 1, median 1, max 416075, avg 12.90
```
## 5) Remove Singeltons
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sortbysize uniques_filtered_merged.fastq -fastqout nosigs_uniques_filtered_merged.fastq -minsize 2

### output ###
Sorting 25982 sequences
```
## 6) Precluster Sequences
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -cluster_fast nosigs_uniques_filtered_merged.fastq -centroids_fastq denoised_nosigs_uniques_filtered_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size

### output ###
Seqs  25982 (26.0k)
  Clusters  555
  Max size  712259 (712.3k)
  Avg size  1628.0
  Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
   Max mem  672Mb
      Time  1.00s
Throughput  26.0k seqs/sec.
```
## 7) Closed Reference-based OTU Picking Using SILVA_132 Database
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -usearch_global denoised_nosigs_uniques_filtered_merged.fastq -id 0.97 -db /mnt/research/ShadeLab/WorkingSpace/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -strand plus -uc ref_seqs.uc -dbmatched SILVA_closed_reference.fasta -notmatchedfq failed_closed.fq

### output ###
100.0% Searching, 63.6% matched

Closed reference OTU picking:
Pick OTUs based on the Silva database in the home directory.
Produce some output files - ref_seqs.uc (pre-clustered), SILVA_closed_reference.fasta will be the matched ones, and failed_closed.fq will be used in de novo OTU picking
```
## 8) De novo OTU picking
```
# sort by size
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sortbysize failed_closed.fq -fastaout sorted_failed_closed.fq

# cluster de novo
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -cluster_otus sorted_failed_closed.fq -minsize 2 -otus denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up

### output ###
100.0% 46 OTUs, 48 chimeras
```
## 9) Combine the Rep Sets Between De novo and SILVA Reference-based OTU Picking
```
cat SILVA_closed_reference.fasta denovo_otus.fasta > FULL_REP_SET.fna
```
## 10) Map 'FULL_REP_SET.fna' Back to Pre-dereplicated Sequences and Make OTU Tables
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64  -usearch_global mergedfastq/merged.fq -db FULL_REP_SET.fna -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom

### output ###
961230 / 962109 mapped to OTUs (99.9%) 
```
# Part II: Switch to QIIME 1.9.1

## 1) Assign taxonomy to SILVA_132_QIIME_release with UCLUSTq
```
assign_taxonomy.py -i FULL_REP_SET.fna -o taxonomy -r /mnt/home/bintarti/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -t /mnt/home/bintarti/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt
#-r reference -> path to silva db 
# -t taxonomy
```
## 2) Convert OTU_table.txt. to OTU_table.from_txt_json.biom
```
biom convert -i OTU_table.txt -o OTU_table.biom --table-type="OTU table" --to-json
```
## 3) Add taxonomy to OTU table
```
biom add-metadata -i OTU_table.biom -o OTU_table_tax.biom --observation-metadata-fp=taxonomy/FULL_REP_SET_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy
```
## 4) Filter non-bacteria/archaea
```
filter_taxa_from_otu_table.py -i OTU_table_tax.biom -o OTU_table_tax_filt.biom -n D_4__Mitochondria,D_3__Chloroplast,Unassigned

# summarize OTU table

biom summarize-table -i OTU_table_tax.biom -o OTU_table_tax_sum.txt

biom summarize-table -i OTU_table_tax_filt.biom -o OTU_table_tax_filt_sum.txt
```
## Rarefaction
```
single_rarefaction.py -d 11138 -o single_rare.biom -i OTU_table_tax_filt.biom

biom summarize-table -i single_rare.biom -o single_rare_sum.txt
```
## 5) Convert and add taxonomy
```
biom convert -i single_rare.biom -o single_rare_filt.txt --header-key taxonomy --to-tsv
```









