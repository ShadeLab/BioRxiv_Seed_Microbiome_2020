# Raw Sequences Analysis of ITS gene from bean seed of pat's pilot experiment
# Date: April 2nd 2020
# Author: A. Fina Bintarti

## Analysis of ITS Miseq Data
Here we used the USEARCH pipeline (v10.0.240_i86linux64) for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). Additional analyses were conducted using QIIME

raw sequence data stored on HPCC
/mnt/research/ShadeLab/Sequence/raw_sequence/.....

Moved/copy raw sequences (22 samples: 8 control, 7 water withholding, 7 nutrient addition) to the working space:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/PatSeedData/newseedtest_its

## All analysis results are stored on hpcc: 
"/mnt/research/ShadeLab/WorkingSpace/Bintarti/PatSeedData/newseedtest_its"

### 1. Quality checking and pre-filtering
```
# count read numbers
for fastq in rawreads/*.fastq
do wc -l $fastq
done > reads_raw.counts

# produce reads quality graphs using FastQC
mkdir stats

cat rawreads/*R1_001.fastq > raw_reads_R1.fastq; cat rawreads/*R2_001.fastq > raw_reads_R2.fastq

module load FastQC/0.11.7-Java-1.8.0_162

fastqc raw_reads_R1.fastq raw_reads_R2.fastq -o stats && rm -rf raw_reads_R1.fastq raw_reads_R2.fastq
```
### 2. Merge paired end reads
```
mkdir mergedfastq

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_mergepairs rawreads/*R1*.fastq -relabel @ -tabbedout merged_tabbed.txt -report merged_summary.txt -fastqout mergedfastq/merged.fastq

# output
   1258147  Pairs (1.3M)
    979816  Merged (979.8k, 77.88%)
    547978  Alignments with zero diffs (43.55%)
    268010  Too many diffs (> 5) (21.30%)
     10321  No alignment found (0.82%)
         0  Alignment too short (< 16) (0.00%)
    155721  Staggered pairs (12.38%) merged & trimmed
    192.96  Mean alignment length
    268.94  Mean merged length
      0.87  Mean fwd expected errors
      1.12  Mean rev expected errors
      0.06  Mean merged expected errors
```
### 3. Check sequence quality of merged sequences using Usearch and Vsearch
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 mergedfastq/merged.fastq -output stats_eestats2_USEARCH.txt

# output
979816 reads, max len 466, avg 268.9

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50     938988( 95.8%)     940171( 96.0%)     940254( 96.0%)
   100     887920( 90.6%)     891542( 91.0%)     891796( 91.0%)
   150     884284( 90.3%)     891134( 90.9%)     891773( 91.0%)
   200     881618( 90.0%)     890664( 90.9%)     891701( 91.0%)
   250     822216( 83.9%)     831839( 84.9%)     833020( 85.0%)
   300     421989( 43.1%)     427750( 43.7%)     429001( 43.8%)
   350       6730(  0.7%)       7613(  0.8%)       8245(  0.8%)
   400       2332(  0.2%)       2741(  0.3%)       3073(  0.3%)
   450        984(  0.1%)       1203(  0.1%)       1397(  0.1%)
###############################################################

module load vsearch/2.9.1

vsearch -fastq_stats mergedfastq/merged.fastq -fastq_qmax 42 -log stats_results_VSEARCH.txt
```
### 4. Remove primer and adapters with cutadapt
```
################
CS1-ITS1 (fwd): 5’- CTTGGTCATTTAGAGGAAGTAA – 3’ (EMP/Smith and Peay 2014)
CS2-ITS2 (rev): 5’- GCTGCGTTCTTCATCGATGC – 3’ (EMP/Smith and Peay 2014)
Reverse complement of adapter-reverse primer: GCATCGATGAAGAACGCAGC
################

module load cutadapt (v.2.0)

cutadapt -g CTTGGTCATTTAGAGGAAGTAA -a GCATCGATGAAGAACGCAGC -f fastq -n 2 --discard-untrimmed --match-read-wildcards -o cut_merged.fastq mergedfastq/merged.fastq > cut_adpt_results.txt

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 cut_merged.fastq -output cutdapt_eestats2_USEARCH.txt

# output
979749 reads, max len 444, avg 227.0

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50     890410( 90.9%)     891501( 91.0%)     891538( 91.0%)
   100     886927( 90.5%)     891286( 91.0%)     891534( 91.0%)
   150     883551( 90.2%)     890897( 90.9%)     891521( 91.0%)
   200     823971( 84.1%)     832151( 84.9%)     832983( 85.0%)
   250     448717( 45.8%)     453725( 46.3%)     454568( 46.4%)
   300      13644(  1.4%)      14771(  1.5%)      15443(  1.6%)
   350       2386(  0.2%)       2792(  0.3%)       3112(  0.3%)
   400       1611(  0.2%)       1944(  0.2%)       2203(  0.2%)
```
### 5. Quality filter
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 cut_merged.fastq -output cut_merged.pre_filtered.eestats2.txt -length_cutoffs 100,400,10

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter cut_merged.fastq -fastq_minlen 150 -fastq_maxee 1 -fastaout cut_merged_filtered.fa -fastaout_discarded merged.no_filter.fa -fastqout cut_merged_filtered.fastq

# output
100.0% Filtering, 90.8% passed
    979749  Reads (979.7k)                  
      1974  Discarded reads with expected errs > 1.00
    889564  Filtered reads (889.6k, 90.8%)
```
#################################################################
```

# convert fastq to fasta for ITSx 
module load BBMap/37.93
reformat.sh in=cut_merged.fastq out=cut_merged.fasta

# install Pearl
module load  GCCcore/8.3.0
module load  Perl/5.30.0

# install HMMER
conda install -c bioconda hmmer

# install ITSx
conda install -c bioconda itsx

# To test if ITSx was successfully installed type "ITSx --help" on the command-line. 
# To check for ITS sequences in the test file, type "ITSx -i test.fasta -o test --cpu 2"
```
#################################################################


### 6. Find the set of unique sequences (dereplication)
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_uniques cut_merged_filtered.fastq -fastaout derep_filtered_cut_merged.fasta -sizeout

#output
889564 seqs, 73234 uniques, 45565 singletons (62.2%)
```
### 7. Open reference-based OTU picking (using UNITE_v.8.0 at 97% identity treshhold)
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -usearch_global derep_filtered_cut_merged.fasta -id 0.97 -db /mnt/research/ShadeLab/UNITE_v.8.0/sh_refs_qiime_ver8_97_02.02.2019.fasta  -strand plus -uc ref_seqs.uc -dbmatched UNITE_reference.fasta -notmatched UNITE_failed_closed.fq

# output
100.0% Searching, 74.8% matched
```
### 8. Sorting by size and de novo-based OTU picking on sequences that failed to hit reference
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64  -sortbysize UNITE_failed_closed.fq -fastaout sorted_UNITE_failed_closed.fq

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64  -cluster_otus sorted_UNITE_failed_closed.fq -minsize 2 -otus de_novo_otus.fasta -uparseout uparse_otus.txt -relabel OTU_

# output
100.0% 151 OTUs, 66 chimeras
```
### 9. Combine the seqs of de novo and reference-based OTU picking
```
cat UNITE_reference.fasta de_novo_otus.fasta > REP_seq.fna

# numbering the OTUs for CONSTAX input
/mnt/home/bintarti/python_scripts-master/fasta_number.py REP_seq.fna OTU_ > 04022020_pat.NUMB_REP_seq.fasta
```
### 10. Construct OTU table
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -usearch_global mergedfastq/merged.fastq -db 04022020_pat.NUMB_REP_seq.fasta -strand plus -id 0.97 -uc OTU_map.uc -otutabout 04022020_pat.OpenRef_OTU_table.txt

# output
878716 / 979816 mapped to OTUs (89.7%) 
```
Taxonomic classification using CONSTAX
Please refer to how ‘Running CONSTAX on the MSU HPCC’ on lab guru: https://my.labguru.com/knowledge/documents/330. Here is the publication about CONSTAX tool https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1952-x and how to use it https://github.com/natalie-vandepol/compare_taxonomy. Note: CONSTAX uses python 2.7 to be able to run. HPCC default python is python3 and you can install python 2.7 on hpcc
```
# check python version
python --version

# if you have python3 installed and want to swap to python2
conda create -n python2 python=2.7 anaconda
conda activate python2

# download "CONSTAX_hpcc.tar.gz" from the lab guru link above. Put and extract the file on your home directory on hpcc.
# open "CONSTAX_hpcc" directory and follow the instructions from the lab guru link above.
# use the file 'consensus_taxonomy.txt' in the 'outputs' directory as your taxonomy table.
```
### 11. Convert .txt file to .biom file (if you need .biom file!)
```
biom convert -i 04022020_pat.OpenRef_OTU_table.txt -o 04022020_pat.OpenRef_OTU_table.biom --table-type="OTU table" --to-json

# summarize OTU table
biom summarize-table -i 04022020_pat.OpenRef_OTU_table.biom -o 04022020_pat.OpenRef_OTU_table_sum.txt
```
### 12. Check any eukaryotes contaminant
```
# inspect the 'consensus_taxonomy.txt' generated by CONSTAX tool. 

grep "zoa" consensus_taxonomy.txt

# found two OTUs assigned to Cercozoa (OTU_286 & OTU_253)

grep "Chloro","Mito" consensus_taxonomy.txt  

# transfer 'OpenRef_OTU_table.txt' and 'consensus_taxonomy.txt' to your local working directory.
# remove the contaminant OTUs from taxonomy table and OTU table for further analysis on R.
```
### 13. Rarefy OTU table to the lowest sequencing depth
```
single_rarefaction.py -d 21329 -o 04022020_pat.OpenRef_OTU_table_rare.biom -i 04022020_pat.OpenRef_OTU_table.biom

biom summarize-table -i 04022020_pat.OpenRef_OTU_table_rare.biom -o 04022020_pat.OpenRef_OTU_table_rare_sum.txt

biom convert -i 04022020_pat.OpenRef_OTU_table_rare.biom -o 04022020_pat.OpenRef_OTU_table_rare.txt --to-tsv
```
Here is a useful link https://wiki.hpcc.msu.edu/display/ITH/File+transfer#Filetransfer-UsingFileZilla(forMacandWindows) for file transfer from MSU HPCC to local computer.









