# Mapping bacterial sialic acid metabolism pathways in sequencing data of inflammatory bowel diseases and _Clostridioides difficile_ infection.

**1. Introduction**

This project aims to reprocess and reanalyse public data from inflammatory bowel diseases (IBD), Crohn's Disease (CD) and Ulcerative Colitis (UC), and _Clostridioides difficile_ infection (CDI), searching for associations between those diseases and microbial genes related to sialic acid (Sia) in the gut microbiota.

**2. Methodology**

| ![poster (1)](https://github.com/user-attachments/assets/b5f0aeff-9e5c-48f0-95af-3de81b5fc2af) |
|:------------------------------------------------------------------------------------------:|
| **Pipeline of metagenomic data reprocessing and analysis**                                    |

**2.1. Data acquisition**

The raw data was downloaded from SRA using SRAtool, with the following command:

`fasterq-dump --split-files {accession_number}`

Or from ENA, with the following command:

`wget -nc {data_url}`

The folder ./Metagenomic_analysis/Data_acession_numbers/ in this repository presents the codes of all sequences utilized.

**2.2. Trimming**

Data trimming was performed using Trimmomatic 0.40, using TruSeq3-PE adapters, with the following code:

`java -jar  trimmomatic-0.40-rc1.jar PE {raw_data_1} {raw_data_2} {accession_number}_forward_paired.fastq.gz {accession_number}_forward_unpaired.fastq.gz {accession_number}_reverse_paired.fastq.gz {accession_number}_reverse_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36`

**2.3. Filtering human contaminant sequences**

In order to remove sequences from the human host, bowtie2 2.4.1 was used. The genome index was created using:

`bowtie2-build GRCh38.fa human_genome`

and with the following command, human sequences were removed:

`bowtie2 -x human_genome -1 {accession_number}_forward_paired.fastq.gz -2 {accession_number}_reverse_paired.fastq.gz --un-conc {accession_number}_paired_bt2_filtered\n`

**It's important to note that the data downloaded from public databases had already gone through these processes, since trimming and filtering removed next to no sequence.**

**2.4. Pair-end merging**

Flash 1.2.11 was used to merge the overlapping pair-end sequences, to assist assembly afterwards. The following command was used:

`flash -o {accession_number}  {accession_number}_paired_bt2_filtered.1  {accession_number}_paired_bt2_filtered.2`

**2.5. Assembly**

To assemble the reads, Megahit 1.2.9 was used, with the following command:

`megahit -r {accession_number}.extendedFrags.fastq -1 {accession_number}_paired_bt2_filtered.1 -2 {accession_number}_paired_bt2_filtered.2 -o {accession_number}_output`

**2.6. Coding DNA sequece prediction**

Prodigal 2.6.3 was used to predict CDSs using the following command:

`prodigal -i {accession_number}_output/final.contigs.fa -o {accession_number}_coordinates -a {accession_number}_proteins.faa -d {accession_number}_genes.fasta -f gff`

**
