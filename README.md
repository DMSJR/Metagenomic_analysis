# Mapping bacterial sialic acid metabolism pathways in sequencing data of inflammatory bowel diseases and _Clostridioides difficile_ infection.
<div align="center">
<img width="193" alt="Captura de Tela 2025-06-09 às 09 58 10" src="https://github.com/user-attachments/assets/8818ef69-0de2-4283-8e07-f37f7458a15b" />
<img width="182" alt="Captura de Tela 2025-06-09 às 09 57 56" src="https://github.com/user-attachments/assets/1e085ad8-f184-4e28-a370-d6d9d940490e" />
<img width="182" alt="Captura de Tela 2025-06-09 às 09 57 56" src="https://github.com/user-attachments/assets/9fc3cb8c-c59b-44e3-ad68-3b1d5779e4c9" />
</div>

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

In order to remove sequences from the human host, bowtie2 2.4.1 was used with the GRCh38 genome. The genome index was created using:

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

**2.6. Coding DNA sequeces (CDS) prediction**

Prodigal 2.6.3 was used to predict CDSs using the following command:

`prodigal -i {accession_number}_output/final.contigs.fa -o {accession_number}_coordinates -a {accession_number}_proteins.faa -d {accession_number}_genes.fasta -f gff`

**2.7. Search for Sia-related proteins by alignment**

Diamond 0.9.24 was used to search for genes related to Sia metabolism and transport. The selection of reference sequences is described in section 2.8. To create databases for aligment the following command was used, for each one of the files containing reference sequences:

diamond makedb --in {reference_protein_name}.fasta -d reference_protein_name

The alignment was performed using the following command:

`diamond blastp -q {reference_protein_name}_proteins.faa -d /diamond_db/{protein} -o {accession_number}_{protein}.tsv -f 6 qseqid qlen sallseqid slen evalue length pident`

**2.8. Selection of reference protein sequences**

To gather reference protein sequences, the annotated known DNA sequences were downloaded from NCBI and used on a blastx search at Uniprot. Sequences with at least 40% identity and coverage were used for the reference database, in agreement with Pearson (2013) work (doi: 10.1002/0471250953.bi0301s42). The list of proteins searched can be found in the file protein_list.txt

**2.9. Selection of aligned proteins**

The proteins aligned to the reference, with at least 40% identity and coverage, were used for further analysis. To obtain the sequence of such sequences, it was used the script get_cds_prodigal_diamond.py, that can be found in the repository, with the following command:

`python3 get_cds_prodigal_diamond.py ../diamond/{acession_number}_{protein}.tsv ../prodigal/{acession_number}_{protein}.faa {acession_number}_{protein}_proteins_aligned.faa`

**2.10. Taxonomic identification of detected sequences**

Diamond 0.9.24 was used for aligning the selected CDSs against the NCBI nr database. The database was created with the command:

'diamond makedb --in nr.faa -d nr_diamond --taxonmap  prot.accession2taxid.gz --taxonnodes nodes.dmp'

The alignment was performed with the following command:

`diamond blastp -q ../aligned_proteins/{acession_number}_{protein}_proteins_aligned.faa -d ../nr-faa/nr_diamond -o {acession_number}_{protein}_diamond_tax.tsv --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --max-target-seqs 10 --evalue 1e-5 --threads 4`

Afterwards, blast2lca was used to perform the taxonomic attribution, with the following command:

`blast2lca -i {acession_number}_{protein}_diamond_tax.tsv -o {acession_number}_{protein}_blast2lca.txt --mapDB megan-nucl-Feb2022.db`

**2.11. Calculation of sequence abundance**

To calculate the abundance of reads for each gene, the CDSs predicted on step 2.6. were used. Bowtie2 2.4.1 was used for the alignment with the following command to create the indexes:

`bowtie2-build ../prodigal/{acession_number}_genes.fasta {acession_number}_index`

and the following for alignemnt:

`bowtie2 -x {acession_number}_index -1 {accession_number}_paired_bt2_filtered.1.gz -2 {accession_number}_paired_bt2_filtered.2.gz  -S {accession_number}.sam --very-sensitive`

The sam files were converted to bam with samtools command:

`samtools view -bS {acession_number}.sam | samtools sort -o {acession_number}.sorted.bam`

and the reads were counted with the command:

`samtools idxstats {acession_number}.sorted.bam > {acession_number}_counts.txt`

The normalization was performed using the script python_normalization_V2.py, that can be found in the repository, with the following command:

`python3 python_normalization_v2.py {acession_number}_counts.txt  {acession_number}_normalized_v2.txt`

The counting was carried out with the script python_count.py, with the following command:

`python3 python_count.py {acession_number}_normalized_v2.txt ../diamond/{acession_number}_{protein}.tsv {acession_number}_{protein}_counted.tsv`

**2.12. Other analysis**

The output from step 2.10 and 2.11 (abundance_data.csv) were  used for the analysis in the Metagenomic_analysis.ipynb, together with the conditions.csv file. Please refer to the notebook for further information.

