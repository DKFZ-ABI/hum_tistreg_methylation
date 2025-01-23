#!/bin/bash

# This script quantifies TE expression using TEspeX.
# Author: Niklas Beumer.


# Define a location on /yyy.
b330_space="/yyy"
location=$b330_space/yyy/hm_treg_bs_rgnsbg

# Environment settings.
source activate TEspeX_deps
module load zlib/1.2.12
export CFLAGS='-I/software/zlib/1.2.12/packages/include'
export LDFLAGS='-L/software/zlib/1.2.12/packages/lib'
export PATH=\
"$b330_space/nbeumer/Miniconda/miniconda3/envs/TEspeX_deps/bin:$PATH"

# Specify the RNA sample mapping file (only biological replicates).
rna_sample_mapping_file=$location/sample_mapping_rnaseq_only_biol_rep.txt

# Specify the directory holding RNA-seq Fastq files
fastq_dir=$b330_space/Regensburg/RNA_seq/evalrseq_agnes

# Specify the file containing the TE consensus sequences.
te_cons_fasta=\
$location/external_data/Giri_repbase_data_2023-11-03/\
Human_TEs_identified_by_browsing.fasta

# Specify the file containing full-gene-body-sequences of protein-coding 
# transcripts.
# Note: I'm using full-gene-body sequences instead of only exonic sequences
# to avoid artefacts arising from TEs that are in intronic sequences and might
# appear differentially expressed just because the gene they are located within
# is differentially expressed.
coding_transcr_fasta=\
$location/te_analysis/\
sequences_full_gene_bodies_gencode_v19_hs37d5_plus_strand.fa

# Specify the file containing lncRNA sequences derived from Gencode V19.
lncrna_fasta=\
$location/external_data/2023-12-19_gencode.v19.lncRNA_transcripts.fa.gz

# Identify the relevant RNA-seq sample names.
samples=`cut -f1 $rna_sample_mapping_file | grep -v Sample | grep p024`

# Specify the paths to the relevat RNA-seq Fastqs.
# For one of the samples, there is a typo in the Fastq file name. Thus, the
# sample name needs to be changed.
fastqs=$(\
for sample in $samples; do
    if [[ $sample == "p024_RNAseq_CCR8_Treg_Yellow5_S25_R1" ]]; then
        sample="p024_RNAseq_CCRR_Treg_Yellow5_S25_R1"
    fi
    actual_file_rel=`ls $fastq_dir | grep $sample | grep fastq$`
    actual_file_abs=$fastq_dir/$actual_file_rel
    echo $actual_file_abs
done)

# Generate the file that holds all the Fastq file names.
fastq_overview_file=\
$location/te_analysis/cons_seq_fastq_filenames_for_tespex.txt
rm $fastq_overview_file
for fastq in $fastqs; do
    echo $fastq >> $fastq_overview_file
done

# Specify the output directory that will later hold all TEspeX output.
output_dir=$location/te_analysis/cons_seq_tespex_results_bulk_rna_all

# Run TEspeX.
# Use lncRNA sequences as the non-coding transcripts according to what is 
# suggested in https://github.com/fansalon/TEspeX; 19-Dec-2023.
# Use the "--remove F" option to keep all Bam files.
python3 $b330_space/nbeumer/my_tespex_installation/TEspeX/TEspeX.py \
--TE $te_cons_fasta \
--cdna $coding_transcr_fasta \
--ncrna $lncrna_fasta \
--sample $fastq_overview_file \
--paired F \
--length 76 \
--out $output_dir \
--strand 'no' \
--num_threads 20 \
--remove F