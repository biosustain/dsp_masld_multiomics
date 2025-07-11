# CD-HIT Analysis

## Input Data
- X.

## Renaming: Adding "MASLD" to Gene Names
cp masldmice_genes.fna MASLD
awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' MASLD > MASLD_mice_renamed.fasta

## Check Count of Genes in MASLD
# Total gene count:
10043913 MASLD

## Calculate Gene Lengths
https://github.com/Juassis/mg_catalog_functions/blob/main/count_gene_length.py masldmice_genes.fna > masldmice_genes.fna_length
cat masldmice_genes.fna_length | sort -k2n | awk '$2 >= 90' | cut -f1 > extract_f

## Check Count of Genes After Filtering by Length
# Genes with length >= 90:
9803556 extract_f

## Extracting Sequences by Length
makeblastdb -in MASLD_mice_renamed.fasta -dbtype nucl -parse_seqids -out MASLD
blastdbcmd -db MASLD -entry_batch extract_f -out extracted_sequences.fasta

## EMGC Processing

### Check Count of Genes
# Total gene count in EMGC:
5862027 EMGC.fasta

### Renaming: Adding "EMGC" to Gene Names
awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' EMGC > EMGC_renamed.fasta

## Merging Gene Sets
cat extracted_sequences.fasta EMGC_renamed.fasta > merged_gene_sets_sorted.fna

### Check Count of Merged Genes
# Total gene count after merging:
15665583

## Sort Merged Genes by Length
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' merged_gene_sets.fna | \
    awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' | \
    sort -k1,1n | cut -f 2- | tr "\t" "\n" > merged_gene_sets_sorted.fna

## Running CD-HIT

### Job Submission Script
### General options
#BSUB -J OpenMPjob          # Set the job Name
#BSUB -n 4                  # Number of cores
#BSUB -R "span[hosts=1]"    # Ensure all cores are on a single host
#BSUB -W 24:00              # Set walltime limit
#BSUB -R "rusage[mem=18GB]" # Request memory
#BSUB -u jasge       # Email notifications
#BSUB -B                    # Notify at start
#BSUB -N                    # Notify at completion
#BSUB -o Output_%J.out      # Output file
#BSUB -e Output_%J.err      # Error file

# Set number of threads
OMP_NUM_THREADS=$LSB_DJOB_NUMPROC
export OMP_NUM_THREADS

# Load necessary modules
module load mpi/4.0.7-gcc-9.5.0-binutils-2.38

# Activate Conda environment
source /zhome/2e/3/214187/anaconda3/bin/activate
conda init bash
conda activate dartqc

# Run cd-hit-est tool
mpirun cd-hit-est -T 4 -M 20000 -i /work3/jasge/02_Data/merged_gene_sets_sorted.fna \
  -o /work3/jasge/02_Data/results/clustered_mice_gut_catalog.fa -c 0.95 -n 8 -G 0 \
  -aS 0.9 -g 1 -r 1 -d 0