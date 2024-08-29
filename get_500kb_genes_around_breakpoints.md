# Getting genes from 500kb region on both side of breakpoints

## Got the breakpoints based on genome comparison of Ae. umbellulata with Ae. tauschii

## get the genes 500kb upstream and 500kb downstream of breakpoint in Ae. umbellulata

### Step I: Get the mRNAs 
```bash
awk '$1 == "Chr5U" && $4 >= 594670863 && $5 <= 595670863 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_1a.gff3
awk '$1 == "Chr6U" && $4 >= 29254356 && $5 <= 30254356 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_1b.gff3
awk '$1 == "Chr4U" && $4 >= 144277208 && $5 <= 145277208 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_2.gff3
awk '$1 == "Chr6U" && $4 >= 246059785 && $5 <= 247059785 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_3.gff3
awk '$1 == "Chr4U" && $4 >= 266709752 && $5 <= 267709752 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_4.gff3
awk '$1 == "Chr4U" && $4 >= 89152159 && $5 <= 90152159 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_5.gff3
awk '$1 == "Chr2U" && $4 >= 546858291 && $5 <= 547858291 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_6.gff3
awk '$1 == "Chr4U" && $4 >= 59341477 && $5 <= 60341477 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_7.gff3
awk '$1 == "Chr7U" && $4 >= 439090183 && $5 <= 440090183 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_8.gff3
awk '$1 == "Chr7U" && $4 >= 588579447 && $5 <= 589579447 && $3 == "mRNA"' HC_brakerRNAref_IsoSeqcomp_functional_renamed.gff3 > aumb_event_9.gff3
```
### Step II: Get the mRNA names
```bash
for file in aumb_event_*; do
  cut -f 9 $file | cut -d ";" -f 1 | sed 's/ID=//g' > $(basename $file .gff3)_mRNA_names.txt
done
```

### Step III: Get the number of genes in each event
```bash
for file in aumb_event_*; do
  wc -l $file
done
```

### Step IV: Get the protein sequences of each gene (mRNA)
```bash
for file in *_mRNA_names.txt; do
  seqtk subseq /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/GrainGenes2/HC_brakerRNAref_IsoSeqcomp_functional_renamed.pep $file > $(basename $file _names.txt).pep
done
```

### Step IV: Get the CDS sequences of each gene (mRNA)
```bash
for file in *_mRNA_names.txt; do
  seqtk subseq /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/GrainGenes2/HC_brakerRNAref_IsoSeqcomp_functional_renamed.cds $file > $(basename $file _names.txt).cds
done
```

### Step V: Get genomic region 1Mb
```bash
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr5U:594670863-595670863" > aumb_event_1a_genomic_region.fasta
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr6U:29254356-30254356" > aumb_event_1b_genomic_region.fasta
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr4U:144277208-145277208" > aumb_event_2_genomic_region.fasta
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr6U:246059785-247059785" > aumb_event_3_genomic_region.fasta
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr4U:266709752-267709752" > aumb_event_4_genomic_region.fasta
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr4U:89152159-90152159" > aumb_event_5_genomic_region.fasta
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr2U:546858291-547858291" > aumb_event_6_genomic_region.fasta
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr4U:59341477-60341477" > aumb_event_7_genomic_region.fasta
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr7U:439090183-440090183" > aumb_event_8_genomic_region.fasta
samtools faidx /mmfs1/projects/upinder.gill/AUMB_PROJECT_DATA/FINAL_AeU_asm.fasta "Chr7U:588579447-589579447" > aumb_event_9_genomic_region.fasta
```

### Step VI: Python code to combine the gff3, protein, and cds information in one csv file (get_info_in_csv.py)
```python3
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Process GFF3 and sequence files to add protein and CDS sequences and mark breakpoints.")
    parser.add_argument("-g", "--gff3", required=True, help="GFF3 file")
    parser.add_argument("-p", "--pep", required=True, help="Protein sequence file")
    parser.add_argument("-c", "--cds", required=True, help="CDS sequence file")
    parser.add_argument("-b", "--breakpoint", type=int, required=True, help="Breakpoint position")
    return parser.parse_args()

def load_sequences(file):
    seq_dict = {}
    current_id = None
    with open(file, "r") as f:
        for line in f:
            if line.startswith(">"):
                current_id = line.split()[0][1:]
            else:
                if current_id in seq_dict:
                    seq_dict[current_id] += line.strip()
                else:
                    seq_dict[current_id] = line.strip()
    return seq_dict

def add_sequence(attributes, seq_dict):
    for attribute in attributes.split(";"):
        if attribute.startswith("ID="):
            seq_id = attribute.split("=")[1]
            return seq_dict.get(seq_id, ".")
    return "."

def main():
    args = parse_args()

    # Load the GFF3 file
    gff3_df = pd.read_csv(args.gff3, sep="\t", header=None, comment='#')
    gff3_df.columns = ["Chromosome", "source", "type", "Start_pos", "End_pos", "score", "Strand", "phase", "attributes"]

    # Load the sequences
    pep_dict = load_sequences(args.pep)
    cds_dict = load_sequences(args.cds)

    # Add sequences to the GFF3 DataFrame
    gff3_df['Protein_sequence'] = gff3_df['attributes'].apply(lambda x: add_sequence(x, pep_dict))
    gff3_df['CDS_sequence'] = gff3_df['attributes'].apply(lambda x: add_sequence(x, cds_dict))

    # Calculate protein and CDS sequence lengths
    gff3_df['Protein_seq_length'] = gff3_df['Protein_sequence'].apply(lambda seq: len(seq.replace('*', '')))
    gff3_df['CDS_seq_length'] = gff3_df['CDS_sequence'].apply(lambda seq: len(seq))

    # Extract specific fields from the 'attributes' column
    gff3_df['Gene_description'] = gff3_df['attributes'].str.extract(r'product=([^;]+)')
    gff3_df['Ontology_term'] = gff3_df['attributes'].str.extract(r'Ontology_term=([^;]+)')
    gff3_df['Uniprot_id'] = gff3_df['attributes'].str.extract(r'Uniprot_id=([^;]+)')
    gff3_df['Gene_name'] = gff3_df['attributes'].str.extract(r'Parent=([^;]+)')
    gff3_df['Transcript_name'] = gff3_df['attributes'].str.extract(r'ID=([^;]+)')

    # Determine if the gene is before or after the breakpoint
    gff3_df['Pos_relative_to_breakpoint'] = gff3_df.apply(lambda row: 'Before Breakpoint' if row['End_pos'] < args.breakpoint else ('After Breakpoint' if row['Start_pos'] > args.breakpoint else 'Spans Breakpoint'), axis=1)

    # Reorder and print desired columns
    result_df = gff3_df[['Chromosome', 'Start_pos', 'End_pos', 'Strand', 'Pos_relative_to_breakpoint', 'Gene_name', 'Transcript_name', 'Protein_seq_length', 'CDS_seq_length', 'Gene_description', 'Ontology_term', 'Uniprot_id', 'Protein_sequence', 'CDS_sequence', 'attributes']]

    # Save the result to a CSV file
    output_file = args.gff3.replace('.gff3', '_processed.csv')
    result_df.to_csv(output_file, index=False)
    print(f"Processed data saved to {output_file}")

if __name__ == "__main__":
    main()
```

#### Running python script with breakpoint infromation
```bash
python3 get_info_in_csv.py -g aumb_event_1a.gff3 -p aumb_event_1a_mRNA.pep -c aumb_event_1a_mRNA.cds -b 595170863
python3 get_info_in_csv.py -g aumb_event_1b.gff3 -p aumb_event_1b_mRNA.pep -c aumb_event_1b_mRNA.cds -b 29754356
python3 get_info_in_csv.py -g aumb_event_2.gff3 -p aumb_event_2_mRNA.pep -c aumb_event_2_mRNA.cds -b 144777208
python3 get_info_in_csv.py -g aumb_event_3.gff3 -p aumb_event_3_mRNA.pep -c aumb_event_3_mRNA.cds -b 246559785
python3 get_info_in_csv.py -g aumb_event_4.gff3 -p aumb_event_4_mRNA.pep -c aumb_event_4_mRNA.cds -b 267209752
python3 get_info_in_csv.py -g aumb_event_5.gff3 -p aumb_event_5_mRNA.pep -c aumb_event_5_mRNA.cds -b 89652159
python3 get_info_in_csv.py -g aumb_event_6.gff3 -p aumb_event_6_mRNA.pep -c aumb_event_6_mRNA.cds -b 547358291
python3 get_info_in_csv.py -g aumb_event_7.gff3 -p aumb_event_7_mRNA.pep -c aumb_event_7_mRNA.cds -b 59841477
python3 get_info_in_csv.py -g aumb_event_8.gff3 -p aumb_event_8_mRNA.pep -c aumb_event_8_mRNA.cds -b 439590183
python3 get_info_in_csv.py -g aumb_event_9.gff3 -p aumb_event_9_mRNA.pep -c aumb_event_9_mRNA.cds -b 589079447
```

### Combining all csv files in one excel file as sheets
```python3
import pandas as pd
import os

# List of CSV files
csv_files = [
    'aumb_event_1a_processed.csv',
    'aumb_event_1b_processed.csv',
    'aumb_event_2_processed.csv',
    'aumb_event_3_processed.csv',
    'aumb_event_4_processed.csv',
    'aumb_event_5_processed.csv',
    'aumb_event_6_processed.csv',
    'aumb_event_7_processed.csv',
    'aumb_event_8_processed.csv',
    'aumb_event_9_processed.csv'
]

# Directory containing the CSV files
directory = '/mmfs1/scratch/jatinder.singh/05_Aumb_CRs/01_Aumb_breaking_point_genes'

# Create a new Excel writer object and write all CSV files to it
with pd.ExcelWriter('combined_events.xlsx', engine='openpyxl') as excel_writer:
    for csv_file in csv_files:
        file_path = os.path.join(directory, csv_file)
        sheet_name = os.path.splitext(csv_file)[0]  # Use file name without extension as sheet name
        df = pd.read_csv(file_path)
        df.to_excel(excel_writer, sheet_name=sheet_name, index=False)
```
