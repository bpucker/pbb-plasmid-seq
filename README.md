# pbb-plasmid-seq
## Background
This repository describes the ONT-based plasmid sequencing workflow conducted by the [Plant Biotechnology and Bioinformatics](https://www.pbb.uni-bonn.de) group at Uni Bonn.

Plasmids are pooled at equimolar amounts based on NanoDrop measurements (estimations). The library preparation is done with the rapid library prep kit from ONT. Sequencing is conducted on R10 flow cells on a PromethION after completion of plant genome sequencing projects. The remaining pores offer sufficient capacity for the sequencing of dozens of plasmids. Basecalling is performed with dorado. Expected sequences of all plasmids are collected in a multiple FASTA file. Headers in this file will be cleaned from illegal characters as part of the data analysis process. A read mapping with [minimap2](https://github.com/lh3/minimap2) against all expected sequences is conducted. The resulting mapping is split per plasmid reference using [samtools](https://www.htslib.org/). A variant calling per plasmid reference is performed with [bcftools](https://samtools.github.io/bcftools/bcftools.html). Reads per plasmid are extracted with samtools and converted into FASTQ. Due to very high coverage, a subsampling is possible to reduce the amount of reads subjected to the following assembly step. [Miniasm](https://github.com/lh3/miniasm) is used for a de novo assembly per plasmid. Racon is applied to polish the assembled plasmid sequence. All result files are made accessible to submitting persons via cloud transfer.


## Installation

```
sudo apt install minimap2 && \
sudo apt install samtools && \
sudo apt install seqkit && \
sudo apt install seqtk && \
sudo apt install miniasm && \
sudo apt install racon
```

## Cleaning and merging FASTA files
```
python3 clean_plasmid_seq_input.py \
--in ./folder_with_plasmid_FASTAs/ \
--fasta clean_plasmid_sequences.fasta \
--doc documentation.txt
```

`--in` specifies the input folder containing FASTA files. File extensions fa, fas, fasta, FA, FAS, and FASTA are considered.

`--fasta` specifies the output FASTA file with clean sequences.

`--doc` specifies the documentation file. Mapping of original sequence names to cleaned sequence names is stored here.


## Running data analysis

```
python3 /vol/data/pbb_plasmid_seq/pbb_plasmid_seq.py \
--reads /vol/data/pbb_plasmid_seq/RP080_99f090c1.mod.fastq.gz \
--ref /vol/data/pbb_plasmid_seq/20260307_plasmid_expectations.fasta \
--out /vol/data/pbb_plasmid_seq/test1/ \
--tmp /vol/data/pbb_plasmid_seq/tmp1/ \
--threads 30 \
> /vol/data/pbb_plasmid_seq/test1.txt 2>&1 &
```

`--reads` specifies the input FASTQ file.

`--ref` specifies the reference FASTA file.

`--out` specifies the output folder.

`--tmp` specifies a temporary output folder.

`--threads` specifies the number of threads to use. Default: 4.


## How to cite? Reference
This repository.
