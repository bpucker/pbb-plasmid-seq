# pbb-plasmid-seq
This repository describes the ONT-based plasmid sequencing workflow conducted by the [Plant Biotechnology and Bioinformatics](https://www.pbb.uni-bonn.de) group at Uni Bonn.

Plasmids are pooled at equimolar amounts based on NanoDrop measurements (estimations). The library prepration is done with the rapid library prep kit from ONT. Sequencing is conducted on R10 flow cells on a PromethION after completion of plant genome sequencing projects. The remaining pores offer sufficient capacity for the sequencing of dozens of plasmids. Basecalling is performed with dorado.

Expected sequences of all plasmids are collected in a multiple FASTA file. Headers in this file will be cleaned from illegal characters as part of the data analysis process.

## Installation
sudo apt install minimap2
sudo apt install samtools
sudo apt install seqkit
sudo apt install seqtk
sudo apt install miniasm
sudo apt install racon


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
