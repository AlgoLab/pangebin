# Plasmid database

This database is used to obtain seed thresholds and gene density for new pan-assembly fragments.

The seed threshold pairs are obtained by:

1. assembling short illumina paired-end reads
2. aligning contigs into their associated plasmid reference genomes
3. computing contig gene densities

Two files serve as references:

* `illumina_biosamples.yaml` contains a list of biosamples linked to Illumina paired reads SRAs and plasmid accessions
* `non_illumina_biosamples.yaml` contains a list of biosamples linked to plasmid accessions

These files were obtained by fetching BioSample and SRA NCBI databases from a list of plasmid accessions in `accessions.txt`.

The file `genes.fasta` contains all the gene associated to all the plasmid references in `accessions.txt`.
