# Nextflow tuning pipeline

## Seed parameters

**Inputs:**

* Plasmid BioSamples IDs (see YAML file format below)

**Outputs:**

* Seed thresholds YAML file

**Tasks:**

* [x] For each entry in the plasmid BioSamples IDs:
  * [x] Download FASTQ
  * [x] Assemble FASTQ
  * [x] Get plasmid and non-plasmid file contig with `pangebin sub ground-truth` command
  * [x] Get contig gene densities with `pangebin sub map contig` command
* [x] Create input file for `pangebin sub seed thresholds` command
* [x] Run `pangebin sub seed thresholds` command

### Plasmid BioSamples IDs YAML format

Example:

```yaml
BioSample_ID: SAMN01823701
SRA_ID: SRS452288
Plasmid_GenBank_IDs:
- CP006054
---
BioSample_ID: SAMN01832085
SRA_ID: SRS382234
Plasmid_GenBank_IDs:
- CP004026
- CP004028
# ...
```

### Example

1. Init the BioSample seed thresholds dataset file

   ```sh
   seed_outdir="test/seed"
   mkdir $seed_outdir
   seed_contig_thresholds_dataset="$seed_outdir/thresholds_dataset.tsv"
   touch $seed_contig_thresholds_dataset
   ```

2. Loop on BioSample data in the dataset:
   1. Example

      ```sh
      contigs_fasta="test/SAMN16357463/data/std_asm_graph/unicycler.fasta"
      plasmid_genbank_ids="CP069995 CP004028"
      ground_truth_dir="test/SAMN16357463/result/ground_truth"
      email_address="name.surname@domain.ext"

      gene_fasta="src/pangebin/database/genes.fasta"
      gene_mapping_on_contigs_sam="$ground_truth_dir/gene_mapping_on_contigs.sam"
      filter_config_yaml="test/config/filter_config.yaml"
      filtered_sam="$ground_truth_dir/gene_mapping_on_contigs.filtered.sam"

      contig_gene_densities="$ground_truth_dir/contig_gene_densities.tsv"
      ```

   2. Run `pangebin sub ground-truth`

      ```sh
      pangebin sub ground-truth $contigs_fasta $plasmid_genbank_ids $ground_truth_dir --email $email_address
      ```

   3. Obtain gene density on the contigs

      1. Map the gene on the contigs

         ```sh
         pangebin utils map blast $gene_fasta $contigs_fasta $gene_mapping_on_contigs_sam
         ```

      2. Filter the gene mappings

         ```sh
         pangebin utils map filter $gene_mapping_on_contigs_sam $filtered_sam --query-fasta $gene_fasta --config $filter_config_yaml
         ```

      3. Obtain the gene density on the fragments

         ```sh
         pangebin sub gd fasta $contigs_fasta $filtered_sam $contig_gene_densities
         ```

   4. Add BioSample line to seed thresholds dataset file

      ```sh
      echo -e "$ground_truth_dir/plasmid_contigs.tsv\t$ground_truth_dir/non_plasmid_contigs.tsv\t$contig_gene_densities" >> $seed_contig_thresholds_dataset
      ```

3. Run `pangebin sub seed thresholds`

   ```sh
   pangebin sub seed thresholds $seed_contig_thresholds_dataset -o $seed_outdir
   ```
