# Seed sequences

## Determining seed length and gene density thresholds

**Inputs:**

* File such that each line contains:
  * Plasmid contig list file
  * Non-plasmid contig list file
  * Contig gene densities file
* Ranges of lengths and gene densities to test

**Tasks:**

* [x] Mix all the contig:
  * [x] add length of contigs for both plasmid and non-plasmid contig list files
  * [x] plasmid label if they are contained in the plasmid contig list file
  * [x] non-plasmid label if they are not contained in the plasmid contig list file
  * [x] assign gene density
* [ ] For each threshold pair, compute the chromosone-plasmid difference (see paper)
* [ ] Choose the threshold pair maximizing the chromosone-plasmid difference
* [ ] Output the means and thresholds

### Input file format

`seed_contig_thresholds_dataset.tsv` (tab-separated values)

```html
<path/to/plasmid_contigs.tsv>  <path/to/non_plasmid_contigs.tsv>  <path/to/contig_gene_density.tsv>
  ...
```

### Output file format

`seed_thresholds.yaml`

```yaml
means:
  length: <float>
  gene_density: <float>

best_sp_nps: <int>

thresholds:
  - length: <int>
    gene_density: <float>
  # ...
```

* [ ] #TODO add other score original seed threshold script gives

## Seed fragments

**Inputs:**

* pan-assembly graph
* Fragment gene densities
* Seed thresholds file

**Tasks:**

* [ ] Compute mean of fragments lengths and gene densities
* [ ] Get a list of fragments name passing the normalized thresholds
* [ ] Output to file

### Output file format

(Tabulation separated)

```html
<segment_id>  <length>  <gene_density>
...
```
