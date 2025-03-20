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
* [x] For each threshold pair, compute the chromosone-plasmid difference (see paper)
* [x] Choose the threshold pair maximizing the chromosone-plasmid difference
* [x] Output the means and thresholds

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

* [ ] #TODO add other score original seed threshold script gives(?)

## Seed fragments

**Inputs:**

* pan-assembly graph
* Fragment gene densities
* Seed thresholds file

### From seed thresholds

**Tasks:**

<!-- TODO seed fragments from seed thresholds -->

* [ ] Compute mean of fragments lengths and gene densities
* [ ] Get a list of fragments name passing the normalized thresholds
* [ ] Output to file

### From positive gene densities

**Tasks:**

* [x] Return fragments whose gene density is positive

### From gene density distribution

**Tasks:**

* [ ] Return fragments whose gene density is above a percentage times max gene density

### From plasmidness score

**Tasks:**

* [ ] Compute plasmidness score:

  $$\rho_i = \frac{ gd_i }{1 +  \max_{ j \in \mathcal{F} } \{ |j| \} - |i| }$$

* [ ] Return fragments whose plasmidness score is above a percentage times max plasmidness score

### Output file format

(Tabulation separated)

```html
<segment_id>  <length>  <gene_density>  [<plasmidness_score>]
...
```
