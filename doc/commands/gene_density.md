# Gene density module

## Gene density for fragments

**Input:**

* mapping file genes on contigs

**Tasks:**

* [x] Compute interval of matches on contigs
* [x] Compute interval of matches on fragments from contig intervals
  * [x] Be carefull of the name of the contigs, they were formatted for standartize pggb fasta
* [x] Compute fragment gene density
* [x] Output to file

## Gene density for contigs

**Input:**

* mapping file genes on contigs

**Tasks:**

* [x] Compute interval of matches on contigs
* [x] Compute contig gene density
* [x] Output to file

## File format

(Tabulation separated)

```html
<segment_id>  <gene_density>  <start>:<end>,...
...
```
