# To-do ground-truth

## Determining plasmid sequences

**Inputs:**

* Assembled contigs (FASTA)
* List of Plasmid GenBank IDs

**Outputs:**

* Plasmid contig list file
* Non-plasmid contig list file

**Tasks:**

* [x] Download all the plasmid sequences and merge them (merged FASTA)
* [x] Map the contigs (FASTA) to the plasmid references (merged FASTA)
* [x] For each contig, for each plasmid it maps to, filter according to coverage of the contig on the plasmid
* [x] Output plasmid and non-plasmid contig lists to files

## Plasmid contig list file format

`plasmid_contigs.tsv` (tab-separated values)

```html
<plasmid_id>  <contig_id>  <plasmid_length>  <contig_length>  <contig_coverage>
...
```

## Non-plasmid contig list file format

`non_plasmid_contigs.tsv` (tab separated values)

```html
<contig_id>  <contig_length>
...
```
