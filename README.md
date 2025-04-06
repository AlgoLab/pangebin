# Pangebin

## Setup Python virtual environment

<!-- DOCU condaenv for dev -> change when user's one is ready -->
* [*For dev*] Create the conda environment

  ```sh
  conda env create -n pangebin-dev -f config/condaenv_311-dev.yml
  ```

* [*For dev*] Activate the conda environment

  ```sh
  conda activate pangebin-dev
  ```

## Usage

Example with the dataset `SAMN16357463`:

```sh
dataset_dir="test/SAMN16357463"
```

> **Note:** each script has a command `clean`
>
> ```sh
> ./script.sh clean $dataset_dir
> ```

1. standardize GFA assembly graphs

   ```sh
   ./test/std_asm_graph.sh run $dataset_dir
   ```

2. make pangenome graph with nextflow (make sur you have installed the command for the nextflow profile)

   ```sh
   ./test/pangenome.sh run $dataset_dir
   ```

3. make pan-assembly graph

   ```sh
   ./test/panassembly.sh run $dataset_dir
   ```

4. Obtain the GC probability scores of the fragments

   ```sh
   ./test/gc_prob_scores.sh run $dataset_dir
   ```

5. Obtain gene density on the fragments

   1. Map the gene on the contigs from the two assemblers

      ```sh
      ./test/gene_mapping.sh blast $dataset_dir
      ```

   2. Filter the gene mappings

      ```sh
      ./test/gene_mapping.sh filter $dataset_dir
      ```

   3. Obtain the gene density on the fragments

      ```sh
      ./test/frag_gene_densities.sh run $dataset_dir
      ```

6. Obtain the seed from positive gene densities

   ```sh
   ./test/fragment_seeds.sh run $dataset_dir
   ```

7. Execute PlasBin-Flow modified for pan-assembly

   ```sh
   ./test/plasbin_panasm.sh run $dataset_dir
   ```

## Pangebin-PlasBin-flow conversion

### Use PlasBin-flow inputs to Pangebin

```sh
# pbf : PlasBin-flow
# pg : Pangebin
# Convert plasmidness
pangebin utils pbf-comp plm pbf_plasmidness.tsv pg_plasmidness.tsv
# Convert seeds
pangebin utils pbf-comp seeds pbf_seeds.tsv pg_seeds.tsv
# Recompute GC contents
pangebin sub gc from-gfa graph.gfa pg_gc_scores.tsv
# Run Pangebin on assembly graph
pangebin sub plasbin asm graph.gfa pg_seeds.tsv pg_gc_scores.tsv pg_plasmidness.tsv --outdir pg_outdir
```

### Convert Pangebin outputs to PlasBin-flow output

```sh
pangebin utils pbf-comp bins pg_outdir pbf_bins.tsv
```

## Going further into the details

Understanding GFA tags system:

* [doc/gfa/standardized.md](doc/gfa/standardized.md)
* [doc/gfa/panassembly.md](doc/gfa/panassembly.md)
