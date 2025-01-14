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

1. standardize GFA assembly graphs

   ```sh
   ./test/std_asm_graph.sh run $dataset_dir
   # to clean:
   # ./test/std_asm_graph.sh clean $dataset_dir
   ```

2. make pangenome graph with nextflow (make sur you have installed the command for the nextflow profile)

   ```sh
   ./test/pangenome.sh run $dataset_dir
   # to clean:
   # ./test/pangenome.sh clean $dataset_dir
   ```

3. make pan-assembly graph

   ```sh
   ./test/panassembly.sh run $dataset_dir
   # to clean:
   # ./test/panassembly.sh clean $dataset_dir
   ```

4. Execute PlasBin-Flow modified for pan-assembly

   ```sh
   ./test/plasbin.sh run $dataset_dir
   # to clean:
   # ./test/plasbin.sh clean $dataset_dir
   ```

## Going further into the details

Understanding GFA tags system [`doc/gfa_tags.md`](doc/gfa_tags.md)
