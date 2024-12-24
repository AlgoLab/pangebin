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

Example with the sample `SAMN16357463`:

```sh
dataset_dir="test/SAMN16357463"
```

1. preprocess input files

   ```sh
   ./test/preprocess.sh run $dataset_dir
   # to clean:
   # ./test/preprocess.sh clean $dataset_dir
   ```

2. make pangenome graph using `${mixed_fasta}`

   ```sh
   # DOCU pggb ...
   pangenome="test/data/SAMN16357463/out/SAMN16357463.1.pan.gfa"
   cp ${pggb_out} ${pangenome}
   ```

3. make pan-assembly graph

   ```sh
   ./test/panassembly.sh run $dataset_dir
   # to clean:
   # ./test/panassembly.sh clean $dataset_dir
   ```

4. Execute pangebin

   ```sh
   ./test/pangebin.sh run $dataset_dir
   # to clean:
   # ./test/pangebin.sh clean $dataset_dir
   ```
