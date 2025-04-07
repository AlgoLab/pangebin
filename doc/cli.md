# `pangebin`

**Usage**:

```console
pangebin [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `run`: Run the main PangeBin pipeline.
* `seed-thresholds`: Obtain the seed threshold pairs from...
* `configs`: Write default configuration files for the...
* `asm-pbf`: PangeBin asm-pbf pipeline
* `sub`: Subcommands
* `utils`: Utility commands

## `pangebin run`

Run the main PangeBin pipeline.

**Usage**:

```console
pangebin run [OPTIONS]
```

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

## `pangebin seed-thresholds`

Obtain the seed threshold pairs from paired Illumina BioSamples.

**Usage**:

```console
pangebin seed-thresholds [OPTIONS] DATATEST_YAML GENE_FASTA
```

**Arguments**:

* `DATATEST_YAML`: Paired Illumina datatest YAML file  [required]
* `GENE_FASTA`: Genes FASTA file  [required]

**Options**:

* `--config PATH`: The configuration file path  [default: /media/profchep/linux_work/PROJETS/Pangenome/pangebin/src/pangebin/pipeline/seed_thresholds/pipeline_seed_thresholds_cfg.yaml]
* `--outdir PATH`: Output directory  [default: seed_thresholds]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

## `pangebin configs`

Write default configuration files for the pipelines

**Usage**:

```console
pangebin configs [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `seed-thresholds`: Write the configuration file for the seed...
* `asm-pbf`: Write default configuration files for the...

### `pangebin configs seed-thresholds`

Write the configuration file for the seed thresholds pipeline.

**Usage**:

```console
pangebin configs seed-thresholds [OPTIONS] [CONFIG_PATH]
```

**Arguments**:

* `[CONFIG_PATH]`: Output directory  [default: pipeline_seed_thresholds_cfg.yaml]

**Options**:

* `--help`: Show this message and exit.

### `pangebin configs asm-pbf`

Write default configuration files for the asm-pbf pipelines

**Usage**:

```console
pangebin configs asm-pbf [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `decomp`: Write the configuration files for the...
* `binlab`: Write the configuration files for the...
* `once`: Write the configuration files for the once...

#### `pangebin configs asm-pbf decomp`

Write the configuration files for the decomp approach.

**Usage**:

```console
pangebin configs asm-pbf decomp [OPTIONS] OUTPUT_DIRECTORY
```

**Arguments**:

* `OUTPUT_DIRECTORY`: Output directory  [required]

**Options**:

* `--help`: Show this message and exit.

#### `pangebin configs asm-pbf binlab`

Write the configuration files for the binlab approach.

**Usage**:

```console
pangebin configs asm-pbf binlab [OPTIONS] OUTPUT_DIRECTORY
```

**Arguments**:

* `OUTPUT_DIRECTORY`: Output directory  [required]

**Options**:

* `--help`: Show this message and exit.

#### `pangebin configs asm-pbf once`

Write the configuration files for the once approach.

**Usage**:

```console
pangebin configs asm-pbf once [OPTIONS] OUTPUT_DIRECTORY
```

**Arguments**:

* `OUTPUT_DIRECTORY`: Output directory  [required]

**Options**:

* `--help`: Show this message and exit.

## `pangebin asm-pbf`

PangeBin asm-pbf pipeline

**Usage**:

```console
pangebin asm-pbf [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `decomp`: Run Pangebin decomp approach.
* `binlab`: Run Pangebin binlab approach.
* `once`: Run Pangebin once approach.

### `pangebin asm-pbf decomp`

Run Pangebin decomp approach.

**Usage**:

```console
pangebin asm-pbf decomp [OPTIONS] ASSEMBLY_GFA SEED_CONTIGS_TSV CONTIG_PLASMIDNESS_TSV
```

**Arguments**:

* `ASSEMBLY_GFA`: Assembly GFA file  [required]
* `SEED_CONTIGS_TSV`: TSV file with the seed contigs  [required]
* `CONTIG_PLASMIDNESS_TSV`: TSV file with the contigs and their plasmidness scores  [required]

**Options**:

* `--is-skesa / --no-is-skesa`: SKESA GFA file flag  [default: no-is-skesa]
* `--gc-scores PATH`: TSV file containing the GC scores
* `--bin-cfg PATH`: The general binning configuration YAML path
* `--decomp-cfg PATH`: The decomp approach configuration YAML path
* `--gurobi-cfg PATH`: The Gurobi configuration YAML path
* `--outdir PATH`: Output directory  [default: plasbin]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin asm-pbf binlab`

Run Pangebin binlab approach.

**Usage**:

```console
pangebin asm-pbf binlab [OPTIONS] ASSEMBLY_GFA SEED_CONTIGS_TSV CONTIG_PLASMIDNESS_TSV
```

**Arguments**:

* `ASSEMBLY_GFA`: Assembly GFA file  [required]
* `SEED_CONTIGS_TSV`: TSV file with the seed contigs  [required]
* `CONTIG_PLASMIDNESS_TSV`: TSV file with the contigs and their plasmidness scores  [required]

**Options**:

* `--is-skesa / --no-is-skesa`: SKESA GFA file flag  [default: no-is-skesa]
* `--gc-scores PATH`: TSV file containing the GC scores
* `--bin-cfg PATH`: The general binning configuration YAML path
* `--decomp-cfg PATH`: The decomp approach configuration YAML path
* `--gurobi-cfg PATH`: The Gurobi configuration YAML path
* `--outdir PATH`: Output directory  [default: plasbin]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin asm-pbf once`

Run Pangebin once approach.

**Usage**:

```console
pangebin asm-pbf once [OPTIONS] ASSEMBLY_GFA SEED_CONTIGS_TSV CONTIG_PLASMIDNESS_TSV
```

**Arguments**:

* `ASSEMBLY_GFA`: Assembly GFA file  [required]
* `SEED_CONTIGS_TSV`: TSV file with the seed contigs  [required]
* `CONTIG_PLASMIDNESS_TSV`: TSV file with the contigs and their plasmidness scores  [required]

**Options**:

* `--is-skesa / --no-is-skesa`: SKESA GFA file flag  [default: no-is-skesa]
* `--gc-scores PATH`: TSV file containing the GC scores
* `--bin-cfg PATH`: The general binning configuration YAML path
* `--gurobi-cfg PATH`: The Gurobi configuration YAML path
* `--outdir PATH`: Output directory  [default: plasbin]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

## `pangebin sub`

Subcommands

**Usage**:

```console
pangebin sub [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `std-asm-graph`: Standardize GFA assembly graphs.
* `pangenome`: Produce a pangenome using nf-core/pangenome.
* `panassembly`: Produce a pangenome assembly from the...
* `database`: Create plasmid database (biosamples,...
* `ground-truth`: Create ground truth.
* `plasbin`: Binning plasmids on an assembly or a...
* `gc`: GC content operations.
* `gd`: Gene density operations.
* `seed`: Seed sequences operations.

### `pangebin sub std-asm-graph`

Standardize GFA assembly graphs.

**Usage**:

```console
pangebin sub std-asm-graph [OPTIONS] SKESA_GFA_PATH UNICYCLER_GFA_PATH
```

**Arguments**:

* `SKESA_GFA_PATH`: SKESA GFA assembly graph file  [required]
* `UNICYCLER_GFA_PATH`: Unicycler GFA assembly graph file  [required]

**Options**:

* `--min-contig-length INTEGER`: Minimum contig length threshold  [default: 1]
* `--config PATH`: The configuration file path
* `-o, --outdir PATH`: Output folder  [default: standardize]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin sub pangenome`

Produce a pangenome using nf-core/pangenome.

**Usage**:

```console
pangebin sub pangenome [OPTIONS] MIXED_FASTA
```

**Arguments**:

* `MIXED_FASTA`: Mixed FASTA file (can be gzipped or bgzipped or not)  [required]

**Options**:

* `--release TEXT`: nf-core/pangenome release  [default: 1.1.2]
* `--profile [test|docker|singularity|podman|shifter|charliecloud|apptainer|conda]`: Profile environment  [default: podman]
* `--resume / --no-resume`: Resume workflow  [default: resume]
* `--supplementary-nfcore-pangenome-config-path PATH`: Overiding base nf-core/pangenome configuration file
* `--config PATH`: The configuration file path
* `-o, --outdir PATH`: Output folder  [default: pangenome]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin sub panassembly`

Produce a pangenome assembly from the original assemblers and the pangenome.

**Usage**:

```console
pangebin sub panassembly [OPTIONS] PANGENOME_GFA_PATH STANDARDIZED_SKESA_GFA_PATH STANDARDIZED_UNICYCLER_GFA_PATH
```

**Arguments**:

* `PANGENOME_GFA_PATH`: Pangenome GFA file  [required]
* `STANDARDIZED_SKESA_GFA_PATH`: SKESA GFA standardized graph  [required]
* `STANDARDIZED_UNICYCLER_GFA_PATH`: Unicycler GFA standardized graph  [required]

**Options**:

* `-o, --outdir PATH`: Output folder  [default: panassembly]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin sub database`

Create plasmid database (biosamples, plasmids and paired Illumina SRA ids).

**Usage**:

```console
pangebin sub database [OPTIONS] ACCESSIONS_FILE
```

**Arguments**:

* `ACCESSIONS_FILE`: File with Entrez accessions  [required]

**Options**:

* `-o, --outdir PATH`: Output directory  [default: database]
* `--email TEXT`: Email address to fetch NCBI database
* `--tool TEXT`: Tool name  [default: biopython]
* `--api-key TEXT`: API key
* `--max-tries INTEGER`: Max tries  [default: 3]
* `--sleep-between-tries INTEGER`: Sleep between tries  [default: 15]
* `--entrez-cfg PATH`: The configuration file path
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin sub ground-truth`

Create ground truth.

**Usage**:

```console
pangebin sub ground-truth [OPTIONS] CONTIGS_FASTA_FILE PLASMID_GENBANK_IDS...
```

**Arguments**:

* `CONTIGS_FASTA_FILE`: Contigs FASTA file  [required]
* `PLASMID_GENBANK_IDS...`: GenBank IDs of plasmid sequences  [required]

**Options**:

* `--min-pident FLOAT`: Minimum percent identity threshold (between 0 and 100)  [default: 95]
* `--min-contig-coverage FLOAT`: Minimum contig coverage threshold (between 0 and 1)  [default: 0.95]
* `--config PATH`: The configuration file path
* `--email TEXT`: Email address to fetch NCBI database
* `--tool TEXT`: Tool name  [default: biopython]
* `--api-key TEXT`: API key
* `--max-tries INTEGER`: Max tries  [default: 3]
* `--sleep-between-tries INTEGER`: Sleep between tries  [default: 15]
* `--entrez-cfg PATH`: The configuration file path
* `-o, --outdir PATH`: Output folder  [default: ground_truth]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin sub plasbin`

Binning plasmids on an assembly or a pan-assembly graph.

**Usage**:

```console
pangebin sub plasbin [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `decomp`: Hiearchical decomposition binning method
* `binlab`: Binning-and-labelling method
* `once`: Once method

#### `pangebin sub plasbin decomp`

Hiearchical decomposition binning method

**Usage**:

```console
pangebin sub plasbin decomp [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `asm`: PlasBin an assembly graph with the decomp...
* `panasm`: PlasBin a pan-assembly graph with the...

##### `pangebin sub plasbin decomp asm`

PlasBin an assembly graph with the decomp approach.

**Usage**:

```console
pangebin sub plasbin decomp asm [OPTIONS] ASSEMBLY_GFA SEED_CONTIGS_TSV CONTIG_GC_SCORES_TSV CONTIG_PLASMIDNESS_TSV
```

**Arguments**:

* `ASSEMBLY_GFA`: Assembly GFA file  [required]
* `SEED_CONTIGS_TSV`: TSV file with the seed contigs  [required]
* `CONTIG_GC_SCORES_TSV`: TSV file with the contigs and their GC scores  [required]
* `CONTIG_PLASMIDNESS_TSV`: TSV file with the contigs and their plasmidness scores  [required]

**Options**:

* `--sink-arcs-domain [all|seeds]`: Sink-arcs domain  [default: all]
* `--min-flow FLOAT`: Minimum flow  [default: 0.0001]
* `--min-cumulative-len INTEGER`: Minimum cumulative length  [default: 1000]
* `--circular / --no-circular`: The flow is circular  [default: no-circular]
* `--obj-fun-domain [all|seeds]`: Objective function domain  [default: all]
* `--bin-cfg PATH`: The configuration file path
* `--gamma-mcf FLOAT`: Gamma MCF coefficient  [default: 0.9]
* `--gamma-mgc FLOAT`: Gamma MGC coefficient  [default: 0.9]
* `--decomp-cfg PATH`: The configuration file path
* `--mip-gap FLOAT`: MIP gap
* `--time-limit FLOAT`: Time limit
* `--threads INTEGER`: Threads
* `--gurobi-cfg PATH`: The configuration file path
* `--outdir PATH`: Output directory  [default: plasbin]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

##### `pangebin sub plasbin decomp panasm`

PlasBin a pan-assembly graph with the decomp approach.

**Usage**:

```console
pangebin sub plasbin decomp panasm [OPTIONS] PANASSEMBLY_GFA SEED_FRAGMENTS_TSV FRAGMENT_GC_SCORES_TSV FRAGMENT_PLASMIDNESS_TSV
```

**Arguments**:

* `PANASSEMBLY_GFA`: Pan-assembly GFA file  [required]
* `SEED_FRAGMENTS_TSV`: TSV file with the seed fragments  [required]
* `FRAGMENT_GC_SCORES_TSV`: TSV file with the fragments and their GC scores  [required]
* `FRAGMENT_PLASMIDNESS_TSV`: TSV file with the fragments and their plasmidness scores  [required]

**Options**:

* `--sink-arcs-domain [all|seeds]`: Sink-arcs domain  [default: all]
* `--min-flow FLOAT`: Minimum flow  [default: 0.0001]
* `--min-cumulative-len INTEGER`: Minimum cumulative length  [default: 1000]
* `--circular / --no-circular`: The flow is circular  [default: no-circular]
* `--obj-fun-domain [all|seeds]`: Objective function domain  [default: all]
* `--bin-cfg PATH`: The configuration file path
* `--gamma-mcf FLOAT`: Gamma MCF coefficient  [default: 0.9]
* `--gamma-mgc FLOAT`: Gamma MGC coefficient  [default: 0.9]
* `--decomp-cfg PATH`: The configuration file path
* `--mip-gap FLOAT`: MIP gap
* `--time-limit FLOAT`: Time limit
* `--threads INTEGER`: Threads
* `--gurobi-cfg PATH`: The configuration file path
* `--outdir PATH`: Output directory  [default: plasbin]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin sub plasbin binlab`

Binning-and-labelling method

**Usage**:

```console
pangebin sub plasbin binlab [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `asm`: PlasBin an assembly graph with the binlab...
* `panasm`: PlasBin a pan-assembly graph with the...

##### `pangebin sub plasbin binlab asm`

PlasBin an assembly graph with the binlab approach.

**Usage**:

```console
pangebin sub plasbin binlab asm [OPTIONS] ASSEMBLY_GFA SEED_CONTIGS_TSV CONTIG_GC_SCORES_TSV CONTIG_PLASMIDNESS_TSV
```

**Arguments**:

* `ASSEMBLY_GFA`: Assembly GFA file  [required]
* `SEED_CONTIGS_TSV`: TSV file with the seed contigs  [required]
* `CONTIG_GC_SCORES_TSV`: TSV file with the contigs and their GC scores  [required]
* `CONTIG_PLASMIDNESS_TSV`: TSV file with the contigs and their plasmidness scores  [required]

**Options**:

* `--sink-arcs-domain [all|seeds]`: Sink-arcs domain  [default: all]
* `--min-flow FLOAT`: Minimum flow  [default: 0.0001]
* `--min-cumulative-len INTEGER`: Minimum cumulative length  [default: 1000]
* `--circular / --no-circular`: The flow is circular  [default: no-circular]
* `--obj-fun-domain [all|seeds]`: Objective function domain  [default: all]
* `--bin-cfg PATH`: The configuration file path
* `--gamma-mbs FLOAT`: Gamma MBS coefficient  [default: 0.9]
* `--binlab-cfg PATH`: The configuration file path
* `--mip-gap FLOAT`: MIP gap
* `--time-limit FLOAT`: Time limit
* `--threads INTEGER`: Threads
* `--gurobi-cfg PATH`: The configuration file path
* `--outdir PATH`: Output directory  [default: plasbin]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

##### `pangebin sub plasbin binlab panasm`

PlasBin a pan-assembly graph with the binlab approach.

**Usage**:

```console
pangebin sub plasbin binlab panasm [OPTIONS] PANASSEMBLY_GFA SEED_FRAGMENTS_TSV FRAGMENT_GC_SCORES_TSV FRAGMENT_PLASMIDNESS_TSV
```

**Arguments**:

* `PANASSEMBLY_GFA`: Pan-assembly GFA file  [required]
* `SEED_FRAGMENTS_TSV`: TSV file with the seed fragments  [required]
* `FRAGMENT_GC_SCORES_TSV`: TSV file with the fragments and their GC scores  [required]
* `FRAGMENT_PLASMIDNESS_TSV`: TSV file with the fragments and their plasmidness scores  [required]

**Options**:

* `--sink-arcs-domain [all|seeds]`: Sink-arcs domain  [default: all]
* `--min-flow FLOAT`: Minimum flow  [default: 0.0001]
* `--min-cumulative-len INTEGER`: Minimum cumulative length  [default: 1000]
* `--circular / --no-circular`: The flow is circular  [default: no-circular]
* `--obj-fun-domain [all|seeds]`: Objective function domain  [default: all]
* `--bin-cfg PATH`: The configuration file path
* `--gamma-mbs FLOAT`: Gamma MBS coefficient  [default: 0.9]
* `--binlab-cfg PATH`: The configuration file path
* `--mip-gap FLOAT`: MIP gap
* `--time-limit FLOAT`: Time limit
* `--threads INTEGER`: Threads
* `--gurobi-cfg PATH`: The configuration file path
* `--outdir PATH`: Output directory  [default: plasbin]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin sub plasbin once`

Once method

**Usage**:

```console
pangebin sub plasbin once [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `asm`: PlasBin an assembly graph with the once...
* `panasm`: PlasBin a pan-assembly graph with the once...

##### `pangebin sub plasbin once asm`

PlasBin an assembly graph with the once approach.

**Usage**:

```console
pangebin sub plasbin once asm [OPTIONS] ASSEMBLY_GFA SEED_CONTIGS_TSV CONTIG_GC_SCORES_TSV CONTIG_PLASMIDNESS_TSV
```

**Arguments**:

* `ASSEMBLY_GFA`: Assembly GFA file  [required]
* `SEED_CONTIGS_TSV`: TSV file with the seed contigs  [required]
* `CONTIG_GC_SCORES_TSV`: TSV file with the contigs and their GC scores  [required]
* `CONTIG_PLASMIDNESS_TSV`: TSV file with the contigs and their plasmidness scores  [required]

**Options**:

* `--sink-arcs-domain [all|seeds]`: Sink-arcs domain  [default: all]
* `--min-flow FLOAT`: Minimum flow  [default: 0.0001]
* `--min-cumulative-len INTEGER`: Minimum cumulative length  [default: 1000]
* `--circular / --no-circular`: The flow is circular  [default: no-circular]
* `--obj-fun-domain [all|seeds]`: Objective function domain  [default: all]
* `--bin-cfg PATH`: The configuration file path
* `--mip-gap FLOAT`: MIP gap
* `--time-limit FLOAT`: Time limit
* `--threads INTEGER`: Threads
* `--gurobi-cfg PATH`: The configuration file path
* `--outdir PATH`: Output directory  [default: plasbin]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

##### `pangebin sub plasbin once panasm`

PlasBin a pan-assembly graph with the once approach.

**Usage**:

```console
pangebin sub plasbin once panasm [OPTIONS] PANASSEMBLY_GFA SEED_FRAGMENTS_TSV FRAGMENT_GC_SCORES_TSV FRAGMENT_PLASMIDNESS_TSV
```

**Arguments**:

* `PANASSEMBLY_GFA`: Pan-assembly GFA file  [required]
* `SEED_FRAGMENTS_TSV`: TSV file with the seed fragments  [required]
* `FRAGMENT_GC_SCORES_TSV`: TSV file with the fragments and their GC scores  [required]
* `FRAGMENT_PLASMIDNESS_TSV`: TSV file with the fragments and their plasmidness scores  [required]

**Options**:

* `--sink-arcs-domain [all|seeds]`: Sink-arcs domain  [default: all]
* `--min-flow FLOAT`: Minimum flow  [default: 0.0001]
* `--min-cumulative-len INTEGER`: Minimum cumulative length  [default: 1000]
* `--circular / --no-circular`: The flow is circular  [default: no-circular]
* `--obj-fun-domain [all|seeds]`: Objective function domain  [default: all]
* `--bin-cfg PATH`: The configuration file path
* `--mip-gap FLOAT`: MIP gap
* `--time-limit FLOAT`: Time limit
* `--threads INTEGER`: Threads
* `--gurobi-cfg PATH`: The configuration file path
* `--outdir PATH`: Output directory  [default: plasbin]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin sub gc`

GC content operations.

**Usage**:

```console
pangebin sub gc [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `from-gfa`: Compute GC scores from GFA assembly graph.

#### `pangebin sub gc from-gfa`

Compute GC scores from GFA assembly graph.

**Usage**:

```console
pangebin sub gc from-gfa [OPTIONS] GFA_FILE GC_SCORES_TSV
```

**Arguments**:

* `GFA_FILE`: GFA assembly graph file  [required]
* `GC_SCORES_TSV`: Output TSV file with the sequences and their GC scores  [required]

**Options**:

* `--gc-content-interval-tsv PATH`: GC content intervals TSV file  [default: /media/profchep/linux_work/PROJETS/Pangenome/pangebin/src/pangebin/gc_content/default_gc_content_intervals.tsv]
* `--pseudo-count INTEGER`: Pseudo count for GC probabilities  [default: 10]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin sub gd`

Gene density operations.

**Usage**:

```console
pangebin sub gd [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `fasta`: Compute the gene densities from the...
* `frag`: Compute fragment gene densities from the...

#### `pangebin sub gd fasta`

Compute the gene densities from the mapping of genes against the sequences.

**Usage**:

```console
pangebin sub gd fasta [OPTIONS] FASTA_FILE GENE_MAPPING_SAM OUTPUT_FILE
```

**Arguments**:

* `FASTA_FILE`: FASTA file  [required]
* `GENE_MAPPING_SAM`: Gene mapping SAM file  [required]
* `OUTPUT_FILE`: Output file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin sub gd frag`

Compute fragment gene densities from the mapping of genes against the contigs.

**Usage**:

```console
pangebin sub gd frag [OPTIONS] PANASSEMBLY_GFA GENE_MAPPING_TO_CONTIGS_SAM OUTPUT_FILE
```

**Arguments**:

* `PANASSEMBLY_GFA`: Pan-assembly GFA file  [required]
* `GENE_MAPPING_TO_CONTIGS_SAM`: Gene mapping to contigs SAM file  [required]
* `OUTPUT_FILE`: Output file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin sub seed`

Seed sequences operations.

**Usage**:

```console
pangebin sub seed [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `pos-gd`: Extract seed sequences with positive gene...
* `ctgs-to-frags`: Extract seed fragments from seed contigs.
* `thresholds`: Seed thresholds.

#### `pangebin sub seed pos-gd`

Extract seed sequences with positive gene density.

**Usage**:

```console
pangebin sub seed pos-gd [OPTIONS] GENE_DENSITY_FILE OUTPUT_SEED_TSV
```

**Arguments**:

* `GENE_DENSITY_FILE`: Gene density file  [required]
* `OUTPUT_SEED_TSV`: Seed sequence output file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin sub seed ctgs-to-frags`

Extract seed fragments from seed contigs.

**Usage**:

```console
pangebin sub seed ctgs-to-frags [OPTIONS] SKESA_SEEDS_TSV UNICYCLER_SEEDS_TSV PANASSEMBLY_GFA OUTPUT_SEED_TSV
```

**Arguments**:

* `SKESA_SEEDS_TSV`: SKESA seeds TSV file  [required]
* `UNICYCLER_SEEDS_TSV`: Unicycler seeds TSV file  [required]
* `PANASSEMBLY_GFA`: Pan-assembly GFA file  [required]
* `OUTPUT_SEED_TSV`: Seed sequence output file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin sub seed thresholds`

Seed thresholds.

**Usage**:

```console
pangebin sub seed thresholds [OPTIONS] DATATEST
```

**Arguments**:

* `DATATEST`: Datatest file  [required]

**Options**:

* `--min-length INTEGER`: Minimum length threshold  [default: 50]
* `--max-length INTEGER`: Maximum length threshold  [default: 5000]
* `--step-length INTEGER`: Step length threshold  [default: 50]
* `--min-gene-density FLOAT`: Minimum gene density threshold  [default: 0.01]
* `--max-gene-density FLOAT`: Maximum gene density threshold  [default: 1.0]
* `--step-gene-density FLOAT`: Step gene density threshold  [default: 0.01]
* `--config PATH`: The configuration file path
* `-o, --outdir PATH`: Output folder  [default: seed_thresholds]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

## `pangebin utils`

Utility commands

**Usage**:

```console
pangebin utils [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `gfa`: GFA operations.
* `map`: Mapping operations.
* `pbf-comp`: PlasBin-flow conversion.

### `pangebin utils gfa`

GFA operations.

**Usage**:

```console
pangebin utils gfa [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `check-skesa`: Check a Skeza GFA file.
* `fix-skesa`: Fix a Skeza GFA file.
* `to-fasta`: Convert GFA to FASTA.
* `is-standardized`: Check if a GFA is standardized.
* `stats`: Print GFA stats.

#### `pangebin utils gfa check-skesa`

Check a Skeza GFA file.

**Usage**:

```console
pangebin utils gfa check-skesa [OPTIONS] IN_GFA
```

**Arguments**:

* `IN_GFA`: Input GFA file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin utils gfa fix-skesa`

Fix a Skeza GFA file.

**Usage**:

```console
pangebin utils gfa fix-skesa [OPTIONS] IN_GFA [OUT_GFA]
```

**Arguments**:

* `IN_GFA`: Input GFA file  [required]
* `[OUT_GFA]`: Output GFA file, must be different from input if provided

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin utils gfa to-fasta`

Convert GFA to FASTA.

**Usage**:

```console
pangebin utils gfa to-fasta [OPTIONS] GFA_PATH
```

**Arguments**:

* `GFA_PATH`: Input GFA file  [required]

**Options**:

* `--attribute-string-separator TEXT`: String separator for attributes  [default:  ]
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin utils gfa is-standardized`

Check if a GFA is standardized.

**Usage**:

```console
pangebin utils gfa is-standardized [OPTIONS] GFA_PATH
```

**Arguments**:

* `GFA_PATH`: Input GFA file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin utils gfa stats`

Print GFA stats.

**Usage**:

```console
pangebin utils gfa stats [OPTIONS] GFA_PATH
```

**Arguments**:

* `GFA_PATH`: GFA file path  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin utils map`

Mapping operations.

**Usage**:

```console
pangebin utils map [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `blast`: Map a query FASTA file to a subject FASTA...
* `filter`: Filter a SAM file.

#### `pangebin utils map blast`

Map a query FASTA file to a subject FASTA file.

**Usage**:

```console
pangebin utils map blast [OPTIONS] QUERY_FASTA_FILE SUBJECT_FASTA_FILE OUT_MAPPING_FILE
```

**Arguments**:

* `QUERY_FASTA_FILE`: Query FASTA file  [required]
* `SUBJECT_FASTA_FILE`: Subject FASTA file  [required]
* `OUT_MAPPING_FILE`: Output file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin utils map filter`

Filter a SAM file.

**Usage**:

```console
pangebin utils map filter [OPTIONS] INPUT_SAM [FILTERED_SAM]
```

**Arguments**:

* `INPUT_SAM`: Input mapping file  [required]
* `[FILTERED_SAM]`: Output file

**Options**:

* `--query-fasta PATH`: Query FASTA file
* `--subject-fasta PATH`: Subject FASTA file
* `--min-length INTEGER`: Minimum length threshold  [default: 1]
* `--min-pident FLOAT`: Minimum percent identity threshold (between 0 and 100)  [default: 0]
* `--min-q-cov FLOAT`: Minimum query coverage threshold (between 0 and 1) (only if query FASTA file is provided)  [default: 0]
* `--min-s-cov FLOAT`: Minimum subject coverage threshold (between 0 and 1) (only if subject FASTA file is provided)  [default: 0]
* `--config PATH`: The configuration file path
* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

### `pangebin utils pbf-comp`

PlasBin-flow conversion.

**Usage**:

```console
pangebin utils pbf-comp [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `plm`: Convert PlasBin-flow plasmidness file to...
* `seeds`: Convert PlasBin-flow seed sequences file...
* `bins`: Convert PangeBin-flow bins into...

#### `pangebin utils pbf-comp plm`

Convert PlasBin-flow plasmidness file to PangeBin plasmidness TSV file.

**Usage**:

```console
pangebin utils pbf-comp plm [OPTIONS] PBF_PLASMIDNESS_FILE PG_PLASMIDNESS_TSV
```

**Arguments**:

* `PBF_PLASMIDNESS_FILE`: PlasBin-flow plasmidness file  [required]
* `PG_PLASMIDNESS_TSV`: Pangebin plasmidness TSV file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin utils pbf-comp seeds`

Convert PlasBin-flow seed sequences file to PangeBin seed sequences TSV file.

**Usage**:

```console
pangebin utils pbf-comp seeds [OPTIONS] PBF_SEEDS_FILE PG_SEEDS_TSV
```

**Arguments**:

* `PBF_SEEDS_FILE`: PlasBin-flow seed sequences file  [required]
* `PG_SEEDS_TSV`: Pangebin seed sequences TSV file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.

#### `pangebin utils pbf-comp bins`

Convert PangeBin-flow bins into PlasBin-flow bin info TSV file.

**Usage**:

```console
pangebin utils pbf-comp bins [OPTIONS] PG_BINS_PARENT_DIR PBF_BIN_INFO_TSV
```

**Arguments**:

* `PG_BINS_PARENT_DIR`: Pangebin bins parent directory  [required]
* `PBF_BIN_INFO_TSV`: PlasBin-flow bin info TSV file  [required]

**Options**:

* `--debug / --no-debug`: Debug mode  [default: no-debug]
* `--help`: Show this message and exit.
