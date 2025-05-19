# Change log

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> 'M', 'm' and 'p' are respectively corresponding to major, minor and patch

<!-- The order of keywords:
## [Unreleased] - yyyy-mm-dd

### Added

### Changed

### Deprecated

### Removed

### Fixed

### Security
-->

<!-- next-header -->
## [Unreleased] - yyyy-mm-dd

### Added

* GFA `predecessors` and `successors` functions
* GFA subgraph function according to a a lis of segments to keep and a radius on the neighborhood

### Fixed

* Only keep seeds contained in the graph

## [0.2.2] - 2025-05-13

### Changed

* `to_fasta` gfa utils subcommand now takes a fasta file path as argument

### Fixed

* Non py data files are now listed

## [0.2.1] - 2025-04-10

### Fixed

* Some `pangebin` module names shadowed standard python libraries

## [0.2.0] - 2025-04-07

### Added

* `pangebin configs asm-pbf {decomp|binlab|once} OUTPUT_DIR` commands to write defaults configuration files for the `asm-pbf` pipelines
* Add `doc/cli.md` documentation on CLI

## [0.1.0] - 2025-04-06

First release tag of Pangebin

* Pipeline to launch Pangebin on PlasBin-flow inputs (assembly graph)
  * `pangebin asm-pbf {decomp|binlab|once} ARGUMENTS [OPTIONS]`
* Subcommands available through `pangebin sub`
* Utilities available through `pangebin utils`
