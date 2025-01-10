# GFA tags

## Standardized assembly graphs

Standardize the tags of Skesa and Unicycler assembly graphs.

### Header

| ID   | Type       | Values | Description                     |
| ---- | ---------- | ------ | ------------------------------- |
| `Sd` | `A` (char) | `Y\|N` | The graph is (not) standardized |

### Segments

#### Names

* `ske_<k>` for Skesa
* `uni_<k>` for Unicycler

For $k \in \{1...|S|\}$, with $S$ the set of segments.

#### Tags

| ID   | Type        | Values            | Description                  |
| ---- | ----------- | ----------------- | ---------------------------- |
| `LN` | `i` (int)   | $\mathbb{N}_{>0}$ | Length of the sequence       |
| `dp` | `f` (float) | $\mathbb{R}_{>0}$ | Normalized sequence coverage |

### Links

#### Overlap match

`0M`

## Pan-assembly graph

### Header

#### Tags

| ID   | Type         | Values | Description |
| ---- | ------------ | ------ | ----------- |
| `VN` | `Z` (string) | `M.m`  | GFA version |

### Segments

#### Names

For pangenome segments:

* `panske_<k>` for Skesa subcontigs (with $k \in \{1...|S_{ps}|\}$ where $S_{ps}$ is the set of pangenome segments corresponding to Skesa subcontigs)
* `panuni_<k>` for Unicycler subcontigs (with $k \in \{1...|S_{pu}|\}$ where $S_{pu}$ is the set of pangenome segments corresponding to Unicycler subcontigs)
* `panboth_<k>` for Skesa-Unicycler subcontigs (with $k \in \{1...|S_{pb}|\}$ where $S_{pb}$ is the set of pangenome segments corresponding to Skesa-Unicycler subcontigs)

For whole contigs:

* `ske_<k>` for Skesa (with $k \in \{1...|S|\}$, with $S$ the set of Skesa whole contigs)
* `uni_<k>` for Unicycler (with $k \in \{1...|S|\}$, with $S$ the set of Unicycler whole contigs)

#### Tags

| ID   | Type              | Values                    | Description                                        |
| ---- | ----------------- | :-----------------------: | -------------------------------------------------- |
| `LN` | `i` (int)         | $\mathbb{N}_{>0}$         | Length of the sequence                             |
| `dp` | `f` (float)       | $\mathbb{R}_{>0}$         | Normalized sequence coverage                       |
| `OC` | `i` (int)         | $\mathbb{N}_{>0}$         | From how many assemblies comes from the sequence   |
| `cl` | `Z` (string)      | `name1,name2,...`         | Contigs from which the sequence comes from         |
| `cp` | `B` (float array) | list of $\mathbb{R}_{>0}$ | Part of the contig in the order of the contig list |
| `ns` | `A` (char)        | $\{S, U, s, u, b\}$       | Nature of the segment, see below                   |
| `ap` | `f` (float)       | $\mathbb{R}_{>=0}$        | Pangenome penalty, see below                       |

<!-- REFACTOR change ll by cp -->
<!-- REFACTOR always use dp, remove the use of cv -->
<!-- REFACTOR change aa by ns -->

**Nature of the segment (`ns` tag):**

* `S` whole Skesa contig (from Skesa assembly)
* `U` whole Unicycler contig (from Unicycler assembly)
* `s` Skesa subcontig (from pangenome Skesa-Unicycler)
* `u` Unicycler subcontig (from pangenome Skesa-Unicycler)
* `b` Skesa and Unicycler subcontig (from pangenome Skesa-Unicycler)

**Pangenome penalty (`ap` tag):**

* If `aa:A:b`: the penalty equals $0$
* Otherwise: the penalty equals $max(1, length/1000)$

The idea here is to favour the subcontig shared by the two assemblers by penalizing the others.

### Links

#### Overlap match

`0M`

#### Tags

| ID   | Type         | Values        | Description                   |
| ---- | ------------ | ------------- | ----------------------------- |
| `lo` | `A` (char)   | $\{p, s, u\}$ | Origin of the link, see below |
| `lt` | `Z` (string) | `<x><y>`      | Link type, see below          |

<!-- REFACTOR change aa by lo -->

**Origin of the link (`lo` tag):**

* `p` pangenome link (between a subcontig from Skesa or Unicycler, and a subcontig from both Skesa and Unicycler)
* `s` Skesa link
* `u` Unicycler link

**Link type (`lt` tag):**

Between subcontigs:

* `su` Skesa subcontig to Unicycler subcontig
* `bs` Skesa-Unicycler subcontig to Skesa subcontig
* `bu` Skesa-Unicycler subcontig to Unicycler subcontig
* `bb` Skesa-Unicycler subcontig to Skesa-Unicycler subcontig

Between a whole contig and a subcontig:

* `sS` Skesa subcontig to Skesa whole contig
* `uU` Unicycler subcontig to Unicycler whole contig
* `bS` Skesa-Unicycler whole contig to Skesa whole contig
* `bU` Skesa-Unicycler whole contig to Unicycler whole contig

Between whole contigs:

* `SS` Skesa whole contig to Skesa whole contig
* `UU` Unicycler whole contig to Unicycler whole contig
* `SU` Skesa whole contig to Unicycler whole contig

Redundant combinations (and their reverse order):

* `ss` corresponds to `bs` or `sb`
* `uu` corresponds to `bu` or `ub`
* `sU` corresponds to `bU`
* `uS` corresponds to `bS`

### Paths

<!-- DOCU path tags -->

#### Names

* `ske_<k>` for Skesa (with $k \in \{1...|S|\}$, with $S$ the set of Skesa whole contigs)
* `uni_<k>` for Unicycler (with $k \in \{1...|S|\}$, with $S$ the set of Unicycler whole contigs)

<!-- REFACTOR remove tag for paths -->