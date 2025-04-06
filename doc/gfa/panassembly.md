# Pan-assembly graph

## Header

### Tags

| ID   | Type         | Values        | Description       |
| ---- | ------------ | ------------- | ----------------- |
| `VN` | `Z` (string) | `M.m`         | GFA version       |
| `PA` | `Z` (string) | `panassembly` | Pan-assembly type |

## Segments

### Names

For pangenome segments:

* `panske_<k>` for SKESA subcontigs (with $k \in \{1...|S_{ps}|\}$ where $S_{ps}$ is the set of pangenome segments corresponding to SKESA subcontigs)
* `panuni_<k>` for Unicycler subcontigs (with $k \in \{1...|S_{pu}|\}$ where $S_{pu}$ is the set of pangenome segments corresponding to Unicycler subcontigs)
* `panboth_<k>` for SKESA-Unicycler subcontigs (with $k \in \{1...|S_{pb}|\}$ where $S_{pb}$ is the set of pangenome segments corresponding to SKESA-Unicycler subcontigs)

### Tags

| ID   | Type              |          Values           | Description                                        |
| ---- | ----------------- | :-----------------------: | -------------------------------------------------- |
| `LN` | `i` (int)         |     $\mathbb{N}_{>0}$     | Length of the sequence                             |
| `dp` | `f` (float)       |     $\mathbb{R}_{>0}$     | Normalized sequence coverage                       |
| `OC` | `i` (int)         |     $\mathbb{N}_{>0}$     | From how many assemblies comes from the sequence   |
| `cl` | `Z` (string)      |     `name1,name2,...`     | Contigs from which the sequence comes from         |
| `cp` | `B` (float array) | list of $\mathbb{R}_{>0}$ | Part of the contig in the order of the contig list |
| `sn` | `A` (char)        |       $\{s, u, b\}$       | Nature of the segment, see below                   |
| `pp` | `f` (float)       |    $\mathbb{R}_{>=0}$     | Pangenome penalty, see below                       |

**Segment nature (`sn` tag):**

* `s` SKESA subcontig (from pangenome SKESA-Unicycler)
* `u` Unicycler subcontig (from pangenome SKESA-Unicycler)
* `b` SKESA and Unicycler subcontig (from pangenome SKESA-Unicycler)

**Pangenome penalty (`ap` tag):**

* If `aa:A:b`: the penalty equals $0$
* Otherwise: the penalty equals $min(1, length/1000)$

The idea here is to favour the subcontig shared by the two assemblers by penalizing the others.

## Links

### Overlap match

`0M`

### Tags

| ID   | Type         | Values        | Description                   |
| ---- | ------------ | ------------- | ----------------------------- |
| `lo` | `A` (char)   | $\{p, s, u\}$ | Origin of the link, see below |
| `lt` | `Z` (string) | `<x><y>`      | Link type, see below          |

**Origin of the link (`lo` tag):**

<!-- FIXME p lo value description is not correct, because of multiedges... see example of bb lt -->
* `p` pangenome link (between a subcontig from SKESA or Unicycler, and a subcontig from both SKESA and Unicycler)
* `s` SKESA link
* `u` Unicycler link

**Link type (`lt` tag):**

Between subcontigs:

* `su` SKESA subcontig to Unicycler subcontig
* `bs` SKESA-Unicycler subcontig to SKESA subcontig
* `bu` SKESA-Unicycler subcontig to Unicycler subcontig
* `bb` SKESA-Unicycler subcontig to SKESA-Unicycler subcontig

Redundant combinations (and their reverse order):

* `ss` corresponds to `bs` or `sb`
* `uu` corresponds to `bu` or `ub`

## Paths

### Names

* `ske_<k>` for SKESA (with $k \in \{1...|S|\}$, with $S$ the set of SKESA whole contigs)
* `uni_<k>` for Unicycler (with $k \in \{1...|S|\}$, with $S$ the set of Unicycler whole contigs)

### Tags

| ID   | Type        |      Values       | Description                  |
| ---- | ----------- | :---------------: | ---------------------------- |
| `LN` | `i` (int)   | $\mathbb{N}_{>0}$ | Length of the sequence       |
| `dp` | `f` (float) | $\mathbb{R}_{>0}$ | Normalized sequence coverage |
