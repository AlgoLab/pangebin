# PlasBin-flow adapted to pan-assembly graph

**Inputs:**

* pan-assembly graph
* Seed fragment file
* Plasmidness score file
* GC score file

**Outputs:**

* Bin directories with files

**Tasks:**

* [x] Loop until there is no remaining seed fragments
  * [x] Compute bin

## Bin file formats

In each bin directory `bin_<k>`, there are two files.

`bin_<k>_stats.yaml`

```yaml
cumulative_sequence_length: <int>  # Cumulative fragment length, without considering fragment repetition
coverage_flow: <float>  # Total flow value, use to normalize the coverages
GC_content_interval:
  - <float>  # Min GC content
  - <float>  # Max GC content

milp_stats:
  MCF:
    coverage_score: <float>  # CoverageScore for MCF
  MGC:
    coverage_score: <float>  # CoverageScore for MGC
    GC_probability_score: <float>  # GCProbabilityScore for MGC
  MPS:
    coverage_score: <float>  # CoverageScore for MPS
    GC_probability_score: <float>  # GCProbabilityScore for MPS
    plasmidness_score: <float>  # PlasmidnessScore for MPS

```

`bin_<k>_fragments.tsv` (tab-separated values)

```html
Fragment_ID    Normalized_coverage
<fragment_id>  <fragment_normalized_coverage_in_bin>
...
```
