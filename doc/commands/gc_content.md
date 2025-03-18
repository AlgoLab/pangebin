# GC content module

## GC scores (from panasm graph)

**Inputs:**

* pan-assembly graph GFA file
* GC content interval file

**Outputs:**

* GC score file

**Tasks:**

* [x] Compute GC probability from sequences
* [x] Compute GC scores from GC probabilities
* [x] Out GC score file

### GC score file format

`gc_prob_scores.tsv`

```html
sequence_id  <interval_1>  <interval_2> ... <interval_n>
<id_1>       <score_1>     <score_2>    ... <score_n>
```

## GC content intervals

* [ ] Understand how to obtain them as for seed thresholds

### GC content interval file format

`gc_content_intervals.txt`

```txt
0.0
<step_1>
...
1.0
```

<!-- FIXME find which value is the equilibrum -->