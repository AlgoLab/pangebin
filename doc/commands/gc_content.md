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

**Refactoring:**

* [ ] #REFACTOR GC content score should be an array -> adapt classes
* [ ] #REFACTOR GC content proba should be an array -> adapt classes


### GC score file format

`gc_prob_scores.yaml`

```yaml
intervals:
  - - <min_gc>
    - <max_gc>
  # ...

sequences:
  <fragment_id>:
    - - <gc_prob>
      - <gc_score>
  # ...
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