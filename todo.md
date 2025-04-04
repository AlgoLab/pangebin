# To-do

* [x] Mapping genes to reference sequences
* [x] Gene density `src/pangebin/gene_density/todo.md`
* [ ] Seed sequences `src/pangebin/seed/todo.md`
  * [x] Threshold: Nextflow pipeline
* [x] GC content `src/pangebin/gc_content/todo.md`
* [x] Adapt PlasBin-flow `src/pangebin/plasbin/todo.md`
* [ ] PangeBin-flow on pansm from PBf inputs
  * [ ] Correct the contigs names as in the plasmidness and seed files they are not formatted as in the panassembly graph (from std asm graph)

## Perspectives

### Methods

* [ ] Result filter utilities (on bins...)

### Features

* [ ] Verify parameters given in configs (through typer and configs API)

## Remarks

* Minimum flow:
  * to adress non-0 flow issue
  * how to interpret dp? (why so low?)
  * when the MipGap is to high, the total flow always equal to the minimum flow: big loop

From log gurobi:

```log
Warning: Model contains large matrix coefficient range
         Consider reformulating model or setting NumericFocus parameter
         to avoid numerical issues.
```
