# Standardized assembly graphs

Standardize the tags of SKESA and Unicycler assembly graphs.

## Header

| ID   | Type       | Values | Description                     |
| ---- | ---------- | ------ | ------------------------------- |
| `Sd` | `A` (char) | `Y\|N` | The graph is (not) standardized |

## Segments

### Names

* `ske_<k>` for SKESA
* `uni_<k>` for Unicycler

For $k \in \{1...|S|\}$, with $S$ the set of segments.

### Tags

| ID   | Type        | Values            | Description                  |
| ---- | ----------- | ----------------- | ---------------------------- |
| `LN` | `i` (int)   | $\mathbb{N}_{>0}$ | Length of the sequence       |
| `dp` | `f` (float) | $\mathbb{R}_{>0}$ | Normalized sequence coverage |

## Links

### Overlap match

`0M`
