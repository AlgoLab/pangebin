# Pangebin (WIP)
Check ISMB2025 branch.


### activate env
`source venv/bin/activate`

### 1. preprocess input files
`sample_id="SAMN16357463"`
`unicycler="test/data/SAMN16357463/unicycler.gfa.gz`
`skesa="test/data/SAMN16357463/skesa.gfa.gz"`
`out_dir="test/data/SAMN16357463/out"`
`thr="1"`

`python3 src/pangebin preprocess ${sample_id} ${unicycler} ${skesa} --outdir {out_dir} --thr {thr}`

`mixed_fasta="test/data/SAMN16357463/out/SAMN16357463.1.mix.fasta"`
### 2. make pangenome graph using ${mixed_fasta}
(pggb ...)
`pangenome="test/data/SAMN16357463/out/SAMN16357463.1.pan.gfa"`
`cp ${pggb_out} ${pangenome`
### 3. make pan-assembly graph
`pangenome="test/data/SAMN16357463/out/SAMN16357463.1.pan.gfa"`
`skesa="test/data/SAMN16357463/out/SAMN16357463.1.s.gfa"` 
`unicycler="test/data/SAMN16357463/out/SAMN16357463.1.u.gfa"`
`sample_id="SAMN16357463"`
`out_dir="test/data/SAMN16357463/out/"`

`python3 src/pangebin panassembly $pangenome $skesa $unicycler ${sample_id} ${out_dir}`



