#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/plasbin_pangenome.sh {run|clean} test/<dataset>
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

outdir="$dataset_dir/result/plasbin"

run() {

    local datadir="$dataset_dir/data"

    local panassembly_datadir="$datadir/panassembly"

    local panassembly_gfa="$panassembly_datadir/panassembly.gfa"
    local gc_probabilities="$panassembly_datadir/gc_probabilities.tsv"
    local plasmid_scores="$panassembly_datadir/plasmid_scores.tsv"

    local assembler="pangenome"

    local pangenome_seedlen=1000
    local pangenome_seedscore=0.5
    local pangenome_minplaslen=1000

    local outfile="$outdir/bins.tsv"                    # REFACTOR to not precise
    local log_file="$outdir/plasbin-flow_pangenome.log" # REFACTOR to not precise

    pangebin sub plasbin $panassembly_gfa $gc_probabilities $plasmid_scores $assembler \
        --seed-len-threshold $pangenome_seedlen \
        --seed-score-threshold $pangenome_seedscore \
        --min-plasmid-length $pangenome_minplaslen \
        --out-dir $outdir \
        --out-file $outfile \
        --log-file $log_file
}

case $command in
run)
    run
    ;;
clean)
    echo "Cleaning $out_dir"
    rm -rf "$out_dir"
    ;;
*)
    echo "Unknown command '$command'"
    exit 1
    ;;
esac
