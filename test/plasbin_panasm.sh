#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/plasbin_panasm.sh {run|clean} test/<dataset> [--debug]
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_dir="$dataset_dir/result/plasbin"

debug=""
if [ "$3" == "--debug" ]; then
    debug="--debug"
fi

run() {

    local datadir="$dataset_dir/data"

    local panassembly_datadir="$datadir/panassembly"

    local panassembly_gfa="$panassembly_datadir/panassembly.gfa"
    local gc_scores="$panassembly_datadir/gc_scores.tsv"
    local plasmidness="$panassembly_datadir/plasmidness.tsv"
    local seeds="$panassembly_datadir/seeds.tsv"

    pangebin sub plasbin $panassembly_gfa $gc_scores $plasmidness $seeds \
        --outdir $out_dir \
        $debug
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
