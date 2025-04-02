#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/plasbin_binlab_panasm.sh {run|clean} test/<dataset> [--debug]
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_dir="$dataset_dir/result/plasbin/binlab"

binning_config="test/config/binning_config.yaml"
binlab_config="test/config/plasbin_binlab_config.yaml"
gurobi_config="test/config/gurobi_config.yaml"

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

    pangebin sub plasbin binlab panasm $panassembly_gfa $seeds $gc_scores $plasmidness \
        --bin-cfg $binning_config \
        --binlab-cfg $binlab_config \
        --gurobi-cfg $gurobi_config \
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
