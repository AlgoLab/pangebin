#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/preprocess.sh {run|clean} test/<dataset>
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_dir="$dataset_dir/result/preprocess"

run() {
    local datadir="$dataset_dir/data"
    local assembly_datadir="$datadir/assembly"

    local sample_id=$(basename $dataset_dir) # REFACTOR why need the sample ID?
    local unicycler_gfa="$assembly_datadir/unicycler.gfa.gz"
    local skesa_gfa="$assembly_datadir/skesa.gfa.gz"

    local thr=1

    pangebin preprocess $sample_id $unicycler_gfa $skesa_gfa \
        --outdir $out_dir \
        --thr $thr
}

case $command in
run)
    run
    ;;
clean)
    rm -rf "$out_dir"
    ;;
*)
    echo "Unknown command '$command'"
    exit 1
    ;;
esac
