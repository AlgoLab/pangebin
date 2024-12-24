#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/preprocess.sh {run|clean} test/<dataset> [--debug]
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_dir="$dataset_dir/result/preprocess"

debug=""
if [ "$3" == "--debug" ]; then
    debug="--debug"
fi

run() {
    local datadir="$dataset_dir/data"
    local assembly_datadir="$datadir/assembly"

    local unicycler_gfa="$assembly_datadir/unicycler.gfa.gz"
    local skesa_gfa="$assembly_datadir/skesa.gfa.gz"

    local min_contig_length=1

    pangebin preprocess $unicycler_gfa $skesa_gfa \
        --min-contig-length $min_contig_length \
        --outdir $out_dir \
        $debug
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
