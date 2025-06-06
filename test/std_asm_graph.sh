#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/std_asm_graph.sh {run|clean} test/<dataset> [--debug]
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_dir="$dataset_dir/result/std_asm_graph"

debug=""
if [ "$3" == "--debug" ]; then
    debug="--debug"
fi

run() {
    local datadir="$dataset_dir/data"
    local assembly_datadir="$datadir/assembly"

    local skesa_gfa="$assembly_datadir/skesa.gfa.gz"
    local unicycler_gfa="$assembly_datadir/unicycler.gfa.gz"

    pangebin sub std-asm-graph $skesa_gfa $unicycler_gfa \
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
