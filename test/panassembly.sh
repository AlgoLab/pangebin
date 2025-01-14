#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/panassembly.sh {run|clean} test/<dataset>
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_dir="$dataset_dir/result/panassembly"

run() {
    local datadir="$dataset_dir/data"
    local std_asm_graph_datadir="$datadir/std_asm_graph"
    local pangenome_datadir="$datadir/pangenome"

    local pangenome="$pangenome_datadir/pangenome.gfa"
    local skesa="$std_asm_graph_datadir/skesa.gfa"
    local unicycler="$std_asm_graph_datadir/unicycler.gfa"

    pangebin panassembly $pangenome $skesa $unicycler \
        --outdir $out_dir
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
