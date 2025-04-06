#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/fragment_seeds.sh {run|clean} test/<dataset> [--debug]
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_file="$dataset_dir/result/panassembly/seeds.tsv"

debug=""
if [ "$3" == "--debug" ]; then
    debug="--debug"
fi

run() {
    local datadir="$dataset_dir/data"

    local panassembly_datadir="$datadir/panassembly"
    local gene_density_file="$panassembly_datadir/frag_gene_densities.tsv"

    pangebin sub seed pos-gd $gene_density_file $out_file \
        $debug
}

case $command in
run)
    run
    ;;
clean)
    echo "Cleaning $out_file"
    rm "$out_file"
    rmdir ${out_file%/*} 2>/dev/null
    exit 0
    ;;
*)
    echo "Unknown command '$command'"
    exit 1
    ;;
esac
