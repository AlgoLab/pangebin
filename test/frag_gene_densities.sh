#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/frag_gene_densities.sh {run|clean} test/<dataset> [--debug]
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_file="$dataset_dir/result/panassembly/frag_gene_densities.tsv"

debug=""
if [ "$3" == "--debug" ]; then
    debug="--debug"
fi

run() {
    local datadir="$dataset_dir/data"

    local panassembly_datadir="$datadir/panassembly"
    local panassembly_gfa="$panassembly_datadir/panassembly.gfa"

    local gene_mapping_datadir="$datadir/gene_mapping"
    local filtered_sam="$gene_mapping_datadir/gene_map_on_contigs.filtered.sam"

    pangebin gd frag $panassembly_gfa $filtered_sam $out_file \
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
