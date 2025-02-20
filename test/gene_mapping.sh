#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/pangenome.sh {blast|filter|clean} test/<dataset> [--debug]
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_dir="$dataset_dir/result/gene_mapping"

debug=""
if [ "$3" == "--debug" ]; then
    debug="--debug"
fi

gene_fasta="src/pangebin/database/genes.fasta"

blast() {
    mkdir $out_dir 2>/dev/null

    local datadir="$dataset_dir/data"
    local std_asm_graph_datadir="$datadir/std_asm_graph"
    local mixed_fasta="$std_asm_graph_datadir/mixed.fasta"

    local original_sam="$out_dir/gene_map_on_contigs.sam"

    pangebin utils map blast $gene_fasta $mixed_fasta $original_sam \
        $debug
}

filter() {
    local datadir="$dataset_dir/data"
    local std_asm_graph_datadir="$datadir/std_asm_graph"
    local mixed_fasta="$std_asm_graph_datadir/mixed.fasta"

    local gene_mapping_datadir="$datadir/gene_mapping"
    local original_sam="$gene_mapping_datadir/gene_map_on_contigs.sam"
    local filter_config_yaml="test/config/filter_config.yaml"

    local filtered_sam="$out_dir/gene_map_on_contigs.filtered.sam"

    pangebin utils map filter $original_sam $filtered_sam \
        --query-fasta $gene_fasta \
        --config $filter_config_yaml \
        $debug
}

case $command in
blast)
    blast
    ;;
filter)
    filter
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
