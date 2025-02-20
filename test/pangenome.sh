#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/pangenome.sh {run|clean} test/<dataset> [--debug]
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_dir="$dataset_dir/result/pangenome"

debug=""
if [ "$3" == "--debug" ]; then
    debug="--debug"
fi

supp_nfcore_pangenome_config_path="test/config/supp_nfcore_pangenome_config.yaml"

nfcore_pan_work_dir="work"
nfcore_pan_logs=".nextflow.log*"

run() {
    local datadir="$dataset_dir/data"
    local std_asm_graph_datadir="$datadir/std_asm_graph"

    local mixed_fasta="$std_asm_graph_datadir/mixed.fasta"

    pangebin sub pangenome $mixed_fasta \
        --supplementary-nfcore-pangenome-config-path $supp_nfcore_pangenome_config_path \
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
    echo "Cleaning $nfcore_pan_work_dir"
    rm -rf "$nfcore_pan_work_dir"
    echo "Cleaning $nfcore_pan_logs"
    rm -rf "$nfcore_pan_logs"
    ;;
*)
    echo "Unknown command '$command'"
    exit 1
    ;;
esac
