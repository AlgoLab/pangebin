#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# USAGE
#
# ./test/gc_scores.sh {run|clean} test/<dataset> [--debug]
# ------------------------------------------------------------------------------------ #

command=$1
dataset_dir=$2

out_file="$dataset_dir/result/panassembly/gc_scores.tsv"

debug=""
if [ "$3" == "--debug" ]; then
    debug="--debug"
fi

run() {
    local datadir="$dataset_dir/data"
    local panassembly_datadir="$datadir/panassembly"
    local panassembly_gfa="$panassembly_datadir/panassembly.gfa"

    pangebin sub gc from-gfa $panassembly_gfa $out_file \
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
