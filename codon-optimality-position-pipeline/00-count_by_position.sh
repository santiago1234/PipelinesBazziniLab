#!/bin/bash
set -e
set -u
set -o pipefail

# codon position pipeline

if [ "$#" -lt 3 ]
then
    echo "error: too few arguments, you provided $#, 3 requiered"
    echo "usage: script.sh <input_library> <size> <out_prefix>"
    exit
fi

#initialize variables
input_library=$1
size=$2
out_prefix=$3

echo "geting start and end positions >>>"
Rscript 01-select_concatenate.R $input_library $size tmp_01

echo "obtaining codons by position >>>"
sh 02-codons-by-positions.sh tmp_01 $size

echo "counting codons >>>"
python 03-count-codons-by.py

echo "saving results"
Rscript 04-write-results.R $out_prefix


# remove tmp files
rm tmp_counts
rm seqs_tmp
rm tmp_01
rm positions
