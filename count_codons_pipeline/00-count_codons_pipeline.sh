#!/bin/bash
set -e
set -u
set -o pipefail


# check for command-line arguments

if [ "$#" -lt 4 ]
then
    echo "error: to few arguments, you $#, 4 requiered"
    echo "usage: script.sh <input_file> <minimum_codon_length> <results_dir> <out_prefix>"
    exit 1
fi


# initialize variables
input_file=$1
minimum_codon_length=$2
results_dir=$3
out_prefix=$4

echo "creating output dir: $results_dir"
mkdir $results_dir


# detec premature stop codon
echo "detecting sequences with premature stop >>>"
python 01-detect_premature_stop.py $input_file $out_prefix.01 30


# make visualization plots
echo "making quality check plots >>>"
Rscript 02-visualization.R $out_prefix.01 $out_prefix.pdf


# retrive no premature stop transcrits
echo "getting no premature stop sequences >>>"
Rscript 03-retrive_no_premature_stop.R $out_prefix.01 $out_prefix.02 $minimum_codon_length

# count codons
echo "counting codons >>>"
python 04-count_codons.py $out_prefix.02 $out_prefix.03 30


mv $out_prefix* $results_dir/
echo "task completed!!!!"

