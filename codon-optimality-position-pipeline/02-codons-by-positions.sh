#!/bin/bash
set -e
set -u
set -o pipefail

# get codons by position and makes sequences
if [ "$#" -lt 2 ]
then
    echo "error: too few arguments, 2 needed"
    echo "usage: script.sh <input_file> <position_size>"
    exit 1
fi

# initialize variable with user args
input_file="$1"
number_of_codons="$2"

# run python script
python codons_by_position.py $input_file $number_of_codons

# put codons by possition together

if test -e seqs_tmp
then
    rm seqs_tmp
fi

touch seqs_tmp
positions=($(cut -f 1 "positions"))

for i in ${positions[@]}
do
    echo "$i of $number_of_codons*2"
    seq=($(cut -f $i tmp_seqs | tr -d '\n'))
    echo "$i $seq" >>seqs_tmp
done

rm tmp_seqs
echo "job completed"

