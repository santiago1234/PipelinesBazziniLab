Count codons by position pipeline

Description:
    This pipeline counts the codons by positions.

This is the protocol

1. Get begining and end positions: 01-select_concatenate.R
    the user gives and input library and a position size, this size is in nucleotides,
    the first k codons and the last k codons are concatenated to produce new sequences.
    Note: the sequences that are smaller than the minimum size k*3*2 are note consider,
        no Ns in the sequences only nucletides are accepted {A, T, G, C}

2. Create unique sequences representing each position: 02-codons-by-positions.sh
    For instance all the codons in position j
    of each sequence in the library will be represented by a unique sequence.

3. Count codons by possition: 03-count-codons-by.py
    counts the codons by position

4. Visualize results and saves output of pipeline:  04-write-results.R

The output of the program is a table like this:

Position Codon Frequency Density
1        ATG   10        0.1
1        GGG   12        0.3
...
...

Finaly to run the pipeline use the bash script: 00-count_by_position.sh
