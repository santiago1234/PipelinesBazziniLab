Pilpeline for counting codons in random transcriptome libraries
Bazzini lab


usage: sh 00-count_codons_pipeline.sh <input_file> <minimum_codon_length> <results_dir> <output_prefix>


Args:
    input_file: a tabular input file with 2 columns: id transcrit &  dna sequence transcrit
    minimum_codon_length: a positive integer for filtering transcripts that are less than minimum_codon_length
                          usually this number should be 10
    
    results_dir: the name of the output directory to write the files genereted during the run
    
    output_prefix: the ouput prefix for the files generated in the run.


Returns:
    this pipeline wrrite the following files:
    
    output_prefix.01: the same as input_file but with 2 extra column indication the position of the
                      premature stop codon & the length in codons of the sequence, if there is not
                      prematures stop a -1 is pinted.

    output_prefix.02: the same as output_prefix.01 but only transcripts with no premature stop are display
                      & with length grater than minimum_codon_length.

    output_prefix.03: tabular file with a row representing each transcript and columns indicating the density
                      of each codon & transcript length

    
    output_prefix.pdf: plots to visualize the results of the pipeline.


NOTE: all the transcrit sequences are assumed to be coding dna sequences and only the first open reading
      frame is used. python3 version is used.
