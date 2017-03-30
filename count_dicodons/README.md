## Count di-codones pipeline

This set of python scripts counts the dicodon composition of a dna sequence.

Usage: `python3 count_dicodones.py <input_library> <output_counts>

#### Description of arguments

**input_library |** tabular files with seq id followed by sequence of nucleotides in uppercase, example:

```
id	sequence
AD_1988		ATCGCGTAGCTGAGTTCAGCGTCGAGTCGATTCGTTGAT
AD_0203		CGATGCTACATCGACTGATCGATGCACGACTACGTACGT
GKJG_283	GGCTGACTCAGCTAGCTACGTACGTAGCTAGCTAGCTGACGACGAGCG
```

**output_counts |** name of output tabular file with the counts

#### Files in this directory

**helper.py |** script containing the main function to count the codons.
