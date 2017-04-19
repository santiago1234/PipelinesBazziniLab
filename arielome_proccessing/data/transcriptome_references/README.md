
# Transcriptome References

I got these references from Ariel Paulson (bioinformatics core), I copied these files (from: /core/Bioinformatics/analysis/Bazzini/arb/cbio.arb.101/code/refs)into this directory
and renamed:

```
mv danRer10.Ens_84.longest.fa fish.fa
mv dm6.Ens_84.longest.fa fly.fa
mv pombe_ASM294v2.EnsGen_31.longest.fa pombe.fa
mv smedSxl_v31.MAKER_20150401.longest.fa planaria.fa
```

These reference are the longest transcript sequences if multiple different transcript exist.

To create the indexes with bowtie I run the following command:

```
ls *fa | parallel "bowtie-build  {} {.}"
```

Notes:

*owtie-align version 1.2* 
