## Files to run pipeline

**barcodes_r1.txt |** these barcodes are only for the spplitting proccess, some barcodes where trimmed
to make all of theme equall lenght

**minigenes_barcodes.txt |** these are artificial barcodes to quantify the minigenes, these barcodes where generated from
the minigenes libraries *Mini_genes_library.txt* using the following command:
```
sed 's/CTACACGACGCTCTTCCGATCTCCTGCAGGcaccATG//g' Mini_genes_library.txt  | perl -e 'while($id = <>){$seq = <>; $seq = substr($seq,0,15); print $id.$seq."\n"}' | fasta_formatter -t > minigenes_barcodes.txt
```
