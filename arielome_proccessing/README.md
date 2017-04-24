## Processing Arielome Pipeline

##### Santiago Medina *santiagogmm@hotmail.com*

:horse:

This is a pipeline written in python using the [luigi](https://github.com/spotify/luigi) module to process the **arielome multi-species library**


#### arielome multi-species library description


#### Input Parameters

The pipeline should be run inside this directory, the input parameters are:

1. `--GlobalConfig-r1 r1.fq`: the fastq mate r1 file
2. `--GlobalConfig-r2 r2.fq`: the fastq mate r2 file (second mate)
3. `--GlobalConfig-outdir out/`: the output directory to output the pipeline files, always append the `/` to the end.

#### Runing the pipeline

To run the pipeline:

`python arielome_pipe.py ExtractSeqs --GlobalConfig-r1 r1.fq --GlobalConfig-r2 r2.fq \`
`--GlobalConfig-outdir out/ --local-scheduler`

see the file *requeriments.txt* for the software requirements, see the file *doc.md* to understand what the pipeline is doing.





