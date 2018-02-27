# Shiny app: Salmon and kallisto differences on Geuvadis and SEQC data

## How to load the Shiny app

1) Clone or Download this repo from GitHub (either clone or Download ZIP and unzip)

2) Open `geuvadis.R` or `seqc.R` in RStudio. For `seqc.R` you can change the `dataset` variable to `"A"` or `"B"` for the two different types of samples processed by SEQC.

3) Click **Run App**

4) Clicking on a point on the left side will open a plot on the right side showing the normalized counts for all the samples for the transcript you clicked (green and orange colored points) and for the other transcripts of the gene that had differences. This example shows kallisto with very low counts for the samples from the CNL sequencing center. Green indicates sequencing centers with more fragment sequence bias, while orange indicates less fragment sequence bias. The diagonal lines are `y = 2x` and `x = 2y`.

![Example of clicking a transcript](example.png)

The SEQC quantifications are those from the [Salmon publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/). Quantifications were made against the RefSeq gene annotation file and the genome FASTA contained within the hg19 Illumina iGenome using Salmon 0.8.0 with `--gcBias --seqBias` and kallisto 0.43.0 with `--bias`. The LRT in DESeq2 was used to generate p-values for each transcript, comparing a model with sequencing center to a model with just the intercept.

The Geuvadis quantifications were made against Gencode v26 CHR transcripts using Salmon 0.8.2 with `--gcBias --seqBias` and kallisto 0.43.1 with `--bias`. Limma-voom was used to generate F statistics for each transcript, comparing a full model with sequencing center, sex, and population to a reduced model with only sex and population.
