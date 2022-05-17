[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6558500.svg)](https://doi.org/10.5281/zenodo.6558500)

## 2021 ISB Virtual Microbiome Symposium   

# Day 2 course data

This repository contains cached data and processing steps for day 2 of the symposium. This is split into two major pipelines: [1] Obtaining data from the BioML data set and processing it into assemblies and [2] gnerating carveME reconstructions for all assemblies with decent GTDB assignments.

**Required compute:** About 1000 CPU hours

### Obtaining data and processing it

This is all wrapped into a nextflow pipeline which is provided along with this repository: [assemly.nf](assembly.nf). There is [conda environment file](conda.yml) to set up all required dependencies. it covers the following steps.

1. Downloading the first 1000 isolate genomes from the [BioML paper](https://doi.org/10.1038/s41591-019-0559-3).
2. Quality filtering and trimming with FASTP.
3. Assembly with MEGAHIT.
4. Taxonomic placement with the GTDB toolkit.

After that the data is curated by hand to remove isolates with no clear GTDB bacterial assignment. This is contained in an [Rstudio notebook](curation.rmd). This will leave a little less than 980 assemblies.

### Model reconstruction

This done using the [Gibbons Lab model builder pipeline](https://github.com/Gibbons-Lab/pipelines/tree/master/model_builder). The required media are provided in the repository as well.
