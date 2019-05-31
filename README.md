# BSGatlas workflow

The folder structure in this project is:

* `analysis`: RMarkdown files describing each step in detail
* `data-raw`: Utilized external data
* `scripts`:  longer scripts and helper functions that are used

The individual raw data is described in a separate step.
The rmarkdown notebooks were compiled via rstudio.

For improved reproducbility, an exact description of the used conda environment
is stored in `conda.yml` with

    conda env export > conda.yml

If you wish to recreate this environment, install it and activate it wiht

    conda env create --file cona.yml
    conda activate bsgatlas

