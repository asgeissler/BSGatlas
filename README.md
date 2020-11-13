# BSGatlas workflow

These are the scripts that were used to create the [BSGatlas](https://rth.dk/resources/bsgatlas/)
annotation for *Bacillus subtilis* strain 168.

**BSGatlas: An enhanced annotation integration of the Bacillus subtilis genome**  
Geissler AS, Anthon C, González-Tortuero E, Poulsen LD, Kallehauge TK, Seemann SE, Vinther J, Gorodkin J . Under review 

The folder structure in this project is:

* `analysis`: R scripts that created the annotation and associated analyses
* `data-raw`: Utilized external data in their raw form
* `data`: Parsed form of the external data. Contains both scripts and Rdata objects 
* `data-hub`: Browser hub folder for the final BSGatlas annotation
* `data-gff`: GFF version of the BSGatlas
* `scripts`:  longer scripts and helper functions that were used


For improved reproducbility, an exact description of the used conda environment
is stored in `conda.yml` with

    conda env export > conda.yml

If you wish to recreate this environment, install it and activate it with

    conda env create --file cona.yml
    conda activate bsgatlas


## License

Due to the external resources, this repository's content has *mixed licenses*.

* The **BSGatlas** and the scripts that generated it are under the **Apache License Version 2.0**
* [**SubtiWiki**](http://subtiwiki.uni-goettingen.de) is provided under the **CC BY 4.0*
* [**DBTBS**](http://dbtbs.hgc.jp) was provided by the authors Yuko Makita and Kenta Nakai and is subjected to their copyright.
* The use of the [**BsubCyc**](https://bsubcyc.org) annotation is granted by the royality-free and re-distribution allowing **proprietary academic** license.
  (Includes BioCycTM pathway/genome databases under license from SRI International)
* The results of the [**Rfam**](http://rfam.org) scan were created as part of the BSGatlas and is thus available under the **Apache License Version 2.0**
* The [**RefSeq**](https://www.ncbi.nlm.nih.gov/refseq/) annotation is in the **public domain**.
* The data tables from the supplementary materials by [Dar *et al.*](https://science.sciencemag.org/content/352/6282/aad9822) and [*Nicolas et al.*](https://science.sciencemag.org/content/335/6072/1103) are used for the academic useage


## Note on Nicolas et al. tiling-array signals

With the BSGatlas v1.1, the signals of Nicolas *et al.*'s tiling-array experiments
were added to the browser hub. 
Due their massive size, the raw-files are not included in this repository,
however, they can easily be reloaded from GEO. If you intend to re-run
the script that process the raw data for the browser hub, download the files 
and recreate the following file hierarchy under `data-raw/GEO`.


        data-raw/GEO
        ├── GPL8486_family.soft.gz
        ├── GSE27219_RAW
        │   ├── GPL13149_070910_BaSysBio_expr.ndf.gz
        │   ├── GPL13149_070910_BaSysBio_expr.pos.gz
        │   ├── GSM672549.pair.gz
        │   ├── GSM672550.pair.gz
        │   ├── GSM672551.pair.gz
        │   ├── ...
        │   ├── GSM672813.pair.gz
        │   ├── GSM672814.pair.gz
        │   ├── GSM672815.pair.gz
        │   ├── GSM672816.pair.gz
        │   └── GSM672817.pair.gz
        ├── GSE27219_RAW.tar
        └── GSE27219_family.soft.gz


