
## Set the directory

In the `script` folder, search for *everything* that contains `/Users/jefft/Desktop/p53_project` and replace with your own `xxx/p53_project`.

If you want to do analysis other than eQTL, correct everything in the `p53_project`.


## Part A. eQTL analysis

Get eQTL analysis output for TCGA Pan-cancer datasets

1. Check the master config file in `scripts/eQTL/config/default_vsWT.yaml`
 
    Pay attention to these fields: `experiment_home`: the output folder.

    Make `output` directory under `experiment_home`

2. Check the config file for each cancer type, e.g. `tcga_brca.yaml` (annotation of configs available in this BRCA config)
   
    Pay attention to the covariates.

    Set the "meta mutant" to test: meta mutant refers to a group of mutants stred in `dataset`

3. Run main.R, check using cache or not.
   
   This will generate a series of folders under `output`, for each cancer type.