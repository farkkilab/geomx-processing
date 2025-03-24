# **geomx-processing**

pipeline for pre-processing and basic analysis of geomx data

## Pipeline steps
This pipeline starts with (almost) raw geomx dsp data - DCC files (obtained from fastq and typically delivered by sequencing centre) and performs the following pre-processing and analysis steps:

1. QC
    * filtering AOIs based on basic QC parameters (nr of aligned reads, hight negative control counts)
    * filtering AOIs based on negative probes modelling
    * filtering probes based on geometric mean and grubbs test
    * filtering segments based on LOQ (limit of quantification)
    * calculating gene detection rate  
3. Normalisation
    * quantile (Q3) normalisation
    * Deseq2 normalisation
    * vst (variance stabilising transformation) on deseq2 normalised data
    * UMAP and t-SNE projections on all types of normalisation
4. Batch effect correction
    * PVCA variation assessment on deseq2 normalised data to find variables responsible for batch effect
    * batch effect correction with limma
    * batch effect correction with harmony
    * PVCA on batch-effect corrected data to assess the correction results
    * UMAP and t-SNE od batch-effect corrected data

5. Deconvolution
   * preparing reference scRNAseq dataset (incl removing of low complexity genes)
   * computing cell fractions and ct-specific expression profiles with BayesPrism
   * post-processing of ct-specific expression profiles
       * removing non-reliable predictions
       * vst
       * removing genes with variance = 0 (the same value across all AOIs - artifact of BP + vst)
       * harmony and limma batch effect correction
   * computing cell fractions with SpatialDecon
   
6. Pathway analysis (on full and/or deconvoluted signal)
   *  removing low complexity genes
   *  calculating ssGSEA/GSVA for signatures from selected categories of msigdb database
   *  calculating PROGENY scores (!! currently disabled due to incompatibility issues)
9. differential gene expression (on full and/or deconvoluted signal)
    * calculating differentially expressed genes between specified group of AOIs


## Requirements

The pipeline was set up and tested under **UBUNTU 22.04** on the desktop machine with **20 cores and 64G RAM (32G + 32G swap)** 
and it took around **5h to complete** for ~300 DCC files. 

**R 4.2.2** version was used. All the required packages and versions can be found in **renv** file. 

The weight of the ~300 DCC dataset is around 2G, so it's possible to run the pipeline on the smaller machine (around 8G RAM) 
but the required tests haven't been performed.
The most heavy computational steps are (3) Batch effect correction (PVCA plots), (4) Deconvolution and (6) Differential Gene Expression


## Input and output data TODO

Input files:
* **DCC files** - raw expression data for each AOI
* **.pkc file** - provided by GeoMx - information about sequencing probes and library used
* **annotation .xls file** - spreadheed containing all metadata information
    * have to contain metadata in the spreadsheet named 'Sheet1'
    * HAVE TO CONTAIN the following columns (with the exact same names): 'Sample_ID', 'Slide_Name',  'Aoi', 'Roi' and 'Panel'
    * have to contain all other columns specified in the "set up metadata variables names" section
    * It's best to avoid any whitespace (" ") in the column names
    * if during the experiment NTC AOIs were messed up and doesn't follow 1 NTC/batch_nr - 'NTC_ID' column have to be added manually
      with correct dcc_filename of the NTC AOI corresponding to each other AOI
* **scRNAseq reference file** in .RDS format - Seurat object with reference scRNAseq dataset, used for deconvolution
    * have to contain 'cell_type' column in metadata with cell type label (optionally other grouping variable specified while running deconvolution script)


## How to use the pipeline

This is an academia-grade software. It's not fully optimised for the usage of resources (computation power and memory), usage friendliness and clean code.
Despite the fact that it can be run at one shot, it requires the user to carefully check the intermediate results and fine-tune the parameters 
(eg thr for QC, batch effect variables etc) if the results are not satisfactory. For the first-time run on the new data it's advisable to run it step-by-step and check the results.

Read the comments in the code to learn about different variables and carefully examine the terminal output to learn more about the run. 
If you encounter an error: read an error message - often it's an easy-to-fix problem such as wrong variable/pathway name.

How to use the pipeline:
1. Import the environment using renv (https://rstudio.github.io/renv/articles/renv.html) **TODO**
2. Change the pathways to in/out directory and all input files in the "define variables and paths" section
3. Carefully check your metadata (annotation) xls file and change the column names within "set up metadata variables names" section. 
4. Just run the next 2 sections - "load util functions and create dirs" and "define intermediate output paths"
5. Run the actual pipeline steps one-by-one:
   * adjust important parameters in each step section
   * each step is run by run_unless_exists() function which check for the existence of given file (passed as a 2nd argument of the function - typically intermediate RDS file or logs txt file)
     and if cannot find it - run the script for a specific step. That means that the user have to be careful with the file paths.
     Every time one of the main parameters is changed, the output file path changes and the code is re-run
     (eg for rerunning GSEA or DGE with multiple different signatures/groups). All the parameters are saved in the log file for further examination and comparison.
     If the user wish to re-run the step with the same parameters (or changed the advanced parameters, which doesn't change the intermediate filename),
     they can simply re-name/move/delete the step-output filename and re-run the step.
    * more advanced parameters which typically don't need further adjustments are listed on top of each step-script within the "define variables" section.
      User can experiment with them if needed (as concluded from the inspection of intermediate data)
7. Run the whole pipeline with desired parameters - once the parameters for the crucial steps (1-4) have been adjusted,
   the pipeline may be re-run as a whole with different set of parameters for downstream analysis
   (steps 5-6 - GSEA and DGE) - the already-computed steps (for which the output file exists) will be ommited.

## Final Remarks

The author of the pipeline is **Iga Niemiec** (https://github.com/igbiga). Contact me in case of issues (iga.niemiec@helsinki.fi).

If you find a bug in the code or flaw in the analysis - please contact me directly or report a bug on github. I'd highly appreciate it!

Further steps:

We're planning to add 2 more downstream analysis steps:
1. Ligand-Receptor analysis
2. WGCNA

There are as well multiple 'TODO' in pipeline code - these additional analysis may be added in the future. 
The pipeline is not intended to do any type of more specific downstream analysis since it depends heavily on the data and scientific question

Have fun!
