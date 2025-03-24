# **geomx-processing**

pipeline for pre-processing and analysis of geomx data

## Pipeline steps TODO
This pipeline starts with (almost) raw geomx dsp data - DCC files (obtained from fastq and typically delivered by sequencing centre) and perform the following pre-processing and analysis steps:

1. QC
2. Normalisation
3. Batch effect correction

4. Deconvolution
   
5. Pathway analysis with ssGSEA/GSVA
6. differential gene expression


## Requirements TODO

The pipeline was set up and tested under UBUNTU 22.04 on the desktop machine with 20 cores and 64G RAM (32G + 32G swap) and it took around 5h to complete.
The most heavy computational steps are 3.Batch effect correction (PVCA plots), 4.Deconvolution and 6.Differential Gene Expression


## Input and output data TODO

anno - It's best to avoid any whitespace (" ") in the column names.

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
   5.1 adjust important parameters in each step section
   5.2 each step is run by run_unless_exists() function which check for the existence of given file (passed as a 2nd argument of the function - typically intermediate RDS file or logs txt file)
   and if cannot find it - run the script for a specific step. That means that the user have to be careful with the file paths. Every time one of the main parameters is changed, the output file path changes and the code is re-run
   (eg for rerunning GSEA or DGE with multiple different signatures/groups). All the parameters are being saved in the log file for further examination and comparison. If the user wish to re-run the step with the same parameters
   (or changed the advanced parameters, which doesn't change the intermediate filename), they can simply re-name/move/delete the step-output filename and re-run the step.
   5.3 more advanced parameters which typically don't need further adjustments are listed on top of each step-script within the "define variables" section. User can experiment with them if needed (as concluded from the inspection of 
   intermediate data)
7. Run the whole pipeline with desired parameters - once the parameters for the crucial steps (1-4) have been adjusted, the pipeline may be re-run as a whole with different set of parameters for downstream analysis (steps 5-6 - GSEA and DGE) - the already-computed steps (for which the output file exists) will be ommited.

## Final Remarks

The author of the pipeline is **Iga Niemiec** (https://github.com/igbiga). Contact me in case of issues (iga.niemiec@helsinki.fi).

If you find a bug in the code or flaw in the analysis - please contact me directly or report a bug on github. I'd highly appreciate it!

Further steps:

We're planning to add 2 more downstream analysis steps:
7. Ligand-Receptor analysis
8. WGCNA

There are as well multiple 'TODO' in pipeline code - these additional analysis may be added in the future. 
The pipeline is not intended to do any type of more specific downstream analysis since it depends heavily on the data and scientific question


