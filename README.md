# DDD_chrX
Code used to generate results for the DDD chrX paper


This repository contains a series of R scripts and input files used to run the analysis in "The contribution of X-linked coding variation to severe developmental disorders" https://www.medrxiv.org/content/10.1101/2020.03.18.20037960v1.

Please copy this reposity to a new location and then, in each R script, alter the setwd() command to set the directory to this new location.

##### IMPORTANT: These R scipts have been written to run with R version 3.4.0.
The following R packages are required in at least one of the R scripts: epiDisplay, data.table, ratesci, plyr, metap

Note that these scripts make use of files of post-QC data which have been deidentified by removing unnecessary detail on individuals' and assigning randomised IDs.

####
1. script1.compare_phenotypes_between_sexes_in_DDD.R      : This script runs the phenotype comparisons presented in the paper, and produces Supplementary Figure 1, Supplementary Table 1 and  Supplementary Table 2.

Note that you cannot run this script fully because we are not releasing the full phenotype data due to patient confidentiality reasons.
However, the HPO terms for probands can be obtained by applying for dataset EGAD00001004388 on EGA.
We have included a short file of fake phenotypes to illustrate the format required: fake_phenotype_data.example_file_to_illustrate_format.txt.
####
2. script2.case_control_analysis.R      : This script runs the case/control analysis including boostrapping.

The part of boostrapping takes > 12 hours, so users may wish to submit this R job to a cluster by modifying script2.run_chrX_case_control_analysis_including_bootstrapping.sh.
Alternatively, they can use the output of the boostrapping I ran which is in results_of_bootstrapping_in_case_control_analysis.all.RData
####
3. script3.burden_analysis_and_per_gene_tests_on_de_novos.R  : This script runs the burden analysis and per-gene tests on de novo mutations.

The output table de_novo_enrichment_tests_per_gene.with_sex_bias_test.txt feeds into Supplementary Table 5
.
If you wish to skip the slow boostrapping part, you can make use of the boostrapping I ran which is in results_of_bootstrapping_for_DNM_burden.RData.
####
4. script4.prepare_data_for_TADA.R   : This script munges the case/control variants and de novo mutations to prepare data for input to TADA.
####
5. script5.run_TADA_to_do_per_gene_tests_on_males.R   : This script runs the per-gene tests with TADA  making use of the script TADA.v.1.2.R which was written by Bert Klei.

The output table output_from_TADA.frac_risk_genes_0.05.txt feeds into Supplementary Table 4. Note that the TADA results achieveed should be similar to the paper, but won't be identical, as there is a degree of stochasticity. This script outputs Supp_Table_7.parameters_run_with_TADA.txt too.

####
6. script6.extract_and_plot_results_from_burden_analysis_for_chrX_paper.R    : This script takes the output of scripts 2-5, extracts key summary statistics reported in the paper, and plots the results of the burden analysis and per-gene tests. It makes Figures 1 and 2, Supplementary Figures 2, 3, 5.

Note that the boostrapping results will be slightly different to those reported in the paper, since there is stochasticity in the sampling.

####
7. script7.compare_PRS_between_patient_subsets_for_chrX_paper.R    : This script compares polygenic scores between different groups of male probands. It produces Supplementary Table 6.

