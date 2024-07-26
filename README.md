# NMR Assay Dataframe Curation
## Crystal R. Musser, M.S. (Wet Lab Manager, Katrina Claw Lab @ CU Anschutz, Department of Biomedical Informatics)
## v.05/16/2024

## PURPOSE: 
Combine quantified  LC/MS nicotine metabolite data with phase survey self-reported data, from specified cohort.

## OBJECTS:
"t_allSamplesMetabolites" for curated nicotine metabolite dataset
"t_allSamples_postExclusions" as the master (all variables all samples dataframe, post-QC, after specific exclusions)
"clean_df" as the finalized, shared dataset, 18 variables, n=816 (inc. 28 controls)

## SUPPLEMENTARY FILES to this script: 
Data dictionaries (Supplementary File 1)
"Clean" variable consensus/logic UML (Supplementary File 2)

## REQUIRED INPUTS: 5
1) **v6RawData_combined.csv** : combined experimental LC-MS/MS calculated and QC'd data in [ng/mL] serum concentrations for: , Cotinine, 3-hyroxy-Cotinine, Nicotine, CO, NO, Nornicotine (in this order)
2) **phase2.xlsx** : cohort survey results with associated data dictionary for interpretation  
3) **smoking.xlsx** : dataset with relevant smoking-related phenotypes not included in file 2
4) **phase1.xlsx** : phase I survey results with associated data dictionary for gender, smoking histotry, height, weight
5) **shs677_elig.sas7bdat** : age at phase 1
6) **AlleleActivityScoring.xlsx** : CYP2A6 haplotype variant activity score reference table
7) **pypgx_cyp2a6_cs_20240524.xlsx** : @CarissaSherman's PyPGx output dataset for SHS ph2 cohort

## OUTPUT FILES GENERATED: 
1) A list of all combined results for all 824 samples; t_allSamples_postExclusions, file = "SHSph2_fullDF_05152024CRM.csv"
2) A "clean" dataset with 18 variables for the final 816 participants to use for EWAS/GWAS; clean_df, file = "SHSph2_NMRwSurveyDF_clean05152024CRM.csv"
3) A list of just sampleIDs of final 816 participants; clean_df$idNo, file = "SHSph2_IDs_clean05152024CRM.csv"

## NOTE:
- I'm bringing in calculated metabolic data using Kaja's data output from Teddy's script referencing Laura's library, and am assuming reported values pass QC

### This script: 
Calculated metabolite concentration data > Dataframe curation Script > Data Cleaning Script 
### Connecting it from: 
Raw experimental LC-MS/MS data > metabolite standard curve generation > concentation back-calculations
### Connecting it to: 
Data Analysis Script > Figures & Tables Script

## STILL TO DO:
1) write in a check that there's not data being overwritten for that sample's concentration variable (ie reruns in the same file), handle this how?
2) generate a NMR assay rerun sample list and connect data back in
3) generate a sample list for phase 3 sample IDs
   
