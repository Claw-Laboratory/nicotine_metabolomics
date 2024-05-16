# NMR Assay Dataframe Curation
## Crystal R. Musser, M.S. (Wet Lab Manager, Katrina Claw Lab @ CU Anschutz, Department of Biomedical Informatics)
## v.05/16/2024

## PURPOSE: 
Combine quantified  LC/MS nicotine metabolite data with phase survey self-reported data, from specified cohort.

## OBJECTS:
"t_allSamplesMetabolites" for curated nicotine metabolite dataset
"t_allSamples" as the master (all variables all samples dataframe, post-QC, after specific exclusions)
"clean_df" as the finalized, shared dataset, 18 variables, n=816 (inc. 28 controls)

## SUPPLEMENTARY FILES to this script: 
Data dictionaries (Supplementary File 1)
"Clean" variable consensus/logic UML (Supplementary File 2)

## REQUIRED INPUTS: 4
1) combined experimental LC-MS/MS calculated and QC'd data in [ng/mL] serum concentrations for: , Cotinine, 3-hyroxy-Cotinine, Nicotine, CO, NO, Nornicotine (in this order)      --->   "v6RawData_combined.csv"
2, 3) cohort survey results with associated data dictionary for interpretation                                                                                                 ---> "phase2.xlsx" and "smoking.xlsx"
4) phase I survey results with associated data dictionary for interpretating age, gender, smoking histotry, height, weight                                                  ---> "phase1.xlsx"

## OUTPUT FILES GENERATED: 
1) A list of all combined results for all 824 samples; t_allSamples_postExclusions, file = "SHSph2_fullDF_05152024CRM.csv"
2) A "clean" dataset with 18 variables for the final 816 participants to use for EWAS/GWAS; clean_df, file = "SHSph2_NMRwSurveyDF_clean05152024CRM.csv"
3) A list of just sampleIDs of final 816 participants; clean_df$idNo, file = "SHSph2_IDs_clean05152024CRM.csv"

## NOTE:
- I'm bringing in calculated metabolic data using Kaja's data output from Teddy's script referencing Laura's library, and am assuming reported values pass QC
- I intend on chaning the hardcoded bits into variables for scalability in a larger bioinformatic pipeline (dev in progress, 05/15/2024 CRM)
### This script: 
Calculated metabolite concentration data > Dataframe curation Script > Data Cleaning Script 
### Connecting it from: 
Raw experimental LC-MS/MS data > metabolite standard curve generation > concentation back-calculations
### Connecting it to: 
Data Analysis Script > Figures & Tables Script

## STILL TO DO:
1) write in a check that there's not data being overwritten for that sample's concentration variable (ie reruns in the same file), handle this how?
2) join with samples on my rerun list .xlsx
