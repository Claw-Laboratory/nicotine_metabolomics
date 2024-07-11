## NMR Assay Dataframe Curation
# Crystal R. Musser, M.S. (Wet Lab Manager, Katrina Claw Lab @ CU Anschutz, Department of Biomedical Informatics)
# v.05/16/2024

# PURPOSE: 
# Combine quantified  LC/MS nicotine metabolite data with phase survey self-reported data, from specified cohort.

# OBJECTS:
# "t_allSamplesMetabolites" for curated nicotine metabolite dataset
# "t_allSamples" as the master (all variables all samples dataframe, post-QC, after specific exclusions)
#  "clean_df" as the finalized, shared dataset, 18 variables, n=816 (inc. 28 controls)

# SUPPLEMENTARY FILES to this script: 
# Data dictionaries (Supplementary File 1)
# "Clean" variable consensus/logic UML (Supplementary File 2)

# REQUIRED INPUTS: 4
# 1) combined experimental LC-MS/MS calculated and QC'd data in [ng/mL] serum concentrations for: , Cotinine, 3-hyroxy-Cotinine, Nicotine, CO, NO, Nornicotine (in this order)      --->   "v6RawData_combined.csv"
# 2, 3) cohort survey results with associated data dictionary for interpretation                                                                                                 ---> "phase2.xlsx" and "smoking.xlsx"
# 4) phase I survey results with associated data dictionary for interpretating age, gender, smoking histotry, height, weight                                                  ---> "phase1.xlsx"

# OUTPUT FILES GENERATED: 
# 1) A list of all the sampleIDs excluded from the final 816; allExcluded_list, file = "listIDs_excludedSHSph2"
# 2) A list of all combined results for all 824 samples; t_allSamples_postExclusions, file = "SHSph2_fullDF"
# 3) A "clean" dataset with 18 variables for the final 816 participants to use for EWAS/GWAS; clean_df, file = "SHSph2_NMRwSurveyDF_clean"
# 4) A list of just sampleIDs of final 816 participants; clean_df$idNo, file = "SHSph2_IDs_clean"

##########       Setup       ##########
# Set working directory (the location with the necessary input files)
setwd("/home/clawlab/Projects/SHS")

# Install packages if not installed
if(!"dplyr" %in% installed.packages()) install.packages("dplyr")
if(!"lubridate" %in% installed.packages()) install.packages("lubridate")
if(!"readr" %in% installed.packages()) install.packages("readr")
if(!"stringr" %in% installed.packages()) install.packages("stringr")
if(!"haven" %in% installed.packages()) install.packages("haven")
# Load libraries
library(readr) # Used Tidyverse:readr library to read multiple raw data .csv files into a single data frame
library(readxl) # Used Tidyverse:readxl library to raw survey data and data dictionary .xlsx files
library(stringr) #Used the stringr package to extract the middle 4-6 integer values from a string using regular expressions in R
library(dplyr) #removed duplicate rows from a data frame in R using the dplyr package, which provides the distinct function.
library(lubridate) # Used to calculate date-formatted time deltas
library(haven)

##########       File Input       ##########
t_calcConcData <- read.csv("v6RawData_combined.csv", header=TRUE) #per 'SampID' grab 'Analyte', 'Valid', 'CalcConc',
t_SHSph1_rawSurveyData  <- read_excel("phase1.xlsx") ## contains gender, birthdate, smoked >=100 in lifetime, ade started smoking, age stopped smoking, # years smoked
t_SHSph2_rawSurveyData  <- read_excel("phase2.xlsx") ##date plasma drawn,smoked now, smoked within 4 hours
t_ph1and2_rawSmokingSurveyData <- read_excel("smoking.xlsx") ## S1SMOKE, S1SMKD, S1PPY, S2SMOKE, S2SMKD, S2PPY
SHSph1_age <- read_sas("shs677_elig.sas7bdat") # S1AGE
# this is Teddy's script for survey df: #t_otherSurveyData  <- read_excel("SHSII_HPlasma_ConfirmInventory_Out_2022.03.29_FULL_INVENTORY_withDemographics.xlsx")

##########       File Output       ##########
# Set variables that will be used to automate filename generation
name_initials <- "KRF"
today <- format(Sys.Date(), "%d%m%Y")
# List of the base name for the output files
outfiles <- c("listIDs_excludedSHSph2", "SHSph2_fullDF", 
              "SHSph2_NMRwSurveyDF_clean", "SHSph2_IDs_clean")
# Set the output filenames 
excl_filename = paste0(outfiles[1], "_", today, name_initials, ".csv")
full_filename = paste0(outfiles[2], "_", today, name_initials, ".csv")
clean_filename = paste0(outfiles[3], today, name_initials, ".csv")
cleanIDs_filename = paste0(outfiles[4], today, name_initials, ".csv")

##########    Curate NMR Dataset     ##########
t_allSamplesMetabolites <- data.frame( # new df, indexed by by idNo, with columns:
  idNo = numeric(),
  logNMR = numeric(),
  NMRcalc = numeric(),
  ExpSet = character(),
  COTconc = numeric(),
  tHCconc = numeric(),
  NICconc = numeric(),
  COconc = numeric(),
  NOconc = numeric(),
  NORNconc = numeric()
)

# Normalize naming convention
t_funkNames <- data.frame()
for (each in 1:nrow(t_calcConcData)) { #5327
  extracted_digits <- str_extract(t_calcConcData$SampleID[each], "\\d{4,6}")  # extract sample id using regex:str_extract to extract the middle 4-6 digits
  if (is.na(extracted_digits)) {
    this_id <- each
    t_funkNames <- rbind(t_funkNames,this_id)
    #print(" not named within 6-digit SHS idenitifier parameters")
  }
  else {                      
    t_allSamplesMetabolites <- rbind(t_allSamplesMetabolites, list(
      idNo = as.numeric(extracted_digits),
      ExpSet = NA,
      NMRcalc = NA,
      COTconc = NA,
      tHCconc = NA,
      NICconc = NA,
      COconc = NA,
      NOconc = NA,
      NORNconc = NA)    )
  }
# print(paste(t_calcConcData$SampleID[each], "extracted: ", extracted_digits)) # Print the result. 
}
# Remove duplicates by keeping the first occurrence of each set of duplicates in the data frame, and remove the subsequent ones.
t_allSamplesMetabolites <- t_allSamplesMetabolites[!duplicated(t_allSamplesMetabolites), ]

# Ensure BW controls are the only difference between survey data and LC/MS data? -> survey doesn't have BW controls: #3955,5829,2720,3348
diff_dfs <- anti_join(t_allSamplesMetabolites, t_SHSph2_rawSurveyData, by = c("idNo" = "IDNO"))

# Add experimental date sets by labeling each "Set" appropriately in the raw LC/MS-MS dataset with exp. dates associated to origin file naming convention.
t_calcConcData[1:324,"Set"] <- "A"      #"OUTPUT_12.13.22.csv"
t_calcConcData[325:1386,"Set"] <- "B"   #"OUTPUT_01.27.23.csv"
t_calcConcData[1387:1908,"Set"] <- "C"  #"OUTPUT_02.14.23.csv"
t_calcConcData[1909:2430,"Set"] <- "D"  #"OUTPUT_2.28.23.csv"
t_calcConcData[2431:2970,"Set"] <- "E"  #"OUTPUT_3.06.23.csv"
t_calcConcData[2971:3492,"Set"] <- "F"  # "OUTPUT_3.14.23.csv"
t_calcConcData[3493:4014,"Set"] <- "G"  # "OUTPUT_3.17.23.csv"
t_calcConcData[4015:4536,"Set"] <- "H"  # "OUTPUT_3.23.23.csv"
t_calcConcData[4537:5058,"Set"] <- "I"  # "OUTPUT_4.11.23.csv"
t_calcConcData[5059:5532,"Set"] <- "J"  # "OUTPUT_4.18.23.csv"

#TODO: write in a check that there's not data being overwritten for that sample's concentration variable (ie reruns in the same file), handle how?

# Add data from LC/MS-MS raw analyte concentration data
# For each of the valid samples to be anlayzed, find that sample's object, set it's validity flag to true, set the sample's experimentally determined  plasma metabolite concentration values, and add to a working object list for analysis
#  Check if there is any index in the vector of names from l_validSamples that matches the name of the current sample (foundSamp_inList[[1]]@name). If there is a match, it means the sample is already in the list, and the condition will be TRUE. Otherwise, it will be FALSE.
for (each in 1:nrow(t_calcConcData)) { #5327
  extracted_digits <- str_extract(t_calcConcData$SampleID[each],"\\d{4,6}")     #extract sample id using regex:str_extract to extract the middle 4-6 digits
  this_Analyte <- t_calcConcData[each, "Analyte"]
  this_AnalyteValue <- t_calcConcData[each,"CalcConc"]
  this_ValidAnalyte <- t_calcConcData[each,"Valid"]
  this_Set <- t_calcConcData[each,"Set"]
  if (this_ValidAnalyte==FALSE){
    print(paste("extracted",extracted_digits,"\t",this_Analyte,"\t",this_AnalyteValue,"\t",this_ValidAnalyte)) # Print the result.
  }else{
    if (extracted_digits %in% t_allSamplesMetabolites$idNo) { # Check if sample already exists in list 
      t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"ExpSet"] <- this_Set
      if (!this_ValidAnalyte == TRUE){  # Update analyte value with CalcConc
        if (this_Analyte == "Cotinine") {                      # set [COT]calc
          cotinineValue <- t_calcConcData[each,"CalcConc"]
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"COTconc"] = as.numeric(cotinineValue)
        } else if (this_Analyte == "3-Hydroxycotinine") {      # set [3HC]calc
          tHCValue <- t_calcConcData[each,"CalcConc"]
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"tHCconc"] = as.numeric(tHCValue)
        } else if (this_Analyte == "Nicotine") {               # set [NIC]calc
          nicValue <- t_calcConcData[each,"CalcConc"]
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"NICconc"] = as.numeric(nicValue)
        } else if (this_Analyte == "Nicotine-1-oxide") {               # set [NO]calc
          noValue <- t_calcConcData[each,"CalcConc"]
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"NOconc"] = as.numeric(noValue)
        } else if (this_Analyte == "Cotinine-N-oxide") {               # set [CO]calc
          coValue <- t_calcConcData[each,"CalcConc"]
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"COconc"] = as.numeric(coValue)
        } else if (this_Analyte == "Nornicotine") {               # set [NORN]calc
          nornValue <- t_calcConcData[each,"CalcConc"]
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"NORNconc"] = as.numeric(nornValue)
        }else{ print("Valid analyte not accounted for")}
      }else{ # if analyte value is false, then set pertinent variable to zero.
        if (this_Analyte == "Cotinine") {                      # set [COT]calc
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"COTconc"] = 0
        } else if (this_Analyte == "3-Hydroxycotinine") {      # set [3HC]calc
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"3HCconc"] = 0
        } else if (this_Analyte == "Nicotine") {               # set [NIC]calc
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"NICconc"] = 0
        } else if (this_Analyte == "Nicotine-1-oxide") {               # set [NO]calc
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"NOconc"] = 0
        } else if (this_Analyte == "Cotinine-N-oxide") {               # set [CO]calc
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"COconc"] = 0
        } else if (this_Analyte == "Nornicotine") {               # set [NORN]calc
          t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == extracted_digits,"NORNconc"] = 0
        }else{ print("Valid analyte not accounted for")}
        print("analyte is invalid")}      # or? Update anlayte value with "0"}
    }else{
      #print(paste(extracted_digits,"IF",extracted_digits %in% t_allSamplesMetabolites$idNo))
    }
  }
}

# Calculate NMR
for (s in 1:nrow(t_allSamplesMetabolites)) {
  cotConc <- t_allSamplesMetabolites[s,"COTconc"]
  thcConc <- t_allSamplesMetabolites[s,"tHCconc"]
  calculatedNMR <-thcConc/cotConc
  temp_logNMR = log10(calculatedNMR) #logNMR= log(base 10) of calculated NMR value
  tempID <- t_allSamplesMetabolites$idNo[s]
  t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == tempID,"NMRcalc"] = as.numeric(calculatedNMR)
  t_allSamplesMetabolites[t_allSamplesMetabolites$idNo == tempID,"logNMR"] = as.numeric(temp_logNMR)
} #### use t_allSamplesMetabolites for curated nicotine metabolite dataset ###

t_allSamples <- t_allSamplesMetabolites # n=824 #instantiate master data frame with metabolic data as base architecture

# List of all control (non-smoker) samples
control_samples <- c("103006", "102041", "102167", "103672", "202748",
                    "202056", "202808", "303153", "102155", "202348",
                    "203304", "202542", "202642", "302971", "202699",
                    "102507", "102192", "303222", "102565", "103388",
                    "202870", "302063", "202978", "303527", "302648",
                    "303139", "303151", "103696", "303127", "303013"
                    )

# Recode NAs for control sample measurements to 0
for (ctr in control_samples){
  t_allSamples[t_allSamples$idNo==ctr, "NMRcalc"] <- as.numeric(0)
}

##########       Join Data Sets       ##########
# Join subsetted raw survey data as columns to sample ID rows
t_allSamples <- left_join(t_allSamples, t_SHSph1_rawSurveyData, by = c("idNo" = "IDNO"))  
t_allSamples <- left_join(t_allSamples, t_SHSph2_rawSurveyData, by = c("idNo" = "IDNO"))  
t_allSamples <- left_join(t_allSamples, t_ph1and2_rawSmokingSurveyData, by = c("idNo" = "IDNO"))  
t_allSamples <- left_join(t_allSamples, SHSph1_age, by = c("idNo"="IDNO"))
### t_allSamples ### use t_allSamples as the master (all variables all samples dataframe, post-QC, after specific exclusions)

##########       Clean Variable Determination       ##########
# add supplementary columns for inferred data, interpreted using data dictionaries (Supplementary 1)
t_allSamples["b_SmokerStatus"] <- ""
t_allSamples["b_Gender"] <- ""
t_allSamples["BMIcalc"] <- ""
t_allSamples["BMIcategory"] <- ""
t_allSamples["b_isControl"] <- ""
t_allSamples["NMR_category"] <- ""
t_allSamples["S2AGE"] <- ""

## "b_isControl"
###############################################################
t_allSamples <- t_allSamples  %>% 
                  mutate(b_isControl = ifelse(idNo %in% control_samples, 1, 0))

## BMI variables: "BMIcalc" and "BMIcategory"
## According to CDC, BMI = weight [kg] / (height [m])^2 *10000[cm/m]
# Function to calculate BMI
calculate_bmi <- function(weight_kg, height_cm) {
  bmi <- as.numeric(weight_kg) / ((as.numeric(height_cm)^2)*10000) # Calculate BMI
  return(bmi) # Return the result
}

# Calculate BMI for all samples
t_allSamples <- t_allSamples  %>% 
                  mutate(BMIcalc = calculate_bmi(EX2_7, EX2_8))

# Categorize continuous BMI values
t_allSamples <- t_allSamples %>%
  mutate(
    BMIcategory = case_when(
      is.na(BMIcalc) ~ NA_character_,
      BMIcalc < 25.0 ~ "Underweight/Healthy",
      BMIcalc >= 25 & BMIcalc < 30 ~ "Overweight",
      BMIcalc >= 30 ~ "Obese"
    )
  )

#"b_SmokerStatus"
###############################################################
## Determining smoker status per logic table made upon consensus (w/ KF, CS)
t_allSamples <- t_allSamples  %>% 
                  mutate(selfDeclaredSmoker = INT22_4,
                        cpd = INT22_5,
                        b_SmokerStatus = case_when(
                          is.na(selfDeclaredSmoker) & is.na(cpd) ~ "exclude",
                          is.na(selfDeclaredSmoker) & cpd > 0 ~ "Smoker",
                          is.na(selfDeclaredSmoker) & cpd == 0 ~ "Non-smoker",
                          selfDeclaredSmoker==1 ~ "Smoker",
                          selfDeclaredSmoker==2 & is.na(cpd) ~ "Non-smoker",
                          selfDeclaredSmoker==2 & cpd == 0 ~ "Non-smoker",
                          selfDeclaredSmoker==2 & cpd == 99 ~ "Non-smoker",
                          selfDeclaredSmoker==2 & cpd > 0 ~ "Smoker"
                        ))

for (s in 1:nrow(t_allSamples)) {
  this_selfDeclared_smokerStatus <- t_allSamples[s,"INT22_4"]
  this_CPD <- t_allSamples[s,"INT22_5"]
#  this_CPYr <- t_allSamples[s,"Smoke.More.Than.100.Cigs.P2"]
  if (any(is.na(this_selfDeclared_smokerStatus))) {
    if (is.na(this_CPD)){
      t_allSamples[s,"b_SmokerStatus"] <- "exclude"
    }else if (this_CPD>0) {
      t_allSamples[s,"b_SmokerStatus"] <- "Smoker"
    } else {t_allSamples[s,"b_SmokerStatus"] <- "Non-smoker"}
  }
  else if (this_selfDeclared_smokerStatus==1) {t_allSamples[s,"b_SmokerStatus"] <- "Smoker"} 
  else if (this_selfDeclared_smokerStatus==2) {
    if (is.na(this_CPD)){
      t_allSamples[s,"b_SmokerStatus"] <- "Non-smoker"
    }else if (this_CPD>0) {
        t_allSamples[s,"b_SmokerStatus"] <- "Smoker"
    } else {t_allSamples[s,"b_SmokerStatus"] <- "Non-smoker"} #this_CPD == 0, this_CPD == 99, ,  this_CPYr=="Never"
  }
}

## "S2age"
###############################################################
t_allSamples <- t_allSamples  %>% 
                  mutate(S1EXDATE = as.Date(S1EXDATE.x),
                        S2EXDATE = as.Date(S2EXDATE.x),
                        delta_time_yrs = as.numeric(difftime(S2EXDATE, S1EXDATE, units="days"))/ 365.25,
                        S2AGE = round(S1AGE + delta_time_yrs, digits=1))



## "b_Gender"
###############################################################
## from phase 1 survey, gender (1=male, 2=female) is "INT2_1"
ex_NAgender <- data.frame()
for (s in 1:nrow(t_allSamples)) {
  this_gender <- t_allSamples$INT2_1[s]
  if (is.na(this_gender)) rbind(ex_NAgender,s) # running list of participants with undisclosed gender to be excluded from analysis data set
  else if (this_gender==1) {t_allSamples$b_Gender[s] = "male"} #1=male
  else if (this_gender==2) {t_allSamples$b_Gender[s] = "female"} #2=female
  else print("Error processing gender for ",t_allSamples$idNo[s])
}

##Calculated summary of gender:mean NMRs
t_allSamples %>%
  group_by(b_Gender) %>%
  summarise(gender_avgNMR = mean(NMRcalc,na.rm = TRUE),
    gender_stdDedvNMR = sd(NMRcalc, na.rm=TRUE))

#NMR metabolizer categories:NMRcategory can be high or low,
# high: "logNMR" >= -0.31
# low: "logNMR" < -0.31
ex_NAnmr <- data.frame()
for (s in 1:nrow(t_allSamples)) {
  this_logNMR <- t_allSamples[s,"logNMR"]; print(this_logNMR)
  this_id <- str(t_allSamples$idNo[s])
  if (is.na(this_logNMR)) {
    ex_NAnmr <- rbind(ex_NAnmr,s)
    t_allSamples$NMR_category[s]="na?"}
    # seperate as possible controls, otherwise potentially exclude smoker sample?
  else if (this_logNMR>= -0.31) {t_allSamples$NMR_category[s]="high"}
  else t_allSamples$NMR_category[s]="low"
}

##########       Specific Exclusions      ##########
## Exclude samples with missing data for BMI, age, gender, or CPD information since complete data set for these variables is needed for GWAS/EWAS analysis
# missing BMI
excluded_smokers_naBMI <- subset(t_allSamples,is.na(t_allSamples$BMIcalc)) #BW controls + 2 smokers: 103154,103554
# missing age
excluded_smokers_naAge <- subset(t_allSamples,is.na(t_allSamples$S2AGE)) # 2 smokers: 203292, 303493
# missing gender
excluded_smokers_naGender <- subset(t_allSamples,is.na(t_allSamples$b_Gender)) #none found

# Define QC NMR values and smoking behaviors to check;
sub_smokers <- subset(t_allSamples,t_allSamples$b_SmokerStatus=="Smoker") # smokers subset
sub_nonSmokers <- subset(t_allSamples, t_allSamples$b_SmokerStatus=="Non-smoker")  # non-smokers subset; pulls n=30 (all the controls)
# did NOT exclude smokers with na NMR (they all had [3HC]<LLOQ)
sub_smokers_naNMR <- subset(sub_smokers,is.na(sub_smokers$NMRcalc)) # there are 57 smokers with NMRcalc==na
# excluded non-smoker controls with NMR>0
excluded_nonSmokers_wNMR <- subset(sub_nonSmokers, (sub_nonSmokers$logNMR<0)) ##idNo="202642", ""202699"
# did NOT exclude smokers with missing CPD
sub_smokers_naCPD <- subset(sub_smokers,is.na(sub_smokers$INT22_5)) # 13 found, not excluding for now

# Specified Exclusions Made
sample_exclusions <- c("202642", "202699", "103154", "103554",
                        "2720", "3348", "3955", "5829")
t_allSamples_postExclusions <- t_allSamples
t_allSamples_postExclusions <- t_allSamples_postExclusions[!t_allSamples_postExclusions$idNo %in% sample_exclusions,] 

##########       Output       ##########
# Outputs: NMR dataframe, list of excluded ids for requesting more plasma to re-run,A list of excluded sample ids to request SHS for more plasma to re-run assay
allExcluded_list <- anti_join(t_allSamples,t_allSamples_postExclusions,by = c("idNo" = "idNo")) # results in 8
# TODO: join with samples on my rerun list .xlsx
write.csv(allExcluded_list, file = excl_filename, row.names = TRUE)
cat("Successfully output all excluded sample IDs list as: ", excl_filename)
# Output: A list of all combined results for all 824 samples; t_allSamples_postExclusions, file = "SHSph2_fullDF_05152024CRM.csv"
write.csv(t_allSamples_postExclusions, file = full_filename, row.names = TRUE)
cat("Successfully output t_allSamples_postExclusions as: " full_filename)
# Output: A "clean" dataset for EWAS/GWAS; clean_df, file = "SHSph2_NMRwSurveyDF_clean05152024CRM.csv"        # interpreted data from survey, see logic in UML diagram (Supplementary Figure 2)
clean_df <-t_allSamples_postExclusions  %>% 
            dplyr::select(idNo, NMRcalc, logNMR, COTconc, tHCconc, NICconc,  # data from LC-MS and NMR calcs
            CENTER.y, S2EXDATE, S2SMOKE, S2SMKD, S2PPY, b_SmokerStatus, # data from raw survey variables
            S2AGE, b_Gender, BMIcalc, BMIcategory, b_isControl)  %>% 
            rename("CENTER.y" = "CENTER") # remove the '.y' for convenience
write.csv(clean_df, file = clean_filename, row.names = TRUE)
cat("Successfully output clean_df as: ", clean_filename)
# Output: just sampleIDs of 816 participants; clean_df$idNo, file = "SHSph2_IDs_clean05152024CRM.csv"
write.csv(clean_df$idNo, file = cleanIDs_filename, row.names = TRUE)
cat("Successfully output sample IDs for clean_df as: ", cleanIDs_filename)