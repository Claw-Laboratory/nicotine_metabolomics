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
# 1) A list of all combined results for all 824 samples; t_allSamples_postExclusions, file = "SHSph2_fullDF_05152024CRM.csv"
# 2) A "clean" dataset with 18 variables for the final 816 participants to use for EWAS/GWAS; clean_df, file = "SHSph2_NMRwSurveyDF_clean05152024CRM.csv"
# 3) A list of just sampleIDs of final 816 participants; clean_df$idNo, file = "SHSph2_IDs_clean05152024CRM.csv"

##########       Setup       ##########
setwd("/Users/musserc/NicotineMetabolomics_project")
if(!"dplyr" %in% installed.packages()) install.packages("dplyr")
if(!"lubridate" %in% installed.packages()) install.packages("lubridate")
if(!"readr" %in% installed.packages()) install.packages("readr")
if(!"stringr" %in% installed.packages()) install.packages("stringr")
if(!"lubridate" %in% installed.packages()) install.packages("lubridate")
library(readr) # Used Tidyverse:readr library to read multiple raw data .csv files into a single data frame
library(readxl) # Used Tidyverse:readxl library to raw survey data and data dictionary .xlsx files
library(stringr) #Used the stringr package to extract the middle 4-6 integer values from a string using regular expressions in R
library(dplyr) #removed duplicate rows from a data frame in R using the dplyr package, which provides the distinct function.
library(lubridate) # Used to calculate date-formatted time deltas

##########       File Input       ##########
t_calcConcData <- read.csv("v6RawData_combined.csv", header=TRUE) #per 'SampID' grab 'Analyte', 'Valid', 'CalcConc',
t_SHSph1_rawSurveyData  <- read_excel("phase1.xlsx") ## contains gender, birthdate, smoked >=100 in lifetime, ade started smoking, age stopped smoking, # years smoked
t_SHSph2_rawSurveyData  <- read_excel("phase2.xlsx") ##date plasma drawn,smoked now, smoked within 4 hours
t_ph1and2_rawSmokingSurveyData <- read_excel("smoking.xlsx") ## S1SMOKE, S1SMKD, S1PPY, S2SMOKE, S2SMKD, S2PPY
# this is Teddy's script for survey df: #t_otherSurveyData  <- read_excel("SHSII_HPlasma_ConfirmInventory_Out_2022.03.29_FULL_INVENTORY_withDemographics.xlsx")

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
    print(" not named within 6-digit SHS idenitifier parameters")
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
#be sure to include non-smoker controls that have NMR=0, set to 0 from na so not excluded in join methods
t_allSamples[t_allSamples$idNo=="103006","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="102041","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="102167","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="103672","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="202748","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="202056","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="202808","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="303153","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="102155","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="202348","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="203304","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="202542","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="202642","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="302971","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="202699","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="102507","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="102192","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="303222","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="102565","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="103388","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="202870","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="302063","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="202978","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="303527","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="302648","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="303139","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="303151","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="103696","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="303127","NMRcalc"] <- as.numeric(0) 
t_allSamples[t_allSamples$idNo=="303013","NMRcalc"] <- as.numeric(0) 

##########       Join Data Sets       ##########
# Join subsetted raw survey data as columns to sample ID rows
t_allSamples <- left_join(t_allSamples, t_SHSph1_rawSurveyData, by = c("idNo" = "IDNO"))  
t_allSamples <- left_join(t_allSamples, t_SHSph2_rawSurveyData, by = c("idNo" = "IDNO"))  
t_allSamples <- left_join(t_allSamples, t_ph1and2_rawSmokingSurveyData, by = c("idNo" = "IDNO"))  
### t_allSamples ### use t_allSamples as the master (all variables all samples dataframe, post-QC, after specific exclusions)

##########       Clean Variable Determination       ##########
# add supplementary columns for inferred data, interpreted using data dictionaries (Supplementary 1)
t_allSamples["b_SmokerStatus"] <- ""
t_allSamples["b_Gender"] <- ""
t_allSamples["AgeAtPhase2"] <-""
t_allSamples["AgeCategory"] <- ""
t_allSamples["BMIcalc"] <- ""
t_allSamples["BMIcategory"] <- ""
t_allSamples["b_isControl"] <- ""
t_allSamples["NMR_category"] <- ""
## "b_isControl"
t_allSamples[t_allSamples$idNo=="103006","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="102041","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="102167","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="103672","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="202748","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="202056","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="202808","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="303153","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="102155","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="202348","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="203304","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="202542","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="202642","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="302971","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="202699","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="102507","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="102192","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="303222","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="102565","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="103388","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="202870","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="302063","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="202978","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="303527","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="302648","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="303139","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="303151","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="103696","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="303127","b_isControl"] <- as.numeric(1) 
t_allSamples[t_allSamples$idNo=="303013","b_isControl"] <- as.numeric(1) 
## BMI variables: "BMIcalc" and "BMIcategory"
## According to CDC, BMI = weight [kg] / (height [m])^2 *10000[cm/m]
calculate_bmi <- function(weight_kg, height_cm) {
  bmi <- weight_kg / (height_cm^2)*10000 # Calculate BMI
  return(bmi) # Return the result
}
for (s in 1:nrow(t_allSamples)) {
  this_height_cm <- t_allSamples[s,"EX2_7"]
  this_weight_kg <- t_allSamples[s,"EX2_8"]
  calcBMI <- calculate_bmi(this_weight_kg,this_height_cm)
  t_allSamples[s,"BMIcalc"]= as.numeric(calcBMI)
}
for (s in 1:nrow(t_allSamples)) {
  this_BMI <- t_allSamples[s,"BMIcalc"]
  if (is.na(this_BMI)) next
  else if (this_BMI < 25.0) { t_allSamples[s,"BMIcategory"] <- "Underweight/Healthy"}
  else if (this_BMI >= 25.0 & this_BMI < 30.0 ) { t_allSamples[s,"BMIcategory"] <- "Overweight"}
  else if (this_BMI >= 30.0 ) { t_allSamples[s,"BMIcategory"] <- "Obese"}
  else {print("Error fitting BMI into categories")}
}
#"b_SmokerStatus"
## Determining smoker status per logic table made upon consensus (w/ KF, CS)
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
## Age variables: "AgeAtPhase2" and "AgeCategory
## "AgeAtPhase2": calculated age from Birthday INT2_3 (MMDDYY format) to survey exam date (Phase 2) (for calculating age on the day plasma drawn)
ex_NAage <- data.frame()
for (s in 1:nrow(t_allSamples)) {
  this_examDate <- t_allSamples[s,"S2EXDATE.x"]
  this_birthday <- t_allSamples[s, "INT2_3"]
  if (is.na(this_examDate)) {ex_NAage <- rbind(ex_NAage,s)}
  else {
    yrs_delta <- abs(interval(this_birthday,this_examDate)/years(1)) #uses lubridate package to calculate time interval in years
    t_allSamples[s,"AgeAtPhase2"] <- as.double(yrs_delta)
  }
}
##"AgeCategory": defined younger or older stratifications compared to overall average age
mean_age <- mean(as.double(t_allSamples$AgeAtPhase2), na.rm=TRUE)
for (s in 1:nrow(t_allSamples)) {
  this_age <- as.double(t_allSamples[s,"AgeAtPhase2"])
  this_ageBoolean <- (this_age <= mean_age)
  if (is.na(this_age)) {ex_NAage <- rbind(ex_NAage,s)} ## this doesn't make sense
  else if (this_ageBoolean == TRUE) {
    print(this_age)
    t_allSamples[s,"AgeCategory"] <- "Younger"
    }
  else t_allSamples[s,"AgeCategory"] <- "Older"
}
## "b_Gender"
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
excluded_smokers_naBMI <- subset(t_allSamples,is.na(t_allSamples$BMI_calc)) #BW controls + 2 smokers: 103154,103554
# missing age
excluded_smokers_naAge <- subset(t_allSamples,is.na(t_allSamples$AgeAtPhase2)) #none found
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
t_allSamples_postExclusions <- t_allSamples
t_allSamples_postExclusions <- subset(t_allSamples_postExclusions, !(t_allSamples_postExclusions$idNo=="202642"))
t_allSamples_postExclusions <- subset(t_allSamples_postExclusions, !(t_allSamples_postExclusions$idNo=="202699"))
t_allSamples_postExclusions <- subset(t_allSamples_postExclusions, !(t_allSamples_postExclusions$idNo=="103154"))
t_allSamples_postExclusions <- subset(t_allSamples_postExclusions, !(t_allSamples_postExclusions$idNo=="103554")) #824 filtered down to 820
t_allSamples_postExclusions <- subset(t_allSamples_postExclusions, !(t_allSamples_postExclusions$idNo=="2720"))
t_allSamples_postExclusions <- subset(t_allSamples_postExclusions, !(t_allSamples_postExclusions$idNo=="3348"))
t_allSamples_postExclusions <- subset(t_allSamples_postExclusions, !(t_allSamples_postExclusions$idNo=="3955"))
t_allSamples_postExclusions <- subset(t_allSamples_postExclusions, !(t_allSamples_postExclusions$idNo=="5829")) #down to 816 without BW controls

##########       Output       ##########
# Outputs: NMR dataframe, list of excluded ids for requesting more plasma to re-run,A list of excluded sample ids to request SHS for more plasma to re-run assay
allExcluded_list <- anti_join(t_allSamples,t_allSamples_postExclusions,by = c("idNo" = "idNo")) # results in 8
# TODO: join with samples on my rerun list .xlsx
write.csv(allExcluded_list, file = "listIDs_excludedSHSph2_05152024CRM.csv", row.names = TRUE)
print("Successfully output all excluded sample IDs list as listIDs_excludedSHSph2_05152024CRM.csv")
# Output: A list of all combined results for all 824 samples; t_allSamples_postExclusions, file = "SHSph2_fullDF_05152024CRM.csv"
write.csv(t_allSamples_postExclusions, file = "SHSph2_fullDF_05152024CRM.csv", row.names = TRUE)
print("Successfully output t_allSamples_postExclusions as SHSph2_fullDF_05152024CRM.csv")
# Output: A "clean" dataset for EWAS/GWAS; clean_df, file = "SHSph2_NMRwSurveyDF_clean05152024CRM.csv"
cleanVars <- c("idNo","NMRcalc","logNMR","COTconc","tHCconc","NICconc",    # data from LC-MS and NMR calcs
               "CENTER.y","S2EXDATE.x","S2SMOKE","S2SMKD","S2PPY",         # data from raw survey variables
               "b_SmokerStatus","b_Gender","AgeAtPhase2","AgeCategory","BMIcalc","BMIcategory","b_isControl")          # interpreted data from survey, see logic in UML diagram (Supplementary Figure 2)
clean_df <-t_allSamples_postExclusions[,cleanVars]
write.csv(clean_df, file = "SHSph2_NMRwSurveyDF_clean05152024CRM.csv", row.names = TRUE)
print("Successfully output clean_df as SHSph2_NMRwSurveyDF_clean05152024CRM.csv")
# Output: just sampleIDs of 816 participants; clean_df$idNo, file = "SHSph2_IDs_clean05152024CRM.csv"
write.csv(clean_df$idNo, file = "SHSph2_IDs_clean05152024CRM.csv", row.names = TRUE)
print("Successfully output sample IDs for clean_df as SHSph2_IDs_clean05152024CRM.csv")