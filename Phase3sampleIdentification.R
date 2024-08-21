# Phase3sampleIdentification
# Goal: identify which SHS phase 2 cohort "smoker" participant plasma samples had expected NMR and [nicotine metabolite] concentrations. 

# shs_deceased_phase3.xlsx # variable "PID" are ids listed have deceased before phase 3 cohort
# the clean NMR dataset the group is working with

##########       Setup       ##########
# Set working directory (the location with the necessary input files)
setwd("/Users/musserc/NicotineMetabolomics_project")

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
library(ggplot2)


##########       File Input       ##########
l_deceasedIDs <- read_excel("shs_deceased_phase3.xlsx")
df <- read.csv("SHSph2_NMRwSurveyDF_clean20240729CRM.csv",header=TRUE)

##########       File Output       ##########
# Set variables that will be used to automate filename generation
name_initials <- "CRM"
today <- format(Sys.Date(), "%Y%m%d")
# List of the base name for the output files
outfiles <- c("listIDs_samplesSHSph3")
# # Set the output filenames 
ph2samples_filename = paste0(outfiles[1], "_", today, name_initials, ".csv")
# full_filename = paste0(outfiles[2], "_", today, name_initials, ".csv")
# clean_filename = paste0(outfiles[3], today, name_initials, ".csv")
# cleanIDs_filename = paste0(outfiles[4], today, name_initials, ".csv")
# 
# ##########    Curate  Dataset     ##########
l_aliveIDs <- data.frame(idNo = integer())
# rename l_deceasedIDs$PID to be l_deceasedIDs$idNo
l_deceasedIDs <- l_deceasedIDs %>% rename(idNo = PID) %>% mutate(idNo = as.integer(idNo))
df <- df %>% mutate(idNo = as.integer(idNo))
# join datasets to discover who in my dataset is alive in phase 3
l_aliveIDs <- anti_join(df, l_deceasedIDs, by = "idNo") #select all those alive (n=731)

# now identify all (alive) smokers within 1 z-score of the mean NMR
selected_rows <- l_aliveIDs %>% filter(b_isControl == "samp") #select smokers (516)

# Calculate the mean and standard deviation of logNMR
mean_NMRcalc <- mean(selected_rows$NMRcalc, na.rm = TRUE)
mean_logNMR <- mean(selected_rows$logNMR, na.rm = TRUE)
sd_logNMR <- sd(selected_rows$logNMR, na.rm = TRUE)
half_sd_logNMR <- sd_logNMR/2

# Select rows within 1 z-score of the mean
final_selected_rows <- selected_rows %>%
  filter(logNMR >= (mean_logNMR - half_sd_logNMR) & logNMR <= (mean_logNMR + half_sd_logNMR)) # gives 200 smoker sample with log NMR between -0.4069 - -0.1798

# now verify that all smoker survey variables make sense
# years smoked
# Calculate the mean and standard deviation of number of years smoked
mean_s2smkd <- mean(final_selected_rows$S2SMKD, na.rm = TRUE) #mean=30 years
sd_s2smkd <- sd(final_selected_rows$S2SMKD, na.rm = TRUE) #sd= +/-14 years

final_selected_rows <- final_selected_rows %>%
  filter(S2SMKD >= (mean_s2smkd - sd_s2smkd) & S2SMKD <= (mean_s2smkd + sd_s2smkd)) #removed 72, to 128 obs.

# CPD
mean_cpd <- mean(final_selected_rows$cpd, na.rm = TRUE) #mean=15.4 CPD
sd_cpd <- sd(final_selected_rows$cpd, na.rm = TRUE) #sd= 11.5 CPD

final_selected_rows <- final_selected_rows %>%
  filter(cpd >= (mean_cpd - sd_cpd) & cpd <= (mean_cpd + sd_cpd)) #removed 29, still 99 obs.


# PPY
mean_ppy <- mean(final_selected_rows$S2PPY, na.rm = TRUE) #mean=24.4 packs
sd_ppy <- sd(final_selected_rows$S2PPY, na.rm = TRUE) #sd= +/-14.9 packs

final_selected_rows <- final_selected_rows %>%
  filter(S2PPY >= (mean_ppy - sd_ppy) & S2PPY <= (mean_ppy + sd_ppy)) #removed 29, to 70 obs.


#BMIcalc
mean_bmi <- mean(final_selected_rows$BMIcalc, na.rm = TRUE) #mean_bmi=28.98
sd_bmi <- sd(final_selected_rows$BMIcalc, na.rm = TRUE) #sd_bmi=4.72

final_selected_rows <- final_selected_rows %>%
  filter(BMIcalc >= (mean_bmi - sd_bmi) & BMIcalc <= (mean_bmi + sd_bmi)) # removed 19, to 51 obs.

# select ~50 IDs
print(final_selected_rows$idNo)

# Selected group Visualizations
#  final_selected_rows is data frame and NMRcalc is the column of interest
# NMR calc
ggplot(final_selected_rows, aes(x = NMRcalc)) +
  geom_histogram(binwidth = 0.025, fill = "blue", color = "black") +
  labs(title = "Histogram of NMRcalc", x = "NMRcalc", y = "Frequency")

#logNMR
ggplot(final_selected_rows, aes(x = logNMR)) +
  geom_histogram(binwidth = 0.025, fill = "red", color = "black") +
  labs(title = "Histogram of logNMR", x = "logNMR", y = "Frequency")

#CPD
ggplot(final_selected_rows, aes(x = cpd)) +
  geom_histogram(binwidth = 2, fill = "green", color = "black") +
  labs(title = "Histogram of CPD", x = "CPD", y = "Frequency")

#PPY
ggplot(final_selected_rows, aes(x = S2PPY)) +
  geom_histogram(binwidth = 2, fill = "orange", color = "black") +
  labs(title = "Histogram of PPY", x = "PPY", y = "Frequency")



######## Output ########
write.csv(final_selected_rows$idNo, file = ph2samples_filename, row.names = TRUE)
cat(paste0("Successfully output sample IDs in: ", ph2samples_filename, "\n"))
