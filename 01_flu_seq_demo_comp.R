#################################################################################
# data compilation influenza - lauring lab
# Leigh Papalambros
# 2025-03-02

library(tidyverse)
library(lubridate)
library(stringr)
library(dplyr)
library(readxl)
library(janitor)


place1 <- "/Users/leighbak/University of Michigan Dropbox/MED-LauringLab/"
start_loc <- paste0(place1, "External_Projects_DataRequests/INFLUENZA/")

###############################################
#################################
### full compiled flu data (Lauring Lab data)
iav_list <- read.csv(paste0(place1, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character", na.strings = c("NA", "", " "))
#as of 05/08/2026 21399
# as of 07/2026 21495 (end of '25 season)


ibv_list <- read.csv(paste0(place1, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/FinalSummary/full_compiled_data.csv"), colClasses = "character", na.strings = c("NA", "", " "))
# as of 05/08/2026 3801
# as of 07/2026 3989 (end of '25 season)


# filtering lab data for columns of interest
iav_list_small <- iav_list %>% select(subject_id, sample_id, coll_date, sample_per_subject, nextclade_HA_clade, nextclade_HA_completeness, nextclade_HA_totalMissing,
                                      nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs,
                                      nextclade_HA_runDate, nextclade_HA_type, genbank_SubmissionID,
                                      genbank_HAH1, genbank_HAH3, genbank_MP, genbank_NAN1, genbank_NAN2, genbank_NP, genbank_NS, genbank_PA, genbank_PB1, genbank_PB2,
                                      Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, 
                                      NS.Segment_Id, HE.Segment_Id, P3.Segment_Id, received_source)


# filtering lab data for columns of interest
ibv_list_small <- ibv_list %>% select(subject_id, sample_id, coll_date, sample_per_subject, nextclade_HA_clade, nextclade_HA_completeness, nextclade_HA_totalMissing,
                                      nextclade_HA_qcOverallScore, nextclade_HA_qcOverallStatus, nextclade_HA_totalMutations, nextclade_HA_totalNonACGTNs,
                                      nextclade_HA_runDate, nextclade_HA_type, genbank_SubmissionID,
                                      genbank_HA, genbank_MP, genbank_NA, genbank_NP, genbank_NS, genbank_PA, genbank_PB1, genbank_PB2,
                                      Isolate_Id, PB2.Segment_Id, PB1.Segment_Id, PA.Segment_Id, HA.Segment_Id, NP.Segment_Id, NA.Segment_Id, MP.Segment_Id, 
                                      NS.Segment_Id, received_source)




# combining the Lauring Lab IAV and IBV data sets into one compiled file for creating sub data sets
# and the final compiled data sets

iab_small <- bind_rows(iav_list_small, ibv_list_small)
# as of 07/2026 25484

write.csv(iab_small, paste0(start_loc, "Big_Data_Set/processed_data/lauring_lab_sequencing_data_ab_20261007.csv"), row.names = FALSE, na = "")

################################################################################
### ascension, henry ford, michigan medicine, trinity patient data that is sent to us
## file paths
asj_fp_a <- paste0(place1, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/Manifests/ASCENSION/ASCENSION_MDSS/")
hfhs_fp_a <- paste0(place1, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/Manifests/HENRYFORD/HFHS_MDSS/")
mm_fp_a <- paste0(start_loc, "MISAPPHIRE_FLU_RETURN/")
asj_fp_b <- paste0(place1, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/Manifests/ASCENSION/ASJ_MDSS/")
hfhs_fp_b <- paste0(place1, "SEQUENCING/INFLUENZA_B/4_SequenceSampleMetadata/Manifests/HENRYFORD/HENRYFORD_MDSS/")



## reading in the patient data files and formatting them to be able to create a master file of patient data
## INFLUENZA A

asj_storage_a <- data.frame()
hfhs_storage_a <- data.frame()
mm_storage_a <- data.frame()
full_a <- data.frame()
full_a_geo <- data.frame()
mm_geo_a <- data.frame()
asj_geo_a <- data.frame()
hfhs_geo_a <- data.frame()
mm_name_dob_id_a <- data.frame()
asj_name_dob_id_a <- data.frame()
hfhs_name_dob_id_a <- data.frame()

asj_file_lis_a <- list.files(pattern = "*.csv", path = asj_fp_a)
hfhs_file_list_a <- list.files(pattern = "*.csv", path = hfhs_fp_a)
mm_file_list_a <- list.files(pattern = "*.csv", path = mm_fp_a)

all_a_file_list <- c(asj_file_lis_a, hfhs_file_list_a, mm_file_list_a)

for (each_file in all_a_file_list){
  if (grepl("ASJ", each_file)){
    fileone <- read.csv(paste0(asj_fp_a, "/", each_file), colClasses = "character")
    # calculate age at time of specimen collection
    file_read_in <- fileone %>% mutate(SPECIMEN_COLLECTION_DATE_format = paste0(substr(SPECIMEN_COLLECTION_DATE, 1, 4), "-", substr(SPECIMEN_COLLECTION_DATE, 5, 6), "-", substr(SPECIMEN_COLLECTION_DATE, 7, 8)), 
                                       PATIENT_DOB_format = paste0(substr(PATIENT_DOB, 1, 4), "-", substr(PATIENT_DOB, 5, 6), "-", substr(PATIENT_DOB, 7, 8)))
    file_read_in$dob_orig <- file_read_in$PATIENT_DOB
    
    #fixing the wrong years in a coll_date
    file_read_in$SPECIMEN_COLLECTION_DATE_format <- sub("^2525", "2025", file_read_in$SPECIMEN_COLLECTION_DATE_format)
    
    file_read_in <- file_read_in %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SPECIMEN_COLLECTION_DATE_format), as_date(PATIENT_DOB_format), units = "days"))/365.25))
    
    name_dob_id <- file_read_in %>% select(PATIENT_MRN, ACCESSION_ID, SID, PATIENT_LNAME, PATIENT_FNAME, dob_orig)
    
    
    asj_name_dob_id_a <- rbind(name_dob_id, asj_name_dob_id_a)
    
    file_read_in_set <- file_read_in %>% select(PATIENT_ZIP, PATIENT_COUNTY, SPECIMEN_COLLECTION_DATE_format, 
                                                age_years_at_time_of_specimen_collection, 
                                                PATIENT_RACE, PATIENT_ETHNICITY, PATIENT_SEX, PATIENT_MRN, ACCESSION_ID, SID, PATIENT_LNAME, PATIENT_FNAME, PATIENT_DOB_format) %>% distinct()
    asj_storage_a <- rbind(file_read_in_set, asj_storage_a) 
    # as of 08/2026 966
    geo_readin_asj <- file_read_in %>% select(PATIENT_MRN, ACCESSION_ID, SID, PATIENT_ADDR1, PATIENT_ADDRE2, PATIENT_CITY, PATIENT_STATE, PATIENT_ZIP, PATIENT_COUNTY)
    #asj_geo_a <- data.frame()
    colnames(geo_readin_asj)[colnames(geo_readin_asj) == 'PATIENT_ADDRE2'] <- 'PATIENT_ADDR2'
    asj_geo_a <- rbind(asj_geo_a, geo_readin_asj)
    asj_geo_a <- asj_geo_a[!apply(asj_geo_a == "", 1, all),]
    # as of 08/2026 946
  } else if (grepl("HFHS", each_file)){
    fileone <- read.csv(paste0(hfhs_fp_a, "/", each_file), colClasses = "character")
    # calculate age at time of specimen collection
    file_read_in <- fileone %>% mutate(SPECIMEN_COLLECTION_DATE_format = paste0(substr(SPECIMEN_COLLECTION_DATE, 1, 4), "-", substr(SPECIMEN_COLLECTION_DATE, 5, 6), "-", substr(SPECIMEN_COLLECTION_DATE, 7, 8)), 
                                       PATIENT_DOB_format = paste0(substr(PATIENT_DOB, 1, 4), "-", substr(PATIENT_DOB, 5, 6), "-", substr(PATIENT_DOB, 7, 8)))
    file_read_in$dob_orig <- file_read_in$PATIENT_DOB
    
    file_read_in <- file_read_in %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SPECIMEN_COLLECTION_DATE_format), as_date(PATIENT_DOB_format), units = "days"))/365.25))
    
    name_dob_id <- file_read_in %>% select(PATIENT_MRN, ACCESSION_ID, SID, PATIENT_LNAME, PATIENT_FNAME, dob_orig)
    
    
    hfhs_name_dob_id_a <- rbind(name_dob_id, hfhs_name_dob_id_a)
    
    file_read_in_set <- file_read_in %>% select(PATIENT_ZIP, PATIENT_COUNTY, SPECIMEN_COLLECTION_DATE_format, 
                                                age_years_at_time_of_specimen_collection, 
                                                PATIENT_RACE, PATIENT_ETHNICITY, PATIENT_SEX, PATIENT_MRN, ACCESSION_ID, SID, PATIENT_LNAME, PATIENT_FNAME, PATIENT_DOB_format) %>% distinct()
    
    hfhs_storage_a <- rbind(file_read_in_set, hfhs_storage_a)
    # as of 08/2026 5720
    geo_readin_hfhs <- file_read_in %>% select(PATIENT_MRN, ACCESSION_ID, SID, PATIENT_ADDR1, PATIENT_ADDR2, PATIENT_CITY, PATIENT_STATE, PATIENT_ZIP, PATIENT_COUNTY)
    #hfhs_geo_a <- data.frame()
    hfhs_geo_a <- rbind(hfhs_geo_a, geo_readin_hfhs)
    # as of 08/2026 5744
  } else {
    fileone <- read.csv(paste0(mm_fp_a, "/", each_file), colClasses = "character", stringsAsFactors = FALSE)
    
    #fileone$PATIENT_DOB <- as.character(as.Date(mdy_hms(fileone$PATIENT_DOB)))
    fileone <- separate(fileone, PATIENT_DOB, c("PATIENT_DOB_date", "PATIENT_DOB_time"), sep = " ")
    fileone$dob_orig <- fileone$PATIENT_DOB_date
    #fileone$PATIENT_DOB_date <- as.POSIXct(fileone$PATIENT_DOB_date, format = "%m/%d/%Y")
    fileone$PATIENT_DOB_date <- as.Date(fileone$PATIENT_DOB_date, format = "%m/%d/%y")
    #fileone$PATIENT_DOB_date <- format(fileone$PATIENT_DOB_date, "%Y-%m-%d")
    
    fileone <- separate(fileone, SampleCollectionDate, c("SampleCollectionDate_date", "SampleCollectionDate_time"), sep = " ")
    fileone$SampleCollectionDate_date <- as.POSIXct(fileone$SampleCollectionDate_date, format = "%m/%d/%Y")
    #fileone$SampleCollectionDate_date <- as.Date(fileone$SampleCollectionDate_date, format = "%m/%d/%y")
    #fileone$SampleCollectionDate_date <- format(fileone$PSampleCollectionDate_date, "%Y-%m-%d")
    
    ## need to correct DOB and to calculate the correct age
    fileone$PATIENT_DOB_date_corrected <- as.Date(ifelse(fileone$PATIENT_DOB_date > Sys.Date(), format(fileone$PATIENT_DOB_date, "19%y-%m-%d"), format(fileone$PATIENT_DOB_date)))
    
    # calculate age at time of specimen collection
    file_read_in <- fileone %>% mutate(SPECIMEN_COLLECTION_DATE_format = paste0(substr(SampleCollectionDate_date, 1, 4), "-", substr(SampleCollectionDate_date, 6, 7), "-", substr(SampleCollectionDate_date, 9, 10)), 
                                       PATIENT_DOB_format = paste0(substr(PATIENT_DOB_date_corrected, 1, 4), "-", substr(PATIENT_DOB_date_corrected, 6, 7), "-", substr(PATIENT_DOB_date_corrected, 9, 10)))
    #correcting the collection dates
    file_read_in$SampleCollectionDate_date <- sub("^00", "20", file_read_in$SampleCollectionDate_date)
    needs_fix <- grepl("^\\d{2}-\\d{2}-\\d{2}$", file_read_in$SampleCollectionDate_date)
    
    # Fix only those rows
    file_read_in$SampleCollectionDate_date[needs_fix] <- sub("^(\\d{2})-", "20\\1-", file_read_in$SampleCollectionDate_date[needs_fix])
    
    
    file_read_in <- file_read_in %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SampleCollectionDate_date), as_date(PATIENT_DOB_format), units = "days"))/365.25))
    
    name_dob_id <- file_read_in %>% select(MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, dob_orig)
    
    
    mm_name_dob_id_a <- rbind(name_dob_id, mm_name_dob_id_a)
    
    
    
    file_read_in_set_2 <- file_read_in %>% select(PATIENT_ZIP, PATIENT_COUNTY, 
                                                  age_years_at_time_of_specimen_collection, 
                                                  PATIENT_RACE, ETHNIC_GROUP_DESC, PATIENT_SEX, MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, PATIENT_DOB_format, SampleCollectionDate_date) %>% distinct()
    
    colnames(file_read_in_set_2)[colnames(file_read_in_set_2) == 'ETHNIC_GROUP_DESC'] <- 'PATIENT_ETHNICITY'
    colnames(file_read_in_set_2)[colnames(file_read_in_set_2) == 'MRN'] <- 'PATIENT_MRN'
    colnames(file_read_in_set_2)[colnames(file_read_in_set_2) == 'SampleID'] <- 'ACCESSION_ID'
    colnames(file_read_in_set_2)[colnames(file_read_in_set_2) == 'SampleCollectionDate_date'] <- 'SPECIMEN_COLLECTION_DATE_format'
    
    mm_storage_a <- rbind(file_read_in_set_2, mm_storage_a)
    # as of 08/2026 1163
    geo_readin_mm1 <- file_read_in %>% select(SampleID, MRN, PATIENT_ADDR1, PATIENT_ADDR2, PATIENT_CITY, PATIENT_STATE, PATIENT_ZIP, PATIENT_COUNTY)
    colnames(geo_readin_mm1)[colnames(geo_readin_mm1) == 'SampleID'] <- 'ACCESSION_ID'
    colnames(geo_readin_mm1)[colnames(geo_readin_mm1) == 'MRN'] <- 'PATIENT_MRN'
    #mm_geo_a <- data.frame()
    mm_geo_a <- rbind(mm_geo_a, geo_readin_mm1)
    # as of 08/2026 1163
  } 
  
  full_a <- bind_rows(asj_storage_a, hfhs_storage_a, mm_storage_a)
  # as of 08/2026 7849
  full_a_geo <- bind_rows(mm_geo_a, hfhs_geo_a, asj_geo_a)
  # as of 08/2026 7853
  
}


colnames(mm_name_dob_id_a)[colnames(mm_name_dob_id_a) == 'MRN'] <- 'PATIENT_MRN'
colnames(mm_name_dob_id_a)[colnames(mm_name_dob_id_a) == 'SampleID'] <- 'ACCESSION_ID'

## processing the old23/24 data in the old format separately because it is being weird for michigan medicine
### This is also data for IAV and IBV combined (not split out into respective sub folders 06/16/2026)

flu_23 <- read.csv(paste0(mm_fp_a, "/Old/Flu_Patients_10192023.csv"), colClasses = "character", stringsAsFactors = FALSE)

#formatting dates
flu_23 <- separate(flu_23, PATIENT_DOB, c("PATIENT_DOB_date", "PATIENT_DOB_time"), sep = " ")
flu_23$dob_orig <- flu_23$PATIENT_DOB_date
flu_23$PATIENT_DOB_date <- as.Date(flu_23$PATIENT_DOB_date, format = "%m/%d/%y")
flu_23$PATIENT_DOB_date_corrected <- as.Date(ifelse(flu_23$PATIENT_DOB_date > Sys.Date(), format(flu_23$PATIENT_DOB_date, "19%y-%m-%d"), format(flu_23$PATIENT_DOB_date)))

#flu_23 <- separate(flu_23, SampleCollectionDate, c("SampleCollectionDate_date", "SampleCollectionDate_time"), sep = " ")
flu_23$SampleCollectionDate_date <- as.Date(flu_23$SampleCollectionDate, format = "%m/%d/%y")
#flu_23$SampleCollectionDate_date <- format(flu_23$SampleCollectionDate, "%Y-%m-%d")
#str(flu_2425_season)

flu_23 <- flu_23 %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SampleCollectionDate_date), as_date(PATIENT_DOB_date_corrected), units = "days"))/365.25))

colnames(flu_23)[colnames(flu_23) == 'SampleCollectionDate_date'] <- 'SPECIMEN_COLLECTION_DATE_format' 
colnames(flu_23)[colnames(flu_23) == 'PATIENT_DOB_date_corrected'] <- 'PATIENT_DOB_format' 


flu_23_final_demo <- flu_23 %>% select(PATIENT_ZIP, PATIENT_COUNTY, SPECIMEN_COLLECTION_DATE_format, 
                                       age_years_at_time_of_specimen_collection, 
                                       PATIENT_RACE, ETHNIC_GROUP_DESC, PATIENT_SEX, PATIENT_MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, PATIENT_DOB_format) %>% distinct()

colnames(flu_23_final_demo)[colnames(flu_23_final_demo) == 'SampleID'] <- 'ACCESSION_ID'
colnames(flu_23_final_demo)[colnames(flu_23_final_demo) == 'ETHNIC_GROUP_DESC'] <- 'PATIENT_ETHNICITY'

#creating a data set for a compiled key that is Names, ID's, and DOB column(s) ######
flu_23_names_dob_id <- flu_23 %>% select(PATIENT_MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, dob_orig)
colnames(flu_23_names_dob_id)[colnames(flu_23_names_dob_id) == 'SampleID'] <- 'ACCESSION_ID'

s23_geo <- flu_23 %>% select(SampleID, PATIENT_MRN, PATIENT_ADDR1, PATIENT_ADDR2, PATIENT_CITY, PATIENT_STATE, PATIENT_ZIP, PATIENT_COUNTY)
colnames(s23_geo)[colnames(s23_geo) == 'SampleID'] <- 'ACCESSION_ID'

## Data for the start of the 23/24 season
### This is also data for IAV and IBV combined (not split out into respective sub folders 06/16/2026)

flu_24 <- read.csv(paste0(mm_fp_a, "/Old/Flu_Patients_20240118.csv"), colClasses = "character", stringsAsFactors = FALSE)

#formatting dates
flu_24 <- separate(flu_24, PATIENT_DOB, c("PATIENT_DOB_date", "PATIENT_DOB_time"), sep = " ")
flu_24$dob_orig <- flu_24$PATIENT_DOB_date
flu_24$PATIENT_DOB_date <- as.Date(flu_24$PATIENT_DOB_date, format = "%m/%d/%y")
flu_24$PATIENT_DOB_date_corrected <- as.Date(ifelse(flu_24$PATIENT_DOB_date > Sys.Date(), format(flu_24$PATIENT_DOB_date, "19%y-%m-%d"), format(flu_24$PATIENT_DOB_date)))

flu_24 <- separate(flu_24, SampleCollectionDate, c("SampleCollectionDate_date", "SampleCollectionDate_time"), sep = " ")
flu_24$SampleCollectionDate_date <- as.Date(flu_24$SampleCollectionDate_date, format = "%m/%d/%y")
flu_24$SampleCollectionDate_date_c <- format(flu_24$SampleCollectionDate_date, "%Y-%m-%d")
#str(flu_24)

flu_24 <- flu_24 %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SampleCollectionDate_date), as_date(PATIENT_DOB_date_corrected), units = "days"))/365.25))

colnames(flu_24)[colnames(flu_24) == 'SampleCollectionDate_date'] <- 'SPECIMEN_COLLECTION_DATE_format' 
colnames(flu_24)[colnames(flu_24) == 'PATIENT_DOB_date_corrected'] <- 'PATIENT_DOB_format' 


flu_24_final_demo <- flu_24 %>% select(PATIENT_ZIP, PATIENT_COUNTY, SPECIMEN_COLLECTION_DATE_format, 
                                       age_years_at_time_of_specimen_collection, 
                                       PATIENT_RACE, ETHNIC_GROUP_DESC, PATIENT_SEX, PATIENT_MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, PATIENT_DOB_format) %>% distinct()

colnames(flu_24_final_demo)[colnames(flu_24_final_demo) == 'SampleID'] <- 'ACCESSION_ID'
colnames(flu_24_final_demo)[colnames(flu_24_final_demo) == 'ETHNIC_GROUP_DESC'] <- 'PATIENT_ETHNICITY'

#creating a data set for a compiled key that is Names, ID's, and DOB column(s) ######
flu_24_names_dob_id <- flu_24 %>% select(PATIENT_MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, dob_orig)
colnames(flu_24_names_dob_id)[colnames(flu_24_names_dob_id) == 'SampleID'] <- 'ACCESSION_ID'

s24_geo <- flu_24 %>% select(SampleID, PATIENT_MRN, PATIENT_ADDR1, PATIENT_ADDR2, PATIENT_CITY, PATIENT_STATE, PATIENT_ZIP, PATIENT_COUNTY)
colnames(s24_geo)[colnames(s24_geo) == 'SampleID'] <- 'ACCESSION_ID'

## processing the new 2025 data seperately because it is being weird for michigan medicine
### This is also data for IAV and IBV combined (not split out into respective sub folders 06/16/2026)

flu_2425_season <- read.csv(paste0(mm_fp_a, "/New/data_20250905_returned.csv"), colClasses = "character", stringsAsFactors = FALSE)


#formatting dates
flu_2425_season <- separate(flu_2425_season, PATIENT_DOB, c("PATIENT_DOB_date", "PATIENT_DOB_time"), sep = " ")
flu_2425_season$dob_orig <- flu_2425_season$PATIENT_DOB_date
flu_2425_season$PATIENT_DOB_date <- as.Date(flu_2425_season$PATIENT_DOB_date, format = "%m/%d/%y")
flu_2425_season$PATIENT_DOB_date_corrected <- as.Date(ifelse(flu_2425_season$PATIENT_DOB_date > Sys.Date(), format(flu_2425_season$PATIENT_DOB_date, "19%y-%m-%d"), format(flu_2425_season$PATIENT_DOB_date)))

flu_2425_season <- separate(flu_2425_season, SampleCollectionDate, c("SampleCollectionDate_date", "SampleCollectionDate_time"), sep = " ")
flu_2425_season$SampleCollectionDate_date <- as.Date(flu_2425_season$SampleCollectionDate_date, format = "%m/%d/%y")
#flu_2425_season$SampleCollectionDate_date_c <- format(flu_2425_season$SampleCollectionDate_date, "%Y-%m-%d")
#str(flu_2425_season)

flu_2425_season <- flu_2425_season %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SampleCollectionDate_date), as_date(PATIENT_DOB_date_corrected), units = "days"))/365.25))

colnames(flu_2425_season)[colnames(flu_2425_season) == 'SampleCollectionDate_date'] <- 'SPECIMEN_COLLECTION_DATE_format' 
colnames(flu_2425_season)[colnames(flu_2425_season) == 'PAT_MRN'] <- 'PATIENT_MRN' 
colnames(flu_2425_season)[colnames(flu_2425_season) == 'PATIENT_DOB_date_corrected'] <- 'PATIENT_DOB_format' 
colnames(flu_2425_season)[colnames(flu_2425_season) == 'ETHNIC_GROUP_DESC'] <- 'PATIENT_ETHNICITY' 
#colnames(flu_2425_season)[colnames(flu_2425_season) == 'SampleCollectionDate_date'] <- 'SampleCollectionDate'


s2425_final_demo <- flu_2425_season %>% select(PATIENT_ZIP, PATIENT_COUNTY, SPECIMEN_COLLECTION_DATE_format, 
                                               age_years_at_time_of_specimen_collection, 
                                               PATIENT_RACE, PATIENT_ETHNICITY, PATIENT_SEX, PATIENT_MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, PATIENT_DOB_format) %>% distinct()

colnames(s2425_final_demo)[colnames(s2425_final_demo) == 'SampleID'] <- 'ACCESSION_ID'

#creating a data set for a compiled key that is Names, ID's, and DOB column(s) ######
s2425_names_dob_id <- flu_2425_season %>% select(PATIENT_MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, dob_orig)
colnames(s2425_names_dob_id)[colnames(s2425_names_dob_id) == 'SampleID'] <- 'ACCESSION_ID'

s2425_geo <- flu_2425_season %>% select(SampleID, PATIENT_MRN, PATIENT_ADDR1, PATIENT_ADDR2, PATIENT_CITY, PATIENT_STATE, PATIENT_ZIP, PATIENT_COUNTY)
colnames(s2425_geo)[colnames(s2425_geo) == 'SampleID'] <- 'ACCESSION_ID'

str(flu_24_final_demo)

########## the 25/26 season data returned from MM as well as old data from the outbreak ##########
flu_2526_season <- read.csv(paste0(start_loc, "/jules_flu_compilation/25_26_data/Returned_data_influenza_20260526.csv"), colClasses = "character", stringsAsFactors = FALSE)
ob_flu <- read.csv(paste0(start_loc, "/jules_flu_compilation/25_26_data/returned_add_data_influenza_06012026.csv"), colClasses = "character", stringsAsFactors = FALSE)

ob_2526_flu <- bind_rows(flu_2526_season, ob_flu)
# as of 07/2025 5076

## formatting and calculating age_at_specimen_collection
ob_2526_flu <- separate(ob_2526_flu, PATIENT_DOB, c("PATIENT_DOB_format", "PATIENT_DOB_time"), sep = " ")
ob_2526_flu$dob_orig <- ob_2526_flu$PATIENT_DOB_format
ob_2526_flu$SPECIMEN_COLLECTION_DATE_format <- as.Date(ob_2526_flu$SampleCollectionDate, format = "%m/%d/%Y")
ob_2526_flu <- ob_2526_flu %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SPECIMEN_COLLECTION_DATE_format), as_date(PATIENT_DOB_format), units = "days"))/365.25))

ff_ob_2526_flu <- ob_2526_flu %>% select(MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, PATIENT_ZIP, PATIENT_COUNTY, PATIENT_DOB_format, PATIENT_SEX, PATIENT_RACE, ETHNIC_GROUP_DESC, age_years_at_time_of_specimen_collection, SPECIMEN_COLLECTION_DATE_format)
colnames(ff_ob_2526_flu)[colnames(ff_ob_2526_flu) == 'SampleID'] <- 'ACCESSION_ID'
colnames(ff_ob_2526_flu)[colnames(ff_ob_2526_flu) == 'MRN'] <- 'PATIENT_MRN'
colnames(ff_ob_2526_flu)[colnames(ff_ob_2526_flu) == 'ETHNIC_GROUP_DESC'] <- 'PATIENT_ETHNICITY'

#creating a data set for a compiled key that is Names, ID's, and DOB column(s) ######
ob_2526_names_dob_id <- ob_2526_flu %>% select(MRN, SampleID, PATIENT_LNAME, PATIENT_FNAME, dob_orig)
colnames(ob_2526_names_dob_id)[colnames(ob_2526_names_dob_id) == 'SampleID'] <- 'ACCESSION_ID'
colnames(ob_2526_names_dob_id)[colnames(ob_2526_names_dob_id) == 'MRN'] <- 'PATIENT_MRN'


ff_geo_ob_25_26 <- ob_2526_flu %>% select(MRN, SampleID, PATIENT_ADDR1, PATIENT_ADDR2, PATIENT_CITY, PATIENT_STATE, PATIENT_ZIP, PATIENT_COUNTY)
colnames(ff_geo_ob_25_26)[colnames(ff_geo_ob_25_26) == 'SampleID'] <- 'ACCESSION_ID'
colnames(ff_geo_ob_25_26)[colnames(ff_geo_ob_25_26) == 'MRN'] <- 'PATIENT_MRN'

## this section is for processing the manifests from MDHHS that have some of the patient data that we need
### This is also data for IAV and IBV combined (not split out into respective sub folders 06/16/2026)

mdhhs_manifest_fp <- paste0(place1, "SEQUENCING/INFLUENZA_A/4_SequenceSampleMetadata/Manifests/MDHHS/Original")
mdhhs_files_24 <- list.files(mdhhs_manifest_fp, pattern = "*.xlsx")

flu_mdhhs_manifest <- data.frame()
manifest_mdhhs_bd <- data.frame()

for (each_mdhhs in mdhhs_files_24){
  m_in <- read_xlsx(paste0(mdhhs_manifest_fp, "/", each_mdhhs), sheet = 1)#, detectDates = TRUE) #
  
  # put full manifest listing from mdhhs into full list
  full_mdhhs <- rbind(flu_mdhhs_manifest, m_in)
  # we'll write this out in full later, we'll need it for MDSS submissions
  
  # select only columns we care about
  m_in_bd <- m_in %>% select(`Lab ID`, `Last Name`, `First Name`,  DOB, `Patient Cnty`, Collected)
  
  # we need to create a set with: 
  # position, sample_id, subject_id, coll_date, flag, received_date, received_source
  m_in_bd <- m_in_bd %>% mutate(SPECIMEN_COLLECTION_DATE_format = paste0(substr(Collected, 1, 4), "-", substr(Collected, 6, 7), "-", substr(Collected, 9, 10)))
  m_in_bd <- m_in_bd %>% mutate(PATIENT_DOB_format = paste0(substr(DOB, 1, 4), "-", substr(DOB, 6, 7), "-", substr(DOB, 9, 10)))
  
  # add in 2 new columns: received_date and received_source (from file name)
  #m_in$received_date <- date_from_file(each_mdhhs)
  rec_date <- trimws(as.character(strsplit(each_mdhhs, "_")[[1]][1]))
  rec_date <- as.Date(rec_date, "%m%d%Y")
  m_in_bd$received_date <- rec_date
  
  
  rec_source <- trimws(as.character(strsplit(each_mdhhs, "_")[[1]][2]))
  m_in_bd$received_source <- rec_source
  
  m_in_bd <- m_in_bd %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SPECIMEN_COLLECTION_DATE_format), as_date(PATIENT_DOB_format), units = "days"))/365.25))
  
  manifest_mdhhs_bd <- rbind(manifest_mdhhs_bd, m_in_bd)
  
  ff_mdhhs_demo <- manifest_mdhhs_bd %>% select(`Lab ID`, `Last Name`, `First Name`, age_years_at_time_of_specimen_collection, SPECIMEN_COLLECTION_DATE_format, PATIENT_DOB_format)
  ff_mdhhs_geo <- manifest_mdhhs_bd %>% select(`Lab ID`, `Patient Cnty`)
}

# renaming columns for binding with the other demo and geo file sources
colnames(ff_mdhhs_demo)[colnames(ff_mdhhs_demo) == 'Lab ID'] <- 'ACCESSION_ID'
colnames(ff_mdhhs_demo)[colnames(ff_mdhhs_demo) == 'Last Name'] <- 'PATIENT_LNAME'
colnames(ff_mdhhs_demo)[colnames(ff_mdhhs_demo) == 'First Name'] <- 'PATIENT_FNAME'

colnames(ff_mdhhs_geo)[colnames(ff_mdhhs_geo) == 'Lab ID'] <- 'ACCESSION_ID'
colnames(ff_mdhhs_geo)[colnames(ff_mdhhs_geo) == 'Patient Cnty'] <- 'PATIENT_COUNTY'


s2425_final_demo <- s2425_final_demo %>%
  mutate(
    SPECIMEN_COLLECTION_DATE_format = as.character(SPECIMEN_COLLECTION_DATE_format),
    PATIENT_DOB_format = as.character(PATIENT_DOB_format)
  )

flu_24_final_demo <- flu_24_final_demo %>%
  mutate(
    SPECIMEN_COLLECTION_DATE_format = as.character(SPECIMEN_COLLECTION_DATE_format),
    PATIENT_DOB_format = as.character(PATIENT_DOB_format)
  )

flu_23_final_demo <- flu_23_final_demo %>%
  mutate(
    SPECIMEN_COLLECTION_DATE_format = as.character(SPECIMEN_COLLECTION_DATE_format),
    PATIENT_DOB_format = as.character(PATIENT_DOB_format)
  )

ff_ob_2526_flu <- ff_ob_2526_flu %>%
  mutate(
    SPECIMEN_COLLECTION_DATE_format = as.character(SPECIMEN_COLLECTION_DATE_format),
    PATIENT_DOB_format = as.character(PATIENT_DOB_format)
  )
### Creating the big final data sets
full_a_demo_all <- bind_rows(asj_storage_a, hfhs_storage_a, mm_storage_a, s2425_final_demo, flu_24_final_demo, flu_23_final_demo, ff_ob_2526_flu)
# as of 05/05/2026 11833
# as of 07/2026 16909

### double check all of the DOBS in the formatted column against the orginal formatted column ###
## creating a full compiled list of names, ids, and original DOB from the patient sample information
# asj_name_dob_id_a hfhs_name_dob_id_a mm_name_dob_id_a flu_23_names_dob_id flu_24_names_dob_id s2425_names_dob_id ob_2526_names_dob_id

all_a_names_dob_id <- bind_rows(asj_name_dob_id_a, hfhs_name_dob_id_a, mm_name_dob_id_a, flu_23_names_dob_id, flu_24_names_dob_id, s2425_names_dob_id, ob_2526_names_dob_id)
# as of 08/2026 35595

# need to remove the rows that are all blank in every column
clean_all_a_names_dob_id <- all_a_names_dob_id %>% 
  mutate(across(everything(), ~ na_if(., ""))) %>% 
  remove_empty(which = "rows")
# as of 08/2026 16913

# now creating a new column that is the DOB all formatted the same way
# this code converts dates that are like this in one column 1/1/44, 1/1/1944, 19440101, and 1941-01-01 (character column to start)
# and converts them all to 1941-01-01 (character column)
clean_dates_all_a_names_dob_id <- clean_all_a_names_dob_id %>%
  mutate(
    # 1. Add "yjmd" to the orders to catch unseparated formats like 20000411
    parsed_dob = as.Date(parse_date_time(dob_orig, orders = c("mdy", "ymd", "yjmd"))), # date format column
    
    # 2. Fix the 2-digit year century roll-over issue
    fixed_dob = if_else(parsed_dob > Sys.Date(), parsed_dob - years(100), parsed_dob), # date format column
    
    # 3. Format the final output as a character string (YYYY-MM-DD)
    clean_dob = format(fixed_dob, "%Y-%m-%d") # character format column
  )

# as of 08/2026 16913

# adding in mdhhs that have the missing subject_id
full_a_demo <- bind_rows(full_a_demo_all, ff_mdhhs_demo)
# as of 07/2026 17095

full_a_geo_all <- bind_rows(asj_geo_a, hfhs_geo_a, mm_geo_a, s2425_geo, s24_geo, s23_geo, ff_geo_ob_25_26)
# as of 05/05/2026 11837
# as of 07/2026 16913

## get rid of the trailing and leading whitespaces in the SID column
full_a_geo$SID <- trimws(full_a_geo$SID)
## seeing if there are any odd characters in the SID column (TRUE = Odd character)
full_a_geo$has_odd_chars <- str_detect(full_a_geo$SID, "[^a-zA-Z0-9\\s]")

## problem solving date formatting
wrong_form_date <- full_a_demo %>% filter(grepl("--", SPECIMEN_COLLECTION_DATE_format)) #as of 08/2026 21
wrong_form_date_2 <- full_a_demo %>% filter(grepl("--", PATIENT_DOB_format)) #as of 08/2026 45

## need to filter out blank rows (the date rows have "--" in it that I can remove based on)

full_a_demo_final <- filter(full_a_demo, !grepl("--", SPECIMEN_COLLECTION_DATE_format))
full_a_demo_final <- filter(full_a_demo_final, !grepl("--", PATIENT_DOB_format))

# pulling the one sample that has the wrong DOB that can't be corrected (their partial phone number is listed as the DOB)
full_a_demo_final <- full_a_demo_final %>%
  filter(PATIENT_MRN != "L008317696")
# as of 05/05/2026 11809 after the filtering
# as of 07/2026 16885 after the filtering

## data quality checks for ages to make sure that dates are correct
old_age_p <- full_a_demo_final %>% filter(age_years_at_time_of_specimen_collection >= 125)
neg_age_p <- full_a_demo_final %>% filter(age_years_at_time_of_specimen_collection < 0)

# one person (MRN:009978835) has the wrong year of DOB it is suppose to be 1925 not 2025 (from CBR returned data 09/05/2025)
full_a_demo_final <- full_a_demo_final %>%
  mutate(PATIENT_DOB_format = if_else(PATIENT_MRN == "009978835", "1925-12-02", PATIENT_DOB_format))


full_a_demo_final <- full_a_demo_final %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SPECIMEN_COLLECTION_DATE_format), as_date(PATIENT_DOB_format), units = "days"))/365.25))

# looking at how many people have multiple samples because they will have potentially different ages at date of sample collection
num_of_samples <- full_a_demo_final %>% group_by(PATIENT_MRN) %>% summarize(count = n())

only_multi_samps <- merge(full_a_demo_final, num_of_samples, by = "PATIENT_MRN", all.x = TRUE) %>% filter(count >= 2)
# there are 858 instances of a person have 2 or more collection dates as of 5/8/2026
# there are 1663 instances of a person have 2 or more collection dates as of 5/8/2026

# looking at how many people have multiple samples based on sample_id
num_of_samples_2 <- full_a_demo_final %>% group_by(ACCESSION_ID) %>% summarize(count = n())

only_multi_samps_2 <- merge(full_a_demo_final, num_of_samples_2, by = "ACCESSION_ID", all.x = TRUE) %>% filter(count >= 2)
# as of 08/2026 826

################################################################################
## INFLUEZNA B ##
asj_storage_b <- data.frame()
hfhs_storage_b <- data.frame()
full_b <- data.frame()
full_b_geo <- data.frame()
#mm_geo_b <- data.frame()
asj_geo_b <- data.frame()
hfhs_geo_b <- data.frame()
asj_name_dob_id_b <- data.frame()
hfhs_name_dob_id_b <- data.frame()


asj_file_lis_b <- list.files(pattern = "*.csv", path = asj_fp_b)
hfhs_file_list_b <- list.files(pattern = "*.csv", path = hfhs_fp_b)

all_b_file_list <- c(asj_file_lis_b, hfhs_file_list_b)

for (each_file in all_b_file_list){
  if (grepl("ASJ", each_file)){
    fileone_2 <- read.csv(paste0(asj_fp_b, "/", each_file), colClasses = "character")
    # calculate age at time of specimen collection
    file_read_in_2 <- fileone_2 %>% mutate(SPECIMEN_COLLECTION_DATE_format = paste0(substr(SPECIMEN_COLLECTION_DATE, 1, 4), "-", substr(SPECIMEN_COLLECTION_DATE, 5, 6), "-", substr(SPECIMEN_COLLECTION_DATE, 7, 8)), 
                                           PATIENT_DOB_format = paste0(substr(PATIENT_DOB, 1, 4), "-", substr(PATIENT_DOB, 5, 6), "-", substr(PATIENT_DOB, 7, 8)))
    file_read_in_2$dob_orig <- file_read_in_2$PATIENT_DOB
    #fixing the wrong years in a coll_date
    file_read_in_2$SPECIMEN_COLLECTION_DATE_format <- sub("^2525", "2025", file_read_in_2$SPECIMEN_COLLECTION_DATE_format)
    
    file_read_in_2 <- file_read_in_2 %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SPECIMEN_COLLECTION_DATE_format), as_date(PATIENT_DOB_format), units = "days"))/365.25))
    
    name_dob_id <- file_read_in_2 %>% select(PATIENT_MRN, ACCESSION_ID, SID, PATIENT_LNAME, PATIENT_FNAME, dob_orig)
    
    asj_name_dob_id_b <- rbind(name_dob_id, asj_name_dob_id_b)
    
    file_read_in_set_2 <- file_read_in_2 %>% select(PATIENT_ZIP, PATIENT_COUNTY, SPECIMEN_COLLECTION_DATE_format, 
                                                    age_years_at_time_of_specimen_collection, 
                                                    PATIENT_RACE, PATIENT_ETHNICITY, PATIENT_SEX, PATIENT_MRN, ACCESSION_ID, SID, PATIENT_LNAME, PATIENT_FNAME, PATIENT_DOB_format, SID) %>% distinct()
    asj_storage_b <- rbind(file_read_in_set_2, asj_storage_b) 
    
    
    geo_readin_asj_b <- file_read_in_2 %>% select(PATIENT_MRN, ACCESSION_ID, SID, PATIENT_ADDR1, PATIENT_ADDRE2, PATIENT_CITY, PATIENT_STATE, PATIENT_ZIP, PATIENT_COUNTY)
    #asj_geo_a <- data.frame()
    colnames(geo_readin_asj_b)[colnames(geo_readin_asj_b) == 'PATIENT_ADDRE2'] <- 'PATIENT_ADDR2'
    asj_geo_b <- rbind(asj_geo_b, geo_readin_asj_b)
    asj_geo_b <- asj_geo_b[!apply(asj_geo_b == "", 1, all),]
    
  }else{
    fileone_2 <- read.csv(paste0(hfhs_fp_b, "/", each_file), colClasses = "character")
    # calculate age at time of specimen collection
    file_read_in_2 <- fileone_2 %>% mutate(SPECIMEN_COLLECTION_DATE_format = paste0(substr(SPECIMEN_COLLECTION_DATE, 1, 4), "-", substr(SPECIMEN_COLLECTION_DATE, 5, 6), "-", substr(SPECIMEN_COLLECTION_DATE, 7, 8)), 
                                           PATIENT_DOB_format = paste0(substr(PATIENT_DOB, 1, 4), "-", substr(PATIENT_DOB, 5, 6), "-", substr(PATIENT_DOB, 7, 8)))
    file_read_in_2$dob_orig <- file_read_in_2$PATIENT_DOB
    
    file_read_in_2 <- file_read_in_2 %>% mutate(age_years_at_time_of_specimen_collection = floor(as.numeric(difftime(as_date(SPECIMEN_COLLECTION_DATE_format), as_date(PATIENT_DOB_format), units = "days"))/365.25))
    
    name_dob_id <- file_read_in_2 %>% select(PATIENT_MRN, ACCESSION_ID, SID, PATIENT_LNAME, PATIENT_FNAME, dob_orig)
    
    hfhs_name_dob_id_b <- rbind(name_dob_id, hfhs_name_dob_id_b)
    
    file_read_in_set_2 <- file_read_in_2 %>% select(PATIENT_ZIP, PATIENT_COUNTY, SPECIMEN_COLLECTION_DATE_format, 
                                                    age_years_at_time_of_specimen_collection, 
                                                    PATIENT_RACE, PATIENT_ETHNICITY, PATIENT_SEX, PATIENT_MRN, ACCESSION_ID, SID, PATIENT_LNAME, PATIENT_FNAME, PATIENT_DOB_format, SID) %>% distinct()
    
    hfhs_storage_b <- rbind(file_read_in_set_2, hfhs_storage_b)
    
    geo_readin_hfhs_b <- file_read_in_2 %>% select(PATIENT_MRN, ACCESSION_ID, SID, PATIENT_ADDR1, PATIENT_ADDR2, PATIENT_CITY, PATIENT_STATE, PATIENT_ZIP, PATIENT_COUNTY)
    #hfhs_geo_a <- data.frame()
    hfhs_geo_b <- rbind(hfhs_geo_b, geo_readin_hfhs_b)
  }
  
  full_b <- bind_rows(asj_storage_b, hfhs_storage_b)
  full_b_geo <- bind_rows(hfhs_geo_b, asj_geo_b)
}



str(asj_storage_b)
str(hfhs_storage_b)
## need to filter out blank rows (the date rows have "--" in it that I can remove based on)

full_b <- filter(full_b, !grepl("--", SPECIMEN_COLLECTION_DATE_format))

#b_2425_season <- full_b %>% filter(SPECIMEN_COLLECTION_DATE_format >= as.Date ("2024-07-01"))

## get rid of the trailing and leading whitespaces in the SID column
full_b_geo$SID <- trimws(full_b_geo$SID)
full_b$SID <- trimws(full_b$SID)


## seeing if there are any odd characters in the SID column (TRUE = Odd character)
full_b_geo$has_odd_chars <- str_detect(full_b_geo$SID, "[^a-zA-Z0-9\\s]")

## data quality checks for ages to make sure that dates are correct
old_age_p_b <- full_b %>% filter(age_years_at_time_of_specimen_collection >= 125)
neg_age_p_b <- full_b %>% filter(age_years_at_time_of_specimen_collection < 0)

### double check all of the DOBS in the formatted column against the orginal formatted column ###
## creating a full compiled list of names, ids, and original DOB from the patient sample information
# asj_name_dob_id_a hfhs_name_dob_id_a mm_name_dob_id_a flu_23_names_dob_id flu_24_names_dob_id s2425_names_dob_id ob_2526_names_dob_id

all_b_names_dob_id <- bind_rows(asj_name_dob_id_b, hfhs_name_dob_id_a)
# as of 08/2026 43708

# need to remove the rows that are all blank in every column
clean_all_b_names_dob_id <- all_b_names_dob_id %>% 
  mutate(across(everything(), ~ na_if(., ""))) %>% 
  remove_empty(which = "rows")
# as of 08/2026 5897

# now creating a new column that is the DOB all formatted the same way
# this code converts dates that are like this in one column 1/1/44, 1/1/1944, 19440101, and 1941-01-01 (character column to start)
# and converts them all to 1941-01-01 (character column)
clean_dates_all_b_names_dob_id <- clean_all_b_names_dob_id %>%
  mutate(
    # 1. Add "yjmd" to the orders to catch unseparated formats like 20000411
    parsed_dob = as.Date(parse_date_time(dob_orig, orders = c("mdy", "ymd", "yjmd"))), # date format column
    
    # 2. Fix the 2-digit year century roll-over issue
    fixed_dob = if_else(parsed_dob > Sys.Date(), parsed_dob - years(100), parsed_dob), # date format column
    
    # 3. Format the final output as a character string (YYYY-MM-DD)
    clean_dob = format(fixed_dob, "%Y-%m-%d") # character format column
  )



###########################################################
## combining full IAV and IBV into 1 file


#un_full_a <- full_a %>% distinct()
#un_full_b <- full_b %>% distinct()


all_ab_people <- rbind(full_a_demo_final, full_b) # as of 05/08/2026 13158 as of 07/2026 18234

all_geo_ab <- rbind(full_a_geo, full_b_geo) # as of 05/08/2026 13192 as of 07/2026 18454

all_ab_names_dob_id <- rbind(clean_dates_all_a_names_dob_id, clean_dates_all_b_names_dob_id)
### writing out the working geo and demographic data for ab combined

write.csv(all_ab_people, paste0(start_loc, "jules_flu_compilation/final_data/working_files/compiled_demographic_ab_20260818.csv"), row.names = FALSE, na = "")
write.csv(all_geo_ab, paste0(start_loc, "jules_flu_compilation/final_data/working_files/compiled_geolocation_ab_20260818.csv"), row.names = FALSE, na = "")
write.csv(all_ab_names_dob_id, paste0(start_loc, "jules_flu_compilation/final_data/working_files/compiled_names_dob_id_ab_20260818.csv"), row.names = FALSE, na = "")
#write.csv(all_geo_ab, paste0(start_loc, "jules_flu_compilation/geo_data/working_files/full_compiled_address_flu_07312026.csv"), row.names = FALSE, na = "")

### writing out the deidentifed demographic file
di_people_flu <- all_ab_people %>% select(-PATIENT_LNAME, -PATIENT_FNAME, -PATIENT_DOB_format)

write.csv(di_people_flu, paste0(start_loc, "jules_flu_compilation/final_deidentified_files/full_compiled_flu_deidentified_demo_data_20260818.csv"), row.names = FALSE, na = "")























