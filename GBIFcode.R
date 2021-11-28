

library(readxl)
library(tidyverse)
library(here)
####Saving GBIF as one####

#   
# 
# Sheet1 = read_excel("C:/Users/Lump/OneDrive/Desktop/GBIF_TissPres_Yes_V1.xlsx", sheet = "Sheet1")
# 
# warnings()
# 
# Sheet2 = read_excel("C:/Users/Lump/OneDrive/Desktop/GBIF_TissPres_Yes_V2.xlsx", sheet = "Sheet2")
# dx
# 
# Sheet2.1 = Sheet2 %>%
#   mutate(eventDate = as.character(eventDate)) %>%
#   mutate(verbatimEventDate = as.character(verbatimEventDate))
# 
# Sheet1.1 = Sheet1 %>%
#   mutate(level3Name = as.character(level3Name))
# 
# dataset = bind_rows(Sheet1.1,Sheet2.1)
# 
# 

#loading in the combined dataset
load(here("GBIFdataset.rdata"))

#view a random sample
#sample_n(dataset, 1000) %>%
#  view()

#create frozen and tissue variables
dataset = dataset %>%
  mutate(frozen = case_when(
    str_detect(tolower(preparations), "froz") ~ 1,
    str_detect(tolower(preparations), "congelad") ~ 1,
    str_detect(tolower(preparations), "freez") ~ 1,
    str_detect(tolower(dynamicProperties), "froz") ~ 1,
    str_detect(tolower(dynamicProperties), "congelad") ~ 1,
    str_detect(tolower(dynamicProperties), "freez") ~ 1,
    str_detect(tolower(occurrenceRemarks), "froz") ~ 1,
    str_detect(tolower(occurrenceRemarks), "congelad") ~ 1,
    str_detect(tolower(occurrenceRemarks), "freez") ~ 1,
    TRUE ~ 0
  )) %>%
    mutate(tissue= case_when(
      str_detect(tolower(preparations), "tiss") ~ 1,
      str_detect(tolower(preparations), "tejido") ~ 1,
      str_detect(tolower(dynamicProperties), "tiss") ~ 1,
      str_detect(tolower(dynamicProperties), "tejido") ~ 1,
      str_detect(tolower(occurrenceRemarks), "tiss") ~ 1,
      str_detect(tolower(occurrenceRemarks), "tejido") ~ 1,
      TRUE ~ 0
    ))

#view a breakdown of data in frozen and tissue by order
# dataset %>% 
#   filter(order != "order") %>%
#   group_by(order) %>%
#   summarize(frozen_m = mean(frozen, na.rm = T),
#             frozen_c = sum(frozen, na.rm = T),
#             tissue_m = mean(tissue, na.rm = T),
#             tissue_c = sum(tissue, na.rm = T),
#             n = n())

#keep only the variables of interest
dataset2 = dataset %>%
  select(institutionCode,
         year,
         verbatimEventDate,
         countryCode,
         stateProvince,
         decimalLatitude,
         decimalLongitude,
         order,
         family,
         genus,
         speciesKey,
         frozen,
         tissue)

#replacing year with year from verbatimEventDate when year is NA
dataset2 = dataset2 %>%
  mutate(year = case_when(
    is.na(year) ~ as.numeric(str_extract(verbatimEventDate, "^\\d{4}(?=-)")),
    TRUE ~ year
  ))

#removing years listed as 1800 and replaced with NA
dataset2$year[dataset2$year <1800] <-NA

#making a list of unique species keys
dataset3 = dataset2 %>%
  filter(!duplicated(speciesKey)) %>%
  filter(!is.na(speciesKey)) %>%
  select(speciesKey)

#querying GBIF database to match species names with species key
test = taxize::classification(dataset3$speciesKey, db = "gbif")

#intializing empty dataset
dataset_match = data.frame()

#looping through API response to obtain species names from 
for(n in 1:length(test)){
  temp1 = names(test[n])
  if(length(test[[n]]$name[test[[n]]$rank=="species"])>0){
    temp2 = test[[n]]$name[test[[n]]$rank=="species"]
  } else {
    temp2 = NA
  }
  temp3 = data.frame(speciesKey = temp1,
                     species = temp2)
  dataset_match = bind_rows(dataset_match, temp3)
}

dataset_match = dataset_match %>%
  mutate(speciesKey = as.numeric(speciesKey))

dataset2 = dataset2 %>%
  left_join(dataset_match, by = c("speciesKey"= "speciesKey"))

#checking to see if we missed muscle instead of tissue - only ~1000
#dataset %>%
 # filter(str_detect(tolower(preparations), "muscle") | str_detect(tolower(dynamicProperties), "muscle") | str_detect(tolower(occurrenceRemarks), "muscle"))

dataset2 = dataset2 %>%
  mutate(dataset = "GBIF")
