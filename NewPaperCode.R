require(tidyverse)
require(sqldf)
require(here)
require(readxl)
require(foreign)
# 
# Dataset = read_delim("C:/Users/Lump/OneDrive/Desktop/MCMSpring/CaseStudy/FINAL_PAPER_STUFF/0124064-210914110416597/occurrence.txt")
# 
# Dataset = Dataset %>% select("gbifID", 
#                                "accessRights", 
#                                "bibliographicCitation", 
#                                "identifier",
#                                "language",
#                                "license",
#                                "publisher",
#                                "references",
#                                "type",
#                                "institutionCode",
#                                "basisOfRecord",
#                                "dynamicProperties",
#                                "occurrenceID",
#                                "sex",
#                                "preparations",
#                                "eventDate",
#                                "year",
#                                "verbatimEventDate",
#                                "continent",
#                                "countryCode",
#                                "stateProvince",
#                                "decimalLatitude",
#                                "decimalLongitude",
#                                "order",
#                                "family",
#                                "genus",
#                                "genericName",
#                                "specificEpithet",
#                                "taxonRank",
#                                "datasetKey",
#                                "taxonKey",
#                                "species",
#                                "verbatimScientificName",
#                              "occurrenceRemarks"
# )

load(here("NewPaperEnviro_wRemarks.RData"))

Dataset %>% group_by(language) %>% summarize(n=n())

Dataset %>% group_by(order) %>% summarize(n=n())

Dataset = Dataset %>% filter(order == "Chiroptera" | order == "Rodentia")

Dataset %>% filter(is.na(year) & is.na(eventDate))

tissue2 = c("blood", "sangre", "tiss", "tejid", "serum", "suero", "plasma", "biopsy", "biopsia", "swab", "torunda") %>%
  paste(collapse = "|")

buffer2 = c("ethanol", "etanol", "EtOH", "VTM", "RNAlater","shield", "alcohol", "lysis", "DMSO","EDTA") %>%
  paste(collapse = "|")

frozen2 = c("froze", "freez", "congelad") %>%
  paste(collapse = "|")

#create frozen and tissue variables
Dataset = Dataset %>%
  mutate(frozen = case_when(
    str_detect(tolower(preparations), frozen2) ~ 1,
    str_detect(tolower(dynamicProperties), frozen2) ~ 1,
    str_detect(tolower(occurrenceRemarks), frozen2) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(tissue= case_when(
    str_detect(tolower(preparations), tissue2) ~ 1,
    str_detect(tolower(dynamicProperties), tissue2) ~ 1,
    str_detect(tolower(occurrenceRemarks), tissue2) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(buffer= case_when(
    str_detect(tolower(preparations), buffer2) ~ 1,
    str_detect(tolower(dynamicProperties), buffer2) ~ 1,
    str_detect(tolower(occurrenceRemarks), buffer2) ~ 1,
    TRUE ~ 0
  ))

Dataset %>%
  summarize(frozen = sum(as.numeric(frozen), na.rm = T), 
            tissue = sum(as.numeric(tissue), na.rm = T),  
            buffer = sum(as.numeric(buffer), na.rm = T))

Dataset %>%
  filter(str_detect(occurrenceRemarks, "congelad") | str_detect(dynamicProperties, "congelad") | str_detect(preparations, "congelad")) %>%
  view()

extra.countries = data.frame("Country Code" = c("SS", "GG"),
                                         "Country" = c("South Sudan", "Guernsey"),
                                         "Latitude" = c(6.877, 49.4482),
                                         "Longitude" = c(31.307, 2.5895)) %>%
  rename("Country Code" = Country.Code)

countries = read_csv(here("LatLong.csv")) %>%
  bind_rows(extra.countries)


Dataset$countryCode[Dataset$countryCode == "BQ"] <- "AN"

Dataset = Dataset %>%
  left_join(countries, by = c("countryCode" = "Country Code"))


Dataset %>%
  filter(is.na(decimalLatitude)) %>%
  group_by(countryCode) %>%
  summarize(n = n()) 



Dataset = Dataset %>%
  mutate(decimalLatitude = case_when(
    !is.na(decimalLatitude) ~ decimalLatitude,
    is.na(decimalLatitude) ~ Latitude
  ), decimalLongitude = case_when(
    !is.na(decimalLongitude) ~ decimalLongitude,
    is.na(decimalLongitude) ~ Longitude
  ))

Dataset = Dataset %>%
  mutate(tissUbuff = case_when(
    tissue == 1 ~ 1,
    buffer == 1 ~ 1,
    TRUE ~ 0
  ))

Dataset = Dataset %>%
  mutate(ethanolPrep = case_when(
    str_detect(tolower(preparations), "ethanol|etanol|EtOH|alcohol") ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(ethanolAll= case_when(
    str_detect(tolower(preparations), "ethanol|etanol|EtOH|alcohol") ~ 1,
    str_detect(tolower(dynamicProperties), "ethanol|etanol|EtOH|alcohol") ~ 1,
    str_detect(tolower(occurrenceRemarks), "ethanol|etanol|EtOH|alcohol") ~ 1,
    TRUE ~ 0
  ))
Dataset %>%
  summarize(ethanolPrep = sum(ethanolPrep, na.rm = TRUE),
            ethanolAll = sum(ethanolAll, na.rm = TRUE)) %>%
  mutate(proportion = ethanolPrep/ethanolAll)

Dataset %>%
  group_by(order, species)%>%
  summarize(n = n())%>%
  write_csv(here("SpeciesFrequency.csv"))

Dataset %>%
  group_by(order, species)%>%
  summarize(frozen = sum(frozen, na.rm = T),
            tissue = sum(tissue, na.rm = T),
            buffer = sum(buffer, na.rm = T),
            tissUbuff = sum(tissUbuff, na.rm = T))%>%
  write_csv(here("SpeciesFTBfrequency3.4.22.csv"))

Dataset2 = Dataset %>%
  select(-eventDate)

#exporting data as .csv
Dataset %>%
  write_csv(here("Dataset_3.4.22.csv"))

#exporting data as a .dbf
write.dbf(as.data.frame(Dataset2), here("Dataset_2.21.dbf"))
