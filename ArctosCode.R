
## libraries
library(tidyverse)
library(here)

## load in Arctos data
data = read_csv(here("Arctos_ALL.csv")) ## takes 1-2 minutes to load

## subset to Chiroptera and Rodentia
data$order2 = plyr::revalue(data$order,
                            c("Rodentia; Rodentia"="Rodentia",
                              "Rodentia; Rodentia; Rodentia; Rodentia"="Rodentia",
                              "Chiroptera; Chiroptera"="Chiroptera"))

data = data[data$order2 %in% c("Rodentia","Chiroptera"),]

## fix order
data$order=data$order2
data$order2=NULL

## parse prep type with grepl
# data$frozen=ifelse(grepl("frozen|freeze",data$preparationProcess),1,0)

#create frozen and tissue variables
data = data %>%
  mutate(frozen = case_when(
    str_detect(tolower(preparationProcess), "froz") ~ 1,
    str_detect(tolower(preparationProcess), "congelad") ~ 1,
    str_detect(tolower(preparationProcess), "freez") ~ 1,
    TRUE ~ 0
  )) 

## fix year
data$year=ifelse(data$year==207,2017,data$year)

data2 = data %>%
  select(institutionCode,
         year,
         verbatimEventDate,
         countryCode = country,
         stateProvince,
         decimalLatitude,
         decimalLongitude,
         order,
         family,
         genus,
         species,
         frozen)

data2 = data2 %>%
  mutate(dataset = "Arctos")
