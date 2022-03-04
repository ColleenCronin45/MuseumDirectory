library(foreign)
library(readxl)
library(tidyverse)
library(here)
library(ggthemes)
library(scales)
library(rgbif)


#source will run the 2 code files
source(here("GBIFcode.R"))
source(here("ArctosCode.R"))

countries = read_excel(here("country_codes.xls"))

view(head(data2))

data2 = data2 %>%
  left_join(countries, by = c("countryCode" = "Definition"))

#combining GBIF and Arctos datasets
BatRodentData = bind_rows(data2, dataset2) %>%
  mutate(propercodes = case_when(
    dataset == "Arctos" ~ `Code Value`, 
    dataset == "GBIF" ~ countryCode
  ))

view(head(BatRodentData))

#exporting data as a .dbf
write.dbf(as.data.frame(BatRodentData), here("BatRodentData.dbf"))

graph1 = BatRodentData %>% 
  filter(order!="order" & year>=1850) %>%
  ggplot(aes(x=year, color=as.logical(frozen)))+
  geom_histogram(binwidth=5, fill="white")+
  scale_y_continuous(labels = comma)+
  facet_wrap(~order)+
  theme_fivethirtyeight()+
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"))+
  labs(color = "Frozen Samples")

graph2 = BatRodentData %>% 
  filter(order!="order" & year>=1850) %>%
  ggplot(aes(x=year, color=as.logical(frozen)))+
  geom_histogram(binwidth=5, fill="white")+
  scale_y_continuous(labels = comma)+
  #facet_wrap(~order)+
  theme_fivethirtyeight()+
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"))+
  labs(color = "Frozen Samples")

ggsave(plot=graph1, filename=here("order_histograms.png"), device="png", height=5, width=10)
ggsave(plot=graph2, filename=here("all_histograms.png"), device="png", height=5, width=10)


#GBIF citations

dataset2 = dataset %>%
  filter(datasetKey!="datasetKey")

keys = unique(dataset2$datasetKey)

save(keys, file = here("keys.rdata"))

source(here("GBIF_Citations.Rmd"))
  
#Arctos Citations

ArctosFiles = list.files(here("Arctos", "Arctos"))

ArctosGUIDs=data.frame()

for(n in ArctosFiles){
  temp=read_csv(here("Arctos", "Arctos", n), col_select = "GUID")
  ArctosGUIDs=bind_rows(ArctosGUIDs, temp)
}

ArctosGUIDs=ArctosGUIDs %>%
  filter(!duplicated(GUID))



