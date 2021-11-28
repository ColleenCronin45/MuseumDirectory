library(foreign)

#source will run the 2 code files
source(here("GBIFcode.R"))
source(here("ArctosCode.R"))

#combining GBIF and Arctos datasets
BatRodentData = bind_rows(data2, dataset2)

#exporting data as a .dbf
write.dbf(as.data.frame(BatRodentData), here("BatRodentData.dbf"))

