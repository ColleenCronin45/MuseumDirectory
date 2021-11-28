library(foreign)

source(here("GBIFcode.R"))
source(here("ArctosCode.R"))

BatRodentData = bind_rows(data2, dataset2)

write.dbf(as.data.frame(BatRodentData), here("BatRodentData.dbf"))

