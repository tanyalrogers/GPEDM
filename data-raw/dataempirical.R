## Empirical data

library(dplyr)
library(tidyr)

#brown shrimp from gulf of mexico
sh=read.csv("E:/Documents/NOAA/Shrimp/BrownShrimpGulfMeanCPUEbyYear_cht.csv")
#there are 9 zones
shrimp=gather(sh,zone,cpue,zone11:zone21)

usethis::use_data(shrimp, overwrite = TRUE)

#bottom trawl data
d0=read.csv("E:/Documents/POSTDOC/NRC/matlab/bethany_bottomtrawl_fall.csv",stringsAsFactors = F)
colnames(d0)

trawl=d0 %>% group_by(YEAR,name) %>% 
  summarise(shortfin_squid=mean(Illex_illecebrosus_northern_shortfin_squid),
            longfin_squid=mean(Loligo_pealeii_longfin_squid),
            silver_hake=mean(Merluccius_bilinearis_silver_hake)) %>% 
  rename(REGION = name) %>% 
  mutate(REGION = factor(REGION, levels = c("GME","GEO","SNE","MAB"))) %>% 
  arrange(REGION) %>% 
  ungroup() %>% 
  as.data.frame()

usethis::use_data(trawl, overwrite = TRUE)