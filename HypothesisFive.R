rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(broom)

set.seed(1996)

pheno17<- fread("./PhenoDatabase/PhenoLong_vi17.txt")
glimpse(pheno17)
pheno18<- fread("./PhenoDatabase/PhenoVI_18.txt")
glimpse(pheno18)

nested17<- pheno17 %>% 
  filter(Location != "RP") %>% 
  tidylog::select(-Variety, -range, 
                  -column, -Year, -Treated,-OverallTrial) %>% 
  spread(trait_id, phenotypic_value) %>% 
  gather(key = "trait_id", value = "phenotypic_value", 
         GRYLD:`RedEdge_2017-06-09`) %>% 
  group_by(Location, Trial, trait_id) %>% 
  summarise(mean = mean(phenotypic_value),
            variance = var(phenotypic_value),
            stdeviation = sqrt(variance)) %>% 
  drop_na(mean) %>% 
  filter(Location == "MP") %>% 
  separate(trait_id,c("trait_id","phenotype_date"), sep = "_") %>% 
  mutate(phenotype_date = as.Date(phenotype_date),
         phenotype_date = replace_na(phenotype_date,"2017-06-13"))

nested18<- pheno18 %>% 
  filter(trait_id != "PTHT") %>% 
  unite("trait_id",trait_id,phenotype_date, sep = "_") %>% 
  spread(key = trait_id, value = phenotype_value) %>% 
  gather(key = "trait_id",value = "phenotype_value", 
         `CanopyHeight_2017-12-01`:`RE_2018-06-11`) %>% 
  drop_na(phenotype_value) %>% 
  group_by(Location, Trial, trait_id) %>% 
  summarise(mean = mean(phenotype_value),
            variance = var(phenotype_value),
            stdeviation = sqrt(variance)) %>% 
  drop_na(mean) %>% 
  filter(Location == "MP" |
           Location == "RN" |
           Location == "RP") %>% 
  separate(trait_id,c("trait_id","phenotype_date"), sep = "_") %>% 
  mutate(phenotype_date = as.Date(phenotype_date))

nested17 %>% 
  ggplot(aes(x = phenotype_date, y = mean, colour = Trial)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - stdeviation, ymax = mean + stdeviation)) +
  facet_wrap(~trait_id, scales = "free") +
  theme_bw()

