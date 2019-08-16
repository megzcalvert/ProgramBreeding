rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(broom)

set.seed(1996)

pheno17<- fread("./PhenoDatabase/PhenoLong_vi17.txt")
pheno18<- fread("./PhenoDatabase/PhenoVI_18.txt")

nested17<- pheno17 %>% 
  tidylog::select(-entity_id, -Variety, -range, 
                  -column, -Year, -Treated,-OverallTrial) %>% 
  group_by(Location, Trial, trait_id) %>% 
  nest() %>% 
  mutate(correlation = map(data, ~ lm(GRYLD ~ phenotypic_value, data = .x)),
         tidyCor = map(correlation,glance)) %>% 
  unnest(tidyCor) %>% 
  tidylog::select(-data,-correlation) %>% 
  separate(trait_id,c("trait_id","Date"), sep = "_") %>% 
  mutate(Date = as.Date(Date)) %>% 
  glimpse()

nested18<- pheno18 %>% 
  filter(trait_id != "PTHT") %>% 
  unite("trait_id",trait_id,phenotype_date, sep = "_") %>% 
  spread(key = trait_id, value = phenotype_value) %>% 
  gather(key = "trait_id",value = "phenotype_value", 
         `CanopyHeight_2017-12-01`:`GNDVI_2018-06-11`,
         `NDRE_2017-12-01`:`RE_2018-06-11`) %>% 
  drop_na(phenotype_value) %>% 
  group_by(Location, Trial, trait_id) %>% 
  nest() %>% 
  mutate(correlation = map(data, ~ lm(`GRYLD_2018-06-13` ~ phenotype_value,
                                      data = .x)),
         tidyCor = map(correlation,glance)) %>% 
  unnest(tidyCor) %>% 
  tidylog::select(-data,-correlation) %>% 
  separate(trait_id,c("trait_id","Date"), sep = "_") %>% 
  mutate(Date = as.Date(Date)) %>% 
  glimpse()

nested17 %>% 
  ggplot(aes(x = Date,
             y = adj.r.squared, 
             colour = Trial)) +
  geom_point(alpha = 0.5) +
  facet_wrap(trait_id ~ Location, ncol = 4) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  labs(title = "Linear regression for GRYLD 2016-2017",
       y = "adjusted r.squared")

nested18 %>% 
  ggplot(aes(x = Date,
             y = adj.r.squared, 
             colour = Trial)) +
  geom_point(alpha = 0.5) +
  facet_wrap(trait_id ~ Location, ncol = 4) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  labs(title = "Linear regression for GRYLD 2017-2018",
       y = "adjusted r.squared")
