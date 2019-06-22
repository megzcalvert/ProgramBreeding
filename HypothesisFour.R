rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(asreml)
library(lme4)

##### Read in HTP data, Use separate files at first, unsure if the same
pheno17<- fread("./PhenoDatabase/PhenoLong_vi17.txt")
str(pheno17)

pheno18<- fread("./PhenoDatabase/PhenoVI_18.txt")
str(pheno18)

# Use only Varities with more than one rep per season

pheno17<- pheno17 %>% 
  spread(key = trait_id, value = phenotypic_value) %>% 
  distinct() %>% 
  glimpse() %>% 
  filter(!str_detect(Trial,"PYN")) %>%
  gather(key = trait_id, value = phenotype_value, 
         GRYLD:`RedEdge_2017-06-16`) %>% 
  drop_na(phenotype_value) %>% 
  group_by(Year, Location, Trial, trait_id, Variety) %>% 
  mutate(Count = n()) %>% 
  filter(Count > 1) %>% 
  ungroup() %>% 
  group_by(Year,Location,Trial,trait_id) %>% 
  arrange(Year,Location,Trial,Variety,trait_id) %>% 
  filter(!(Location == "RP")) 

phenoMaxRep18<- pheno18 %>% 
  group_by(Year,Location,Trial,Variety) %>% 
  summarise(max = max(rep)) 

pheno18<- pheno18 %>% 
  unite("trait_id", trait_id, phenotype_date, sep = "_") %>% 
  left_join(phenoMaxRep18, by = c("Year","Location","Trial","Variety")) %>% 
  glimpse() %>% 
  filter(max > 1) %>% 
  filter(Trial != "PYNA") %>%
  filter(Trial != "PYNB") %>%
  filter(Trial != "GSPYN") %>% 
  filter(Trial != "DHPYN") %>% 
  group_by(Year, Location, Trial, trait_id) 

rm(phenoMaxRep18)

###### Calculating heritability

pheno17_H2<- pheno17 %>% 
  nest() %>% 
  mutate(model = map(data,
                     ~ lmer(
                       phenotype_value ~ (1|Variety) + (1|rep) + (1|column) + 
                         (1|range),
                       data = .)),
         VarianceComp = map(model, as.data.frame(VarCorr))) %>% 
  unnest(VarianceComp) 

pheno17_H2<- pheno17_H2 %>% 
  group_by(Year,Location,Trial,trait_id) %>% 
  select(-value.var1,-value.var2,-value.sdcor) %>% 
  spread(key = value.grp, value = value.vcov) %>% 
  mutate(H2 = Variety / (Variety + Residual/2)) %>% 
  separate(trait_id,c("trait_id","Date"), sep = "_") %>% 
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>% 
  mutate(Date = replace_na(Date, "2017-06-13"))
 
  
pheno17_H2 %>% 
  ggplot(aes(x = Date, y = H2, colour = Trial)) +
  geom_point(size = 2) +
  facet_wrap(trait_id ~ Location, scales = "fixed", ncol = 3) +
  scale_x_date(breaks = "1 week",date_labels = "%m/%d") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 16),
        strip.text = element_text(size = 18),
        title = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  labs(title = "Broad-sense Heritability", 
       subtitle = "2016/2017 season")

pheno18_H2<- pheno18 %>% 
  nest() %>% 
  mutate(model = map(data, 
                     ~lmer(phenotype_value ~ 
                             (1|Variety) + (1|rep) + (1|column) + (1|range), 
                            data = .)),
         VarianceComp = map(model, as.data.frame(VarCorr))) %>% 
  unnest(VarianceComp) 

pheno18_H2<- pheno18_H2 %>% 
  group_by(Year,Location,Trial,trait_id) %>% 
  select(-value.var1,-value.var2,-value.sdcor) %>% 
  spread(key = value.grp, value = value.vcov) %>% 
  mutate(H2 = Variety / (Variety + Residual/2))

pheno18_H2 %>% 
  separate(trait_id,c("trait_id","Date"), sep = "_") %>% 
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>% 
  ggplot(aes(x = Date, y = H2, colour = Trial)) +
  geom_point(size = 2) +
  facet_wrap(trait_id ~ Location, scales = "fixed", ncol = 4) +
  scale_x_date(breaks = "3 weeks",date_labels = "%m/%d") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 16),
        strip.text = element_text(size = 18),
        title = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  labs(title = "Broad-sense Heritability", 
       subtitle = "2017/2018 season")
  


















