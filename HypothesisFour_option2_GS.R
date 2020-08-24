rm(list = objects())
ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(broom)
library(rrBLUP)

#### Theme set ####
custom_theme <- theme_minimal() %+replace%
  theme(
    axis.title = element_text(
      colour = "black",
      size = rel(2)
    ),
    axis.title.x = element_text(
      vjust = 0,
      margin = margin(
        t = 0, r = 0.25,
        b = 0, l = 0,
        unit = "cm"
      )
    ),
    axis.title.y = element_text(
      vjust = 1,
      angle = 90,
      margin = margin(
        t = 0, r = 0.25,
        b = 0, l = 0.1,
        unit = "cm"
      )
    ),
    axis.text = element_text(
      colour = "black",
      size = rel(1.5)
    ),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(3, "pt"),
    axis.line = element_line(
      color = "black",
      size = 0.5
    ),
    legend.key.size = unit(2, "lines"),
    # legend.background = element_rect(fill = NULL, colour = NULL),
    # legend.box = NULL,
    legend.margin = margin(
      t = 0, r = 0.75,
      b = 0, l = 0.75,
      unit = "cm"
    ),
    legend.text = element_text(size = rel(2)),
    legend.title = element_text(size = rel(1.5)),
    panel.grid.major = element_line(
      colour = "#969696",
      linetype = 3
    ),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(
      size = rel(2),
      margin = margin(
        t = 0.1, r = 0.1,
        b = 0.1, l = 0.1,
        unit = "cm"
      )
    ),
    plot.margin = margin(
      t = 0.5, r = 0.5,
      b = 0.5, l = 0,
      unit = "cm"
    ),
    plot.title = element_text(
      colour = "black",
      size = rel(3),
      vjust = 0,
      hjust = 0,
      margin = margin(
        t = 0.25, r = 0.25,
        b = 0.5, l = 0.25,
        unit = "cm"
      )
    ),
    plot.subtitle = element_text(
      colour = "black",
      size = rel(1.5)
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      size = 1
    ),
    strip.text = element_text(
      colour = "black",
      size = rel(1.5)
    ),
    complete = F
  )

theme_set(custom_theme)

getwd()
set.seed(1964)

# useful infos **reproducible research**
sessionInfo()
##############################################################################
#### Read in data ####

accuracy17_ayn<- 
  fread("./BeocatScripts/Results/GenomicSelection_80_100_accuracy17_ayn.txt")
group17_ayn_key<- fread("./BeocatScripts/Results/group17_ayn_key.txt")

accuracy17_pyn<- 
  fread("./BeocatScripts/Results/GenomicSelection_80_100_accuracy17_pyn.txt")
group17_pyn_key<- fread("./BeocatScripts/Results/group17_pyn_key.txt")

accuracy18_ayn<- 
  fread("./BeocatScripts/Results/GenomicSelection_80_100_accuracy18_ayn.txt")
group18_ayn_key<- fread("./BeocatScripts/Results/group18_ayn_key.txt")

accuracy18_pyn<- 
  fread("./BeocatScripts/Results/GenomicSelection_80_100_accuracy18_pyn.txt")
group18_pyn_key<- fread("./BeocatScripts/Results/group18_pyn_key.txt")

accuracy19_ayn<- 
  fread("./BeocatScripts/Results/GenomicSelection_80_100_accuracy19_ayn.txt")
group19_ayn_key<- fread("./BeocatScripts/Results/group19_ayn_key.txt")

accuracy19_pyn<- 
  fread("./BeocatScripts/Results/GenomicSelection_80_100_accuracy19_pyn.txt")
group19_pyn_key<- fread("./BeocatScripts/Results/group19_pyn_key.txt")

## 17 
accuracy17_ayn<- accuracy17_ayn %>% 
  mutate(ID = row_number()) %>% 
  pivot_longer(cols = -ID,
    names_to = "group",
               values_to = "Accuracies") %>% 
  left_join(group17_ayn_key) %>% 
  select(-group) %>% 
  separate(trait_id,c("trait_id","phenotype_date","location"),
           sep = "_") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"))

ayn17<- accuracy17_ayn %>% 
  filter(trait_id != "GRYLD") %>% 
  group_by(phenotype_date, location, trait_id) %>% 
  mutate(group = group_indices()) %>% 
  ggplot(aes(x = phenotype_date, y = Accuracies, colour = location)) +
  geom_boxplot(aes(group = group), fill = NA) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  labs(title = "Genomic Prediction accuracy for Vegetation Indices",
       subtitle = "AYN 2016-2017 Cross-Validation")
ayn17
ggpubr::ggexport(ayn17, filename = "./Figures/GP_vi_ayn17.png",
                 width = 1500,
                 height = 750)

accuracy17_pyn<- accuracy17_pyn %>% 
  mutate(ID = row_number()) %>% 
  pivot_longer(cols = -ID,
               names_to = "group",
               values_to = "Accuracies") %>% 
  left_join(group17_pyn_key) %>% 
  select(-group) %>% 
  separate(trait_id,c("trait_id","phenotype_date","location"),
           sep = "_") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"))

pyn17<- accuracy17_pyn %>% 
  filter(trait_id != "GRYLD") %>% 
  group_by(phenotype_date, location, trait_id) %>% 
  mutate(group = group_indices()) %>% 
  ggplot(aes(x = phenotype_date, y = Accuracies, colour = location)) +
  geom_boxplot(aes(group = group), fill = NA) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  labs(title = "Genomic Prediction accuracy for Vegetation Indices",
       subtitle = "PYN 2016-2017 Cross-Validation")
pyn17
ggpubr::ggexport(pyn17, filename = "./Figures/GP_vi_pyn17.png",
                 width = 1500,
                 height = 750)

## 18 
accuracy18_ayn<- accuracy18_ayn %>% 
  mutate(ID = row_number()) %>% 
  pivot_longer(cols = -ID,
               names_to = "group",
               values_to = "Accuracies") %>% 
  left_join(group18_ayn_key) %>% 
  select(-group) %>% 
  separate(trait_id,c("trait_id","phenotype_date","location"),
           sep = "_") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"))

ayn18<- accuracy18_ayn %>% 
  filter(trait_id != "GRYLD") %>% 
  group_by(phenotype_date, location, trait_id) %>% 
  mutate(group = group_indices()) %>% 
  ggplot(aes(x = phenotype_date, y = Accuracies, colour = location)) +
  geom_boxplot(aes(group = group), fill = NA) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  labs(title = "Genomic Prediction accuracy for Vegetation Indices",
       subtitle = "AYN 2017-2018 Cross-Validation")
ayn18
ggpubr::ggexport(ayn18, filename = "./Figures/GP_vi_ayn18.png",
                 width = 1500,
                 height = 750)

accuracy18_pyn<- accuracy18_pyn %>% 
  mutate(ID = row_number()) %>% 
  pivot_longer(cols = -ID,
               names_to = "group",
               values_to = "Accuracies") %>% 
  left_join(group18_pyn_key) %>% 
  select(-group) %>% 
  separate(trait_id,c("trait_id","phenotype_date","location"),
           sep = "_") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"))

pyn18<- accuracy18_pyn %>% 
  filter(trait_id != "GRYLD") %>% 
  group_by(phenotype_date, location, trait_id) %>% 
  mutate(group = group_indices()) %>% 
  ggplot(aes(x = phenotype_date, y = Accuracies, colour = location)) +
  geom_boxplot(aes(group = group), fill = NA) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  labs(title = "Genomic Prediction accuracy for Vegetation Indices",
       subtitle = "PYN 2017-2018 Cross-Validation")
pyn18
ggpubr::ggexport(pyn18, filename = "./Figures/GP_vi_pyn18.png",
                 width = 1500,
                 height = 750)

## 19 
accuracy19_ayn<- accuracy19_ayn %>% 
  mutate(ID = row_number()) %>% 
  pivot_longer(cols = -ID,
               names_to = "group",
               values_to = "Accuracies") %>% 
  left_join(group19_ayn_key) %>% 
  select(-group) %>% 
  separate(trait_id,c("trait_id","phenotype_date","location"),
           sep = "_") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"))

ayn19<- accuracy19_ayn %>% 
  filter(trait_id != "GRYLD") %>% 
  group_by(phenotype_date, location, trait_id) %>% 
  mutate(group = group_indices()) %>% 
  ggplot(aes(x = phenotype_date, y = Accuracies, colour = location)) +
  geom_boxplot(aes(group = group), fill = NA) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  labs(title = "Genomic Prediction accuracy for Vegetation Indices",
       subtitle = "AYN 2018-2019 Cross-Validation")
ayn19
ggpubr::ggexport(ayn19, filename = "./Figures/GP_vi_ayn19.png",
                 width = 1500,
                 height = 750)

accuracy19_pyn<- accuracy19_pyn %>% 
  mutate(ID = row_number()) %>% 
  pivot_longer(cols = -ID,
               names_to = "group",
               values_to = "Accuracies") %>% 
  left_join(group19_pyn_key) %>% 
  select(-group) %>% 
  separate(trait_id,c("trait_id","phenotype_date","location"),
           sep = "_") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"))

pyn19<- accuracy19_pyn %>% 
  filter(trait_id != "GRYLD") %>% 
  group_by(phenotype_date, location, trait_id) %>% 
  mutate(group = group_indices()) %>% 
  ggplot(aes(x = phenotype_date, y = Accuracies, colour = location)) +
  geom_boxplot(aes(group = group), fill = NA) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  labs(title = "Genomic Prediction accuracy for Vegetation Indices",
       subtitle = "PYN 2018-2019 Cross-Validation")
pyn19
ggpubr::ggexport(pyn19, filename = "./Figures/GP_vi_pyn19.png",
                 width = 1500,
                 height = 750)
