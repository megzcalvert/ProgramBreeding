rm(list = objects())
ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(broom)
library(reshape2)
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
      size = rel(1.5),
      hjust = 0
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
# options(scipen = 999)
# useful infos **reproducible research**
sessionInfo()

##############################################################################
#### Read in data ####
# PYN

fileNames <- list.files(
  path = "./BeocatScripts/Results",
  full.names = T,
  pattern = "accuracy_80_100_viYld_pyn"
)
traitNames <- basename(fileNames) %>%
  str_remove_all(c("accuracy_80_100_viYld_|.txt"))

load.file <- function(filename) {
  d <- fread(file = filename, header = TRUE, check.names = F, data.table = F)
  d
}

pheno_pyn <- lapply(fileNames, load.file)

names(pheno_pyn) <- traitNames
pheno_pyn <- plyr::ldply(pheno_pyn, data.frame, .id = "location")
print(colnames(pheno_pyn))

pheno_pyn <- pheno_pyn %>%
  separate(location, c("trial", "location"), sep = "_") %>%
  mutate(year = str_remove(trial, "pyn")) %>%
  pivot_longer(
    cols = GRYLD_20170620_MP:NDVI_20190626_RN,
    names_to = "trait_id",
    values_to = "accuracy"
  ) %>%
  separate(trait_id,
           c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y%m%d")) %>% 
  drop_na(accuracy)

## AYN

fileNames <- list.files(
  path = "./BeocatScripts/Results",
  full.names = T,
  pattern = "accuracy_80_100_viYld_ayn"
)
traitNames <- basename(fileNames) %>%
  str_remove_all(c("accuracy_80_100_viYld_|.txt"))

pheno_ayn <- lapply(fileNames, load.file)

names(pheno_ayn) <- traitNames
pheno_ayn <- plyr::ldply(pheno_ayn, data.frame, .id = "location")
print(colnames(pheno_ayn))

pheno_ayn <- pheno_ayn %>%
  separate(location, c("trial", "location"), sep = "_") %>%
  mutate(year = str_remove(trial, "ayn")) %>%
  pivot_longer(
    cols = GRYLD_20170620_MP:NDVI_20190626_RN,
    names_to = "trait_id",
    values_to = "accuracy"
  ) %>%
  separate(trait_id,
           c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y%m%d")) %>% 
  drop_na(accuracy)

###############################################################################
###Figures ####
ayn17_summary<- pheno_ayn %>% 
  filter(year == "17",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2017-02-01")
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(sd = sd(accuracy),
            accuracy = mean(accuracy))

pheno_ayn %>% 
  filter(year == "17",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2017-02-01")
  ) %>% 
  ggplot(aes(x = phenotype_date,
             y = accuracy,
             colour = location,
             group = phenotype_date)) +
  geom_jitter(alpha = 0.15) +
facet_wrap(~trait_id, ncol = 4, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(date_labels = "%m/%d") +
  scale_colour_manual(values = c('#1b9e77','#d95f02','#7570b3','#e7298a')) +
  #theme(legend.position = "none") +
  labs(title = "Using VI as a co-factor to predict GRYLD",
       subtitle = "AYN 2016-2017")

ggsave("./Figures/GS_cofactor_17ayn.png",
       height = 30,
       width = 55, units = "cm", dpi = 350
)

pyn17_summary<- pheno_pyn %>% 
  filter(year == "17",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2017-02-01")
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(sd = sd(accuracy),
            accuracy = mean(accuracy)) %>% 
  mutate(year = "17")

pheno_pyn %>% 
  filter(year == "17",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2017-02-01")
  ) %>% 
  ggplot(aes(x = phenotype_date,
             y = accuracy,
             colour = location,
             group = phenotype_date)) +
  geom_jitter(alpha = 0.15) +
  geom_point(data = pyn17_summary,
             aes(x = phenotype_date, 
                 y = accuracy,
                 colour = location),
             size = 3) +
  geom_errorbar(data = pyn17_summary,
                aes(x = phenotype_date,
                    ymin = accuracy - sd,
                    ymax = accuracy + sd,
                    colour = location),
                size = 1) +
  facet_wrap(~trait_id, ncol = 4, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(date_labels = "%m/%d") +
  scale_colour_manual(values =  c('#1b9e77','#d95f02','#7570b3','#e7298a')) +
  #theme(legend.position = "none") +
  labs(title = "Using VI as a co-factor to predict GRYLD",
       subtitle = "PYN 2016-2017")

ggsave("./Figures/GS_cofactor_17pyn.png",
       height = 30,
       width = 55, units = "cm", dpi = 350
)

ayn18_summary<- pheno_ayn %>% 
  filter(year == "18",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2018-02-01")
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(sd = sd(accuracy),
            accuracy = mean(accuracy))

pheno_ayn %>% 
  filter(year == "18",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2018-02-01")
  ) %>% 
  ggplot(aes(x = phenotype_date,
             y = accuracy,
             colour = location,
             group = phenotype_date)) +
  geom_jitter(alpha = 0.15) +
  geom_point(data = ayn18_summary,
             aes(x = phenotype_date, 
                 y = accuracy,
                 colour = location),
             size = 3) +
  geom_errorbar(data = ayn18_summary,
                aes(x = phenotype_date,
                    ymin = accuracy - sd,
                    ymax = accuracy + sd,
                    colour = location),
                size = 1) +
  facet_wrap(~trait_id, ncol = 4, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(date_labels = "%m/%d") +
  scale_colour_manual(values = c('#1b9e77','#7570b3','#e7298a','#66a61e')) +
  #theme(legend.position = "none") +
  labs(title = "Using VI as a co-factor to predict GRYLD",
       subtitle = "AYN 2017-2018")

ggsave("./Figures/GS_cofactor_18ayn.png",
       height = 30,
       width = 55, units = "cm", dpi = 350
)

pyn18_summary<- pheno_pyn %>% 
  filter(year == "18",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2018-02-01")
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(sd = sd(accuracy),
            accuracy = mean(accuracy))  %>% 
  mutate(year = "18")

pheno_pyn %>% 
  filter(year == "18",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2018-02-01")
  ) %>% 
  ggplot(aes(x = phenotype_date,
             y = accuracy,
             colour = location,
             group = phenotype_date)) +
  geom_jitter(alpha = 0.15) +
  geom_point(data = pyn18_summary,
             aes(x = phenotype_date, 
                 y = accuracy,
                 colour = location),
             size = 3) +
  geom_errorbar(data = pyn18_summary,
                aes(x = phenotype_date,
                    ymin = accuracy - sd,
                    ymax = accuracy + sd,
                    colour = location),
                size = 1) +
  facet_wrap(~trait_id, ncol = 4, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(date_labels = "%m/%d") +
  scale_colour_manual(values = c('#1b9e77','#7570b3','#e7298a','#66a61e')) +
  #theme(legend.position = "none") +
  labs(title = "Using VI as a co-factor to predict GRYLD",
       subtitle = "PYN 2017-2018")

ggsave("./Figures/GS_cofactor_18pyn.png",
       height = 30,
       width = 55, units = "cm", dpi = 350
)

ayn19_summary<- pheno_ayn %>% 
  filter(year == "19",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2019-02-01")
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(sd = sd(accuracy),
            accuracy = mean(accuracy))

pheno_ayn %>% 
  filter(year == "19",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2019-02-01")
  ) %>% 
  ggplot(aes(x = phenotype_date,
             y = accuracy,
             colour = location,
             group = phenotype_date)) +
  geom_jitter(alpha = 0.15) +
  geom_point(data = ayn19_summary,
             aes(x = phenotype_date, 
                 y = accuracy,
                 colour = location),
             size = 3) +
  geom_errorbar(data = ayn19_summary,
                aes(x = phenotype_date,
                    ymin = accuracy - sd,
                    ymax = accuracy + sd,
                    colour = location),
                size = 1) +
  facet_wrap(~trait_id, ncol = 4, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(date_labels = "%m/%d") +
  scale_colour_manual(values =  c('#d95f02','#7570b3')) +
  #theme(legend.position = "none") +
  labs(title = "Using VI as a co-factor to predict GRYLD",
       subtitle = "AYN 2018-2019")

ggsave("./Figures/GS_cofactor_19ayn.png",
       height = 30,
       width = 55, units = "cm", dpi = 350
)

pyn19_summary<- pheno_pyn %>% 
  filter(year == "19",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2019-02-01")
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(sd = sd(accuracy),
            accuracy = mean(accuracy)) %>% 
  mutate(year = "19")

pheno_pyn %>% 
  filter(year == "19",
         trait_id != "RedEdge",
         trait_id != "Nir",
         trait_id != "RE",
         phenotype_date > as.Date("2019-02-01")
  ) %>% 
  ggplot(aes(x = phenotype_date,
             y = accuracy,
             colour = location,
             group = phenotype_date)) +
  geom_jitter(alpha = 0.15) +
  geom_point(data = pyn19_summary,
             aes(x = phenotype_date, 
                 y = accuracy,
                 colour = location),
             size = 3) +
  geom_errorbar(data = pyn19_summary,
                aes(x = phenotype_date,
                    ymin = accuracy - sd,
                    ymax = accuracy + sd,
                    colour = location),
                size = 1) +
  facet_wrap(~trait_id, ncol = 4, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(date_labels = "%m/%d") +
  scale_colour_manual(values =  c('#7570b3')) +
  #theme(legend.position = "none") +
  labs(title = "Using VI as a co-factor to predict GRYLD",
       subtitle = "PYN 2018-2019")

ggsave("./Figures/GS_cofactor_19pyn.png",
       height = 30,
       width = 55, units = "cm", dpi = 350
)

test_ayn18<- fread("./Results/GS_accuracy_ayn18_test.txt")

test_ayn18 %>% 
  #select(-GRYLD_NA_NA) %>% 
  pivot_longer(cols = everything(),
               names_to = "GSPrediction",
               values_to = "accuracy") %>% 
  mutate(GSPrediction = str_remove(GSPrediction, "_SA")) %>% 
  ggplot(aes(x = GSPrediction,
             y = accuracy)) +
  geom_boxplot() +
  geom_jitter() +
  labs(title = "Genomic prediction accuracies for GRYLD",
       subtitle = "AYN 2017-2018 Salina",
       x = "Factors")

pyn_summary<- bind_rows(pyn17_summary,pyn18_summary,pyn19_summary)

pheno_pyn %>% 
  filter(phenotype_date != "2017-12-01",
         phenotype_date != "2017-12-07",
         phenotype_date != "2017-12-08",
         phenotype_date != "2017-12-19") %>% 
  ggplot(aes(x = phenotype_date,
             y = accuracy,
             colour = factor(location))) +
  geom_jitter(alpha = 0.05,
              size = 1) +
geom_point(data = pyn_summary, 
           aes(x = phenotype_date,
               y = accuracy),
           size = 2) +
  geom_errorbar(data = pyn_summary,
                aes(ymin = accuracy - sd,
                    ymax = accuracy + sd)) +
  scale_colour_manual(values = c("#009E73","#490092","#006DDB",
                                 "#920000","#DB6D00"),
                      name = "Location") +
  scale_x_date(date_labels = "%m/%d",
               date_breaks = "3 weeks") +
  coord_cartesian(ylim = c(-1,1)) +
  facet_wrap(year~trait_id, scales = "free") +
  labs(x = "Date",
       y = "Accuracy")

ggsave(filename = "~/OneDrive - Kansas State University/Dissertation_Calvert/BreedingProgram/Figures/Figure6.png",
       width = 40,
       height = 20,
       units = "cm",
       dpi = 320)

