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
geno <- fread("./GenoDatabase/SelectedGeno_Numeric.txt")
geno[1:10, 1:10]
hapgeno <- geno[, 4:ncol(geno)]
hapgeno[hapgeno == 0] <- -1
hapgeno[hapgeno == 0.5] <- 0

geno <- bind_cols(geno[, 1:3], hapgeno)
geno[1:10, 1:10]
rm(hapgeno)

genoVarieties <- as.data.frame(colnames(geno[, 4:ncol(geno)]))
names(genoVarieties) <- "Variety"

blups17_ayn <- fread("./Results/Blups_ayn_17.txt")
blups17_pyn <- fread("./Results/Blups_pyn_17.txt")

blups18_ayn <- fread("./Results/Blups_ayn_18.txt")
blups18_pyn <- fread("./Results/Blups_pyn_18.txt")

blups17_ayn <- blups17_ayn %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value) %>%
  semi_join(genoVarieties, by = "Variety")

blups17_pyn <- blups17_pyn %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value) %>%
  semi_join(genoVarieties, by = "Variety")

blups18_ayn <- blups18_ayn %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value) %>%
  semi_join(genoVarieties, by = "Variety")

blups18_pyn <- blups18_pyn %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value) %>%
  semi_join(genoVarieties, by = "Variety")

summaryBlup17_ayn <- blups17_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

summaryBlup17_pyn <- blups17_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

summaryBlup17_ayn %>% 
  ggplot(aes(x = phenotype_date, y = average, colour = location)) +
  geom_point() +
  geom_point(aes(x = phenotype_date, y = minimum, colour = location)) +
  geom_point(aes(x = phenotype_date, y = maximum, colour = location)) +
  geom_errorbar(aes(ymin = average - sd, ymax = average + sd)) +
  facet_wrap(~trait_id, scales = "free") +
  labs(title = "BLUPs of VI and Reflectance Values 2017",
       subtitle = "Advanced Yield Trials")

summaryBlup17_pyn %>% 
  ggplot(aes(x = phenotype_date, y = average, colour = location)) +
  geom_point() +
  geom_point(aes(x = phenotype_date, y = minimum, colour = location)) +
  geom_point(aes(x = phenotype_date, y = maximum, colour = location)) +
  geom_errorbar(aes(ymin = average - sd, ymax = average + sd)) +
  facet_wrap(~trait_id, scales = "free") +
  labs(title = "BLUPs of VI and Reflectance Values 2017",
       subtitle = "Preliminary Yield Trials")

summaryBlup18_ayn <- blups18_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

summaryBlup18_pyn <- blups18_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

summaryBlup18_ayn %>% 
  ggplot(aes(x = phenotype_date, y = average, colour = location)) +
  geom_point() +
  geom_point(aes(x = phenotype_date, y = minimum, colour = location)) +
  geom_point(aes(x = phenotype_date, y = maximum, colour = location)) +
  geom_errorbar(aes(ymin = average - sd, ymax = average + sd)) +
  facet_wrap(~trait_id, scales = "free") +
  labs(title = "BLUPs of VI and Reflectance Values 2018",
       subtitle = "Advanced Yield Trials")

summaryBlup18_pyn %>% 
  ggplot(aes(x = phenotype_date, y = average, colour = location)) +
  geom_point() +
  geom_point(aes(x = phenotype_date, y = minimum, colour = location)) +
  geom_point(aes(x = phenotype_date, y = maximum, colour = location)) +
  geom_errorbar(aes(ymin = average - sd, ymax = average + sd)) +
  facet_wrap(~trait_id, scales = "free") +
  labs(title = "BLUPs of VI and Reflectance Values 2018",
       subtitle = "Preliminary Yield Trials")

### Kinship matrix
snpMatrix <- t(geno[, 4:ncol(geno)])

relMat <- A.mat(snpMatrix)
