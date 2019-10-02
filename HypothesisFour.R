rm(list = objects())
ls()

library(tidyverse)
library(tidylog)
library(janitor)
library(lme4)

##### Set up work space ####
#### Theme set
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
    legend.key.size = unit(1.75, "lines"),
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
      size = rel(2)
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      size = 1
    ),
    strip.text = element_text(
      colour = "black",
      size = rel(1)
    ),
    complete = F
  )

theme_set(custom_theme)

set.seed(1935)

##### Read in HTP data, Use separate files at first, unsure if the same
pheno17 <- data.table::fread("./PhenoDatabase/PhenoLong_vi17.txt")
str(pheno17)

pheno18 <- data.table::fread("./PhenoDatabase/PhenoVI_18.txt")
str(pheno18)

# Use only Varities with more than one rep per season

pheno17 <- pheno17 %>%
  spread(key = trait_id, value = phenotypic_value) %>%
  distinct() %>%
  glimpse() %>%
  filter(!str_detect(Trial, "PYN")) %>%
  gather(
    key = trait_id, value = phenotype_value,
    GRYLD:`RedEdge_2017-06-16`
  ) %>%
  drop_na(phenotype_value) %>%
  group_by(Year, Location, Trial, trait_id, Variety) %>%
  mutate(Count = n()) %>%
  ungroup() %>%
  group_by(Year, Location, Trial, trait_id) %>%
  arrange(Year, Location, Trial, Variety, trait_id) %>%
  filter(!(Location == "RP"))

harmonic17<- 1/mean(1/pheno17$rep) 

phenoMaxRep18 <- pheno18 %>%
  group_by(Year, Location, Trial, Variety) %>%
  summarise(max = max(rep))

pheno18 <- pheno18 %>%
  unite("trait_id", trait_id, phenotype_date, sep = "_") %>%
  left_join(phenoMaxRep18, by = c("Year", "Location", "Trial", "Variety")) %>%
  glimpse() %>%
  filter(Trial != "PYNA") %>%
  filter(Trial != "PYNB") %>%
  filter(Trial != "GSPYN") %>%
  filter(Trial != "DHPYN") %>%
  group_by(Year, Location, Trial, trait_id)

harmonic18<- 1/mean(1/pheno18$rep) 

rm(phenoMaxRep18)

###### Calculating heritability

pheno17_H2 <- pheno17 %>%
  nest() %>%
  mutate(
    model = map(
      data,
      ~ lmer(
        phenotype_value ~ (1 | Variety) + (1 | rep) + (1 | column) +
          (1 | range),
        data = .
      )
    ),
    VarianceComp = map(model, ~ as.data.frame(lme4::VarCorr(.x)))
  ) %>%
  unnest(VarianceComp)

pheno17_H2 <- pheno17_H2 %>%
  group_by(Year, Location, Trial, trait_id) %>%
  select(-var1, -var2, -sdcor) %>%
  spread(key = grp, value = vcov) %>%
  mutate(H2 = Variety / (Variety + Residual / harmonic17)) %>%
  separate(trait_id, c("trait_id", "Date"), sep = "_") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  mutate(Date = replace_na(Date, "2017-06-13")) %>%
  select(-data, -model)

pheno17_H2 %>%
  ggplot(aes(x = Date, y = H2, colour = Trial)) +
  geom_point(size = 2) +
  facet_wrap(trait_id ~ Location, scales = "free", ncol = 3) +
  scale_x_date(breaks = "1 week", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Broad-sense Heritability lme4",
    subtitle = "2016/2017 season"
  )

pheno18_H2 <- pheno18 %>%
  nest() %>%
  mutate(
    model = map(
      data,
      ~ lmer(phenotype_value ~
      (1 | Variety) + (1 | rep) + (1 | column) + (1 | range),
      data = .
      )
    ),
    VarianceComp = map(model, ~ as.data.frame(lme4::VarCorr(.x)))
  ) %>%
  unnest(VarianceComp)

pheno18_H2 <- pheno18_H2 %>%
  group_by(Year, Location, Trial, trait_id) %>%
  select(-var1, -var2, -sdcor) %>%
  spread(key = grp, value = vcov) %>%
  mutate(H2 = Variety / (Variety + Residual / harmonic18)) %>%
  select(-data, -model)

pheno18_H2 %>%
  separate(trait_id, c("trait_id", "Date"), sep = "_") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  ggplot(aes(x = Date, y = H2, colour = Trial)) +
  geom_point(size = 2) +
  facet_wrap(trait_id ~ Location, scales = "free", ncol = 4) +
  scale_x_date(breaks = "3 weeks", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Broad-sense Heritability lme4",
    subtitle = "2017/2018 season"
  )

write_delim(pheno17_H2, "./PhenoDatabase/Heritabilities_lme4_2017.txt",
  delim = "\t"
)
write_delim(pheno18_H2, "./PhenoDatabase/Heritabilities_lme4_2018.txt",
  delim = "\t"
)

##### ASREML version ####
library(asreml)

asreml.license.status()

pheno17asreml <- pheno17 %>%
  ungroup() %>%
  separate(trait_id, c("trait_id", "phenotype_date"), sep = "_") %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = format(phenotype_date, "%Y%m%d")
  ) %>%
  unite("trait_id", c("trait_id", "phenotype_date"), sep = "_") %>%
  unite("trait_id", c("trait_id", "Trial"), sep = "_") %>%
  unite("trait_id", c("trait_id", "Location"), sep = "_") %>%
  spread(key = trait_id, value = phenotype_value) %>%
  mutate(
    Variety = as.factor(Variety),
    rep = as.factor(rep),
    range = as.factor(range),
    column = as.factor(column)
  )

pheno18asreml <- pheno18 %>%
  ungroup() %>%
  separate(trait_id, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(
    trait_id != "CanopyHeight",
    trait_id != "PTHT"
  ) %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = format(phenotype_date, "%Y%m%d")
  ) %>%
  unite("trait_id", c("trait_id", "phenotype_date"), sep = "_") %>%
  unite("trait_id", c("trait_id", "Trial"), sep = "_") %>%
  unite("trait_id", c("trait_id", "Location"), sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value) %>%
  mutate(
    Variety = as.factor(Variety),
    rep = as.factor(rep),
    range = as.factor(range),
    column = as.factor(column)
  )

## 2017
t17 <- asreml(
  fixed = GNDVI_20170506_AYN2_MP ~ 1,
  random = ~ Variety + rep + range + column,
  data = pheno17asreml
)
plot(t17)
coef(t17)$random
fitted(t17)
summary(t17)
resid(t17)

h <- as.data.frame(summary(t17)$varcomp)
h

h2 <- as.data.frame(h[4, 1] / (h[4, 1] + (h[5, 1] / harmonic17)))
h2

## 2017 all
effectvars <- names(pheno17asreml) %in% c(
  "block", "rep", "Variety", "year", "column",
  "range", "Plot_ID", "entity_id", "Year",
  "Treated", "OverallTrial", "Count"
)
traits <- colnames(pheno17asreml[, !effectvars])
H2_2017 <- data.frame(traits)
H2_2017$Heritability <- NA
fieldInfo <- pheno17asreml %>%
  tidylog::select(Variety, rep, column, range)
ntraits <- 1:nrow(H2_2017)

for (i in ntraits) {
  print(paste("Working on trait", H2_2017[i, 1]))
  j <- H2_2017[i, 1]

  data <- cbind(fieldInfo, pheno17asreml[, paste(j)])
  names(data) <- c("Variety", "rep", "column", "range", "Trait")

  t17 <- asreml(
    fixed = Trait ~ 1,
    random = ~ Variety + rep + range + column,
    data = data
  )
  pdf(paste0(
    "./Figures/AsremlPlots/ASREML_repRangeColumn17_", H2_2017[i, 1],
    ".pdf"
  ))
  plot(t17)
  dev.off()
  h <- as.data.frame(summary(t17)$varcomp)
  print(paste("Creating Data Frame", j))
  print(h)

  h2 <- (h[4, 1] / (h[4, 1] + (h[5, 1] / harmonic17)))
  h2
  H2_2017[i, 2] <- h2
}

dev.off()

## 2018
t17 <- asreml(
  fixed = GNDVI_20170506_AYN2_MP ~ 1,
  random = ~ Variety + rep + range + column,
  data = pheno17asreml
)
plot(t17)
coef(t17)$random
fitted(t17)
summary(t17)
resid(t17)

h <- as.data.frame(summary(t17)$varcomp)
h

h2 <- as.data.frame(h[4, 1] / (h[4, 1] + (h[5, 1] / harmonic18)))
h2

## 2018 all
effectvars <- names(pheno18asreml) %in% c(
  "block", "rep", "Variety", "year", "column",
  "range", "Plot_ID", "entity_id", "Year",
  "Treated", "OverallTrial", "Count", "max"
)
traits <- colnames(pheno18asreml[, !effectvars])
H2_2018 <- data.frame(traits)
H2_2018$Heritability <- NA
fieldInfo <- pheno18asreml %>%
  tidylog::select(Variety, rep, column, range)
ntraits <- 1:nrow(H2_2018)

for (i in ntraits) {
  print(paste("Working on trait", H2_2018[i, 1]))
  j <- H2_2018[i, 1]

  data <- cbind(fieldInfo, pheno18asreml[, paste(j)])
  names(data) <- c("Variety", "rep", "column", "range", "Trait")

  t18 <- asreml(
    fixed = Trait ~ 1,
    random = ~ Variety + rep + range + column,
    data = data
  )
  pdf(paste0(
    "./Figures/AsremlPlots/ASREML_repRangeColumn18_", H2_2018[i, 1],
    ".pdf"
  ))
  plot(t17)
  dev.off()
  h <- as.data.frame(summary(t18)$varcomp)
  print(paste("Creating Data Frame", j))
  print(h)

  h2 <- (h[4, 1] / (h[4, 1] + (h[5, 1] / harmonic18)))
  h2
  H2_2018[i, 2] <- h2
}

dev.off()

H2_2017 <- H2_2017 %>%
  separate(
    col = traits,
    into = c("trait_id", "phenotype_date", "trial", "location"),
    sep = "_"
  ) %>%
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"))

H2_2018 <- H2_2018 %>%
  separate(
    col = traits,
    into = c("trait_id", "phenotype_date", "trial", "location"),
    sep = "_"
  ) %>%
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"))

H2_2017 %>%
  filter(trait_id != "GRYLD") %>%
  ggplot(aes(x = phenotype_date, y = Heritability, colour = trial)) +
  geom_point(size = 2) +
  facet_wrap(c("trait_id", "location"), scales = "free", ncol = 3) +
  scale_x_date(date_breaks = "14 days", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Broad-sense Heritability Asreml",
    subtitle = "2016/2017 season"
  )

H2_2018 %>%
  filter(trait_id != "GRYLD") %>%
  ggplot(aes(x = phenotype_date, y = Heritability, colour = trial)) +
  geom_point(size = 2) +
  facet_wrap(c("trait_id", "location"), scales = "free", ncol = 4) +
  scale_x_date(date_breaks = "20 days", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Broad-sense Heritability Asreml",
    subtitle = "2017/2018 season"
  )

write_delim(H2_2017, "./PhenoDatabase/Heritabilities_asreml_2017.txt",
  delim = "\t"
)
write_delim(H2_2018, "./PhenoDatabase/Heritabilities_asreml_2018.txt",
  delim = "\t"
)

