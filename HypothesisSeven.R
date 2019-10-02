rm(list = objects())
ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(broom)
library(asreml)
library(ggrepel)
library(Polychrome)

set.seed(1990)
asreml.license.status()

#### Set up work space ####

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
    legend.key.size = unit(1, "lines"),
    # legend.background = element_rect(fill = NULL, colour = NULL),
    # legend.box = NULL,
    legend.margin = margin(
      t = 0, r = 0.75,
      b = 0, l = 0.75,
      unit = "cm"
    ),
    legend.text = element_text(size = rel(1)),
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
      colour = "#512888",
      size = rel(3),
      vjust = 0,
      hjust = 0,
      margin = margin(
        t = 0.25, r = 0.25,
        b = 0.5, l = 0.25,
        unit = "cm"
      )
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

#### Read in data
pheno17 <- fread("./PhenoDatabase/PhenoLong_vi17.txt")

pheno18 <- fread("./PhenoDatabase/PhenoVI_18.txt")

#### Move from long to wide
pheno17 <- pheno17 %>%
  spread(key = trait_id, value = phenotypic_value) %>%
  mutate(
    Location = as.factor(Location),
    Variety = as.factor(Variety),
    Trial = as.factor(Trial),
    rep = as.factor(rep),
    range = as.factor(range),
    column = as.factor(column)
  ) %>%
  glimpse()

pheno18 <- pheno18 %>%
  filter(trait_id != "PTHT") %>%
  unite("trait_id", trait_id, phenotype_date) %>%
  spread(key = trait_id, value = phenotype_value) %>%
  mutate(
    Location = as.factor(Location),
    Variety = as.factor(Variety),
    Trial = as.factor(Trial),
    rep = as.factor(rep),
    range = as.factor(range),
    column = as.factor(column)
  ) %>%
  glimpse()

#### Asreml trial

t17 <- asreml(
  fixed = GRYLD ~ 1,
  random = ~ Variety
  + Location
    + Location:Variety
    + Location:Trial
    + Location:Trial:range
    + Location:Trial:column,
  data = pheno17
)

plot(t17)
summary.asreml(t17)

blues <- setDT(as.data.frame(coef(t17)$random), keep.rownames = T)
blues$rn <- str_remove(blues$rn, "Variety_")
dat17 <- blues %>%
  rename(GRYLD = effect) %>%
  glimpse()

colnames(pheno17)
# pheno17<- pheno17 %>%
#   select(-`Nir_2017-06-16`,-`RedEdge_2017-05-24`,-`RedEdge_2017-06-16`)

##### Generating all 2017 VI BLUEs

effectvars <- names(pheno17) %in% c(
  "entity_id", "Variety", "Year", "Trial",
  "Location", "Treated", "rep", "range", "column",
  "OverallTrial", "GRYLD"
)
traits <- colnames(pheno17[, !effectvars])
traits
fieldInfo <- pheno17 %>%
  tidylog::select(Variety, Location, Trial, range, column)



for (i in traits) {
  print(paste("Working on trait", i))

  data <- cbind(fieldInfo, pheno17[, paste(i)])
  names(data) <- c("Variety", "Location", "Trial", "range", "column", "Trait")
  print(colnames(data))

  t17 <- asreml(
    fixed = Trait ~ 1,
    random = ~ Variety + Trial:column + Trial:range,
    data = data
  )

  pdf(paste0("./Figures/AsremlPlots/ASREML_TrialColumnRange_17_", i, ".pdf"))
  plot(t17)
  dev.off()
  print(summary(t17))
  blues <- setDT(as.data.frame(coef(t17)$random), keep.rownames = T)
  blues$rn <- str_remove(blues$rn, "Variety_")
  colnames(blues)[colnames(blues) == "effect"] <- paste(i)
  dat17 <- blues %>%
    dplyr::inner_join(dat17, by = "rn")
}

warnings()

#### Asreml trial
pheno18 <- pheno18 %>%
  clean_names()

t18 <- asreml(
  fixed = gryld_2018_06_13 ~ 1,
  random = ~ variety
  + location
    + location:variety
    + location:trial
    + location:trial:range
    + location:trial:column,
  data = pheno18
)

plot(t18)
summary.asreml(t18)

blues <- setDT(as.data.frame(coef(t18)$random), keep.rownames = T)
blues$rn <- str_remove(blues$rn, "variety_")
dat18 <- blues %>%
  rename(gryld_2018_06_13 = effect) %>%
  glimpse()

colnames(pheno18)

##### Generating all 2017 VI BLUEs

effectvars <- names(pheno18) %in% c(
  "entity_id", "variety", "year", "trial",
  "location", "treated", "rep", "range", "column",
  "overall_trial", "gryld_2018_06_13"
)
traits <- colnames(pheno18[, !effectvars])
traits
fieldInfo <- pheno18 %>%
  tidylog::select(variety, location, trial, range, column)

for (i in traits) {
  print(paste("Working on trait", i))

  data <- cbind(fieldInfo, pheno18[, paste(i)])
  names(data) <- c("Variety", "Location", "Trial", "range", "column", "Trait")
  print(colnames(data))

  t18 <- asreml(
    fixed = Trait ~ 1,
    random = ~ Variety + Trial:column + Trial:range,
    data = data
  )

  pdf(paste0("./Figures/AsremlPlots/ASREML_TrialColumnRange_18_", i, ".pdf"))
  plot(t18)
  dev.off()
  print(summary(t18))
  blues <- setDT(as.data.frame(coef(t18)$random), keep.rownames = T)
  blues$rn <- str_remove(blues$rn, "Variety_")
  colnames(blues)[colnames(blues) == "effect"] <- paste(i)
  dat18 <- blues %>%
    inner_join(dat18)
}

warnings()

write.table(dat17, "./PhenoDatabase/BLUPs_TrialColumnRange_2017.txt",
  quote = F, sep = "\t",
  row.names = F, col.names = T
)
write.table(dat18, "./PhenoDatabase/BLUPs_TrialColumnRange_2018.txt",
  quote = F, sep = "\t",
  row.names = F, col.names = T
)

#### Calculating Blups for individual trait date trial location ####

pheno17 <- data.table::fread("./PhenoDatabase/PhenoLong_vi17.txt")
str(pheno17)

pheno18 <- data.table::fread("./PhenoDatabase/PhenoVI_18.txt")
str(pheno18)

# Use only Varities with more than one rep per season

pheno17 <- pheno17 %>%
  spread(key = trait_id, value = phenotypic_value) %>%
  distinct() %>%
  glimpse() %>% 
  gather(
    key = trait_id, value = phenotype_value,
    GRYLD:`RedEdge_2017-06-16`
  ) %>%
  drop_na(phenotype_value) %>%
  group_by(Year, Location, Trial, trait_id) %>%
  arrange(Year, Location, Trial, Variety, trait_id) 

pheno18 <- pheno18 %>%
  unite("trait_id", trait_id, phenotype_date, sep = "_") %>% 
  group_by(Year, Location, Trial, trait_id)

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
blues <- setDT(as.data.frame(coef(t17)$random), keep.rownames = T)
blues$rn <- str_remove(blues$rn, "Variety_")
dat17 <- blues %>%
  rename(GNDVI_20170506_AYN2_MP = effect) %>%
  glimpse()

## 2017 all
effectvars <- names(pheno17asreml) %in% c(
  "block", "rep", "Variety", "year", "column",
  "range", "Plot_ID", "entity_id", "Year",
  "Treated", "OverallTrial", "Count"
)
traits <- colnames(pheno17asreml[, !effectvars])
traits
fieldInfo <- pheno17asreml %>%
  tidylog::select(Variety, range, column, rep)

for (i in traits) {
  print(paste("Working on trait", i))
  
  data <- cbind(fieldInfo, pheno17asreml[, paste(i)])
  names(data) <- c("Variety", "range", "column", "rep", "Trait")
  print(colnames(data))
  
  t17 <- asreml(
    fixed = Trait ~ 1,
    random = ~ Variety + rep + range + column,
    data = data
  )
  pdf(paste0(
    "./Figures/AsremlPlots/ASREML_repRangeColumn17_", i,
    ".pdf"
  ))
  plot(t17)
  blues <- setDT(as.data.frame(coef(t17)$random), keep.rownames = T)
  blues$rn <- str_remove(blues$rn, "Variety_")
  colnames(blues)[colnames(blues) == "effect"] <- paste(i)
  dat17 <- blues %>%
    dplyr::inner_join(dat17, by = "rn")
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
blues <- setDT(as.data.frame(coef(t17)$random), keep.rownames = T)
blues$rn <- str_remove(blues$rn, "Variety_")
dat18 <- blues %>%
  rename(GRYLD = effect) %>%
  glimpse()

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
  plot(t18)
  dev.off()
  plot(t18)
  blues <- setDT(as.data.frame(coef(t18)$random), keep.rownames = T)
  blues$rn <- str_remove(blues$rn, "Variety_")
  colnames(blues)[colnames(blues) == "effect"] <- paste(i)
  dat18 <- blues %>%
    dplyr::inner_join(dat18, by = "rn")
}

dev.off()


#### PCA of phenotypes
dat17pca <- as.matrix(t(dat17[, -1]))

pcaMethods::checkData(data = dat17pca)

pcaDat17 <- prcomp(dat17pca)
summary(pcaDat17)

scores <- setDT(as.data.frame(pcaDat17$x), keep.rownames = T)
scores <- scores %>%
  separate(rn, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(trait_id != "Height") %>%
  mutate(phenotype_date = replace_na(phenotype_date, "2017-06-20"))

ggplot(data = scores, aes(
  x = PC1, y = PC2, colour = phenotype_date,
  shape = trait_id
)) +
  geom_point() +
  theme(aspect.ratio = 1:1) +
  scale_shape_manual(values = c(0, 15, 1, 2, 8, 11, 9)) +
  scale_colour_manual(values = c(
    "#875692", "#f38400", "#a1caf1", "#be0032",
    "#848482", "#008856", "#e68fac", "#0067a5",
    "#604e97", "#f6a600", "#222222"
  )) +
  labs(
    title = "PCA of Phenotypic BLUPs",
    x = "PC1 = 52.2%",
    y = "PC2 = 14.44%"
  )

ggplot(data = scores, aes(
  x = PC1, y = PC3, colour = phenotype_date,
  shape = trait_id
)) +
  geom_point() +
  theme(aspect.ratio = 1:1) +
  scale_shape_manual(values = c(0, 15, 1, 2, 8, 11, 9)) +
  scale_colour_manual(values = c(
    "#875692", "#f38400", "#a1caf1", "#be0032",
    "#848482", "#008856", "#e68fac", "#0067a5",
    "#604e97", "#f6a600", "#222222"
  )) +
  labs(
    title = "PCA of Phenotypic BLUPs",
    x = "PC1 = 52.2%",
    y = "PC3 = 10.96%"
  )

ggplot(data = scores, aes(
  x = PC2, y = PC3, colour = phenotype_date,
  shape = trait_id
)) +
  geom_point() +
  theme(aspect.ratio = 1:1) +
  scale_shape_manual(values = c(0, 15, 1, 2, 8, 11, 9)) +
  scale_colour_manual(values = c(
    "#875692", "#f38400", "#a1caf1", "#be0032",
    "#848482", "#008856", "#e68fac", "#0067a5",
    "#604e97", "#f6a600", "#222222"
  )) +
  labs(
    title = "PCA of Phenotypic BLUPs",
    x = "PC2 = 14.44%",
    y = "PC3 = 10.96%"
  )

biplot(pcaDat17)

# ggbiplot2(pcaDat17)

dat18pca <- as.matrix(t(dat18[, -1]))

pcaMethods::checkData(data = dat18pca)

pcaDat18 <- prcomp(dat18pca)
summary(pcaDat18)

scores <- setDT(as.data.frame(pcaDat18$x), keep.rownames = T)
scores <- scores %>%
  separate(rn, c(
    "trait_id", "phenotype_year", "phenotype_month",
    "phenotype_day"
  ), sep = "_") %>%
  filter(trait_id != "canopy") %>%
  unite("phenotype_date", phenotype_year, phenotype_month, phenotype_day,
    sep = "-"
  )

ggplot(data = scores, aes(
  x = PC1, y = PC2, colour = phenotype_date,
  shape = trait_id
)) +
  geom_point() +
  theme(aspect.ratio = 1:1) +
  scale_shape_manual(values = c(0, 15, 1, 2, 8, 11, 9)) +
  scale_colour_manual(values = c(
    "#f3c300", "#875692", "#f38400", "#a1caf1",
    "#be0032", "#848482", "#008856", "#e68fac",
    "#0067a5", "#604e97", "#f6a600", "#b3446c",
    "#882d17", "#8db600", "#e25822", "#222222"
  )) +
  labs(
    title = "PCA of Phenotypic BLUPs",
    x = "PC1 = 38.91%",
    y = "PC2 = 14.25%"
  )

ggplot(data = scores, aes(
  x = PC1, y = PC3, colour = phenotype_date,
  shape = trait_id
)) +
  geom_point() +
  theme(aspect.ratio = 1:1) +
  scale_shape_manual(values = c(0, 15, 1, 2, 8, 11, 9)) +
  scale_colour_manual(values = c(
    "#f3c300", "#875692", "#f38400", "#a1caf1",
    "#be0032", "#848482", "#008856", "#e68fac",
    "#0067a5", "#604e97", "#f6a600", "#b3446c",
    "#882d17", "#8db600", "#e25822", "#222222"
  )) +
  labs(
    title = "PCA of Phenotypic BLUPs",
    x = "PC1 = 38.91%",
    y = "PC3 = 10.73%"
  )

ggplot(data = scores, aes(
  x = PC2, y = PC3, colour = phenotype_date,
  shape = trait_id
)) +
  geom_point() +
  theme(aspect.ratio = 1:1) +
  scale_shape_manual(values = c(0, 15, 1, 2, 8, 11, 9)) +
  scale_colour_manual(values = c(
    "#f3c300", "#875692", "#f38400", "#a1caf1",
    "#be0032", "#848482", "#008856", "#e68fac",
    "#0067a5", "#604e97", "#f6a600", "#b3446c",
    "#882d17", "#8db600", "#e25822", "#222222"
  )) +
  labs(
    title = "PCA of Phenotypic BLUPs",
    x = "PC2 = 14.25%",
    y = "PC3 = 10.73%"
  )

biplot(pcaDat18)

# ggbiplot2(pcaDat18)

##### GRYLD all years

gryld17 <- pheno17 %>%
  select(entity_id:GRYLD) %>%
  clean_names()
gryld18 <- pheno18 %>%
  select(location:overall_trial, gryld_2018_06_13) %>%
  dplyr::rename(gryld = gryld_2018_06_13)

gryld <- gryld17 %>%
  bind_rows(gryld18) %>%
  mutate(
    variety = as.factor(variety),
    year = as.factor(year),
    location = as.factor(location),
    trial = as.factor(trial),
    range = as.factor(range),
    column = as.factor(column)
  )

gryldBlup <- asreml(
  fixed = gryld ~ 1,
  random = ~ variety
  + year
    + year:location
    + year:location:trial,
  data = gryld
)
plot(gryldBlup)
summary(gryldBlup)

gryldBLUP <- setDT(as.data.frame(coef(gryldBlup)$random), keep.rownames = T)
gryldBLUP <- gryldBLUP %>%
  filter(str_detect(rn, "variety_")) %>%
  tidylog::mutate(rn = str_remove(rn, "variety_"))

write.table(gryldBLUP, "./PhenoDatabase/GRYLD_Blups_allyrs.txt",
  quote = F,
  sep = "\t", row.names = F, col.names = T
)
