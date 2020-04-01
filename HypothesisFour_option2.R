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

pheno <- fread("./PhenoDatabase/Phenotype_Master.txt")
glimpse(pheno)

phenoVI_17 <- fread("./PhenoDatabase/PhenoLong_vi17.txt")
glimpse(phenoVI_17)
varieties17_VI<- tibble(varietyVI_17 = unique(phenoVI_17$Variety))

phenoVI_18 <- fread("./PhenoDatabase/PhenoVI_18.txt")
glimpse(phenoVI_18)
varieties18_VI<- tibble(varietyVI_18 = unique(phenoVI_18$Variety))

pheno17 <- pheno %>%
  filter(year == "17")
varieties17_db<- tibble(varietydb_17 = unique(pheno17$Variety))

pheno18 <- pheno %>%
  filter(year == "18") 
varieties18_db<- tibble(varietydb_18 = unique(pheno18$Variety))

## Name comparison between db and VI information

varieties17<- varieties17_db %>% 
  left_join(varieties17_VI, by = c("varietydb_17" = "varietyVI_17"))
varieties18<- varieties18_db %>% 
  left_join(varieties18_VI, by = c("varietydb_18" = "varietyVI_18"))

# Variety names seem to match. Problem may be coming from different entity_id
# or possibly wrong field information. 

phenoVI_17 <- phenoVI_17 %>%
  filter(trait_id != "GRYLD") %>%
  select(entity_id, trait_id, phenotypic_value) %>%
  pivot_wider(
    names_from = trait_id,
    values_from = phenotypic_value
  )

missing17_pheno <- pheno17 %>%
  anti_join(phenoVI_17)
missing17_vi <- phenoVI_17 %>%
  anti_join(pheno17)

pheno17 <- pheno17 %>%
  inner_join(phenoVI_17, by = c("entity_id")) %>%
  pivot_wider(
    names_from = trait_id,
    values_from = phenotype_value
  ) %>%
  select(
    entity_id, Variety, range, column, year, trial, location, treated,
    plot, Trial, block, rep, everything()
  ) %>%
  pivot_longer(
    cols = `NDVI_2017-05-25`:GRYLD, names_to = "trait_ID",
    values_to = "phenotype_value"
  ) %>%
  mutate(
    Variety = as.factor(Variety),
    range = as.factor(range),
    column = as.factor(column),
    trial = as.factor(trial),
    location = as.factor(location),
    treated = as.factor(treated),
    Trial = as.factor(Trial),
    block = as.factor(block),
    rep = as.factor(rep),
    phenotype_date = as.factor(phenotype_date)
  )

phenoVI_18 <- phenoVI_18 %>%
  filter(trait_id != "GRYLD") %>%
  filter(trait_id != "PTHT") %>%
  unite("trait_id", trait_id, phenotype_date) %>%
  select(entity_id, trait_id, phenotype_value) %>%
  pivot_wider(
    names_from = trait_id,
    values_from = phenotype_value
  )

missing18_pheno <- pheno18 %>%
  anti_join(phenoVI_18)
missing18_vi <- phenoVI_18 %>%
  anti_join(pheno18)

pheno18 <- pheno18 %>%
  inner_join(phenoVI_18, by = c("entity_id")) %>%
  pivot_wider(
    names_from = trait_id,
    values_from = phenotype_value
  ) %>%
  select(
    entity_id, Variety, range, column, year, trial, location, treated,
    plot, Trial, block, rep, everything()
  ) %>%
  pivot_longer(
    cols = `CanopyHeight_2018-05-22`:GRYLD, names_to = "trait_ID",
    values_to = "phenotype_value"
  ) %>%
  mutate(
    Variety = as.factor(Variety),
    range = as.factor(range),
    column = as.factor(column),
    trial = as.factor(trial),
    location = as.factor(location),
    treated = as.factor(treated),
    Trial = as.factor(Trial),
    block = as.factor(block),
    rep = as.factor(rep)
  )



#### Trait Summaries ####
traitSummary17 <- pheno17 %>%
  group_by(location, Trial, treated, trait_ID) %>%
  mutate(Count = row_number()) %>%
  summarise(
    Total = max(Count),
    minimum = min(phenotype_value),
    average = mean(phenotype_value),
    maximum = max(phenotype_value),
    variance = var(phenotype_value)
  ) %>%
  drop_na(minimum)

traitSummary18 <- pheno18 %>%
  group_by(location, Trial, treated, trait_ID) %>%
  mutate(Count = row_number()) %>%
  summarise(
    Total = max(Count),
    minimum = min(phenotype_value),
    average = mean(phenotype_value),
    maximum = max(phenotype_value),
    variance = var(phenotype_value)
  ) %>%
  drop_na(minimum)

varietySummary17 <- pheno17 %>%
  group_by(Variety, Trial, treated, location) %>%
  mutate(Count = row_number()) %>%
  summarise(Total = max(Count))

varietySummary18 <- pheno18 %>%
  group_by(Variety, Trial, treated, location) %>%
  mutate(Count = row_number()) %>%
  summarise(Total = max(Count))

pheno17 <- pheno17 %>%
  separate(col = trait_ID, into = c("trait_ID", "Date"), sep = "_") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  filter(trait_ID != "Height") %>%
  drop_na(phenotype_value)

dateLocationSummary_17 <- distinct(pheno17, Date, location)

pheno17_ayn <- pheno17 %>%
  filter(Trial == "AYN")

pheno17_pyn <- pheno17 %>%
  filter(Trial == "PYN")

pheno18 <- pheno18 %>%
  separate(col = trait_ID, into = c("trait_ID", "Date"), sep = "_") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  filter(trait_ID != "CanopyHeight") %>%
  drop_na(phenotype_value)
dateLocationSummary_18 <- distinct(pheno18, Date, location)

pheno18_ayn <- pheno18 %>%
  filter(Trial == "AYN")

pheno18_pyn <- pheno18 %>%
  filter(Trial == "PYN")

pheno18_pyn %>%
  filter(trait_ID != "GRYLD") %>%
  filter(trait_ID != "CanopyHeight") %>%
  ggplot(aes(x = Date, y = phenotype_value, colour = location)) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 1, varwidth = TRUE,
    alpha = 0.5, size = 0.5
  ) +
  facet_wrap(~trait_ID, scales = "free_y") +
  labs(
    title = "Distribution of Vegetation Indices over season",
    subtitle = "2017-2018 Season, Preliminary Yield Nursery"
  )

write.table(pheno17, "./PhenoDatabase/AllPheno_17.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(missing17_pheno, "./PhenoDatabase/missingVI_17.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(pheno18, "./PhenoDatabase/AllPheno_18.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

###############################################################################
##### Heritability #####

library(lmerTest)

pheno17_ayn <- pheno17_ayn %>%
  unite("trait_ID", trait_ID, Date, location, sep = "_") %>%
  group_by(trait_ID) %>%
  mutate(
    group = "group",
    ID = group_indices()
  ) %>%
  unite("ID", group, ID, sep = "_")

keyData_17 <- pheno17_ayn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_17

pheno17_ayn <- pheno17_ayn %>%
  pivot_wider(names_from = ID, values_from = phenotype_value)

calcH2r <- function(dat, fill = NA, ...) {
  r <- length(table(dat$rep))

  effectvars <- names(dat) %in% c(
    "entity_id", "Variety", "range", "column",
    "year", "trial", "trait_ID", "treated", "plot",
    "Trial", "block", "rep", "phenotype_date", "location",
    "phenotype_person"
  )

  t <- colnames(dat[, !effectvars])

  res <- tibble(
    Trait = t, H = NA, varVariety = NA, varRep = NA,
    varRepBlock = NA, varRes = NA
  )
  ntraits <- 1:nrow(res)

  for (i in ntraits) {
    j <- res[i, 1]

    h <- as.data.frame(VarCorr(
      lmer(paste0(j, "~ (1|Variety) + (1|rep) + (1|rep:block)"), data = dat)
    ))
    H2 <- h[1, 4] / (h[1, 4] + (h[4, 4] / r))

    res[i, 2] <- H2
    res[i, 3] <- round(h[1, 4], digits = 5)
    res[i, 4] <- round(h[3, 4], digits = 5)
    res[i, 5] <- round(h[2, 4], digits = 5)
    res[i, 6] <- round(h[4, 4], digits = 5)
  }
  return(res)
}

ayn_17 <- calcH2r(dat = pheno17_ayn)

ayn_17 <- ayn_17 %>%
  inner_join(keyData_17, by = c("Trait" = "ID")) %>%
  select(trait_ID, -Trait, everything()) %>%
  separate(trait_ID, c("trait_ID", "Date", "location"), sep = "_") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

pheno18_ayn <- pheno18_ayn %>%
  unite("trait_ID", trait_ID, Date, location, sep = "_") %>%
  group_by(trait_ID) %>%
  mutate(
    group = "group",
    ID = group_indices()
  ) %>%
  unite("ID", group, ID, sep = "_")

keyData_18 <- pheno18_ayn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_18

pheno18_ayn <- pheno18_ayn %>%
  pivot_wider(names_from = ID, values_from = phenotype_value)

ayn_18 <- calcH2r(dat = pheno18_ayn)

ayn_18 <- ayn_18 %>%
  inner_join(keyData_18, by = c("Trait" = "ID")) %>%
  select(trait_ID, -Trait, everything()) %>%
  separate(trait_ID, c("trait_ID", "Date", "location"), sep = "_") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

ayn_18 %>%
  filter(trait_ID != "GRYLD") %>%
  ggplot(aes(x = Date, y = H, colour = as.factor(location))) +
  geom_point() +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  facet_wrap(~trait_ID, scales = "free") +
  labs(title = "Heritability for 2018 AYN ")

write.table(ayn_17, "./Results/H2_ayn_2017.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(ayn_18, "./Results/H2_ayn_2018.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

###############################################################################
##### BLUPs ######
library(asreml)
asreml.license.status()

## 2017
# Pyn
pheno17_pyn <- pheno17_pyn %>%
  filter(location != "RP") %>%
  unite("trait_ID", trait_ID, Date, location, sep = "_") %>%
  group_by(trait_ID) %>%
  mutate(
    group = "group",
    ID = dplyr::group_indices()
  ) %>%
  unite("ID", group, ID, sep = "_") %>%
  select(-phenotype_person, -range)

glimpse(pheno17_pyn)

keyData_17 <- pheno17_pyn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_17

pheno17_pyn <- pheno17_pyn %>%
  pivot_wider(names_from = ID, values_from = phenotype_value) %>%
  ungroup()

fit <- asreml(
  fixed = group_1 ~ 1,
  random = ~ Variety + block,
  data = pheno17_pyn
)
plot(fit)

blups <- setDT(as.data.frame(coef(fit)$random), keep.rownames = T)
blups$rn <- str_remove(blups$rn, "Variety_")
colnames(blups)[colnames(blups) == "rn"] <- "Variety"
blups <- blups %>%
  select(-effect)

# As a function
pyn_blups_2017 <- function(dat, joinFile, saveFile, ...) {
  effectvars <- names(dat) %in% c(
    "entity_id", "Variety", "column", "year", "trial",
    "trait_ID", "treated", "plot", "Trial", "block",
    "rep", "phenotype_date", "entry", "new", "entryc"
  )
  traits <- colnames(dat[, !effectvars])
  traits
  fieldInfo <- dat %>%
    dplyr::select("Variety", "rep", "block", "column")

  for (i in traits) {
    print(paste("Working on trait", i))

    data <- cbind(fieldInfo, dat[, paste(i)])
    names(data) <- c("Variety", "rep", "block", "column", "Trait")
    data <- data %>% 
      drop_na(Trait)

    t_fit <- asreml(
      fixed = Trait ~ 1,
      random = ~ Variety + block,
      data = data
    )
    pdf(paste0(
      saveFile,
      i, "_pyn.pdf"
    ))
    plot(t_fit)

    blups <- setDT(as.data.frame(coef(t_fit)$random), keep.rownames = T)
    blups$rn <- str_remove(blups$rn, "Variety_")
    colnames(blups)[colnames(blups) == "rn"] <- "Variety"
    colnames(blups)[colnames(blups) == "effect"] <- paste(i)

    joinFile <- blups %>%
      inner_join(joinFile, by = "Variety")
    graphics.off()
  }
  return(joinFile)
}

blups17_pyn <- pyn_blups_2017(
  dat = pheno17_pyn,
  joinFile = blups,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2017/"
)

blups17_pyn <- blups17_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  inner_join(keyData_17, by = "ID") %>%
  select(Variety, trait_ID, phenotype_value) %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "block"))

# Ayn

glimpse(pheno17_ayn)

fit <- asreml(
  fixed = group_1 ~ 1,
  random = ~ Variety + rep + rep:block,
  data = pheno17_ayn
)
plot(fit)

blups <- setDT(as.data.frame(coef(fit)$random), keep.rownames = T)
blups$rn <- str_remove(blups$rn, "Variety_")
colnames(blups)[colnames(blups) == "rn"] <- "Variety"
blups <- blups %>%
  select(-effect)

# As a function
ayn_blups_2017 <- function(dat, joinFile, saveFile, ...) {
  effectvars <- names(dat) %in% c(
    "entity_id", "Variety", "range", "column", "year",
    "trial", "trait_ID", "treated", "plot", "Trial",
    "block", "rep", "phenotype_date", "phenotype_person"
  )
  traits <- colnames(dat[, !effectvars])
  traits
  fieldInfo <- dat %>%
    ungroup() %>%
    dplyr::select("Variety", "rep", "block", "column")

  for (i in traits) {
    print(paste("Working on trait", i))

    data <- cbind(fieldInfo, dat[, paste(i)])
    names(data) <- c("Variety", "rep", "block", "column", "Trait")
    data <- data %>% 
      drop_na(Trait)

    t_fit <- asreml(
      fixed = Trait ~ 1,
      random = ~ Variety + rep + rep:block,
      data = data
    )
    pdf(paste0(
      saveFile,
      i, "_ayn.pdf"
    ))
    plot(t_fit)

    blups <- setDT(as.data.frame(coef(t_fit)$random), keep.rownames = T)
    blups$rn <- str_remove(blups$rn, "Variety_")
    colnames(blups)[colnames(blups) == "rn"] <- "Variety"
    colnames(blups)[colnames(blups) == "effect"] <- paste(i)

    joinFile <- blups %>%
      inner_join(joinFile, by = "Variety")
    graphics.off()
  }
  return(joinFile)
}

blups17_ayn <- ayn_blups_2017(
  dat = pheno17_ayn,
  joinFile = blups,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2017/"
)

blups17_ayn <- blups17_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  inner_join(keyData_17, by = "ID") %>%
  select(Variety, trait_ID, phenotype_value) %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "rep"))

blups17_pyn <- blups17_pyn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2017-06-10")
    )
  )

blups17_ayn <- blups17_ayn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2017-06-10")
    )
  )

write.table(blups17_pyn, "./Results/Blups_pyn_17.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(blups17_ayn, "./Results/Blups_ayn_17.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

## 2018

pheno18_pyn <- pheno18_pyn %>%
  unite("trait_ID", trait_ID, Date, location, sep = "_") %>%
  group_by(trait_ID) %>%
  mutate(
    group = "group",
    ID = dplyr::group_indices()
  ) %>%
  unite("ID", group, ID, sep = "_") %>%
  select(-phenotype_person, -range)

glimpse(pheno18_pyn)

keyData_18 <- pheno18_pyn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_18

pheno18_pyn <- pheno18_pyn %>%
  pivot_wider(names_from = ID, values_from = phenotype_value) %>%
  ungroup()

fit <- asreml(
  fixed = group_1 ~ 1,
  random = ~ Variety + block,
  data = pheno18_pyn
)
plot(fit)

blups <- setDT(as.data.frame(coef(fit)$random), keep.rownames = T)
blups$rn <- str_remove(blups$rn, "Variety_")
colnames(blups)[colnames(blups) == "rn"] <- "Variety"
blups <- blups %>%
  select(-effect)

# As a function
pyn_blups_2018 <- function(dat, joinFile, saveFile, ...) {
  effectvars <- names(dat) %in% c(
    "entity_id", "Variety", "column", "year", "trial",
    "trait_ID", "treated", "plot", "Trial", "block",
    "rep", "phenotype_date", "entry", "new", "entryc"
  )
  traits <- colnames(dat[, !effectvars])
  traits
  fieldInfo <- dat %>%
    dplyr::select("Variety", "rep", "block", "column")

  for (i in traits) {
    print(paste("Working on trait", i))

    data <- cbind(fieldInfo, dat[, paste(i)])
    names(data) <- c("Variety", "rep", "block", "column", "Trait")
    data <- data %>% 
      drop_na(Trait)

    t_fit <- asreml(
      fixed = Trait ~ 1,
      random = ~ Variety + block,
      data = data
    )
    pdf(paste0(
      saveFile,
      i, "_pyn.pdf"
    ))
    plot(t_fit)

    blups <- setDT(as.data.frame(coef(t_fit)$random), keep.rownames = T)
    blups$rn <- str_remove(blups$rn, "Variety_")
    colnames(blups)[colnames(blups) == "rn"] <- "Variety"
    colnames(blups)[colnames(blups) == "effect"] <- paste(i)

    joinFile <- blups %>%
      inner_join(joinFile, by = "Variety")
    graphics.off()
  }
  return(joinFile)
}

blups18_pyn <- pyn_blups_2018(
  dat = pheno18_pyn,
  joinFile = blups,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2018/"
)

blups18_pyn <- blups18_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  inner_join(keyData_18, by = "ID") %>%
  select(Variety, trait_ID, phenotype_value) %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "block"))

# Ayn

glimpse(pheno18_ayn)

fit <- asreml(
  fixed = group_1 ~ 1,
  random = ~ Variety + rep + rep:block,
  data = pheno18_ayn
)
plot(fit)

blups <- setDT(as.data.frame(coef(fit)$random), keep.rownames = T)
blups$rn <- str_remove(blups$rn, "Variety_")
colnames(blups)[colnames(blups) == "rn"] <- "Variety"
blups <- blups %>%
  select(-effect)

# As a function
ayn_blups_2018 <- function(dat, joinFile, saveFile, ...) {
  effectvars <- names(dat) %in% c(
    "entity_id", "Variety", "range", "column", "year",
    "trial", "trait_ID", "treated", "plot", "Trial",
    "block", "rep", "phenotype_date", "phenotype_person"
  )
  traits <- colnames(dat[, !effectvars])
  traits
  fieldInfo <- dat %>%
    ungroup() %>%
    dplyr::select("Variety", "rep", "block", "column")

  for (i in traits) {
    print(paste("Working on trait", i))

    data <- cbind(fieldInfo, dat[, paste(i)])
    names(data) <- c("Variety", "rep", "block", "column", "Trait")
    data <- data %>% 
      drop_na(Trait)

    t_fit <- asreml(
      fixed = Trait ~ 1,
      random = ~ Variety + rep + rep:block,
      data = data
    )
    pdf(paste0(
      saveFile,
      i, "_ayn.pdf"
    ))
    plot(t_fit)

    blups <- setDT(as.data.frame(coef(t_fit)$random), keep.rownames = T)
    blups$rn <- str_remove(blups$rn, "Variety_")
    colnames(blups)[colnames(blups) == "rn"] <- "Variety"
    colnames(blups)[colnames(blups) == "effect"] <- paste(i)

    joinFile <- blups %>%
      inner_join(joinFile, by = "Variety")
    graphics.off()
  }
  return(joinFile)
}

blups18_ayn <- ayn_blups_2018(
  dat = pheno17_ayn,
  joinFile = blups,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2018/"
)

blups18_ayn <- blups18_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  inner_join(keyData_18, by = "ID") %>%
  select(Variety, trait_ID, phenotype_value) %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "rep"))

blups18_pyn <- blups18_pyn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2018-06-15")
    )
  )

blups18_ayn <- blups18_ayn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2018-06-15")
    )
  )

write.table(blups18_pyn, "./Results/Blups_pyn_18.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(blups18_ayn, "./Results/Blups_ayn_18.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
###############################################################################
