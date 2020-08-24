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
allVarieties <- tibble(Variety = unique(pheno$Variety))

pheno17 <- pheno %>%
  filter(year == "17")
varieties17_db <- tibble(varietydb_17 = unique(pheno17$Variety), DB = 1)
entityID17_db <- tibble(entityIDdb_17 = unique(pheno17$entity_id), DB = 1)

pheno18 <- pheno %>%
  filter(year == "18")
varieties18_db <- tibble(
  varietydb_18 = unique(pheno18$Variety),
  DB = 1
)
entityID18_db <- tibble(
  entityIDdb_18 = unique(pheno18$entity_id),
  DB = 1
)

pheno19 <- pheno %>%
  filter(year == "19") # %>%
# mutate(phenotype_date = str_replace(phenotype_date,
#                                     "2019-06-30","2019-07-01"))

varieties19_db <- tibble(varietydb_19 = unique(pheno19$Variety), DB = 1)
entityID19_db <- tibble(entityIDdb_19 = unique(pheno19$entity_id), DB = 1)

# VI

## 2016-2017

fileNames <- list.files(
  path = "./PhenoDatabase/2016_2017",
  full.names = T
)

traitNames <- basename(fileNames) %>%
  str_remove_all(c("traits_2016-2017_|_no_fills.xlsx|.xlsx"))

load.file <- function(filename) {
  d <- filename %>%
    readxl::excel_sheets() %>%
    purrr::set_names() %>%
    map(readxl::read_excel, path = filename)
  d
}

phenoVI_17 <- lapply(fileNames, load.file)

names(phenoVI_17) <- traitNames
phenoVI_17 <- plyr::ldply(phenoVI_17, data.frame, .id = "location")
print(colnames(phenoVI_17))

write.table(phenoVI_17, "./PhenoDatabase/raw_phenoVI_17.txt",
  quote = F,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
phenoVI_17 <- fread("./PhenoDatabase/raw_phenoVI_17.txt")

phenoVI_17 <- phenoVI_17 %>%
  rename(entity_id = NDVI.Plot_ID) %>%
  select(-ends_with("Plot_ID")) %>%
  pivot_longer(
    col = `NDVI.42880`:`Height.42894`,
    names_to = "trait_id",
    values_to = "phenotype_value"
  )

phenoVI_17 <- phenoVI_17 %>%
  separate(trait_id, c("trait_id", "phenotype_date")) %>%
  mutate(
    phenotype_date = as.numeric(phenotype_date),
    phenotype_date = as.Date(phenotype_date, origin = "1899-12-30"),
    phenotype_date = as.Date(phenotype_date, format = "%Y%m%d")
  ) %>%
  drop_na(phenotype_value) %>%
  filter(trait_id != "Height") %>%
  filter(trait_id != "Green") %>%
  filter(entity_id != "Fill") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RedEdge") %>%
  filter(str_detect(entity_id, ".*YN")) %>%
  mutate(phenotype_date = str_remove_all(phenotype_date, "-")) %>%
  unite("trait_id", trait_id, phenotype_date, sep = "_") %>%
  pivot_wider(
    names_from = trait_id,
    values_from = phenotype_value
  ) %>%
  select(-location)

str(phenoVI_17)

write.table(phenoVI_17, "./PhenoDatabase/Pheno_vi17.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = T
)
phenoVI_17 <- fread("./PhenoDatabase/Pheno_vi17.txt")

entityID17_VI <- tibble(entityIDVI_17 = unique(phenoVI_17$entity_id), VI = 1)

## 2017-2018
fileNames <- list.files(
  path = "./PhenoDatabase/2017_2018",
  full.names = T,
  pattern = ".csv$"
)
traitNames <- basename(fileNames) %>%
  str_remove_all(".csv")

load.file <- function(filename) {
  d <- fread(file = filename, header = TRUE, check.names = F)
  d
}

phenoVI_18 <- lapply(fileNames, load.file)

names(phenoVI_18) <- traitNames

names(phenoVI_18) <- traitNames

phenoVI_18 <- plyr::ldply(phenoVI_18, data.frame, .id = "trait_id")

phenoVI_18 <- phenoVI_18 %>%
  select(trait_id, Plot_ID, GNDVI, NDRE, NDVI) %>%
  pivot_longer(
    cols = GNDVI:NDVI,
    names_to = "VI",
    values_to = "phenotype_value"
  ) %>%
  drop_na(phenotype_value) %>%
  mutate(separate_1 = Plot_ID) %>%
  separate(separate_1, c("entity_id", "value_1", "value_2"),
    sep = ":"
  ) %>%
  separate(trait_id, c("phenotype_date", "location", "trait_id"),
    sep = "_"
  ) %>%
  mutate(
    location = as.character(location),
    location = str_replace(location, "HERN", "RN"),
    location = str_replace(location, "MCP", "MP"),
    value_1 = as.numeric(value_1),
    value_2 = as.numeric(value_2),
    column_mc = if_else(location == "MP", value_1, 0),
    column_sa = if_else(location == "SA", 48 - value_1, 0),
    column_hern = if_else(location == "RN", value_2, 0),
    column_rp = if_else(location == "RP", 58 - value_2, 0),
    row_mc = if_else(location == "MP", value_2, 0),
    row_sa = if_else(location == "SA", 26 - value_2, 0),
    row_hern = if_else(location == "RN", 28 - value_1, 0),
    row_rp = if_else(location == "RP", value_1, 0)
  ) %>%
  pivot_longer(
    cols = column_mc:column_rp,
    names_to = "locationColumn",
    values_to = "column"
  ) %>%
  pivot_longer(
    cols = row_mc:row_rp,
    names_to = "locationRow",
    values_to = "range"
  ) %>%
  filter(
    column > 0,
    range > 0
  ) %>%
  select(-value_1, -value_2, -locationRow, -VI, -locationColumn) %>%
  select(trait_id, phenotype_date, everything()) %>%
  unite("trait_id", trait_id:phenotype_date, sep = "_") %>%
  pivot_wider(
    names_from = trait_id,
    values_from = phenotype_value
  )

entityID18_VI <- tibble(entityIDVI_18 = unique(phenoVI_18$entity_id), VI = 1)

## 2018-2019

fileNames <- list.files(
  path = "./PhenoDatabase/2018_2019",
  full.names = T,
  pattern = ".csv$"
)
traitNames <- basename(fileNames) %>%
  str_remove_all(".csv")

load.file <- function(filename) {
  d <- fread(file = filename, header = TRUE, check.names = F, data.table = F)
  d
}

phenoVI_19 <- lapply(fileNames, load.file)

names(phenoVI_19) <- traitNames

phenoVI_19 <- lapply(traitNames, function(x) {
  names(phenoVI_19[[x]])[2] <- x
  phenoVI_19[[x]]
})

names(phenoVI_19) <- traitNames

phenoVI_19 <- plyr::ldply(phenoVI_19, data.frame, .id = "trait_id")

phenoVI_19 <- phenoVI_19 %>%
  filter(str_detect(Plot_ID, ".*AYN.*") | str_detect(Plot_ID, ".*PYN.*")) %>%
  select(-trait_id) %>%
  pivot_longer(
    cols = -Plot_ID,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  drop_na(phenotype_value) %>%
  separate(trait_id, c("phenotype_date", "location", "trial", "trait_id"),
    sep = "_"
  ) %>%
  mutate(
    location = str_replace(location, "ASH", "RL"),
    phenotype_date = str_remove(phenotype_date, "X")
  ) %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  select(-trial) %>%
  unite("trait_id", trait_id, location, phenotype_date,
    sep = "_"
  )

write.table(phenoVI_19, "./PhenoDatabase/raw_phenoVI_19.txt",
  quote = F,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
phenoVI_19 <- fread("./PhenoDatabase/raw_phenoVI_19.txt")

entityID19_VI <- tibble(entityIDVI_19 = unique(phenoVI_19$Plot_ID), VI = 1)

## Name comparison between db and VI information

entityID17 <- entityID17_db %>%
  left_join(entityID17_VI, by = c("entityIDdb_17" = "entityIDVI_17")) %>%
  mutate(
    Total = DB + VI,
    Total = replace_na(Total, 99)
  )

entityID18 <- entityID18_db %>%
  left_join(entityID18_VI, by = c("entityIDdb_18" = "entityIDVI_18")) %>%
  mutate(
    Total = DB + VI,
    Total = replace_na(Total, 99)
  )

entityID19 <- entityID19_db %>%
  left_join(entityID19_VI, by = c("entityIDdb_19" = "entityIDVI_19")) %>%
  mutate(
    Total = DB + VI,
    Total = replace_na(Total, 99)
  )

# Variety names seem to match. Problem may be coming from different entity_id
# or possibly wrong field information.

missing17_pheno <- pheno17 %>%
  anti_join(phenoVI_17) %>%
  select(-phenotype_person, -phenotype_date) %>%
  distinct()

missing17_vi <- phenoVI_17 %>%
  anti_join(pheno17)

colnames(pheno17)
colnames(phenoVI_17)

missing18_pheno <- pheno18 %>%
  anti_join(phenoVI_18) %>%
  select(-trait_id, -phenotype_value, -phenotype_person, -phenotype_date) %>%
  distinct()

write.table(missing18_pheno,
  file = "./missingVI_pheno18.txt", quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(missing17_pheno,
  file = "./missingVI_pheno17.txt", quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

pheno17 <- pheno17 %>%
  filter(trait_id != "PTHT") %>%
  mutate(phenotype_date = if_else(phenotype_date == "2017-06-21",
    "2017-06-22",
    phenotype_date,
    missing = NULL
  )) %>%
  unite("trait_id", trait_id, phenotype_date, sep = "_") %>%
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
    cols = `NDVI_20170525`:`GRYLD_2017-06-22`, names_to = "trait_ID",
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
    trait_ID = str_remove_all(trait_ID, "-")
  ) %>%
  drop_na(phenotype_value)

pheno18 <- pheno18 %>%
  filter(trait_id != "PTHT") %>%
  unite("trait_id", trait_id, phenotype_date, sep = "_") %>%
  inner_join(phenoVI_18, by = c("location", "range", "column")) %>%
  pivot_wider(
    names_from = trait_id,
    values_from = phenotype_value
  ) %>%
  select(
    entity_id.x, Plot_ID, Variety, range, column, year, trial,
    location, treated, plot, Trial, block, rep, everything()
  ) %>%
  select(
    -entity_id.y
  ) %>%
  pivot_longer(
    cols = GNDVI_20171201:`GRYLD_2018-06-29`, names_to = "trait_ID",
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
    trait_ID = str_remove_all(trait_ID, "-")
  ) %>%
  rename(
    entity_id = entity_id.x
  ) %>%
  drop_na(phenotype_value)

keyFile18 <- pheno18 %>%
  select(entity_id, Plot_ID, Variety, range, column, trial, location, rep) %>%
  distinct()

write.table(keyFile18, "./PhenoDatabase/2017_2018/KeyFile.txt",
  quote = FALSE, sep = "\t", row.names = FALSE,
  col.names = TRUE
)

pheno18 <- pheno18 %>%
  select(-Plot_ID)

pheno19 <- pheno19 %>%
  filter(trait_id == "GRYLD")

duplicated(pheno19$entity_id)

pheno19 <- pheno19[!duplicated(pheno19$entity_id), ]

pheno19 <- pheno19 %>%
  unite("trait", trait_id, location, phenotype_date,
    sep = "_"
  ) %>%
  mutate(trait = str_remove_all(trait, "-")) %>%
  pivot_wider(
    names_from = trait,
    values_from = phenotype_value
  ) %>%
  left_join(phenoVI_19, by = c("entity_id" = "Plot_ID")) %>%
  pivot_wider(
    names_from = trait_id,
    values_from = phenotype_value
  ) %>%
  pivot_longer(
    cols = GRYLD_RP_20190717:NDVI_RL_20190617,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  drop_na(phenotype_value) %>%
  separate(trait_id, c("trait_id", "location", "phenotype_date"),
    sep = "_"
  ) %>%
  mutate(
    phenotype_date = str_remove_all(phenotype_date, "-"),
    phenotype_date = replace_na(phenotype_date, "20190720")
  ) %>%
  unite("trait_id", trait_id, phenotype_date, sep = "_") %>%
  drop_na(phenotype_value)

pheno19 <- pheno19 %>%
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

traitSummary19 <- pheno19 %>%
  group_by(location, Trial, treated, trait_id) %>%
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

varietySummary19 <- pheno19 %>%
  group_by(Variety, Trial, treated, location) %>%
  mutate(Count = row_number()) %>%
  summarise(Total = max(Count))

pheno17 <- pheno17 %>%
  separate(col = trait_ID, into = c("trait_ID", "Date"), sep = "_") %>%
  mutate(
    Date = as.Date(Date, format = "%Y%m%d"),
    Date = replace_na(Date, "2017-06-20"),
    Date = as.Date(Date, format = "%Y-%m-%d")
  ) %>%
  filter(trait_ID != "Height") %>%
  drop_na(phenotype_value)

dateLocationSummary_17 <- distinct(pheno17, Date, location)

pheno17_ayn <- pheno17 %>%
  filter(Trial == "AYN")
varieties17_ayn <- tibble(Variety = unique(pheno17_ayn$Variety))

pheno17_pyn <- pheno17 %>%
  filter(Trial == "PYN")
varieties17_pyn <- tibble(Variety = unique(pheno17_pyn$Variety))

pheno18 <- pheno18 %>%
  separate(col = trait_ID, into = c("trait_ID", "Date"), sep = "_") %>%
  mutate(
    Date = replace_na(Date, "20180611"),
    Date = as.Date(Date, format = "%Y%m%d")
  ) %>%
  filter(trait_ID != "CanopyHeight") %>%
  drop_na(phenotype_value)
dateLocationSummary_18 <- distinct(pheno18, Date, location)

pheno18_ayn <- pheno18 %>%
  filter(Trial == "AYN")
varieties18_ayn <- tibble(Variety = unique(pheno18_ayn$Variety))

pheno18_pyn <- pheno18 %>%
  filter(Trial == "PYN")
varieties18_pyn <- tibble(Variety = unique(pheno18_pyn$Variety))

pheno19 <- pheno19 %>%
  separate(
    col = trait_id, into = c("trait_id", "phenotype_date"),
    sep = "_"
  ) %>%
  mutate(
    phenotype_date = str_replace(phenotype_date, "NA", "20190701"),
    phenotype_date = as.Date(phenotype_date, format = "%Y%m%d")
  ) %>%
  drop_na(phenotype_value)

dateLocationSummary_19 <- distinct(pheno19, phenotype_date, location)

pheno19_ayn <- pheno19 %>%
  filter(Trial == "AYN")
varieties19_ayn <- tibble(Variety = unique(pheno19_ayn$Variety))

pheno19_pyn <- pheno19 %>%
  filter(Trial == "PYN")
varieties19_pyn <- tibble(Variety = unique(pheno19_pyn$Variety))

pheno17_pyn %>%
  filter(trait_ID != "GRYLD") %>%
  # filter(trait_ID != "CanopyHeight") %>%
  ggplot(aes(x = Date, y = phenotype_value, colour = location)) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 1, varwidth = TRUE,
    alpha = 0.5, size = 0.5
  ) +
  facet_wrap(~trait_ID,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Distribution of Vegetation Indices over season",
    subtitle = "2016-2017 Season, Preliminary Yield Nursery"
  )
ggsave("./Figures/distribution_VI_17pyn.png",
  height = 30,
  width = 30, units = "cm", dpi = 350
)

pheno17_ayn %>%
  filter(trait_ID != "GRYLD") %>%
  # filter(trait_ID != "CanopyHeight") %>%
  ggplot(aes(x = Date, y = phenotype_value, colour = location)) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 1, varwidth = TRUE,
    alpha = 0.5, size = 0.5
  ) +
  facet_wrap(~trait_ID,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Distribution of Vegetation Indices over season",
    subtitle = "2016-2017 Season, Advanced Yield Nursery"
  )
ggsave("./Figures/distribution_VI_17ayn.png",
  height = 30,
  width = 30, units = "cm", dpi = 350
)

pheno18_pyn %>%
  filter(trait_ID != "GRYLD") %>%
  # filter(trait_ID != "CanopyHeight") %>%
  ggplot(aes(x = Date, y = phenotype_value, colour = location)) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 1, varwidth = TRUE,
    alpha = 0.5, size = 0.5
  ) +
  facet_wrap(~trait_ID,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Distribution of Vegetation Indices over season",
    subtitle = "2017-2018 Season, Preliminary Yield Nursery"
  )
ggsave("./Figures/distribution_VI_18pyn.png",
  height = 30,
  width = 30, units = "cm", dpi = 350
)

pheno18_ayn %>%
  filter(trait_ID != "GRYLD") %>%
  # filter(trait_ID != "CanopyHeight") %>%
  ggplot(aes(x = Date, y = phenotype_value, colour = location)) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 1, varwidth = TRUE,
    alpha = 0.5, size = 0.5
  ) +
  facet_wrap(~trait_ID,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Distribution of Vegetation Indices over season",
    subtitle = "2017-2018 Season, Advanced Yield Nursery"
  )
ggsave("./Figures/distribution_VI_18ayn.png",
  height = 30,
  width = 30, units = "cm", dpi = 350
)

pheno19_pyn %>%
  filter(trait_id != "GRYLD") %>%
  # filter(trait_ID != "CanopyHeight") %>%
  ggplot(aes(x = phenotype_date, y = phenotype_value, colour = location)) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 1, varwidth = TRUE,
    alpha = 0.5, size = 0.5
  ) +
  facet_wrap(~trait_id,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Distribution of Vegetation Indices over season",
    subtitle = "2018-2019 Season, Preliminary Yield Nursery"
  )
ggsave("./Figures/distribution_VI_19pyn.png",
  height = 30,
  width = 30, units = "cm", dpi = 350
)

pheno19_ayn %>%
  filter(trait_id != "GRYLD") %>%
  # filter(trait_ID != "CanopyHeight") %>%
  ggplot(aes(x = phenotype_date, y = phenotype_value, colour = location)) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 1, varwidth = TRUE,
    alpha = 0.5, size = 0.5
  ) +
  facet_wrap(~trait_id,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Distribution of Vegetation Indices over season",
    subtitle = "2018-2019 Season, Advanced Yield Nursery"
  )
ggsave("./Figures/distribution_VI_19ayn.png",
  height = 30,
  width = 30, units = "cm", dpi = 350
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
write.table(pheno19, "./PhenoDatabase/AllPheno_19.txt",
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
    ID = cur_group_id()
  ) %>%
  unite("ID", group, ID, sep = "_")

keyData_ayn_17 <- pheno17_ayn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_ayn_17

pheno17_ayn <- pheno17_ayn %>%
  ungroup() %>%
  select(-trait_ID) %>%
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
  res <- res %>%
    mutate(
      Trait = as.character(Trait),
      H = as.numeric(H),
      varVariety = as.numeric(varVariety),
      varRep = as.numeric(varRep),
      varRepBlock = as.numeric(varRepBlock),
      varRes = as.numeric(varRes)
    )
  ntraits <- 1:nrow(res)

  for (i in ntraits) {
    j <- res[i, 1]

    h <- as.data.frame(VarCorr(
      lmer(paste0(j, "~ (1|Variety) + (1|rep) + (1|rep:block)"), data = dat)
    ))

    res[i, 2] <- h[1, 4] / (h[1, 4] + (h[4, 4] / r))
    res[i, 3] <- round(h[1, 4], digits = 5)
    res[i, 4] <- round(h[3, 4], digits = 5)
    res[i, 5] <- round(h[2, 4], digits = 5)
    res[i, 6] <- round(h[4, 4], digits = 5)
  }
  return(res)
}

ayn_17 <- calcH2r(dat = pheno17_ayn)

ayn_17 <- ayn_17 %>%
  left_join(keyData_ayn_17, by = c("Trait" = "ID")) %>%
  select(trait_ID, -Trait, everything()) %>%
  separate(trait_ID, c("trait_ID", "Date", "location"), sep = "_") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

pheno18_ayn <- pheno18_ayn %>%
  unite("trait_ID", trait_ID, Date, location, sep = "_") %>%
  group_by(trait_ID) %>%
  mutate(
    group = "group",
    ID = cur_group_id()
  ) %>%
  unite("ID", group, ID, sep = "_") %>%
  ungroup()

keyData_ayn_18 <- pheno18_ayn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_ayn_18

pheno18_ayn <- pheno18_ayn %>%
  select(-trait_ID) %>%
  pivot_wider(names_from = ID, values_from = phenotype_value)

ayn_18 <- calcH2r(dat = pheno18_ayn)

ayn_18 <- ayn_18 %>%
  left_join(keyData_ayn_18, by = c("Trait" = "ID")) %>%
  select(trait_ID, -Trait, everything()) %>%
  separate(trait_ID, c("trait_ID", "Date", "location"), sep = "_") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

pheno19_ayn <- pheno19_ayn %>%
  filter(location != "RP") %>%
  unite("trait_ID", trait_id, phenotype_date, location, sep = "_") %>%
  group_by(trait_ID) %>%
  mutate(
    group = "group",
    ID = cur_group_id()
  ) %>%
  unite("ID", group, ID, sep = "_") %>%
  ungroup()

keyData_ayn_19 <- pheno19_ayn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_ayn_19

pheno19_ayn <- pheno19_ayn %>%
  select(-trait_ID) %>%
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
  res <- res %>%
    mutate(
      Trait = as.character(Trait),
      H = as.numeric(H),
      varVariety = as.numeric(varVariety),
      varRep = as.numeric(varRep),
      varRepBlock = as.numeric(varRepBlock),
      varRes = as.numeric(varRes)
    )
  ntraits <- 1:nrow(res)

  for (i in ntraits) {
    j <- res[i, 1]

    h <- as.data.frame(VarCorr(
      lmer(paste0(j, "~ (1|Variety) + (1|rep) + (1|rep:block)"), data = dat)
    ))

    res[i, 2] <- h[1, 4] / (h[1, 4] + (h[4, 4] / r))
    res[i, 3] <- round(h[1, 4], digits = 5)
    res[i, 4] <- round(h[3, 4], digits = 5)
    res[i, 5] <- round(h[2, 4], digits = 5)
    res[i, 6] <- round(h[4, 4], digits = 5)
  }
  return(res)
}

ayn_19 <- calcH2r(dat = pheno19_ayn)

ayn_19 <- ayn_19 %>%
  left_join(keyData_ayn_19, by = c("Trait" = "ID")) %>%
  select(trait_ID, -Trait, everything()) %>%
  separate(trait_ID, c("trait_ID", "Date", "location"), sep = "_") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

ayn_17 %>%
  filter(trait_ID != "GRYLD") %>%
  # filter(Date > as.Date("2018-01-02")) %>%
  ggplot(aes(x = Date, y = H, colour = as.factor(location))) +
  geom_point() +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  scale_x_date(date_labels = "%m/%d") +
  facet_wrap(~trait_ID, scales = "free") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Heritability for 2017 AYN ",
    y = expression("H"^{
      2
    })
  )

ggsave("./Figures/Heritability_VI_2017_ayn.png",
  height = 10,
  width = 35, units = "cm", dpi = 350
)

ayn_18 %>%
  filter(trait_ID != "GRYLD") %>%
  filter(Date > as.Date("2018-01-02")) %>%
  ggplot(aes(x = Date, y = H, colour = as.factor(location))) +
  geom_point() +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  scale_x_date(date_labels = "%m/%d") +
  facet_wrap(~trait_ID, scales = "free") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Heritability for 2018 AYN ",
    y = expression("H"^{
      2
    })
  )

ggsave("./Figures/Heritability_VI_2018_ayn.png",
  height = 10,
  width = 35, units = "cm", dpi = 350
)

ayn_19 %>%
  filter(trait_ID != "GRYLD") %>%
  # filter(Date > as.Date("2018-01-02")) %>%
  ggplot(aes(x = Date, y = H, colour = as.factor(location))) +
  geom_point() +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  scale_x_date(date_labels = "%m/%d") +
  facet_wrap(~trait_ID, scales = "free") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Heritability for 2019 AYN ",
    y = expression("H"^{
      2
    })
  )

ggsave("./Figures/Heritability_VI_2019_ayn.png",
  height = 10,
  width = 35, units = "cm", dpi = 350
)

write.table(ayn_17, "./Results/H2_ayn_2017.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(ayn_18, "./Results/H2_ayn_2018.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(ayn_19, "./Results/H2_ayn_2019.txt",
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
    ID = cur_group_id()
  ) %>%
  unite("ID", group, ID, sep = "_") %>%
  select(-phenotype_person, -range) %>%
  ungroup()

glimpse(pheno17_pyn)

keyData_17 <- pheno17_pyn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_17

pheno17_pyn <- pheno17_pyn %>%
  select(-trait_ID) %>%
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

    joinFile <- joinFile %>%
      tidylog::left_join(blups, by = "Variety")
    graphics.off()
  }
  return(joinFile)
}

blups17_pyn <- pyn_blups_2017(
  dat = pheno17_pyn,
  joinFile = allVarieties,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2017/"
)

blups17_pyn <- blups17_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  left_join(keyData_17, by = "ID") %>%
  select(Variety, trait_ID, phenotype_value) %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "block")) %>%
  drop_na(phenotype_value)

varieties17_pynBlups <- tibble(Variety = unique(blups17_pyn$Variety))

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

    joinFile <- joinFile %>%
      tidylog::left_join(blups, by = "Variety")
    graphics.off()
  }
  return(joinFile)
}

blups17_ayn <- ayn_blups_2017(
  dat = pheno17_ayn,
  joinFile = allVarieties,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2017/"
)

varieties17_aynBlups <- tibble(Variety = unique(blups17_ayn$Variety))

blups17_ayn <- blups17_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  left_join(keyData_ayn_17, by = "ID") %>%
  select(Variety, trait_ID, phenotype_value) %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "rep")) %>%
  drop_na(phenotype_value)

varieties17_aynBlups <- tibble(Variety = unique(blups17_ayn$Variety))

blups17_pyn <- blups17_pyn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2017-06-10")
    )
  ) %>%
  drop_na(phenotype_value)

blups17_ayn <- blups17_ayn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2017-06-10")
    )
  ) %>%
  drop_na(phenotype_value)

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
    ID = cur_group_id()
  ) %>%
  unite("ID", group, ID, sep = "_") %>%
  select(-phenotype_person, -range) %>%
  ungroup()

glimpse(pheno18_pyn)

keyData_18 <- pheno18_pyn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_18

pheno18_pyn <- pheno18_pyn %>%
  select(-trait_ID) %>%
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

    joinFile <- joinFile %>%
      tidylog::left_join(blups, by = "Variety")
    graphics.off()
  }
  return(joinFile)
}

blups18_pyn <- pyn_blups_2018(
  dat = pheno18_pyn,
  joinFile = allVarieties,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2018/"
)

blups18_pyn <- blups18_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  left_join(keyData_18, by = "ID") %>%
  select(Variety, trait_ID, phenotype_value) %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "block")) %>%
  drop_na(phenotype_value)

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

    joinFile <- joinFile %>%
      tidylog::left_join(blups)
    graphics.off()
  }
  return(joinFile)
}

blups18_ayn <- ayn_blups_2018(
  dat = pheno18_ayn,
  joinFile = allVarieties,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2018/"
)

blups18_ayn <- blups18_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  left_join(keyData_ayn_18, by = "ID") %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "rep")) %>%
  drop_na(phenotype_value)

blups18_pyn <- blups18_pyn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2018-06-15")
    )
  ) %>%
  drop_na(phenotype_value)

blups18_ayn <- blups18_ayn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2018-06-15")
    )
  ) %>%
  drop_na(phenotype_value)

write.table(blups18_pyn, "./Results/Blups_pyn_18.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(blups18_ayn, "./Results/Blups_ayn_18.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

## 2019
# Pyn
pheno19_pyn <- pheno19_pyn %>%
  unite("trait_ID", trait_id, phenotype_date, location, sep = "_") %>%
  group_by(trait_ID) %>%
  mutate(
    group = "group",
    ID = dplyr::group_indices()
  ) %>%
  unite("ID", group, ID, sep = "_") %>%
  select(-range) %>%
  ungroup()

glimpse(pheno19_pyn)

keyData_19 <- pheno19_pyn %>%
  select(trait_ID, ID) %>%
  distinct()
keyData_19

pheno19_pyn <- pheno19_pyn %>%
  select(-trait_ID) %>%
  pivot_wider(names_from = ID, values_from = phenotype_value) %>%
  mutate(Variety = as.factor(Variety)) %>%
  ungroup()

fit <- asreml(
  fixed = group_1 ~ 1,
  random = ~ Variety + block,
  data = pheno19_pyn
)
plot(fit)

blups <- setDT(as.data.frame(coef(fit)$random), keep.rownames = T)
blups$rn <- str_remove(blups$rn, "Variety_")
colnames(blups)[colnames(blups) == "rn"] <- "Variety"

# As a function
pyn_blups_2019 <- function(dat, joinFile, saveFile, ...) {
  effectvars <- names(dat) %in% c(
    "entity_id", "Variety", "column", "year", "trial",
    "trait_ID", "treated", "plot", "Trial", "block",
    "rep", "phenotype_date", "entry", "new", "entryc",
    "phenotype_person"
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

    joinFile <- joinFile %>%
      tidylog::left_join(blups, by = "Variety")
    graphics.off()
  }
  return(joinFile)
}

blups19_pyn <- pyn_blups_2019(
  dat = pheno19_pyn,
  joinFile = allVarieties,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2019/"
)

blups19_pyn <- blups19_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  left_join(keyData_19, by = "ID") %>%
  select(Variety, trait_ID, phenotype_value) %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "block")) %>%
  drop_na(phenotype_value)

varieties19_pynBlups <- tibble(Variety = unique(blups19_pyn$Variety))

# Ayn

glimpse(pheno19_ayn)
pheno19_ayn <- pheno19_ayn %>%
  mutate(Variety = as.factor(Variety))

fit <- asreml(
  fixed = group_1 ~ 1,
  random = ~ Variety + rep + rep:block,
  data = pheno19_ayn
)
plot(fit)

blups <- setDT(as.data.frame(coef(fit)$random), keep.rownames = T)
blups$rn <- str_remove(blups$rn, "Variety_")
colnames(blups)[colnames(blups) == "rn"] <- "Variety"

# As a function
ayn_blups_2019 <- function(dat, joinFile, saveFile, ...) {
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

    joinFile <- joinFile %>%
      tidylog::left_join(blups, by = "Variety")
    graphics.off()
  }
  return(joinFile)
}

blups19_ayn <- ayn_blups_2019(
  dat = pheno19_ayn,
  joinFile = allVarieties,
  saveFile = "./Figures/AsremlPlots/BLUPs/season2019/"
)

varieties19_aynBlups <- tibble(Variety = unique(blups19_ayn$Variety))

blups19_ayn <- blups19_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "ID",
    values_to = "phenotype_value"
  ) %>%
  left_join(keyData_ayn_19, by = "ID") %>%
  select(Variety, trait_ID, phenotype_value) %>%
  separate(trait_ID,
    c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>%
  filter(!str_detect(Variety, "rep")) %>%
  drop_na(phenotype_value)

varieties19_aynBlups <- tibble(Variety = unique(blups19_ayn$Variety))

blups19_pyn <- blups19_pyn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2019-07-20")
    )
  ) %>%
  drop_na(phenotype_value)

blups19_ayn <- blups19_ayn %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = replace_na(
      phenotype_date, as.Date("2019-07-20")
    )
  ) %>%
  drop_na(phenotype_value)

write.table(blups19_pyn, "./Results/Blups_pyn_19.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(blups19_ayn, "./Results/Blups_ayn_19.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

###############################################################################
blups17_ayn <- fread("./Results/Blups_ayn_17.txt")
blups17_pyn <- fread("./Results/Blups_pyn_17.txt")

blups18_ayn <- fread("./Results/Blups_ayn_18.txt")
blups18_pyn <- fread("./Results/Blups_pyn_18.txt")

blups19_ayn <- fread("./Results/Blups_ayn_19.txt")
blups19_pyn <- fread("./Results/Blups_pyn_19.txt")

blups19_pyn %>%
  filter(trait_id != "GRYLD") %>%
  ggplot(aes(
    x = as.Date(phenotype_date, format = "%Y-%m-%d"),
    y = phenotype_value,
    colour = location
  )) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 2.75, varwidth = TRUE,
    alpha = 0.5, size = 0.5
  ) +
  facet_wrap(~trait_id, scales = "free") +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "21 days"
  ) +
  labs(
    title = "Distribution of Vegetation Indices BLUPs over season",
    subtitle = "2018-2019 Season, Preliminary Yield Nursery",
    x = "Date"
  )
###############################################################################
#### Correlations ####

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

gryld17 <- fread("./PhenoDatabase/BlupsGRYLD_17_allVarieties.txt")
gryld18 <- fread("./PhenoDatabase/BlupsGRYLD_18_allVarieties.txt")
gryld19 <- fread("./PhenoDatabase/BlupsGRYLD_19_allVarieties.txt")

blups17_ayn <- blups17_ayn %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value)

blups17_pyn <- blups17_pyn %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value)

blups18_ayn <- blups18_ayn %>%
  select(-ID) %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value)

blups18_pyn <- blups18_pyn %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value)

blups19_ayn <- blups19_ayn %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value)

blups19_pyn <- blups19_pyn %>%
  unite("trait_id", trait_id, phenotype_date, location, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = phenotype_value)

cor_17ayn_pheno <- Hmisc::rcorr(as.matrix(pheno17_ayn[, 14:ncol(pheno17_ayn)]))
cor_17ayn_pheno <- flattenCorrMatrix(
  cormat = cor_17ayn_pheno$r,
  pmat = cor_17ayn_pheno$P
)

cor_17ayn <- Hmisc::rcorr(as.matrix(blups17_ayn[, 2:ncol(blups17_ayn)]))
cor_17ayn <- flattenCorrMatrix(
  cormat = cor_17ayn$r,
  pmat = cor_17ayn$P
)

cor_17pyn_pheno <- Hmisc::rcorr(as.matrix(pheno17_pyn[, 12:ncol(pheno17_pyn)]))
cor_17pyn_pheno <- flattenCorrMatrix(
  cormat = cor_17pyn_pheno$r,
  pmat = cor_17pyn_pheno$P
)

cor_17pyn <- Hmisc::rcorr(as.matrix(blups17_pyn[, 2:ncol(blups17_pyn)]))
cor_17pyn <- flattenCorrMatrix(
  cormat = cor_17pyn$r,
  pmat = cor_17pyn$P
)

cor_18ayn_pheno <- Hmisc::rcorr(as.matrix(pheno18_ayn[, 14:ncol(pheno18_ayn)]))
cor_18ayn_pheno <- flattenCorrMatrix(
  cormat = cor_18ayn_pheno$r,
  pmat = cor_18ayn_pheno$P
)

cor_18ayn <- Hmisc::rcorr(as.matrix(blups18_ayn[, 3:ncol(blups18_ayn)]))
cor_18ayn <- flattenCorrMatrix(
  cormat = cor_18ayn$r,
  pmat = cor_18ayn$P
)

cor_18pyn_pheno <- Hmisc::rcorr(as.matrix(pheno18_pyn[, 12:ncol(pheno18_pyn)]))
cor_18pyn_pheno <- flattenCorrMatrix(
  cormat = cor_18pyn_pheno$r,
  pmat = cor_18pyn_pheno$P
)

cor_18pyn <- Hmisc::rcorr(as.matrix(blups18_pyn[, 2:ncol(blups18_pyn)]))
cor_18pyn <- flattenCorrMatrix(
  cormat = cor_18pyn$r,
  pmat = cor_18pyn$P
)

cor_19ayn_pheno <- Hmisc::rcorr(as.matrix(pheno19_ayn[, 14:ncol(pheno19_ayn)]))
cor_19ayn_pheno <- flattenCorrMatrix(
  cormat = cor_19ayn_pheno$r,
  pmat = cor_19ayn_pheno$P
)

cor_19ayn <- Hmisc::rcorr(as.matrix(blups19_ayn[, 2:ncol(blups19_ayn)]))
cor_19ayn <- flattenCorrMatrix(
  cormat = cor_19ayn$r,
  pmat = cor_19ayn$P
)

cor_19pyn_pheno <- Hmisc::rcorr(as.matrix(pheno19_pyn[, 13:ncol(pheno19_pyn)]))
cor_19pyn_pheno <- flattenCorrMatrix(
  cormat = cor_19pyn_pheno$r,
  pmat = cor_19pyn_pheno$P
)

cor_19pyn <- Hmisc::rcorr(as.matrix(blups19_pyn[, 2:ncol(blups19_pyn)]))
cor_19pyn <- flattenCorrMatrix(
  cormat = cor_19pyn$r,
  pmat = cor_19pyn$P
)

cor_17ayn <- cor_17ayn %>%
  filter(str_detect(row, "GRYLD") | str_detect(column, "GRYLD")) %>%
  separate(row,
    c("row_trait_id", "row_phenotype_date", "row_location"),
    sep = "_"
  ) %>%
  separate(column,
    c("column_trait_id", "column_phenotype_date", "column_location"),
    sep = "_"
  ) %>%
  mutate(
    column_phenotype_date = replace_na(
      column_phenotype_date,
      "2017-06-15"
    ),
    column_location = if_else(is.na(column_location),
      "All",
      column_location
    ),
    facet_trait_id = if_else(
      row_trait_id == "GRYLD", column_trait_id, row_trait_id, missing = NULL
    ),
    facet_date = if_else(
      row_trait_id == "GRYLD", column_phenotype_date, row_phenotype_date,
      missing = NULL
    ),
    facet_date = as.Date(facet_date)
  ) %>%
  filter(row_location == column_location) %>%
  mutate(
    column_location = if_else(column_location == "All",
      "Season",
      "Field"
    ),
    row_phenotype_date = as.Date(row_phenotype_date, format = "%Y-%m-%d")
  )

write.table(cor_17ayn,
  file = "./Results/Correlations_ayn17.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

cor_17ayn %>%
  filter(row_trait_id != "GRYLD") %>%
  ggplot(aes(
    x = facet_date, y = cor,
    colour = row_location
  )) +
  geom_point() +
  facet_wrap(~facet_trait_id,
    scales = "free"
  ) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "2 weeks") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  labs(
    title = "Correlation of Vegetation Indices with GRYLD",
    subtitle = "2016-2017 Advanced Yield Nursery",
    x = "Date",
    y = "Correlation"
  )

ggsave("./Figures/correlation_17_ayn.png",
  height = 12,
  width = 35, units = "cm", dpi = 350
)

cor_17pyn <- cor_17pyn %>%
  filter(str_detect(row, "GRYLD") | str_detect(column, "GRYLD")) %>%
  separate(row,
    c("row_trait_id", "row_phenotype_date", "row_location"),
    sep = "_"
  ) %>%
  separate(column,
    c("column_trait_id", "column_phenotype_date", "column_location"),
    sep = "_"
  ) %>%
  mutate(
    column_phenotype_date = replace_na(
      column_phenotype_date,
      "2017-06-15"
    ),
    column_location = if_else(is.na(column_location),
      "All",
      column_location
    )
  ) %>%
  filter(row_location == column_location) %>%
  mutate(
    column_location = if_else(column_location == "All",
      "Season",
      "Field"
    ),
    row_phenotype_date = as.Date(row_phenotype_date, format = "%Y-%m-%d"),
    column_phenotype_date = as.Date(column_phenotype_date, format = "%Y-%m-%d"),
    facet_trait_id = if_else(
      row_trait_id == "GRYLD", column_trait_id, row_trait_id, missing = NULL
    ),
    facet_date = if_else(
      row_trait_id == "GRYLD", column_phenotype_date, row_phenotype_date,
      missing = NULL
    ),
    facet_date = as.Date(facet_date)
  )

write.table(cor_17pyn,
  file = "./Results/Correlations_pyn17.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

cor_17pyn %>%
  filter(row_trait_id != "GRYLD") %>%
  ggplot(aes(
    x = facet_date, y = cor,
    colour = row_location
  )) +
  geom_point() +
  facet_wrap(~facet_trait_id,
    scales = "free"
  ) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "2 weeks") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  labs(
    title = "Correlation of Vegetation Indices with GRYLD",
    subtitle = "2016-2017 Preliminary Yield Nursery",
    x = "Date",
    y = "Correlation"
  )

ggsave("./Figures/correlation_17_pyn.png",
  height = 12,
  width = 35, units = "cm", dpi = 350
)

cor_18ayn <- cor_18ayn %>%
  filter(str_detect(row, "GRYLD") | str_detect(column, "GRYLD")) %>%
  separate(row,
    c("row_trait_id", "row_phenotype_date", "row_location"),
    sep = "_"
  ) %>%
  separate(column,
    c("column_trait_id", "column_phenotype_date", "column_location"),
    sep = "_"
  ) %>%
  mutate(
    column_phenotype_date = replace_na(
      column_phenotype_date,
      "2018-06-15"
    ),
    column_location = if_else(is.na(column_location),
      "All",
      column_location
    )
  ) %>%
  filter(row_location == column_location) %>%
  mutate(
    column_location = if_else(column_location == "All",
      "Season",
      "Field"
    ),
    row_phenotype_date = as.Date(row_phenotype_date, format = "%Y-%m-%d"),
    column_phenotype_date = as.Date(column_phenotype_date, format = "%Y-%m-%d"),
    facet_trait_id = if_else(
      row_trait_id == "GRYLD", column_trait_id, row_trait_id, missing = NULL
    ),
    facet_date = if_else(
      row_trait_id == "GRYLD", column_phenotype_date, row_phenotype_date,
      missing = NULL
    ),
    facet_date = as.Date(facet_date)
  )

write.table(cor_18ayn,
  file = "./Results/Correlations_ayn18.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

cor_18ayn %>%
  filter(row_trait_id != "GRYLD") %>%
  filter(row_phenotype_date > as.Date("2018-02-01")) %>%
  ggplot(aes(
    x = facet_date, y = cor,
    colour = row_location
  )) +
  geom_point() +
  facet_wrap(~facet_trait_id,
    scales = "free"
  ) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "3 weeks") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  labs(
    title = "Correlation of Vegetation Indices with GRYLD",
    subtitle = "2017-2018 Advanced Yield Nursery",
    x = "Date",
    y = "Correlation"
  )

ggsave("./Figures/correlation_18_ayn.png",
  height = 12,
  width = 35, units = "cm", dpi = 350
)

cor_18pyn <- cor_18pyn %>%
  filter(str_detect(row, "GRYLD") | str_detect(column, "GRYLD")) %>%
  separate(row,
    c("row_trait_id", "row_phenotype_date", "row_location"),
    sep = "_"
  ) %>%
  separate(column,
    c("column_trait_id", "column_phenotype_date", "column_location"),
    sep = "_"
  ) %>%
  mutate(
    column_phenotype_date = replace_na(
      column_phenotype_date,
      "2018-06-20"
    ),
    column_location = if_else(is.na(column_location),
      "All",
      column_location
    )
  ) %>%
  filter(row_location == column_location) %>%
  mutate(
    column_location = if_else(column_location == "All",
      "Season",
      "Field"
    ),
    row_phenotype_date = as.Date(row_phenotype_date, format = "%Y-%m-%d"),
    column_phenotype_date = as.Date(column_phenotype_date, format = "%Y-%m-%d"),
    facet_trait_id = if_else(
      row_trait_id == "GRYLD", column_trait_id, row_trait_id, missing = NULL
    ),
    facet_date = if_else(
      row_trait_id == "GRYLD", column_phenotype_date, row_phenotype_date,
      missing = NULL
    ),
    facet_date = as.Date(facet_date)
  )

write.table(cor_18pyn,
  file = "./Results/Correlations_pyn18.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

cor_18pyn %>%
  filter(row_trait_id != "GRYLD") %>%
  filter(row_phenotype_date > as.Date("2018-02-01")) %>%
  ggplot(aes(
    x = facet_date, y = cor,
    colour = row_location
  )) +
  geom_point() +
  facet_wrap(~facet_trait_id,
    scales = "free"
  ) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "3 weeks") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  labs(
    title = "Correlation of Vegetation Indices with GRYLD",
    subtitle = "2017-2018 Preliminary Yield Nursery",
    x = "Date",
    y = "Correlation"
  )

ggsave("./Figures/correlation_18_pyn.png",
  height = 12,
  width = 35, units = "cm", dpi = 350
)

cor_19ayn <- cor_19ayn %>%
  filter(str_detect(row, "GRYLD") | str_detect(column, "GRYLD")) %>%
  separate(row,
    c("row_trait_id", "row_phenotype_date", "row_location"),
    sep = "_"
  ) %>%
  separate(column,
    c("column_trait_id", "column_phenotype_date", "column_location"),
    sep = "_"
  ) %>%
  mutate(
    column_phenotype_date = replace_na(
      column_phenotype_date,
      "2019-06-15"
    ),
    column_location = if_else(is.na(column_location),
      "All",
      column_location
    ),
    facet_trait_id = if_else(
      row_trait_id == "GRYLD", column_trait_id, row_trait_id, missing = NULL
    ),
    facet_date = if_else(
      row_trait_id == "GRYLD", column_phenotype_date, row_phenotype_date,
      missing = NULL
    ),
    facet_date = as.Date(facet_date)
  ) %>%
  filter(
    facet_trait_id != "GRYLD",
    row_location == column_location
  ) %>%
  mutate(
    column_location = if_else(column_location == "All",
      "Season",
      "Field"
    ),
    row_phenotype_date = as.Date(row_phenotype_date, format = "%Y-%m-%d"),
    column_phenotype_date = as.Date(column_phenotype_date, format = "%Y-%m-%d")
  )

write.table(cor_19ayn,
  file = "./Results/Correlations_ayn19.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

cor_19ayn %>%
  # filter(row_trait_id != "GRYLD") %>%
  ggplot(aes(
    x = facet_date, y = cor,
    colour = row_location
  )) +
  geom_point() +
  facet_wrap(~facet_trait_id,
    scales = "free"
  ) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "3 weeks") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  labs(
    title = "Correlation of Vegetation Indices with GRYLD",
    subtitle = "2018-2019 Advanced Yield Nursery",
    x = "Date",
    y = "Correlation"
  )

ggsave("./Figures/correlation_19_ayn.png",
  height = 12,
  width = 35, units = "cm", dpi = 350
)

cor_19pyn <- cor_19pyn %>%
  filter(str_detect(row, "GRYLD") | str_detect(column, "GRYLD")) %>%
  separate(row,
    c("row_trait_id", "row_phenotype_date", "row_location"),
    sep = "_"
  ) %>%
  separate(column,
    c("column_trait_id", "column_phenotype_date", "column_location"),
    sep = "_"
  ) %>%
  mutate(
    column_phenotype_date = replace_na(
      column_phenotype_date,
      "2019-07-15"
    ),
    column_location = if_else(is.na(column_location),
      "All",
      column_location
    ),
    facet_trait_id = if_else(
      row_trait_id == "GRYLD", column_trait_id, row_trait_id, missing = NULL
    ),
    facet_date = if_else(
      row_trait_id == "GRYLD", column_phenotype_date, row_phenotype_date,
      missing = NULL
    ),
    facet_date = as.Date(facet_date)
  ) %>%
  filter(row_location == column_location) %>%
  mutate(
    column_location = if_else(column_location == "All",
      "Season",
      "Field"
    ),
    row_phenotype_date = as.Date(row_phenotype_date, format = "%Y-%m-%d")
  )


write.table(cor_19pyn,
  file = "./Results/Correlations_pyn19.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

cor_19pyn %>%
  filter(row_trait_id != "GRYLD") %>%
  ggplot(aes(
    x = facet_date, y = cor,
    colour = row_location
  )) +
  geom_point() +
  facet_wrap(~facet_trait_id,
    scales = "free"
  ) +
  scale_x_date(date_labels = "%m/%d", date_breaks = "3 weeks") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_colour_manual(
    values = c("#1b9e77", "#7570b3", "#e7298a", "#66a61e"),
    name = "Location"
  ) +
  labs(
    title = "Correlation of Vegetation Indices with GRYLD",
    subtitle = "2018-2019 Preliminary Yield Nursery",
    x = "Date",
    y = "Correlation"
  )

ggsave("./Figures/correlation_19_pyn.png",
  height = 12,
  width = 35, units = "cm", dpi = 350
)


#### writing Blups ####
write.table(blups17_pyn, "./Results/Blups_pyn_17.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(blups17_ayn, "./Results/Blups_ayn_17.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(blups18_pyn, "./Results/Blups_pyn_18.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(blups18_ayn, "./Results/Blups_ayn_18.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(blups19_pyn, "./Results/Blups_pyn_19.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(blups19_ayn, "./Results/Blups_ayn_19.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
