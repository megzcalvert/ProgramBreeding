rm(list = objects())
ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(RMySQL)
library(janitor)
library(broom)
library(reshape2)
library(lme4)

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
    legend.title = element_text(size = rel(2)),
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

getwd()
set.seed(1964)
options(scipen = 999)
# useful infos **reproducible research**
sessionInfo()

#### Connect to database ####
wheatgenetics <- dbConnect(MySQL(),
  user = rstudioapi::askForPassword("Database user"),
  dbname = "wheatgenetics", host = "beocat.cis.ksu.edu",
  password =
    rstudioapi::askForPassword("Database password"),
  port = 6306
)

# SQL Query to get all AM Panel phenotype data

pheno_query <- "SELECT phenotype.entity_id,
phenotype.trait_id,
phenotype.phenotype_value,
phenotype.phenotype_date,
phenotype.phenotype_person,
plot.plot_name AS 'Variety',
plot.location,
plot.range,
plot.column,
plot.rep,
plot.block,
plot.treatment
FROM wheatgenetics.phenotype LEFT JOIN 
wheatgenetics.plot ON plot.plot_id = phenotype.entity_id 
WHERE wheatgenetics.phenotype.entity_id LIKE '%YN%-MP-%' OR
wheatgenetics.phenotype.entity_id LIKE '%YN%-RN-%' OR 
wheatgenetics.phenotype.entity_id LIKE '%YN%-RP-%' OR
wheatgenetics.phenotype.entity_id LIKE '%YN%-SA-%'OR
wheatgenetics.phenotype.entity_id LIKE '%YN%-RL-%' ;"

# run the query to get plot information

pheno <- dbGetQuery(wheatgenetics, pheno_query)
str(pheno)

# save original data
saveRDS(pheno, "./PhenoDatabase/Pheno.RDS")

# disconnect from database
dbDisconnect(wheatgenetics)

rm(wheatgenetics, pheno_query, pheno)

# Read in original data
pheno_long <- readRDS("./PhenoDatabase/Pheno.RDS")
varietyNames <- tabyl(pheno_long$Variety)
names(varietyNames)[1:2] <- c("Variety", "Count")

write.table(varietyNames, "./PhenoDatabase/VarietyNames_Database.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

#### Need to fix some of the names in a text editor because
# they are str components: ^, ~, +

write.table(pheno_long, "./PhenoDatabase/PhenoLong.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

pheno_long <- fread("./PhenoDatabase/PhenoLong.txt", fill = TRUE)

pheno_long <- pheno_long %>%
  mutate(Variety = str_replace(Variety, "-M-", "M-")) %>%
  mutate(Variety = str_replace(Variety, "-K-", "K-")) %>%
  mutate(Variety = tolower(Variety)) %>%
  mutate(Variety = str_squish(Variety)) %>%
  mutate(Variety = str_replace(Variety, "-", "_")) %>%
  mutate(Variety = str_replace(Variety, "wb_4458", "wb4458")) %>%
  mutate(Variety = str_replace(Variety, "symonument", "monument")) %>%
  mutate(Variety = str_replace(Variety, "smith'sgold", "smithsgold"))

varietyNames_2 <- tabyl(pheno_long$Variety)
names(varietyNames_2)[1:2] <- c("Variety", "Count")

write.table(varietyNames_2, "./PhenoDatabase/VarietyNames_Corrected.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

write.table(pheno_long, "./PhenoDatabase/PhenoLongNamesCorrect.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
pheno_long <- fread("./PhenoDatabase/PhenoLongNamesCorrect.txt", fill = TRUE)

geno <- fread("./GenoDatabase/SelectedGeno_Numeric.txt")
geno <- as.data.frame(colnames(geno[, 4:ncol(geno)]))
names(geno)[1] <- "Variety"

variety_names_wgeno <- varietyNames_2 %>%
  semi_join(geno)
variety_names_wogeno <- varietyNames_2 %>%
  anti_join(geno)

pheno_long <- pheno_long %>%
  filter(
    Variety != "fill",
    Variety != "centralblend",
    trait_id == "GRYLD"
  ) %>%
  select(
    entity_id, trait_id, phenotype_value, phenotype_date, phenotype_person,
    Variety, range, column
  ) %>%
  mutate(Sep = entity_id) %>%
  separate(Sep, c("year", "trial", "location", "treated", "plot"), sep = "-")

trialSummary <- tabyl(pheno_long, trial, location, year)
trialSummary

# Function to add a column based on a portion of text in another column
ff <- function(x, patterns, replacements = patterns, fill = NA, ...) {
  stopifnot(length(patterns) == length(replacements))

  ans <- rep_len(as.character(fill), length(x))
  empty <- seq_along(x)

  for (i in seq_along(patterns)) {
    greps <- grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] <- replacements[[i]]
    empty <- empty[!greps]
  }

  return(ans)
}

# Adding a year column based on the entity_id
pheno_long$Trial <- ff(pheno_long$trial,
  c(
    "AYN3", "PYN", "GSPYN", "AYN2", "DHPYN",
    "AYN1", "AYN4", "PYN1", "PYN2", "DHPYN1",
    "DHPYN2", "DHAYN2", "DHAYN1", "PYNA",
    "PYNB", "PYN3", "PYN4", "GSAYN", "ATAYN1", "ATAYN2"
  ),
  c(
    "AYN", "PYN", "PYN", "AYN", "PYN",
    "AYN", "AYN", "PYN", "PYN", "PYN",
    "PYN", "AYN", "AYN1", "PYN",
    "PYN", "PYN", "PYN", "AYN", "AYN", "AYN"
  ),
  "NA",
  ignore.case = TRUE
)

tabyl(pheno_long, year, location)

pheno_long_line_number <- pheno_long %>%
  select(Variety, year, Trial, location) %>%
  distinct()

tabyl(pheno_long_line_number, Trial, location, year)

pheno_long <- pheno_long %>%
  mutate(block = range) %>%
  mutate(
    phenotype_value = as.numeric(phenotype_value),
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_person = as.factor(phenotype_person),
    Variety = as.factor(Variety),
    range = as.factor(range),
    column = as.factor(column),
    year = as.factor(year),
    location = as.factor(location),
    treated = as.factor(treated),
    plot = as.numeric(plot),
    Trial = as.factor(Trial),
    block = as.factor(block)
  ) %>%
  mutate(
    rep = if_else(Trial == "AYN" & plot >= 200, 2, 1),
    rep = as.factor(rep)
  ) %>%
  filter(case_when(
    trait_id == "PTHT" ~ phenotype_value < 200,
    T ~ phenotype_value == phenotype_value
  )) %>%
  filter(
    trait_id != "GRWT",
    trait_id != "PCTHEAD",
    trait_id != "MOIST",
    trait_id != "TESTWT"
  )

write.table(pheno_long, "./PhenoDatabase/Phenotype_Master.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

varietySummary <- pheno_long %>%
  group_by(year, location, trial, Variety, treated) %>%
  mutate(Count = row_number()) %>%
  summarise(Total = max(Count))

traitSummary <- pheno_long %>%
  group_by(year, location, trial, treated, trait_id) %>%
  mutate(Count = row_number()) %>%
  summarise(
    Total = max(Count),
    minimum = min(phenotype_value),
    average = mean(phenotype_value),
    maximum = max(phenotype_value),
    variance = var(phenotype_value)
  )

traitByTrialSummary <- pheno_long %>%
  group_by(year, location, Trial, treated, trait_id) %>%
  mutate(Count = row_number()) %>%
  summarise(
    Total = max(Count),
    minimum = min(phenotype_value),
    average = mean(phenotype_value),
    maximum = max(phenotype_value),
    variance = var(phenotype_value)
  )

write.table(traitSummary, "./Results/traitSummary.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
write.table(traitByTrialSummary, "./Results/traitByTrialSummary.txt",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

#### Dividing data into sets that make sense to predict later ####

pheno_long <- pheno_long %>%
  filter(trait_id == "GRYLD")
glimpse(pheno_long)

### By year
pheno_16 <- pheno_long %>%
  filter(year == 16)
pheno_17 <- pheno_long %>%
  filter(year == 17)
pheno_18 <- pheno_long %>%
  filter(year == 18)
pheno_19 <- pheno_long %>%
  filter(year == 19)

### Exluding whole years
pheno_excluding_16 <- pheno_long %>%
  filter(year != 16)
pheno_excluding_17 <- pheno_long %>%
  filter(year != 17)
pheno_excluding_18 <- pheno_long %>%
  filter(year != 18)
pheno_excluding_19 <- pheno_long %>%
  filter(year != 19)

### Year-location combinations
pheno_16_MP <- pheno_long %>%
  filter(year == 16) %>%
  filter(location == "MP")
pheno_16_RL <- pheno_long %>%
  filter(year == 16) %>%
  filter(location == "RL")
pheno_16_RN <- pheno_long %>%
  filter(year == 16) %>%
  filter(location == "RN")
pheno_16_SA <- pheno_long %>%
  filter(year == 16) %>%
  filter(location == "SA")

pheno_17_MP <- pheno_long %>%
  filter(year == 17) %>%
  filter(location == "MP")
pheno_17_RL <- pheno_long %>%
  filter(year == 17) %>%
  filter(location == "RL")
pheno_17_RN <- pheno_long %>%
  filter(year == 17) %>%
  filter(location == "RN")
pheno_17_RP <- pheno_long %>%
  filter(year == 17) %>%
  filter(location == "RP")

pheno_18_MP <- pheno_long %>%
  filter(year == 18) %>%
  filter(location == "MP")
pheno_18_RL <- pheno_long %>%
  filter(year == 18) %>%
  filter(location == "RL")
pheno_18_RN <- pheno_long %>%
  filter(year == 18) %>%
  filter(location == "RN")
pheno_18_RP <- pheno_long %>%
  filter(year == 18) %>%
  filter(location == "RP")
pheno_18_SA <- pheno_long %>%
  filter(year == 18) %>%
  filter(location == "SA")

pheno_19_RL <- pheno_long %>%
  filter(year == 19) %>%
  filter(location == "RL")
pheno_19_RN <- pheno_long %>%
  filter(year == 19) %>%
  filter(location == "RN")
pheno_19_RP <- pheno_long %>%
  filter(year == 19) %>%
  filter(location == "RP")

### All AYN for each year for heritability
pheno_16_ayn <- pheno_long %>%
  filter(year == 16) %>%
  filter(Trial == "AYN")

pheno_17_ayn <- pheno_long %>%
  filter(year == 17) %>%
  filter(Trial == "AYN")

pheno_18_ayn <- pheno_long %>%
  filter(year == 18) %>%
  filter(Trial == "AYN")

pheno_19_ayn <- pheno_long %>%
  filter(year == 19) %>%
  filter(Trial == "AYN")

### All PYN for each year
pheno_16_pyn <- pheno_long %>%
  filter(year == 16) %>%
  filter(Trial == "PYN")

pheno_17_pyn <- pheno_long %>%
  filter(year == 17) %>%
  filter(Trial == "PYN")

pheno_18_pyn <- pheno_long %>%
  filter(year == 18) %>%
  filter(Trial == "PYN")

pheno_19_pyn <- pheno_long %>%
  filter(year == 19) %>%
  filter(Trial == "PYN")

## Excluding each years PYN
pheno_excluding_16_pyn <- pheno_long %>%
  anti_join(pheno_16_pyn, by = "entity_id")
pheno_excluding_17_pyn <- pheno_long %>%
  anti_join(pheno_17_pyn, by = "entity_id")
pheno_excluding_18_pyn <- pheno_long %>%
  anti_join(pheno_18_pyn, by = "entity_id")
pheno_excluding_19_pyn <- pheno_long %>%
  anti_join(pheno_19_pyn, by = "entity_id")

## Including all data & from one location in predicted year
pheno_excluding_16_including_16MP <- pheno_long %>%
  filter(year != 16) %>%
  bind_rows(pheno_16_MP)
pheno_excluding_16_including_16RL <- pheno_long %>%
  filter(year != 16) %>%
  bind_rows(pheno_16_RL)
pheno_excluding_16_including_16RN <- pheno_long %>%
  filter(year != 16) %>%
  bind_rows(pheno_16_RN)
pheno_excluding_16_including_16SA <- pheno_long %>%
  filter(year != 16) %>%
  bind_rows(pheno_16_SA)

pheno_excluding_17_including_17MP <- pheno_long %>%
  filter(year != 17) %>%
  bind_rows(pheno_17_MP)
pheno_excluding_17_including_17RL <- pheno_long %>%
  filter(year != 17) %>%
  bind_rows(pheno_17_RL)
pheno_excluding_17_including_17RN <- pheno_long %>%
  filter(year != 17) %>%
  bind_rows(pheno_17_RN)
pheno_excluding_17_including_17RP <- pheno_long %>%
  filter(year != 17) %>%
  bind_rows(pheno_17_RP)

pheno_excluding_18_including_18MP <- pheno_long %>%
  filter(year != 18) %>%
  bind_rows(pheno_18_MP)
pheno_excluding_18_including_18RL <- pheno_long %>%
  filter(year != 18) %>%
  bind_rows(pheno_18_RL)
pheno_excluding_18_including_18RN <- pheno_long %>%
  filter(year != 18) %>%
  bind_rows(pheno_18_RN)
pheno_excluding_18_including_18RP <- pheno_long %>%
  filter(year != 18) %>%
  bind_rows(pheno_18_RP)
pheno_excluding_18_including_18SA <- pheno_long %>%
  filter(year != 18) %>%
  bind_rows(pheno_18_SA)

pheno_excluding_19_including_19RL <- pheno_long %>%
  filter(year != 19) %>%
  bind_rows(pheno_19_RL)
pheno_excluding_19_including_19RN <- pheno_long %>%
  filter(year != 19) %>%
  bind_rows(pheno_19_RN)
pheno_excluding_19_including_19RP <- pheno_long %>%
  filter(year != 19) %>%
  bind_rows(pheno_19_RP)

## Leave one location out by year
pheno_16_MP_RL_RN <- pheno_long %>%
  filter(year == 16) %>%
  filter(location != "SA")
pheno_16_MP_RL_SA <- pheno_long %>%
  filter(year == 16) %>%
  filter(location != "RN")
pheno_16_MP_RN_SA <- pheno_long %>%
  filter(year == 16) %>%
  filter(location != "RL")
pheno_16_RN_RL_SA <- pheno_long %>%
  filter(year == 16) %>%
  filter(location != "MP")

pheno_17_MP_RL_RN <- pheno_long %>%
  filter(year == 17) %>%
  filter(location != "RP")
pheno_17_MP_RL_RP <- pheno_long %>%
  filter(year == 17) %>%
  filter(location != "RN")
pheno_17_MP_RN_RP <- pheno_long %>%
  filter(year == 17) %>%
  filter(location != "RL")
pheno_17_RL_RN_RP <- pheno_long %>%
  filter(year == 17) %>%
  filter(location != "MP")

pheno_18_MP_RL_RN_RP <- pheno_long %>%
  filter(year == 18) %>%
  filter(location != "SA")
pheno_18_MP_RL_RN_SA <- pheno_long %>%
  filter(year == 18) %>%
  filter(location != "RP")
pheno_18_MP_RL_RP_SA <- pheno_long %>%
  filter(year == 18) %>%
  filter(location != "RN")
pheno_18_MP_RN_RP_SA <- pheno_long %>%
  filter(year == 18) %>%
  filter(location != "RL")
pheno_18_RL_RN_RP_SA <- pheno_long %>%
  filter(year == 18) %>%
  filter(location != "MP")

pheno_19_RL_RN <- pheno_long %>%
  filter(year == 19) %>%
  filter(location != "RP")
pheno_19_RL_RP <- pheno_long %>%
  filter(year == 19) %>%
  filter(location != "RN")
pheno_19_RN_RP <- pheno_long %>%
  filter(year == 19) %>%
  filter(location != "RL")

#### Testing regression models ####

fit_1 <- lmer(
  phenotype_value ~ (1 | Variety) + (1 | year) + (1 | location),
  data = pheno_long
)

fit_2 <- lmer(
  phenotype_value ~ (1 | Variety) + (1 | year) + (1 | location) + (1 | treated),
  data = pheno_long
)

fit_3 <- lmer(
  phenotype_value ~ (1 | Variety) + (1 | year) + (1 | location) + (1 | treated)
    + (1 | year:location),
  data = pheno_long
)

fit_4 <- lmer(
  phenotype_value ~ (1 | Variety) + (1 | year) + (1 | location) + (1 | treated)
    + (1 | year:location) + (1 | year:location:Variety),
  data = pheno_long
)

fit_5 <- lmer(
  phenotype_value ~ (1 | Variety) + (1 | year) + (1 | location) + (1 | treated)
    + (1 | year:location) + (1 | year:location:Variety)
    + (1 | year:location:rep),
  data = pheno_long
)

fit_6 <- lmer(
  phenotype_value ~ (1 | Variety) + (1 | year)
    + (1 | year:location) + (1 | year:location:Variety)
    + (1 | year:location:rep) + (1 | year:location:treated),
  data = pheno_long
)

anova(fit_1, fit_2, fit_3, fit_4, fit_5, fit_6)
tidy(anova(fit_1, fit_2, fit_3, fit_4, fit_5, fit_6))

summary(fit_1)
summary(fit_2)
summary(fit_3)
summary(fit_4)
summary(fit_5)
summary(fit_6)

plot(fit_1, type = c("p", "smooth"))
plot(fit_2, type = c("p", "smooth"))
plot(fit_3, type = c("p", "smooth"))
plot(fit_4, type = c("p", "smooth"))
plot(fit_5, type = c("p", "smooth"))
plot(fit_6, type = c("p", "smooth"))

plot(fit_1, sqrt(abs(resid(.))) ~ fitted(.),
  type = c("p", "smooth")
)
plot(fit_2, sqrt(abs(resid(.))) ~ fitted(.),
  type = c("p", "smooth")
)
plot(fit_3, sqrt(abs(resid(.))) ~ fitted(.),
  type = c("p", "smooth")
)
plot(fit_4, sqrt(abs(resid(.))) ~ fitted(.),
  type = c("p", "smooth")
)
plot(fit_5, sqrt(abs(resid(.))) ~ fitted(.),
  type = c("p", "smooth")
)
plot(fit_6, sqrt(abs(resid(.))) ~ fitted(.),
  type = c("p", "smooth")
)

pheno_long %>%
  ggplot(aes(x = phenotype_value, colour = location)) +
  geom_density(fill = NA) +
  facet_wrap(Trial ~ year, ncol = 4)

pheno_long %>%
  ggplot(aes(phenotype_value, colour = location, linetype = treated)) +
  geom_density(fill = NA) +
  facet_wrap(Trial ~ year, ncol = 4)

values <- setDT(as.data.frame(ranef(fit_6)$Variety), keep.rownames = TRUE) %>%
  rename(
    GRYLD = `(Intercept)`,
    Variety = rn
  )

values %>%
  ggplot(aes(x = GRYLD)) +
  geom_density(fill = NA)

write.table(values, "./Results/BlupsGRYLD_allYrs_allVarieties.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

varianceComponents <- as.data.frame(VarCorr(fit_6))
varianceComponents

year_count <- pheno_long %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_long %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_long %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

getVariance <- function(asr, comp) {
  variance <- as.data.frame(VarCorr(asr))
  variance <- column_to_rownames(variance, var = "grp")
  idx <- which(rownames(variance) == comp)
  v <- variance$vcov[idx]
  print(paste("variance component", v))
  return(v)
}

Vg <- getVariance(asr = fit_6, comp = "Variety")
Vgl <- getVariance(asr = fit_6, comp = "year:location:Variety")
Ve <- getVariance(asr = fit_6, comp = "Residual")

h_gryld <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld

####### Different years ######
#### GryldBlups Function ####

GryldBlups <- function(dat, identifier, saveFile, saveFigures,
                       equation) {
  fit <- lmer(equation,
    data = dat
  )

  print(summary(fit))
  values <- setDT(as.data.frame(ranef(fit)$Variety), keep.rownames = TRUE) %>%
    rename(
      GRYLD = `(Intercept)`,
      Variety = rn
    )

  values %>%
    ggplot(aes(x = GRYLD)) +
    geom_density(fill = NA) +
    labs(
      title = paste("Distribution of Blups"),
      subtitle = paste("pheno", identifier)
    )
  ggsave(
    filename = paste0("blups", identifier, ".png"),
    path = saveFigures,
    width = 20,
    height = 15,
    units = "cm"
  )
  write.table(values, paste0(saveFile, "BlupsGRYLD", identifier, ".txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )

  return(fit)
}
#### 2015-2016 ####
blups_excluding_16 <- GryldBlups(
  dat = pheno_excluding_16, identifier = "_excluding_16",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_16 %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_16 %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_16 %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_16, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_16, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_16, comp = "Residual")

h_gryld_excluding_16 <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_16

blups_excluding_16_pyn <- GryldBlups(
  dat = pheno_excluding_16_pyn, identifier = "_excluding_16_pyn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_16_pyn %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_16_pyn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_16_pyn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_16_pyn, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_16_pyn, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_16_pyn, comp = "Residual")

h_gryld_excluding_16_pyn <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_16_pyn

blups_excluding_16_including_16MP <- GryldBlups(
  dat = pheno_excluding_16_including_16MP,
  identifier = "_excluding_16_including_16MP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_16_including_16MP %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_16_including_16MP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_16_including_16MP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_excluding_16_including_16MP,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_excluding_16_including_16MP,
  comp = "year:location:Variety"
)
Ve <- getVariance(
  asr = blups_excluding_16_including_16MP,
  comp = "Residual"
)

h_gryld_excluding_16_including_16MP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_16_including_16MP

blups_excluding_16_including_16RL <- GryldBlups(
  dat = pheno_excluding_16_including_16RL,
  identifier = "_excluding_16_including_16RL",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_16_including_16RL %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_16_including_16RL %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_16_including_16RL %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_excluding_16_including_16RL,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_excluding_16_including_16RL,
  comp = "year:location:Variety"
)
Ve <- getVariance(
  asr = blups_excluding_16_including_16RL,
  comp = "Residual"
)

h_gryld_excluding_16_including_16RL <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_16_including_16RL

blups_excluding_16_including_16RN <- GryldBlups(
  dat = pheno_excluding_16_including_16RN,
  identifier = "_excluding_16_including_16RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_16_including_16RN %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_16_including_16RN %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_16_including_16RN %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_excluding_16_including_16RN,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_excluding_16_including_16RN,
  comp = "year:location:Variety"
)
Ve <- getVariance(
  asr = blups_excluding_16_including_16RN,
  comp = "Residual"
)

h_gryld_excluding_16_including_16RN <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_16_including_16RN

blups_excluding_16_including_16SA <- GryldBlups(
  dat = pheno_excluding_16_including_16SA,
  identifier = "_excluding_16_including_16SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_16_including_16SA %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_16_including_16SA %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_16_including_16SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_excluding_16_including_16SA,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_excluding_16_including_16SA,
  comp = "year:location:Variety"
)
Ve <- getVariance(
  asr = blups_excluding_16_including_16SA,
  comp = "Residual"
)

h_gryld_excluding_16_including_16SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_16_including_16SA

blups_16 <- GryldBlups(
  dat = pheno_16, identifier = "_16",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_16 %>%
  group_by(Variety, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_16 %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16, comp = "Variety")
Vgl <- getVariance(asr = blups_16, comp = "location:Variety")
Ve <- getVariance(asr = blups_16, comp = "Residual")

h_gryld_16 <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_16

blups_16_ayn <- GryldBlups(
  dat = pheno_16_ayn, identifier = "_16_ayn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_16_ayn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_16_ayn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16_ayn, comp = "Variety")
Vgl <- getVariance(asr = blups_16_ayn, comp = "location:Variety")
Ve <- getVariance(asr = blups_16_ayn, comp = "Residual")

h_gryld_16_ayn <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_16_ayn

blups_16_pyn <- GryldBlups(
  dat = pheno_16_pyn, identifier = "_16_pyn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_16_pyn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_16_pyn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16_pyn, comp = "Variety")
Vgl <- getVariance(asr = blups_16_pyn, comp = "location:Variety")
Ve <- getVariance(asr = blups_16_pyn, comp = "Residual")

h_gryld_16_pyn <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_16_pyn

blups_16_MP_RL_RN <- GryldBlups(
  dat = pheno_16_MP_RL_RN, identifier = "_16_MP_RL_RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_16_MP_RL_RN %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_16_MP_RL_RN %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16_MP_RL_RN, comp = "Variety")
Vgl <- getVariance(asr = blups_16_MP_RL_RN, comp = "location:Variety")
Ve <- getVariance(asr = blups_16_MP_RL_RN, comp = "Residual")

h_gryld_16_MP_RL_RN <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_16_MP_RL_RN

blups_16_MP_RL_SA <- GryldBlups(
  dat = pheno_16_MP_RL_SA, identifier = "_16_MP_RL_SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_16_MP_RL_SA %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_16_MP_RL_SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16_MP_RL_SA, comp = "Variety")
Vgl <- getVariance(asr = blups_16_MP_RL_SA, comp = "location:Variety")
Ve <- getVariance(asr = blups_16_MP_RL_SA, comp = "Residual")

h_gryld_16_MP_RL_SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_16_MP_RL_SA

blups_16_MP_RN_SA <- GryldBlups(
  dat = pheno_16_MP_RN_SA, identifier = "_16_MP_RN_SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_16_MP_RN_SA %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_16_MP_RN_SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16_MP_RN_SA, comp = "Variety")
Vgl <- getVariance(asr = blups_16_MP_RN_SA, comp = "location:Variety")
Ve <- getVariance(asr = blups_16_MP_RN_SA, comp = "Residual")

h_gryld_16_MP_RN_SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_16_MP_RN_SA

blups_16_MP <- GryldBlups(
  dat = pheno_16_MP, identifier = "_16_MP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1| rep) +
    (1 | range) +
    (1 | column) +
    (1 | treated)"
)

rep_count <- pheno_16_MP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16_MP, comp = "Variety")
Ve <- getVariance(asr = blups_16_MP, comp = "Residual")

h_gryld_16_MP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_16_MP

blups_16_RL <- GryldBlups(
  dat = pheno_16_RL, identifier = "_16_RL",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1| rep) +
    (1 | range) +
    (1 | column) +
    (1 | treated)"
)

rep_count <- pheno_16_RL %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16_RL, comp = "Variety")
Ve <- getVariance(asr = blups_16_RL, comp = "Residual")

h_gryld_16_RL <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_16_RL

blups_16_RN <- GryldBlups(
  dat = pheno_16_RN, identifier = "_16_RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1| rep) +
    (1 | range) +
    (1 | column) +
    (1 | treated)"
)

rep_count <- pheno_16_RN %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16_RN, comp = "Variety")
Ve <- getVariance(asr = blups_16_RN, comp = "Residual")

h_gryld_16_RN <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_16_RN

blups_16_SA <- GryldBlups(
  dat = pheno_16_SA, identifier = "_16_SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1| rep) +
    (1 | range) +
    (1 | column) +
    (1 | treated)"
)

rep_count <- pheno_16_SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_16_SA, comp = "Variety")
Ve <- getVariance(asr = blups_16_SA, comp = "Residual")

h_gryld_16_SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_16_SA

#### 2016-2017 ####

blups_excluding_17 <- GryldBlups(
  dat = pheno_excluding_17, identifier = "_excluding_17",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_17 %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_17 %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_17 %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_17, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_17, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_17, comp = "Residual")

h_gryld_excluding_17 <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_17

blups_excluding_17_pyn <- GryldBlups(
  dat = pheno_excluding_17_pyn, identifier = "_excluding_17_pyn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_17_pyn %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_17_pyn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_17_pyn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_17_pyn, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_17_pyn, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_17_pyn, comp = "Residual")

h_gryld_excluding_17_pyn <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_17_pyn

blups_excluding_17_including_17MP <- GryldBlups(
  dat = pheno_excluding_17_including_17MP,
  identifier = "_excluding_17_including_17MP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_17_including_17MP %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_17_including_17MP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_17_including_17MP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_excluding_17_including_17MP,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_excluding_17_including_17MP,
  comp = "year:location:Variety"
)
Ve <- getVariance(
  asr = blups_excluding_17_including_17MP,
  comp = "Residual"
)

h_gryld_excluding_17_including_17MP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_17_including_17MP

blups_excluding_17_including_17RL <- GryldBlups(
  dat = pheno_excluding_17_including_17RL,
  identifier = "_excluding_17_including_17RL",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_17_including_17RL %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_17_including_17RL %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_17_including_17RL %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_excluding_17_including_17RL,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_excluding_17_including_17RL,
  comp = "year:location:Variety"
)
Ve <- getVariance(
  asr = blups_excluding_17_including_17RL,
  comp = "Residual"
)

h_gryld_excluding_17_including_17RL <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_17_including_17RL

blups_excluding_17_including_17RN <- GryldBlups(
  dat = pheno_excluding_17_including_17RN,
  identifier = "_excluding_17_including_17RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_17_including_17RN %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_17_including_17RN %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_17_including_17RN %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_excluding_17_including_17RN,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_excluding_17_including_17RN,
  comp = "year:location:Variety"
)
Ve <- getVariance(
  asr = blups_excluding_17_including_17RN,
  comp = "Residual"
)

h_gryld_excluding_17_including_17RN <-
  Vg / (Vg +
    (Ve / ((psych::harmonic.mean(year_count$n) *
      psych::harmonic.mean(location_count$n) *
      psych::harmonic.mean(rep_count$n)))) +
    (Vgl / ((psych::harmonic.mean(year_count$n) *
      psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_17_including_17RN

blups_excluding_17_including_17RP <- GryldBlups(
  dat = pheno_excluding_17_including_17RP,
  identifier = "_excluding_17_including_17RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_17_including_17RP %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_17_including_17RP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_17_including_17RP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_excluding_17_including_17RP,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_excluding_17_including_17RP,
  comp = "year:location:Variety"
)
Ve <- getVariance(
  asr = blups_excluding_17_including_17RP,
  comp = "Residual"
)

h_gryld_excluding_17_including_17RP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_17_including_17RP

blups_17 <- GryldBlups(
  dat = pheno_17, identifier = "_17",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_17 %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_17 %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17, comp = "Variety")
Vgl <- getVariance(asr = blups_17, comp = "location:Variety")
Ve <- getVariance(asr = blups_17, comp = "Residual")

h_gryld_17 <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_17

blups_17_ayn <- GryldBlups(
  dat = pheno_17_ayn, identifier = "_17_ayn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_17_ayn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_17_ayn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_ayn, comp = "Variety")
Vgl <- getVariance(asr = blups_17_ayn, comp = "location:Variety")
Ve <- getVariance(asr = blups_17_ayn, comp = "Residual")

h_gryld_17_ayn <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_17_ayn

blups_17_pyn <- GryldBlups(
  dat = pheno_17_pyn, identifier = "_17_pyn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_17_pyn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_17_pyn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_pyn, comp = "Variety")
Vgl <- getVariance(asr = blups_17_pyn, comp = "location:Variety")
Ve <- getVariance(asr = blups_17_pyn, comp = "Residual")

h_gryld_17_pyn <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_17_pyn

blups_17_MP_RL_RN <- GryldBlups(
  dat = pheno_17_MP_RL_RN, identifier = "_17_MP_RL_RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_17_MP_RL_RN %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_17_MP_RL_RN %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_MP_RL_RN, comp = "Variety")
Vgl <- getVariance(asr = blups_17_MP_RL_RN, comp = "location:Variety")
Ve <- getVariance(asr = blups_17_MP_RL_RN, comp = "Residual")

h_gryld_17_MP_RL_RN <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_17_MP_RL_RN

blups_17_MP_RL_RP <- GryldBlups(
  dat = pheno_17_MP_RL_RP, identifier = "_17_MP_RL_RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_17_MP_RL_RP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_17_MP_RL_RP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_MP_RL_RP, comp = "Variety")
Vgl <- getVariance(asr = blups_17_MP_RL_RP, comp = "location:Variety")
Ve <- getVariance(asr = blups_17_MP_RL_RP, comp = "Residual")

h_gryld_17_MP_RL_RP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_17_MP_RL_RP

blups_17_MP_RN_RP <- GryldBlups(
  dat = pheno_17_MP_RN_RP, identifier = "_17_MP_RN_RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_17_MP_RN_RP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_17_MP_RN_RP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_MP_RN_RP, comp = "Variety")
Vgl <- getVariance(asr = blups_17_MP_RN_RP, comp = "location:Variety")
Ve <- getVariance(asr = blups_17_MP_RN_RP, comp = "Residual")

h_gryld_17_MP_RN_RP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_17_MP_RN_RP

blups_17_RL_RN_RP <- GryldBlups(
  dat = pheno_17_RL_RN_RP, identifier = "_17_RL_RN_RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1 | location) +
    (1 | location:Variety) +
    (1 | location:rep) +
    (1 | location:treated)"
)

location_count <- pheno_17_RL_RN_RP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_17_RL_RN_RP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_RL_RN_RP, comp = "Variety")
Vgl <- getVariance(asr = blups_17_RL_RN_RP, comp = "location:Variety")
Ve <- getVariance(asr = blups_17_RL_RN_RP, comp = "Residual")

h_gryld_17_RL_RN_RP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_17_RL_RN_RP

blups_17_MP <- GryldBlups(
  dat = pheno_17_MP, identifier = "_17_MP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1| rep) +
    (1 | range) +
    (1 | column) +
    (1 | treated)"
)

rep_count <- pheno_17_MP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_MP, comp = "Variety")
Ve <- getVariance(asr = blups_17_MP, comp = "Residual")

h_gryld_17_MP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_17_MP

blups_17_RL <- GryldBlups(
  dat = pheno_17_RL, identifier = "_17_RL",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1| rep) +
    (1 | range) +
    (1 | column) +
    (1 | treated)"
)

rep_count <- pheno_17_RL %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_RL, comp = "Variety")
Ve <- getVariance(asr = blups_17_RL, comp = "Residual")

h_gryld_17_RL <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_17_RL

blups_17_RN <- GryldBlups(
  dat = pheno_17_RN, identifier = "_17_RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1| rep) +
    (1 | range) +
    (1 | column) "
)

rep_count <- pheno_17_RN %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_RN, comp = "Variety")
Ve <- getVariance(asr = blups_17_RN, comp = "Residual")

h_gryld_17_RN <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_17_RN

blups_17_RP <- GryldBlups(
  dat = pheno_17_RP, identifier = "_17_RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
    (1| rep) +
    (1 | range) +
    (1 | column) +
    (1 | treated)"
)

rep_count <- pheno_17_RP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_17_RP, comp = "Variety")
Ve <- getVariance(asr = blups_17_RP, comp = "Residual")

h_gryld_17_RP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_17_RP

#### 2017-2018 ####

blups_excluding_18 <- GryldBlups(
  dat = pheno_excluding_18, identifier = "_excluding_18",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_18 %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_18 %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_18 %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_18, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_18, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_18, comp = "Residual")

h_gryld_excluding_18 <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_18

blups_excluding_18_pyn <- GryldBlups(
  dat = pheno_excluding_18_pyn, identifier = "_excluding_18_pyn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_18_pyn %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_18_pyn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_18_pyn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_18_pyn, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_18_pyn, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_18_pyn, comp = "Residual")

h_gryld_excluding_18_pyn <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_18_pyn

blups_excluding_18_including_18MP <- GryldBlups(
  dat = pheno_excluding_18_including_18MP, identifier = "_excluding_18_including_18MP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_18_including_18MP %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_18_including_18MP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_18_including_18MP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_18_including_18MP, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_18_including_18MP, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_18_including_18MP, comp = "Residual")

h_gryld_excluding_18_including_18MP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_18_including_18MP

blups_excluding_18_including_18RL <- GryldBlups(
  dat = pheno_excluding_18_including_18RL, identifier = "_excluding_18_including_18RL",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_18_including_18RL %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_18_including_18RL %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_18_including_18RL %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_excluding_18_including_18RL,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_excluding_18_including_18RL,
  comp = "year:location:Variety"
)
Ve <- getVariance(
  asr = blups_excluding_18_including_18RL,
  comp = "Residual"
)

h_gryld_excluding_18_including_18RL <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_18_including_18RL

blups_excluding_18_including_18RL <- GryldBlups(
  dat = pheno_excluding_18_including_18RL, identifier = "_excluding_18_including_18RL",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_18_including_18RL %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_18_including_18RL %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_18_including_18RL %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_18_including_18RL, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_18_including_18RL, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_18_including_18RL, comp = "Residual")

h_gryld_excluding_18_including_18RL <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_18_including_18RL

blups_excluding_18_including_18RN <- GryldBlups(
  dat = pheno_excluding_18_including_18RN, identifier = "_excluding_18_including_18RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_18_including_18RN %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_18_including_18RN %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_18_including_18RN %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_18_including_18RN, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_18_including_18RN, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_18_including_18RN, comp = "Residual")

h_gryld_excluding_18_including_18RN <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_18_including_18RN

blups_excluding_18_including_18RP <- GryldBlups(
  dat = pheno_excluding_18_including_18RP, identifier = "_excluding_18_including_18RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_18_including_18RP %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_18_including_18RP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_18_including_18RP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_18_including_18RP, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_18_including_18RP, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_18_including_18RP, comp = "Residual")

h_gryld_excluding_18_including_18RP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_18_including_18RP

blups_excluding_18_including_18SA <- GryldBlups(
  dat = pheno_excluding_18_including_18SA, identifier = "_excluding_18_including_18SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_18_including_18SA %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_18_including_18SA %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_18_including_18SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_18_including_18SA, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_18_including_18SA, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_18_including_18SA, comp = "Residual")

h_gryld_excluding_18_including_18SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_18_including_18SA

blups_18 <- GryldBlups(
  dat = pheno_18, identifier = "_18",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_18 %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_18 %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18, comp = "Variety")
Vgl <- getVariance(asr = blups_18, comp = "location:Variety")
Ve <- getVariance(asr = blups_18, comp = "Residual")

h_gryld_18 <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_18

blups_18_ayn <- GryldBlups(
  dat = pheno_18_ayn, identifier = "_18_ayn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_18_ayn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_18_ayn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_ayn, comp = "Variety")
Vgl <- getVariance(asr = blups_18_ayn, comp = "location:Variety")
Ve <- getVariance(asr = blups_18_ayn, comp = "Residual")

h_gryld_18_ayn <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_18_ayn

blups_18_pyn <- GryldBlups(
  dat = pheno_18_pyn, identifier = "_18_pyn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_18_pyn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_18_pyn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_pyn, comp = "Variety")
Vgl <- getVariance(asr = blups_18_pyn, comp = "location:Variety")
Ve <- getVariance(asr = blups_18_pyn, comp = "Residual")

h_gryld_18_pyn <- Vg / (Vg +
                          (Ve / ((psych::harmonic.mean(year_count$n) *
                                    psych::harmonic.mean(location_count$n) *
                                    psych::harmonic.mean(rep_count$n)))) +
                          (Vgl / ((psych::harmonic.mean(year_count$n) *
                                     psych::harmonic.mean(location_count$n)))))
h_gryld_18_pyn

blups_18_MP_RL_RN_RP <- GryldBlups(
  dat = pheno_18_MP_RL_RN_RP, identifier = "_18_MP_RL_RN_RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_18_MP_RL_RN_RP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_18_MP_RL_RN_RP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_MP_RL_RN_RP, comp = "Variety")
Vgl <- getVariance(asr = blups_18_MP_RL_RN_RP, comp = "location:Variety")
Ve <- getVariance(asr = blups_18_MP_RL_RN_RP, comp = "Residual")

h_gryld_18_MP_RL_RN_RP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_18_MP_RL_RN_RP

blups_18_MP_RL_RN_SA <- GryldBlups(
  dat = pheno_18_MP_RL_RN_SA, identifier = "_18_MP_RL_RN_SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_18_MP_RL_RN_SA %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_18_MP_RL_RN_SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_MP_RL_RN_SA, comp = "Variety")
Vgl <- getVariance(asr = blups_18_MP_RL_RN_SA, comp = "location:Variety")
Ve <- getVariance(asr = blups_18_MP_RL_RN_SA, comp = "Residual")

h_gryld_18_MP_RL_RN_SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_18_MP_RL_RN_SA

blups_18_MP_RL_RP_SA <- GryldBlups(
  dat = pheno_18_MP_RL_RP_SA, identifier = "_18_MP_RL_RP_SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_18_MP_RL_RP_SA %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_18_MP_RL_RP_SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_MP_RL_RP_SA, comp = "Variety")
Vgl <- getVariance(asr = blups_18_MP_RL_RP_SA, comp = "location:Variety")
Ve <- getVariance(asr = blups_18_MP_RL_RP_SA, comp = "Residual")

h_gryld_18_MP_RL_RP_SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_18_MP_RL_RP_SA

blups_18_MP_RN_RP_SA <- GryldBlups(
  dat = pheno_18_MP_RN_RP_SA, identifier = "_18_MP_RN_RP_SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_18_MP_RN_RP_SA %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_18_MP_RN_RP_SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_MP_RN_RP_SA, comp = "Variety")
Vgl <- getVariance(asr = blups_18_MP_RN_RP_SA, comp = "location:Variety")
Ve <- getVariance(asr = blups_18_MP_RN_RP_SA, comp = "Residual")

h_gryld_18_MP_RN_RP_SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_18_MP_RN_RP_SA

blups_18_RL_RN_RP_SA <- GryldBlups(
  dat = pheno_18_RL_RN_RP_SA, identifier = "_18_RL_RN_RP_SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_18_RL_RN_RP_SA %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_18_RL_RN_RP_SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_RL_RN_RP_SA, comp = "Variety")
Vgl <- getVariance(asr = blups_18_RL_RN_RP_SA, comp = "location:Variety")
Ve <- getVariance(asr = blups_18_RL_RN_RP_SA, comp = "Residual")

h_gryld_18_RL_RN_RP_SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_18_RL_RN_RP_SA

blups_18_MP <- GryldBlups(
  dat = pheno_18_MP, identifier = "_18_MP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | range) +
                         (1 | column)"
)

rep_count <- pheno_18_MP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_MP, comp = "Variety")
Ve <- getVariance(asr = blups_18_MP, comp = "Residual")

h_gryld_18_MP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_18_MP

blups_18_RL <- GryldBlups(
  dat = pheno_18_RL, identifier = "_18_RL",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | range) +
                         (1 | column)"
)

rep_count <- pheno_18_RL %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_RL, comp = "Variety")
Ve <- getVariance(asr = blups_18_RL, comp = "Residual")

h_gryld_18_RL <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_18_RL

blups_18_RN <- GryldBlups(
  dat = pheno_18_RN, identifier = "_18_RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | range) +
                         (1 | column)"
)

rep_count <- pheno_18_RN %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_RN, comp = "Variety")
Ve <- getVariance(asr = blups_18_RN, comp = "Residual")

h_gryld_18_RN <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_18_RN

blups_18_RP <- GryldBlups(
  dat = pheno_18_RP, identifier = "_18_RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | range) +
                         (1 | column)"
)

rep_count <- pheno_18_RP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_RP, comp = "Variety")
Ve <- getVariance(asr = blups_18_RP, comp = "Residual")

h_gryld_18_RP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_18_RP

blups_18_SA <- GryldBlups(
  dat = pheno_18_SA, identifier = "_18_SA",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | range) +
                         (1 | column)"
)

rep_count <- pheno_18_SA %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_18_SA, comp = "Variety")
Ve <- getVariance(asr = blups_18_SA, comp = "Residual")

h_gryld_18_SA <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_18_SA

#### 2018-2019 ####

blups_excluding_19 <- GryldBlups(
  dat = pheno_excluding_19, identifier = "_excluding_19",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_19 %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_19 %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_19 %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_19, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_19, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_19, comp = "Residual")

h_gryld_excluding_19 <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_19

blups_excluding_19_pyn <- GryldBlups(
  dat = pheno_excluding_19_pyn, identifier = "_excluding_19_pyn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_19_pyn %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_19_pyn %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_19_pyn %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_19_pyn, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_19_pyn, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_19_pyn, comp = "Residual")

h_gryld_excluding_19_pyn <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_19_pyn

blups_excluding_19_including_19RL <- GryldBlups(
  dat = pheno_excluding_19_including_19RL, identifier = "_excluding_19_including_19RL",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_19_including_19RL %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_19_including_19RL %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_19_including_19RL %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_19_including_19RL, 
                  comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_19_including_19RL, 
                   comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_19_including_19RL, 
                  comp = "Residual")

h_gryld_excluding_19_including_19RL <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_19_including_19RL

blups_excluding_19_including_19RN <- GryldBlups(
  dat = pheno_excluding_19_including_19RN, 
  identifier = "_excluding_19_including_19RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_19_including_19RN %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_19_including_19RN %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_19_including_19RN %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_19_including_19RN,
                  comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_19_including_19RN, 
                   comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_19_including_19RN, 
                  comp = "Residual")

h_gryld_excluding_19_including_19RN <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_19_including_19RN

blups_excluding_19_including_19RP <- GryldBlups(
  dat = pheno_excluding_19_including_19RP, identifier = "_excluding_19_including_19RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | year) +
                         (1 | year:location) +
                         (1 | year:location:Variety) +
                         (1 | year:location:rep) +
                         (1 | year:location:treated)"
)

year_count <- pheno_excluding_19_including_19RP %>%
  group_by(Variety, year) %>%
  summarise(n = n()) %>%
  summarise(n = n())

location_count <- pheno_excluding_19_including_19RP %>%
  group_by(Variety, year, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_excluding_19_including_19RP %>%
  group_by(Variety, year, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(asr = blups_excluding_19_including_19RP, comp = "Variety")
Vgl <- getVariance(asr = blups_excluding_19_including_19RP, comp = "year:location:Variety")
Ve <- getVariance(asr = blups_excluding_19_including_19RP, comp = "Residual")

h_gryld_excluding_19_including_19RP <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(year_count$n) *
    psych::harmonic.mean(location_count$n)))))
h_gryld_excluding_19_including_19RP

blups_19 <- GryldBlups(
  dat = pheno_19, identifier = "_19",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_19 %>%
  group_by(Variety, location) %>%
  summarise(n = n())

rep_count <- pheno_19 %>%
  group_by(Variety, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_19,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_19,
  comp = "location:Variety"
)
Ve <- getVariance(asr = blups_19, comp = "Residual")

h_gryld_19 <- Vg / (Vg +
  (Ve / ((psych::harmonic.mean(location_count$n) *
    psych::harmonic.mean(rep_count$n)))) +
  (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_19

blups_19_ayn <- GryldBlups(
  dat = pheno_19_ayn, identifier = "_19_ayn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_19_ayn %>%
  group_by(Variety, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_19_ayn %>%
  group_by(Variety, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_19_ayn,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_19_ayn,
  comp = "location:Variety"
)
Ve <- getVariance(asr = blups_19_ayn, comp = "Residual")

h_gryld_19_ayn <- Vg / (Vg +
                      (Ve / ((psych::harmonic.mean(location_count$n) *
                                psych::harmonic.mean(rep_count$n)))) +
                      (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_19_ayn

blups_19_RL_RN <- GryldBlups(
  dat = pheno_19_RL_RN, identifier = "_19_RL_RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_19_RL_RN %>%
  group_by(Variety, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_19_RL_RN %>%
  group_by(Variety, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_19_RL_RN,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_19_RL_RN,
  comp = "location:Variety"
)
Ve <- getVariance(asr = blups_19_RL_RN, comp = "Residual")

h_gryld_19_RL_RN <- Vg / (Vg +
                          (Ve / ((psych::harmonic.mean(location_count$n) *
                                    psych::harmonic.mean(rep_count$n)))) +
                          (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_19_RL_RN

blups_19_RL_RP <- GryldBlups(
  dat = pheno_19_RL_RP, identifier = "_19_RL_RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_19_RL_RP %>%
  group_by(Variety, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_19_RL_RP %>%
  group_by(Variety, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_19_RL_RP,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_19_RL_RP,
  comp = "location:Variety"
)
Ve <- getVariance(asr = blups_19_RL_RP, comp = "Residual")

h_gryld_19_RL_RP <- Vg / (Vg +
                          (Ve / ((psych::harmonic.mean(location_count$n) *
                                    psych::harmonic.mean(rep_count$n)))) +
                          (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_19_RL_RP

blups_19_RN_RP <- GryldBlups(
  dat = pheno_19_RN_RP, identifier = "_19_RN_RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | location) +
                         (1 | location:Variety) +
                         (1 | location:rep) +
                         (1 | location:treated)"
)

location_count <- pheno_19_RN_RP %>%
  group_by(Variety, location) %>%
  summarise(n = n()) %>%
  summarise(n = n())

rep_count <- pheno_19_RN_RP %>%
  group_by(Variety, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_19_RN_RP,
  comp = "Variety"
)
Vgl <- getVariance(
  asr = blups_19_RN_RP,
  comp = "location:Variety"
)
Ve <- getVariance(asr = blups_19_RN_RP, comp = "Residual")

h_gryld_19_RN_RP <- Vg / (Vg +
                          (Ve / ((psych::harmonic.mean(location_count$n) *
                                    psych::harmonic.mean(rep_count$n)))) +
                          (Vgl / ((psych::harmonic.mean(location_count$n)))))
h_gryld_19_RN_RP

blups_19_pyn <- GryldBlups(
  dat = pheno_19_pyn, identifier = "_19_pyn",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | range) +
                         (1 | column)"
)

rep_count <- pheno_19_pyn %>%
  group_by(Variety, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_19_pyn,
  comp = "Variety"
)

Ve <- getVariance(asr = blups_19_pyn, comp = "Residual")

h_gryld_19_pyn <- Vg / (Vg +
                          (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_19_pyn

blups_19_RL <- GryldBlups(
  dat = pheno_19_RL, identifier = "_19_RL",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | range) +
                         (1 | column)"
)

rep_count <- pheno_19_RL %>%
  group_by(Variety, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_19_RL,
  comp = "Variety"
)

Ve <- getVariance(asr = blups_19_RL, comp = "Residual")

h_gryld_19_RL <- Vg / (Vg +
                          (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_19_RL

blups_19_RN <- GryldBlups(
  dat = pheno_19_RN, identifier = "_19_RN",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | range) +
                         (1 | column)"
)

rep_count <- pheno_19_RN %>%
  group_by(Variety, location, rep) %>%
  summarise(n = n())  %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_19_RN,
  comp = "Variety"
)

Ve <- getVariance(asr = blups_19_RN, comp = "Residual")

h_gryld_19_RN <- Vg / (Vg +
                          (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_19_RN

blups_19_RP <- GryldBlups(
  dat = pheno_19_RP, identifier = "_19_RP",
  saveFile = "./Results/",
  saveFigures = "./Figures/",
  equation = "phenotype_value ~ (1 | Variety) +
                         (1 | range) +
                         (1 | column)"
)

rep_count <- pheno_19_RP %>%
  group_by(Variety, location, rep) %>%
  summarise(n = n()) %>%
  summarise(n = n())

Vg <- getVariance(
  asr = blups_19_RP,
  comp = "Variety"
)

Ve <- getVariance(asr = blups_19_RP, comp = "Residual")

h_gryld_19_RP <- Vg / (Vg +
                          (Ve / ((psych::harmonic.mean(rep_count$n)))))
h_gryld_19_RP

##### All and Plots ####
gryld_heritability <- tibble(
  h_gryld, h_gryld_16, h_gryld_16_ayn, h_gryld_16_pyn,
  h_gryld_17, h_gryld_17_ayn, h_gryld_17_pyn,
  h_gryld_18, h_gryld_18_ayn, h_gryld_18_pyn,
  h_gryld_19, h_gryld_19_ayn)

gryld_heritability <- gryld_heritability %>%
  pivot_longer(
    cols = everything(),
    names_to = "data",
    values_to = "Heritability"
  ) %>%
  separate(col = data, into = c("h", "gryld", "Year", "Trial"), sep = "_") %>%
  mutate(Trial = replace_na(Trial, "All"))

gryld_heritability %>%
  ggplot(aes(
    x = Year,
    y = Heritability,
    colour = Trial
  )) +
  geom_point(size = 2) +
  labs(title = "GRYLD Broad-sense heritabilities")
