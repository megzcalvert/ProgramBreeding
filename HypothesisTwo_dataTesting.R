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
    legend.key.size = unit(4, "lines"),
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

pheno_long <- pheno_long %>%
  filter(
    Variety != "fill",
    Variety != "centralblend"
  ) %>%
  select(
    entity_id, trait_id, phenotype_value, phenotype_date, phenotype_person,
    Variety, range, column
  ) %>%
  mutate(Sep = entity_id) %>%
  separate(Sep, c("year", "trial", "location", "treated", "plot"), sep = "-")

trialSummary <- tabyl(pheno_long, trial)
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

tabyl(pheno_long, Trial)

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

pheno_long <- pheno_long %>%
  filter(trait_id == "GRYLD")
glimpse(pheno_long)

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
  phenotype_value ~ (1 | Variety) + (1 | year) + (1 | location) + (1 | treated)
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

tidy(fit_1)
tidy(fit_2)
tidy(fit_3)
tidy(fit_4)
tidy(fit_5)
tidy(fit_6)

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
  facet_wrap(Trial ~ year)

pheno_long %>%
  ggplot(aes(phenotype_value, colour = location, linetype = treated)) +
  geom_density(fill = NA) +
  facet_wrap(Trial ~ year)

values <- setDT(as.data.frame(ranef(fit_6)$Variety), keep.rownames = TRUE) %>%
  rename(
    GRYLD = `(Intercept)`,
    Variety = rn
  )

values %>%
  ggplot(aes(x = GRYLD)) +
  geom_density(fill = NA)

write.table(values, "./PhenoDatabase/BlupsGRYLD_allYrs_allVarieties.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

varianceComponents <- as.data.frame(VarCorr(fit_6))
varianceComponents

h_gryld <- varianceComponents[2, 4] / (varianceComponents[2, 4] +
  (varianceComponents[9, 4] /
    psych::harmonic.mean(as.numeric(pheno_long$rep)) *
    length(unique(pheno_long$year)) *
    length(unique(pheno_long$location))) +
  (varianceComponents[1, 4] /
    length(unique(pheno_long$year)) *
    length(unique(pheno_long$location))))
h_gryld

