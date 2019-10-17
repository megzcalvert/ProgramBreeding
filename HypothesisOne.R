rm(list = objects())
ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(RMySQL)
library(janitor)
library(broom)
library(reshape2)

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

# Fix all the annoying naming inconsistencies
pheno_long <- pheno_long %>%
  mutate(
    Variety = str_replace(Variety, "~", "-"),
    Variety = str_replace(Variety, "-K-", "K-"),
    Variety = str_replace(Variety, "-M-", "M-"),
    Variety = str_replace(Variety, " ", "-"),
    Variety = str_replace(Variety, "-", "_"),
    Variety = tolower(Variety)
  ) %>%
  mutate(
    treatment = if_else(treatment == "Fungicide", "TREATED", "UNTREATED")
  )
glimpse(pheno_long)

## Divide out all important info from the entity_id
# remove awns and pcthead as they were taken sporadically

pheno_long <- pheno_long %>%
  mutate(
    Sep = entity_id,
    phenotype_value = as.numeric(phenotype_value),
    ID = row_number()
  ) %>%
  separate(Sep,
    c("Year", "Trial", "Location", "Treated", "Plot"),
    sep = "-"
  ) %>%
  filter(trait_id != "AWNS") %>%
  filter(trait_id != "PCTHEAD") %>%
  filter(phenotype_value >= 0) %>%
  filter(phenotype_value <= 6500)

unique(pheno_long$Trial)

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

# Adding an overall trial column
pheno_long$OverallTrial <- ff(pheno_long$Trial,
  c(
    "AYN3", "PYN", "GSPYN", "AYN2", "DHPYN", "AYN1",
    "AYN4", "PYN1", "PYN2", "DHPYN1", "DHPYN2", "DHAYN2",
    "DHAYN1", "PYNA", "PYNB", "PYN3", "PYN4", "GSAYN",
    "ATAYN1", "ATAYN2"
  ),
  c(
    "AYN", "PYN", "PYN", "AYN", "PYN", "AYN", "AYN",
    "PYN", "PYN", "PYN", "PYN", "AYN", "AYN", "PYN",
    "PYN", "PYN", "PYN", "AYN", "AYN", "AYN"
  ),
  "NA",
  ignore.case = TRUE
)

# adjusting factors
pheno_long <- pheno_long %>%
  mutate(
    Location = as.factor(Location),
    Trial = as.factor(Trial),
    Year = as.factor(Year)
  )

# Which phenotypic traits are we working with
traits <- unique(pheno_long$trait_id)

for (i in traits) {
  p <- pheno_long %>%
    tidylog::filter(trait_id == paste(i)) %>%
    ggplot(aes(x = phenotype_value, colour = Location)) +
    geom_density() +
    facet_wrap(Year ~ Trial, scales = "free") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black")) +
    labs(title = paste(i))
  ggsave(paste0("./Figures/Summary_", i, ".png"),
    plot = p, width = 16,
    height = 16,
    units = "cm"
  )
  print(p)
}


write.table(pheno_long, "./PhenoDatabase/PhenoLong.txt",
  quote = F, sep = "\t",
  col.names = T, row.names = F
)

# separating into traits
unique(pheno_long$trait_id)
phenoGryld <- pheno_long %>%
  tidylog::filter(trait_id == "GRYLD")
phenoPtht <- pheno_long %>%
  tidylog::filter(trait_id == "PTHT")
phenoMoist <- pheno_long %>%
  tidylog::filter(trait_id == "MOIST")
phenoTestwt <- pheno_long %>%
  tidylog::filter(trait_id == "TESTWT")

gryld <- phenoGryld %>%
  group_by(Year) %>%
  do(tidy(anova(lm(phenotype_value ~ Trial * Location, data = .))))

ptht <- phenoPtht %>%
  group_by(Year) %>%
  do(tidy(anova(lm(phenotype_value ~ Trial * Location, data = .))))

moist <- phenoMoist %>%
  group_by(Year) %>%
  do(tidy(anova(lm(phenotype_value ~ Trial * Location, data = .))))

testwt <- phenoTestwt %>%
  group_by(Year) %>%
  do(tidy(anova(lm(phenotype_value ~ Trial * Location, data = .))))

pheno_longRep <- pheno_long %>%
  tidylog::select(
    -trait_id, -phenotype_value, -phenotype_date,
    -phenotype_person, -ID, -range, -column, -OverallTrial,
    -location, -block, -rep, -treatment, -Plot
  ) %>%
  group_by(Year, Location, Trial, Treated, Variety) %>%
  distinct() %>%
  mutate(repAssigned = row_number()) %>%
  inner_join(pheno_long, by = c(
    "entity_id", "Variety", "Year", "Trial",
    "Location", "Treated"
  )) %>%
  select(
    ID, entity_id, Variety, Year, Trial,
    Location, Treated, treatment, repAssigned, block, range, column,
    OverallTrial, trait_id, phenotype_value,
    phenotype_date, phenotype_person
  )

write.table(pheno_longRep, "./PhenoDatabase/PhenoLongRep.txt",
  quote = F,
  sep = "\t", col.names = T, row.names = F
)

glimpse(pheno_long)

TrialSummary <- pheno_long %>%
  group_by(Year, Location, Trial) %>%
  summarise(n = n()) %>%
  glimpse()

VarietySummary <- pheno_long %>%
  group_by(Year, Location, Trial, Variety) %>%
  summarise(n = n())

VarSummary <- pheno_long %>%
  group_by(Variety, Year, Location, Trial) %>%
  summarise(n = n())
write.table(VarSummary,
  file = "./PhenoDatabase/VarietySummary.txt", quote = F,
  sep = "\t", row.names = F, col.names = T
)
