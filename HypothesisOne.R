rm(list = objects())
ls()

library(data.table)
library(tidyverse)
library(tidylog)
library(janitor)
library(broom)

#### Setting up workspace ####

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
      size = rel(2)
    ),
    complete = F
  )

theme_set(custom_theme)

getwd()

# useful infos **reproducible research**
sessionInfo()
options(scipen = 999)

################## Initial scenario ##################

budget <- 100000
## per breeding *cycle*
Line_dev_cost <- 10
## the cost of generating an inbred line, applied to either form of evaluation
ayn_number <- 100
## number of AYN to be planted

pheno_cost <- 85
## cost to phenotype in field
pheno_rep <- 6
## number of phenotyping plots planted
pheno_cycle_length <- 1
## number of years phenotypic selection cycle takes

geno_cost <- 5
## cost to genotype
geno_cycle_length <- 1
## number of years prediction selection cycle takes

# Number of plots if only phenotyping
only_pheno <- floor(budget / (Line_dev_cost + (pheno_cost * pheno_rep)))

# Number of plots if only predicting
only_geno <- floor(budget / (Line_dev_cost + geno_cost))

## Visualisng the parameter space
pheno_plots <- c(0, only_pheno)
geno_plots <- c(only_geno, 0)

dat <- tibble::tibble(
  cost = geno_plots, pheno_plots
)

basicPlot <- ggplot(
  data = dat,
  aes(
    x = pheno_plots,
    y = geno_plots
  )
) +
  geom_point() +
  geom_path()
basicPlot

#### Defining it for the entire parameter space ####

SearchingGeneticGainParameters <- function(budget = budget,
                                           # per breeding *cycle*
                                           Line_dev_cost = Line_dev_cost,
                                           # the cost of generating an inbred line, applied to either form of evaluation
                                           ayn_number = ayn_number,
                                           # number of AYN to be planted
                                           pheno_cost = pheno_cost,
                                           # cost to phenotype in field
                                           pheno_rep = pheno_rep,
                                           # number of phenotyping plots planted
                                           pheno_cycle_length = pheno_cycle_length,
                                           # number of years phenotypic selection cycle takes
                                           h_pheno = h_pheno,
                                           # Heritability of yield in pheno
                                           geno_cost = geno_cost,
                                           # cost to genotype
                                           geno_cycle_length = geno_cycle_length,
                                           # number of years prediction selection cycle takes
                                           pred_accuracy = pred_accuracy,
                                           # Prediction accuracy that is synonymous with genetic correlation
                                           h_predTrait = h_predTrait # assumes the genome in an inbred is 100% heritable and all additive
) {
  # Number of plots if only phenotyping
  only_pheno <- floor(budget / (Line_dev_cost + (pheno_cost * pheno_rep)))
  cost_per_pheno <- round(budget / only_pheno, digits = 2)
  number_pheno <- as.integer(0:only_pheno)

  # Number of plots if only predicting
  only_geno <- floor(budget / (Line_dev_cost + geno_cost))
  cost_per_geno <- round(budget / only_geno, digits = 2)
  number_geno <- as.integer(0:only_geno)

  # Creating dataframe of all combinations
  df <- crossing(number_geno, number_pheno) %>%
    mutate(
      geno_spend = number_geno * cost_per_geno, # spent on genotyping
      pheno_spend = number_pheno * cost_per_pheno, # spent on phenotyping
      total_spend = geno_spend + pheno_spend, # total money spent
      proportion_geno = geno_spend / total_spend, # proportion $ in predictin
      proportion_pheno = pheno_spend / total_spend
    ) %>% # proportion $ in phenotyping
    filter(
      total_spend >= 0.5 * budget,
      total_spend <= budget
    ) %>%
    mutate(
      geno_ayn = round(ayn_number * proportion_geno), # AYN selected by prediction
      pheno_ayn = round(ayn_number * proportion_pheno), # AYN selected by phenotyping
      correlated_intensity = geno_ayn / number_geno,
      direct_intensity = pheno_ayn / number_pheno,
      correlated_response =
        (correlated_intensity * pred_accuracy * sqrt(h_predTrait)) /
          geno_cycle_length, # correlated response
      direct_response =
        (direct_intensity * sqrt(h_pheno)) /
          pheno_cycle_length, # direct response
      response_ratio = correlated_response / direct_response
    ) %>%
    filter(
      !is.infinite(response_ratio),
      !is.na(response_ratio),
      response_ratio >= 0
    )

  gainPlot <- ggplot(
    data = df,
    aes(
      x = number_pheno,
      y = number_geno,
      fill = response_ratio
    )
  ) +
    geom_tile() + #interpolate = FALSE) +
    scale_fill_gradient2(
      low = "#762a83",
      mid = "#f7f7f7",
      high = "#1b7837",
      midpoint = quantile(df$response_ratio, 0.75),
      limits = c(min(df$response_ratio), max(df$response_ratio))
    ) +
    guides(fill = guide_colourbar(nbin = 15)) +
    labs(
      title = paste0(
        "Parameter space for Budget $", budget,
        " and AYN number ", ayn_number
      ),
      subtitle = paste0(
        "pheno cost = $", pheno_cost,
        ", geno cost = $", geno_cost,
        ", direct h2 = ", h_pheno,
        ", correlated h2 = ", h_predTrait,
        ", r = ", pred_accuracy,
        ", geno time as proportion pheno time = ",
        geno_cycle_length,
        ", line development cost = $", Line_dev_cost
      )
    )

  print(gainPlot)

  return(df)
}


trial1 <- SearchingGeneticGainParameters(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 10, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 100, # number of AYN to be planted
  pheno_cost = 85, # cost to phenotype in field
  pheno_rep = 6, # number of phenotyping plots planted
  pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  h_pheno = 0.1, # Heritability of yield in pheno
  geno_cost = 5, # cost to genotype
  geno_cycle_length = 0.25, # number of years prediction selection cycle takes
  pred_accuracy = 0.85, # Prediction accuracy that is synonymous with genetic correlation
  h_predTrait = 1
)

trial2 <- SearchingGeneticGainParameters(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 10, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 100, # number of AYN to be planted
  pheno_cost = 85, # cost to phenotype in field
  pheno_rep = 6, # number of phenotyping plots planted
  pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  h_pheno = 0.1, # Heritability of yield in pheno
  geno_cost = 5, # cost to genotype
  geno_cycle_length = 0.2, # number of years prediction selection cycle takes
  pred_accuracy = 0.95, # Prediction accuracy that is synonymous with genetic correlation
  h_predTrait = 1
)

trial3 <- SearchingGeneticGainParameters(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 10, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 200, # number of AYN to be planted
  pheno_cost = 85, # cost to phenotype in field
  pheno_rep = 5, # number of phenotyping plots planted
  pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  h_pheno = 0.1, # Heritability of yield in pheno
  geno_cost = 5, # cost to genotype
  geno_cycle_length = 0.2, # number of years prediction selection cycle takes
  pred_accuracy = 0.9, # Prediction accuracy that is synonymous with genetic correlation
  h_predTrait = 1
)

trial4 <- SearchingGeneticGainParameters(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 20, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 100, # number of AYN to be planted
  pheno_cost = 85, # cost to phenotype in field
  pheno_rep = 6, # number of phenotyping plots planted
  pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  h_pheno = 0.1, # Heritability of yield in pheno
  geno_cost = 5, # cost to genotype
  geno_cycle_length = 0.2, # number of years prediction selection cycle takes
  pred_accuracy = 0.85, # Prediction accuracy that is synonymous with genetic correlation
  h_predTrait = 1
)

trial5 <- SearchingGeneticGainParameters(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 25, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 100, # number of AYN to be planted
  pheno_cost = 85, # cost to phenotype in field
  pheno_rep = 6, # number of phenotyping plots planted
  pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  h_pheno = 0.1, # Heritability of yield in pheno
  geno_cost = 5, # cost to genotype
  geno_cycle_length = 0.2, # number of years prediction selection cycle takes
  pred_accuracy = 0.85, # Prediction accuracy that is synonymous with genetic correlation
  h_predTrait = 1
)

trial6 <- SearchingGeneticGainParameters(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 30, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 100, # number of AYN to be planted
  pheno_cost = 85, # cost to phenotype in field
  pheno_rep = 6, # number of phenotyping plots planted
  pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  h_pheno = 0.1, # Heritability of yield in pheno
  geno_cost = 5, # cost to genotype
  geno_cycle_length = 0.2, # number of years prediction selection cycle takes
  pred_accuracy = 0.85, # Prediction accuracy that is synonymous with genetic correlation
  h_predTrait = 1
)

trial7 <- SearchingGeneticGainParameters(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 10, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 100, # number of AYN to be planted
  pheno_cost = 85, # cost to phenotype in field
  pheno_rep = 6, # number of phenotyping plots planted
  pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  h_pheno = 0.2, # Heritability of yield in pheno
  geno_cost = 5, # cost to genotype
  geno_cycle_length = 0.2, # number of years prediction selection cycle takes
  pred_accuracy = 0.85, # Prediction accuracy that is synonymous with genetic correlation
  h_predTrait = 1
)

trial8 <- SearchingGeneticGainParameters(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 10, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 100, # number of AYN to be planted
  pheno_cost = 85, # cost to phenotype in field
  pheno_rep = 6, # number of phenotyping plots planted
  pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  h_pheno = 0.05, # Heritability of yield in pheno
  geno_cost = 5, # cost to genotype
  geno_cycle_length = 0.2, # number of years prediction selection cycle takes
  pred_accuracy = 0.85, # Prediction accuracy that is synonymous with genetic correlation
  h_predTrait = 1
)

trial9 <- SearchingGeneticGainParameters(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 10, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 100, # number of AYN to be planted
  pheno_cost = 125, # cost to phenotype in field
  pheno_rep = 6, # number of phenotyping plots planted
  pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  h_pheno = 0.2, # Heritability of yield in pheno
  geno_cost = 5, # cost to genotype
  geno_cycle_length = 0.2, # number of years prediction selection cycle takes
  pred_accuracy = 0.85, # Prediction accuracy that is synonymous with genetic correlation
  h_predTrait = 1
)
