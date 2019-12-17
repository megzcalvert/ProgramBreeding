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
pheno_cycle_length <- 7
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
pheno_cycle_length <- 7
## number of years phenotypic selection cycle takes
h_pheno <- 0.8
## Heritability of yield in pheno

geno_cost <- 5
## cost to genotype
geno_cycle_length <- 5
## number of years prediction selection cycle takes
pred_accuracy <- 0.75
## Prediction accuracy that is synonymous with genetic correlation
h_predTrait <- 1
## assumes the genome in an inbred is 100% heritable and all additive

# Number of plots if only phenotyping
only_pheno <- floor(budget / (Line_dev_cost + (pheno_cost * pheno_rep)))
number_pheno <- as.integer(0:only_pheno)

# Number of plots if only predicting
only_geno <- floor(budget / (Line_dev_cost + geno_cost))
number_geno <- as.integer(0:only_geno)

# Creating dataframe of all combinations
df <- crossing(number_geno, number_pheno)

df <- df %>%
  mutate(correlated_response = 
           (((ayn_number / number_geno) * pred_accuracy * sqrt(h_predTrait)) / 
              geno_cycle_length),
         direct_response = 
           (((ayn_number / number_pheno) * sqrt(h_pheno)) / pheno_cycle_length),
         response_ratio = correlated_response / direct_response) %>% 
  filter(!is.infinite(response_ratio)) %>% 
  mutate(log_response = log(response_ratio))


gainPlot <- ggplot(data = df,
                   aes(x = number_pheno, 
                       y = number_geno,
                       fill = log_response)) +
  geom_raster() +
  scale_fill_gradient2(low = "#762a83", 
                       high = "#1b7837")
gainPlot
