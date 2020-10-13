rm(list = objects())
ls()

library(data.table)
library(tidyverse)
library(janitor)
library(broom)
library(msm)

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
      size = rel(2)
    ),
    complete = F
  )

theme_set(custom_theme)

getwd()

# useful infos **reproducible research**
sessionInfo()
# options(scipen = 999)

###############################################################################
###################### Response ratio by numbers ############################
###############################################################################

GainByResponse <- function(
                           budget = budget, # per breeding *cycle*
                           Line_dev_cost = Line_dev_cost, # the cost of generating an inbred line, applied to either form of evaluation
                           ayn_number = ayn_number, # number of AYN to be planted
                           pheno_cost_start = pheno_cost_start, # cost to phenotype in field
                           pheno_cost_end = pheno_cost_end, # cost to phenotype in field
                           pheno_cost_by = pheno_cost_by, # cost to phenotype in field
                           pheno_rep = pheno_rep, # number of phenotyping plots planted
                           #pheno_cycle_length = pheno_cycle_length, # number of years phenotypic selection cycle takes
                           geno_cost_start = geno_cost_start, # cost to genotype
                           geno_cost_end = geno_cost_end, # cost to genotype
                           geno_cost_by = geno_cost_by, # cost to genotype
                           h_pheno_min = h_pheno_min, # Heritability of yield in pheno (min)
                           h_pheno_max = h_pheno_max, # Heritability of yield in pheno (max)
                           h_pheno_by = h_pheno_by,
                           #geno_cycle_length = geno_cycle_length, # number of years prediction selection cycle takes
                           pred_accuracy_min = pred_accuracy_min, # Prediction accuracy that is synonymous with genetic correlation
                           pred_accuracy_max = pred_accuracy_max, # Prediction accuracy that is synonymous with genetic correlation
                           pred_accuracy_by = pred_accuracy_by,
                           h_predTrait = h_predTrait) {

  # Number of plots if only phenotyping
  only_pheno_low <- floor(budget /
    (Line_dev_cost + (pheno_cost_start * pheno_rep)))
  only_pheno_high <- floor(budget /
    (Line_dev_cost + (pheno_cost_end * pheno_rep)))
  only_pheno <- seq(
    from = only_pheno_high,
    to = only_pheno_low,
    by = pheno_cost_by
  )

  # Number of plots if only predicting
  only_geno_low <- floor(budget / (Line_dev_cost + geno_cost_start))
  only_geno_high <- floor(budget / (Line_dev_cost + geno_cost_end))
  only_geno <- seq(
    from = only_geno_high,
    to = only_geno_low,
    by = geno_cost_by
  )

  all_costs <- crossing(only_geno, only_pheno) %>%
    filter(only_pheno > (ayn_number * 1.5) )
  
  print(paste("Minimum phenotyped: ", min(all_costs$only_pheno)))
  print(paste("Maximum phenotyped: ", max(all_costs$only_pheno)))

  costs <- nrow(all_costs)

  dT <- tibble::tibble(
    only_geno_used = NA,
    only_pheno_used = NA,
    ratio_numbers = NA,
    h_pheno = NA,
    pred_accuracy = NA,
    correlated_response = NA,
    direct_response = NA,
    response_ratio = NA,
    correlated_intensity = NA,
    direct_intensity = NA
  )

  for (q in seq(1, costs, by = 1)) {
    only_geno <- as.numeric(all_costs[q, 1])
    only_pheno <- as.numeric(all_costs[q, 2])

    print(paste0("dT = ", q, " / ", costs))
    tVar <- 1 ## total variance

    prop_sel_pheno <- ayn_number / only_pheno
    prop_sel_geno <- ayn_number / only_geno

    h_pheno <- seq(
      from = h_pheno_min, to = h_pheno_max,
      by = h_pheno_by
    )
    pred_accuracy <- seq(
      from = pred_accuracy_min, to = pred_accuracy_max,
      by = pred_accuracy_by
    )

    df <- crossing(h_pheno, pred_accuracy) %>%
      mutate(
        only_geno_used = NA,
        only_pheno_used = NA,
        ratio_numbers = NA,
        correlated_response = NA,
        direct_response = NA,
        response_ratio = NA,
        correlated_intensity = as.numeric(NA),
        direct_intensity = as.numeric(NA)
      )

    ID <- nrow(df)

    for (i in seq(1, ID, by = 1)) {
      h_pheno <- as.numeric(df[i, 1])

      # Only Pheno
      gVar_pheno <- tVar * h_pheno # genetic variance
      eVar_pheno <- tVar * (1 - h_pheno) # enviromental variance

      geno_pheno <- rnorm(only_pheno, mean = 0, sd = sqrt(gVar_pheno)) # breeding value
      env_pheno <- rnorm(only_pheno, mean = 0, sd = sqrt(eVar_pheno)) # enviromental deviation

      pheno_pheno <- geno_pheno + env_pheno

      # Only geno
      gVar_geno <- tVar * h_pheno # genetic variance
      eVar_geno <- tVar * (1 - h_pheno) # enviromental variance

      geno_geno <- rnorm(only_geno, mean = 0, sd = sqrt(gVar_geno)) # breeding value
      env_geno <- rnorm(only_geno, mean = 0, sd = sqrt(eVar_geno)) # enviromental deviation

      pheno_geno <- geno_geno + env_geno

      # direct selection
      threshold_pheno <- sort(pheno_pheno, decreasing = TRUE)[ayn_number]

      rtn_pheno <- rtnorm(
        n = only_pheno,
        mean = mean(pheno_pheno),
        sd = sd(pheno_pheno),
        lower = threshold_pheno
      )
      i_pheno <- mean(rtn_pheno)

      df[i, "direct_intensity"] <- i_pheno

      # Indirect selection
      threshold_geno <- sort(pheno_geno, decreasing = TRUE)[ayn_number]

      rtn_geno <- rtnorm(
        n = only_geno,
        mean = mean(pheno_geno),
        sd = sd(pheno_geno),
        lower = threshold_geno
      )
      i_geno <- mean(rtn_geno)

      df[i, "correlated_intensity"] <- i_geno
    }

    df <- df %>%
      mutate(
        only_geno_used = only_geno,
        only_pheno_used = only_pheno,
        ratio_numbers = round(only_geno_used / only_pheno_used, 2),
        correlated_intensity = correlated_intensity,
        direct_intensity = direct_intensity,
        correlated_response =
          (correlated_intensity * pred_accuracy * sqrt(h_predTrait)), 
        direct_response =
          (direct_intensity * sqrt(h_pheno)),
        response_ratio = correlated_response / direct_response
      ) %>%
      filter(
        !is.infinite(response_ratio)
      ) %>%
      filter(
        !is.na(response_ratio)
      ) %>%
      filter(
        direct_intensity >= 0
      )

    dT <- bind_rows(dT, df)
  }

  return(dT)
}

try1 <- GainByResponse(
  budget = 100000, # per breeding *cycle*
  Line_dev_cost = 40, # the cost of generating an inbred line, applied to either form of evaluation
  ayn_number = 100, # number of AYN to be planted
  pheno_cost_start = 35, # cost to phenotype in field
  pheno_cost_end = 70, # cost to phenotype in field
  pheno_cost_by = 5, # cost to phenotype in field
  pheno_rep = 6, # number of phenotyping plots planted
  #pheno_cycle_length = 1, # number of years phenotypic selection cycle takes
  geno_cost_start = 5, # cost to genotype
  geno_cost_end = 10, # cost to genotype
  geno_cost_by = 1, # cost to genotype
  h_pheno_min = 0, # Heritability of yield in pheno (min)
  h_pheno_max = 1, # Heritability of yield in pheno (max)
  h_pheno_by = 0.1,
  #geno_cycle_length = 0.25, # number of years prediction selection cycle takes
  pred_accuracy_min = 0, # Prediction accuracy that is synonymous with genetic correlation
  pred_accuracy_max = 1, # Prediction accuracy that is synonymous with genetic correlation
  pred_accuracy_by = 0.1,
  h_predTrait = 0.95
)

p <- try1 %>% 
  drop_na(pred_accuracy) %>% 
  ggplot(aes(
  x = ratio_numbers,
  y = response_ratio,
  group = factor(pred_accuracy),
  colour = factor(pred_accuracy)
)) +
  #geom_point(alpha = 0.05, fill = NA) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 1, linetype = 2) +
 # scale_color_manual(
 # values = c(
 #   '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
 #   '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
 #   '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
 #   '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
 # ),
 # name = "Parameters"
 # ) +
  facet_wrap(~h_pheno, ncol = 5, scales = "free_y") +
  coord_cartesian(ylim = c(0,8)) +
  guides(colour = guide_legend(title = "Prediction Accuracy")) +
  labs(
   # title = "The response ratio when using exclusively prediction or phenotyping",
    x = "Phenotyped:Predicted",
    y = "CR / R"
  ) +
theme(axis.text = element_text(size = rel(2)))

p

ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/BreedingProgram/Figures/Figure1.png",
       width = 60,
       height = 30,
       units = "cm",
       dpi = 320)


groupings <- try1_1 %>%
  select(h_pheno, pred_accuracy, groupID) %>%
  distinct()
groupings

ggsave("./Figures/Justification/Response_ratio_Trendlines_predAccuracy.png",
  width = 50, height = 30,
  units = "cm", dpi = 720
)

try1_2<- try1 %>% 
  drop_na() %>%
  filter(h_pheno > 0.88,
         h_pheno < 0.92) %>% 
  group_by(h_pheno, pred_accuracy) %>%
  mutate(groupID = group_indices()) %>%
  # filter(h_pheno > pred_accuracy) %>%
  filter(
    pred_accuracy > 0,
    h_pheno > 0,
    response_ratio < 1000
  )

q <- ggplot(data = try1_2, aes(
  x = ratio_numbers,
  y = response_ratio,
  group = factor(pred_accuracy),
  colour = factor(pred_accuracy)
)) +
  #geom_point(alpha = 0.5) +
  geom_smooth() +
  geom_hline(yintercept = 1, linetype = 2) +
  #scale_color_manual(
  # values = c(
  #   '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
  #   '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
  #   '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
  #   '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
  # ),
  # name = "Parameters"
  # ) +
  #facet_wrap(~h_pheno, ncol = 5, scales = "free_y") +
  labs(
    title = "The response ratio when using exclusively prediction or phenotyping",
    subtitle = paste("Narrow-sense heritability = ", unique(try1_2$h_pheno)),
    x = "phenotyped:predicted",
    y = "correlated_response / direct_response"
  ) 
q
