rm(list = objects())
ls()

library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(ggridges)
library(broom)

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
    legend.key.size = unit(1.5, "lines"),
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
      size = rel(2),
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

set.seed(1642)

#### Load data ####

accuracy_all <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryldA.txt"
)
accuracy_all <- accuracy_all %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "all",
    Exp_trial = "all"
  )

accuracy_all %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "All GRYLD",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/BreedingProgram/Figures/GP_accuracyCV_allGryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_all %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "All years GRYLD",
    x = "Correlation"
  )
ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/BreedingProgram/Figures/GP_accuracyCV_allGryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## 16
accuracy_16 <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_16.txt"
)
accuracy_16 <- accuracy_16 %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "16",
    Exp_trial = "all"
  )

accuracy_16 %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "2015-2016 GRYLD",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_16Gryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_16 %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "2015-2016 GRYLD",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_16Gryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## 16_pyn
accuracy_16_pyn <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_16_pyn.txt"
)
accuracy_16_pyn <- accuracy_16_pyn %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "16",
    Exp_trial = "PYN"
  )

accuracy_16_pyn %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "2015-2016 GRYLD PYN",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_16pynGryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_16_pyn %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "2015-2016 GRYLD PYN",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_16pynGryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## 17
accuracy_17 <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_17.txt"
)
accuracy_17 <- accuracy_17 %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "17",
    Exp_trial = "all"
  )

accuracy_17 %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "2016-2017 GRYLD",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_17Gryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_17 %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "2016-2017 GRYLD",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_17Gryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## 17_pyn
accuracy_17_pyn <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_17_pyn.txt"
)
accuracy_17_pyn <- accuracy_17_pyn %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "17",
    Exp_trial = "PYN"
  )

accuracy_17_pyn %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "2016-2017 GRYLD PYN",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_17pynGryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_17_pyn %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "2016-2017 GRYLD PYN",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_17pynGryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## 18
accuracy_18 <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_18.txt"
)
accuracy_18 <- accuracy_18 %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "18",
    Exp_trial = "all"
  )

accuracy_18 %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "2017-2018 GRYLD",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_18Gryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_18 %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "2017-2018 GRYLD",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_18Gryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## 18_pyn
accuracy_18_pyn <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_18_pyn.txt"
)
accuracy_18_pyn <- accuracy_18_pyn %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "18",
    Exp_trial = "PYN"
  )

accuracy_18_pyn %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "2017-2018 GRYLD PYN",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_18pynGryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_18_pyn %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "2017-2018 GRYLD PYN",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_18pynGryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## 19
accuracy_19 <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_19.txt"
)
accuracy_19 <- accuracy_19 %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "19",
    Exp_trial = "all"
  )

accuracy_19 %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "2018-2019 GRYLD",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_19Gryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_19 %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "2018-2019 GRYLD",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_19Gryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## 19_pyn
accuracy_19_pyn <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_19_pyn.txt"
)
accuracy_19_pyn <- accuracy_19_pyn %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "19",
    Exp_trial = "PYN"
  )

accuracy_19_pyn %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "2018-2019 GRYLD PYN",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_19pynGryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_19_pyn %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "2018-2019 GRYLD PYN",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_19pynGryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## ex16
accuracy_ex16 <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_ex16.txt"
)
accuracy_ex16 <- accuracy_ex16 %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "excluding_16",
    Exp_trial = "all"
  )

accuracy_ex16 %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "All years excluding 2015-2016 GRYLD",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex16Gryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_ex16 %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "All years excluding 2015-2016 GRYLD",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex16Gryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## ex16_pyn
accuracy_ex16_pyn <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_ex16_pyn.txt"
)
accuracy_ex16_pyn <- accuracy_ex16_pyn %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "excluding_16",
    Exp_trial = "PYN"
  )

accuracy_ex16_pyn %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "All years excluding 2015-2016 GRYLD PYN",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex16pynGryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_ex16_pyn %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "All years excluding 2015-2016 GRYLD PYN",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex16pynGryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## ex17
accuracy_ex17 <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_ex17.txt"
)
accuracy_ex17 <- accuracy_ex17 %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "excluding_17",
    Exp_trial = "all"
  )

accuracy_ex17 %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "All years excluding 2016-2017 GRYLD",
    x = "Accuracy",
    y = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex17Gryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_ex17 %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "All years excluding 2016-2017 GRYLD",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex17Gryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## 17 pyn

accuracy_ex17_pyn <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_ex17_pyn.txt"
)
accuracy_ex17_pyn <- accuracy_ex17_pyn %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "excluding_17",
    Exp_trial = "PYN"
  )

accuracy_ex17_pyn %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "All years excluding 2016-2017 GRYLD PYN",
    x = "Accuracy"
  )
ggsave("./Figures/GP_accuracyCV_ex17pynGryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_ex17_pyn %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "All years excluding 2016-2017 GRYLD PYN",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex17pynGryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## ex18
accuracy_ex18 <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_ex18.txt"
)
accuracy_ex18 <- accuracy_ex18 %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "excluding_18",
    Exp_trial = "all"
  )

accuracy_ex18 %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "All years excluding 2017-2018 GRYLD",
    x = "Accuracy"
  )
ggsave("./Figures/GP_accuracyCV_ex18Gryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_ex18 %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "All years excluding 2017-2018 GRYLD",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex18Gryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## ex18_pyn
accuracy_ex18_pyn <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_ex18_pyn.txt"
)
accuracy_ex18_pyn <- accuracy_ex18_pyn %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "excluding_18",
    Exp_trial = "PYN"
  )

accuracy_ex18_pyn %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "All years excluding 2017-2018 GRYLD PYN",
    x = "Accuracy"
  )
ggsave("./Figures/GP_accuracyCV_ex18pynGryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_ex18_pyn %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "All years excluding 2017-2018 GRYLD PYN",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex18pynGryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## ex19
accuracy_ex19 <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_ex19.txt"
)
accuracy_ex19 <- accuracy_ex19 %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "excluding_19",
    Exp_trial = "all"
  )

accuracy_ex19 %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "All years excluding 2018-2019 GRYLD",
    x = "Accuracy"
  )
ggsave("./Figures/GP_accuracyCV_ex19Gryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_ex19 %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "All years excluding 2018-2019 GRYLD",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex19Gryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

## ex19_pyn
accuracy_ex19_pyn <- fread(
  "./Results/GenomicSelection_rrBlup_bglr_80_100_accuracy_gryld_ex19_pyn.txt"
)
accuracy_ex19_pyn <- accuracy_ex19_pyn %>%
  pivot_longer(
    cols = everything(),
    values_to = "phenotype_value",
    names_to = "trial"
  ) %>%
  mutate(
    year = "excluding_19",
    Exp_trial = "PYN"
  )

accuracy_ex19_pyn %>%
  ggplot(aes(x = trial, y = phenotype_value)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  labs(
    title = "Distribution of accuracy of genomic prediction",
    subtitle = "All years excluding 2018-2019 GRYLD PYN",
    x = "Accuracy"
  )
ggsave("./Figures/GP_accuracyCV_ex19pynGryld_boxplot.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

accuracy_ex19_pyn %>%
  ggplot(aes(x = phenotype_value, y = trial)) +
  geom_density_ridges(
    quantile_lines = TRUE,
    scale = 0.9,
    fill = NA
  ) +
  scale_y_discrete(
    breaks = c(
      "rrBlup_modelTraining", "rrBlup_modelTesting",
      "rkhs_modelTraining", "rkhs_modelTesting",
      "bayesC_modelTraining", "bayesC_modelTesting"
    ),
    labels = c(
      "rrBlup_Training", "rrBlup_Testing",
      "rkhs_Training", "rkhs_Testing",
      "bayesC_Training", "bayesC_Testing"
    )
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Cross validation Genomic Prediction",
    subtitle = "All years excluding 2018-2019 GRYLD PYN",
    x = "Correlation"
  )
ggsave("./Figures/GP_accuracyCV_ex19pynGryld_ridges.png",
  width = 30,
  height = 20,
  units = "cm",
  dpi = 320
)

all_accuracies <- bind_rows(
  accuracy_16, accuracy_17, accuracy_18,
  accuracy_19, accuracy_ex16, accuracy_ex17,
  accuracy_ex18, accuracy_ex19, accuracy_all
) %>%
  filter(str_detect(trial, "rrBlup")) %>%
  mutate(
    year = str_replace(year, "excluding", "ex"),
    year = as.factor(year)
  ) %>%
  separate(trial, into = c("rrBlup", "Population"))

all_accuracies %>%
  filter(Exp_trial == "all") %>%
  ggplot(aes(
    x = year,
    y = phenotype_value,
    colour = Population
  )) +
  geom_jitter(alpha = 0.25) +
  geom_boxplot(fill = NA) +
  scale_colour_manual(values = c("#542788", "#b35806")) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    y = "Correlation",
    x = "Season"
  )

ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/BreedingProgram/Figures/Figure5.png",
  height = 25,
  width = 40, units = "cm", dpi = 320
)

all_accuracies %>%
  filter(str_detect(Population, "Testing")) %>%
  group_by(year, Exp_trial, Population) %>%
  summarise(
    average = mean(phenotype_value),
    sd = sd(phenotype_value)
  ) %>%
  ungroup()

#### Blind predictions ####

PredictedValues_ex16toPred16_gryld_rrBlup <- fread(
  "./Results/PredictedValues_ex16_pyntoPred16_pyn_gryld_rrBlup.txt"
)
PredictedValues_ex16_pyntoPred16_pyn_gryld_rrBlup <- fread(
  "./Results/PredictedValues_ex16_pyntoPred16_pyn_gryld_rrBlup.txt"
)
PredictedValues_pred16_MP_MANH_GYP_from_all_else_rrBlup <- fread(
  "./Results/PredictedValues_pred16_MP_MANH_GYP_from_all_else_rrBlup.txt"
)
PredictedValues_ex17toPred17_gryld_rrBlup <- fread(
  "./Results/PredictedValues_ex17toPred17_gryld_rrBlup.txt"
)
PredictedValues_ex17_pyntoPred17_pyn_gryld_rrBlup <- fread(
  "./Results/PredictedValues_ex17_pyntoPred17_pyn_gryld_rrBlup.txt"
)
PredictedValues_pred17_MP_MANH_BEL_from_all_else_rrBlup <- fread(
  "./Results/PredictedValues_pred17_MP_MANH_BEL_from_all_else_rrBlup.txt"
)
PredictedValues_ex18toPred18_gryld_rrBlup <- fread(
  "./Results/PredictedValues_ex18toPred18_gryld_rrBlup.txt"
)
PredictedValues_ex18_pyntoPred18_pyn_gryld_rrBlup <- fread(
  "./Results/PredictedValues_ex18_pyntoPred18_pyn_gryld_rrBlup.txt"
)
PredictedValues_pred18_MP_MANH_BEL_GYP_from_all_else_rrBlup <- fread(
  "./Results/PredictedValues_pred18_MP_MANH_BEL_GYP_from_all_else_rrBlup.txt"
)
PredictedValues_ex19toPred19_gryld_rrBlup <- fread(
  "./Results/PredictedValues_ex19toPred19_gryld_rrBlup.txt"
)
PredictedValues_ex19_pyntoPred19_pyn_gryld_rrBlup <- fread(
  "./Results/PredictedValues_ex19_pyntoPred19_pyn_gryld_rrBlup.txt"
)
PredictedValues_pred19_MANH_BEL_from_all_else_rrBlup <- fread(
  "./Results/PredictedValues_pred19_MANH_BEL_from_all_else_rrBlup.txt"
)

cor.test(
  PredictedValues_ex16_pyntoPred16_pyn_gryld_rrBlup$Obs_rank,
  PredictedValues_ex16_pyntoPred16_pyn_gryld_rrBlup$Pred_rank
)
cor.test(
  PredictedValues_ex16toPred16_gryld_rrBlup$Obs_rank,
  PredictedValues_ex16toPred16_gryld_rrBlup$Pred_rank
)
cor.test(
  PredictedValues_ex17_pyntoPred17_pyn_gryld_rrBlup$Obs_rank,
  PredictedValues_ex17_pyntoPred17_pyn_gryld_rrBlup$Pred_rank
)
cor.test(
  PredictedValues_ex17toPred17_gryld_rrBlup$Obs_rank,
  PredictedValues_ex17toPred17_gryld_rrBlup$Pred_rank
)
cor.test(
  PredictedValues_ex18_pyntoPred18_pyn_gryld_rrBlup$Obs_rank,
  PredictedValues_ex18_pyntoPred18_pyn_gryld_rrBlup$Pred_rank
)
cor.test(
  PredictedValues_ex18toPred18_gryld_rrBlup$Obs_rank,
  PredictedValues_ex18toPred18_gryld_rrBlup$Pred_rank
)

cor.test(PredictedValues_ex19toPred19_gryld_rrBlup$Obs_rank,
         PredictedValues_ex19toPred19_gryld_rrBlup$Pred_rank)

cor.test(PredictedValues_ex19_pyntoPred19_pyn_gryld_rrBlup$Obs_rank,
         PredictedValues_ex19_pyntoPred19_pyn_gryld_rrBlup$Pred_rank)

cor.test(
  PredictedValues_pred16_MP_MANH_GYP_from_all_else_rrBlup$Obs_rank,
  PredictedValues_pred16_MP_MANH_GYP_from_all_else_rrBlup$Pred_rank
)
cor.test(
  PredictedValues_pred17_MP_MANH_BEL_from_all_else_rrBlup$Obs_rank,
  PredictedValues_pred17_MP_MANH_BEL_from_all_else_rrBlup$Pred_rank
)
cor.test(
  PredictedValues_pred18_MP_MANH_BEL_GYP_from_all_else_rrBlup$Obs_rank,
  PredictedValues_pred18_MP_MANH_BEL_GYP_from_all_else_rrBlup$Pred_rank
)
cor.test(
  PredictedValues_pred19_MANH_BEL_from_all_else_rrBlup$Obs_rank,
  PredictedValues_pred19_MANH_BEL_from_all_else_rrBlup$Pred_rank
)

