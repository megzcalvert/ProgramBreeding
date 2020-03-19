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

##### Read in data ####
geno <- fread("./GenoDatabase/SelectedGeno_Numeric.txt")
geno[1:10, 1:10]
hapgeno <- geno[, 4:ncol(geno)]
hapgeno[hapgeno == 0] <- -1
hapgeno[hapgeno == 0.5] <- 0

geno <- bind_cols(geno[, 1:3], hapgeno)
geno[1:10, 1:10]

pheno <- fread("./PhenoDatabase/Blups_allYrs_allVarieties.txt")
pheno[1:10, ]

genoVarieties <- as.data.frame(colnames(geno[, 4:ncol(geno)]))
names(genoVarieties) <- "Variety"
pheno <- pheno %>%
  semi_join(genoVarieties)

snpMatrix <- t(geno[, 4:ncol(geno)])

relMat <- A.mat(snpMatrix)

gwas_Gryld <- GWAS(pheno = pheno, geno = geno, K = relMat, n.PC = 4, 
                   P3D = TRUE, plot = TRUE)

#### Plot ####

don_gwas <- gwas_Gryld %>%
  # Compute chromosome size
  group_by(chrom) %>%
  summarise(chr_len = max(pos)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwas_Gryld, ., by = c("chrom" = "chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum = pos + tot)

axisdf <- don_gwas %>%
  group_by(chrom) %>%
  dplyr::summarize(
    center = (max(BPcum) + min(BPcum)) / 2,
    maximum = max(BPcum)
  )

ggplot(don_gwas, aes(x = BPcum, colour = as.factor(don_gwas$chrom))) +
  # Show all points
  geom_point(aes(y = don_gwas$GRYLD),
    alpha = 0.5, size = 1
  ) +
  scale_color_manual(values = rep(c("#2ca25f", "#8856a7", "#43a2ca"), 21)) +
  # Significance Threshold
  geom_hline(yintercept = -log10(0.05 / nrow(don_gwas)), linetype = 2) +
  # geom_vline(xintercept = )
  # custom X axis:
  scale_x_continuous(label = axisdf$chrom, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0.05)) + # remove space between plot area and x axis
  theme(legend.position = "none") +
  labs(
    title = "GWAS results GRYLD all years",
    subtitle = "Bonferroni corection alpha = 0.05",
    x = "Chromosome",
    y = "-log10(P)"
  )


heterozygosityTest <- geno %>%
  select(-chrom, -pos) %>%
  pivot_longer(cols = -snp, names_to = "Variety") %>%
  group_by(snp, value) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = value, values_from = count) %>%
  mutate(Total = `-1` + `0` + `1`,
         homoRecessive = `-1` / Total,
         heterozygous = `0` / Total,
         homoDominant = `1` / Total,
         Check = homoRecessive + heterozygous + homoDominant)

heterozygosityTest %>% 
  ggplot(aes(x = heterozygous)) +
  geom_histogram() +
  labs(title = "Distribution of SNP states in KSU breeding program",
       #subtitle = "-1 = homozygous recessive, 0 = heterozygous, 1 = homozygous dominant",
       x = "Heterozygous Frequency")

