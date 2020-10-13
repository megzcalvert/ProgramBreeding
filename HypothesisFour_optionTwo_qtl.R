rm(list = objects())
ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(broom)
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
      size = rel(1.5)
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

# useful infos **reproducible research**
sessionInfo()
##############################################################################
#### Read in data ####
geno <- fread("./GenoDatabase/SelectedGeno_Numeric.txt")
geno[1:10, 1:10]
hapgeno <- geno[, 4:ncol(geno)]
hapgeno[hapgeno == 0] <- -1
hapgeno[hapgeno == 0.5] <- 0

geno <- bind_cols(geno[, 1:3], hapgeno)
geno[1:10, 1:10]
rm(hapgeno)

genoVarieties <- as.data.frame(colnames(geno[, 4:ncol(geno)]))
names(genoVarieties) <- "Variety"

blups17_ayn <- fread("./Results/Blups_ayn_17.txt")
blups17_pyn <- fread("./Results/Blups_pyn_17.txt")

blups18_ayn <- fread("./Results/Blups_ayn_18.txt")
blups18_pyn <- fread("./Results/Blups_pyn_18.txt")

blups19_ayn <- fread("./Results/Blups_ayn_19.txt")
blups19_pyn <- fread("./Results/Blups_pyn_19.txt")

blups17_ayn <- blups17_ayn %>%
  semi_join(genoVarieties, by = "Variety") %>% 
  arrange(Variety)

blups17_pyn <- blups17_pyn %>%
  semi_join(genoVarieties, by = "Variety") %>% 
  arrange(Variety)

blups18_ayn <- blups18_ayn %>%
  semi_join(genoVarieties, by = "Variety") %>% 
  arrange(Variety)

blups18_pyn <- blups18_pyn %>%
  semi_join(genoVarieties, by = "Variety") %>% 
  arrange(Variety)

blups19_ayn <- blups19_ayn %>%
  semi_join(genoVarieties, by = "Variety") %>% 
  arrange(Variety)

blups19_pyn <- blups19_pyn %>%
  semi_join(genoVarieties, by = "Variety") %>% 
  arrange(Variety)

summaryBlup17_ayn <- blups17_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
    sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

summaryBlup17_pyn <- blups17_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

summaryBlup18_ayn <- blups18_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

summaryBlup18_pyn <- blups18_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

summaryBlup19_ayn <- blups19_ayn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

summaryBlup19_pyn <- blups19_pyn %>%
  pivot_longer(
    cols = -Variety,
    names_to = "trait_id",
    values_to = "phenotype_value"
  ) %>%
  separate(trait_id, c("trait_id", "phenotype_date", "location"),
           sep = "_"
  ) %>% 
  group_by(location, phenotype_date, trait_id) %>% 
  summarise(n = n(),
            minimum = min(phenotype_value),
            average = mean(phenotype_value),
            sd = sd(phenotype_value),
            maximum = max(phenotype_value))

### Kinship matrix

snpMatrix <- as.data.frame(t(geno[, 4:ncol(geno)]))
snpMatrix <- setDT(snpMatrix, keep.rownames = TRUE)

snpMat_17_ayn<- snpMatrix%>% 
  semi_join(blups17_ayn, by = c("rn" = "Variety")) %>%
  arrange(rn) %>% 
  column_to_rownames(var = "rn")

geno_17_ayn<- as.data.frame(t(snpMat_17_ayn)) %>% 
  bind_cols(geno[,1:3]) %>% 
  select(snp, chrom, pos, everything())

relMat <- A.mat(snpMat_17_ayn)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMat_17_ayn) == blups17_ayn$Variety), sep = ' '))
print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(relMat) == blups17_ayn$Variety), sep = ' '))

ayn17<- GWAS(pheno = blups17_ayn,
             geno = geno_17_ayn,
             n.PC = 3,
             K = relMat,
             plot = FALSE,
             P3D = FALSE)

snpMat_17_pyn<- snpMatrix%>% 
  semi_join(blups17_pyn, by = c("rn" = "Variety")) %>% 
  column_to_rownames(var = "rn")

geno_17_pyn<- as.data.frame(t(snpMat_17_pyn)) %>% 
  bind_cols(geno[,1:3]) %>% 
  select(snp, chrom, pos, everything())

relMat <- A.mat(snpMat_17_pyn)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMat_17_pyn) == blups17_pyn$Variety), sep = ' '))
print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(relMat) == blups17_pyn$Variety), sep = ' '))

pyn17<- GWAS(pheno = blups17_pyn,
             geno = geno_17_pyn,
             n.PC = 3,
             K = relMat,
             plot = FALSE,
             P3D = FALSE)

snpMat_18_ayn<- snpMatrix%>% 
  semi_join(blups18_ayn, by = c("rn" = "Variety")) %>% 
  column_to_rownames(var = "rn")

geno_18_ayn<- as.data.frame(t(snpMat_18_ayn)) %>% 
  bind_cols(geno[,1:3]) %>% 
  select(snp, chrom, pos, everything())

relMat <- A.mat(snpMat_18_ayn)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMat_18_ayn) == blups18_ayn$Variety), sep = ' '))
print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(relMat) == blups18_ayn$Variety), sep = ' '))

blups18_ayn<- blups18_ayn %>% 
  select(-ID)

ayn18<- GWAS(pheno = blups18_ayn,
             geno = geno_18_ayn,
             n.PC = 3,
             K = relMat,
             plot = FALSE,
             P3D = FALSE)

snpMat_18_pyn<- snpMatrix%>% 
  semi_join(blups18_pyn, by = c("rn" = "Variety")) %>% 
  column_to_rownames(var = "rn")

geno_18_pyn<- as.data.frame(t(snpMat_18_pyn)) %>% 
  bind_cols(geno[,1:3]) %>% 
  select(snp, chrom, pos, everything())

relMat <- A.mat(snpMat_18_pyn)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMat_18_pyn) == blups18_pyn$Variety), sep = ' '))
print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(relMat) == blups18_pyn$Variety), sep = ' '))

pyn18<- GWAS(pheno = blups18_pyn,
             geno = geno_18_pyn,
             n.PC = 3,
             K = relMat,
             plot = FALSE,
             P3D = FALSE)

snpMat_19_ayn<- snpMatrix%>% 
  semi_join(blups19_ayn, by = c("rn" = "Variety")) %>% 
  column_to_rownames(var = "rn")

geno_19_ayn<- as.data.frame(t(snpMat_19_ayn)) %>% 
  bind_cols(geno[,1:3]) %>% 
  select(snp, chrom, pos, everything())

relMat <- A.mat(snpMat_19_ayn)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMat_19_ayn) == blups19_ayn$Variety), sep = ' '))
print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(relMat) == blups19_ayn$Variety), sep = ' '))

ayn19<- GWAS(pheno = blups19_ayn,
             geno = geno_19_ayn,
             n.PC = 3,
             K = relMat,
             plot = FALSE,
             P3D = FALSE)

snpMat_19_pyn<- snpMatrix%>% 
  semi_join(blups19_pyn, by = c("rn" = "Variety")) %>% 
  column_to_rownames(var = "rn")

geno_19_pyn<- as.data.frame(t(snpMat_19_pyn)) %>% 
  bind_cols(geno[,1:3]) %>% 
  select(snp, chrom, pos, everything())

relMat <- A.mat(snpMat_19_pyn)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMat_19_pyn) == blups19_pyn$Variety), sep = ' '))
print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(relMat) == blups19_pyn$Variety), sep = ' '))

pyn19<- GWAS(pheno = blups19_pyn,
             geno = geno_19_pyn,
             n.PC = 3,
             K = relMat,
             plot = FALSE,
             P3D = TRUE)

##### Making this a function to plot all ####

myManhattans<- function(dat, traits, saveFileFigures, identifier,
                        colourOne, colourTwo, colourThree, ...){
  
  dat <- dat %>% 
    # Compute chromosome size
    group_by(chrom) %>% 
    summarise(chr_len=max(pos), .groups = "drop_last") %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    tidylog::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(dat, ., by=c("chrom"="chrom")) %>%
    # Add a cumulative position of each SNP
    arrange(chrom, pos) %>%
    mutate( BPcum=pos+tot) %>% 
    select(1:3,tot, BPcum, everything())
  
  axisdf <- dat %>% 
    group_by(chrom) %>% 
    dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2,
                     maximum = max(BPcum),
                     .groups = "drop_last")
  plotList<- list()
  
  for (i in traits) {
    thisPlot<- ggplot(data = dat, 
                      mapping = aes(x=BPcum, 
                                    colour = as.factor(chrom))) +
      # Show all points
      geom_point(aes_string(y = paste(i)),
                 alpha = 0.5, size = 1) +
      scale_color_manual(values = rep(c(colourOne, colourTwo, colourThree), 
                                      21 ),
                         name = "Chromosome") +
      # Significance Threshold
      geom_hline(yintercept = -log10(0.05/nrow(dat)), linetype = 2) +
      # custom X axis:
      scale_x_continuous( label = c("1A","1B","1D",
                                    "2A","2B","2D",
                                    "3A","3B","3D",
                                    "4A","4B","4D",
                                    "5A","5B","5D",
                                    "6A","6B","6D",
                                    "7A","7B","7D"), 
                          breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0.05) ) +
      theme(legend.position = "none") +
      labs(title = paste0("GWAS results ", identifier, " ", i),
           subtitle = "Bonferroni Threshold alpha = 0.05",
           x = "Chromosome",
           y = "-log10(P)") 
    plotList[[i]] <- thisPlot
    ggpubr::ggexport(thisPlot, filename = paste0(saveFileFigures,
                                                 identifier,"_",i,".png"),
                     width = 1000, height = 550)
  }
  return(plotList)
}

traits<- colnames(ayn17[,4:ncol(ayn17)])
ayn17<- ayn17 %>% 
  filter(chrom != "UN")
manhattana17<- myManhattans(dat = ayn17, traits = traits, 
                           saveFileFigures = "./Figures/GWAS/",
                           identifier = "ayn17",
                           colourOne = "#2ca25f", colourTwo = "#43a2ca",
                           colourThree = "#8856a7")

traits<- colnames(pyn17[,4:ncol(pyn17)])
pyn17<- pyn17 %>% 
  filter(chrom != "UN")
manhattanp17<- myManhattans(dat = pyn17, traits = traits, 
                           saveFileFigures = "./Figures/GWAS/",
                           identifier = "pyn17",
                           colourOne = "#2ca25f", colourTwo = "#43a2ca",
                           colourThree = "#8856a7")

traits<- colnames(ayn18[,4:ncol(ayn18)])
ayn18<- ayn18 %>% 
  filter(chrom != "UN")
manhattana18<- myManhattans(dat = ayn18, traits = traits, 
                            saveFileFigures = "./Figures/GWAS/",
                            identifier = "ayn18",
                            colourOne = "#2ca25f", colourTwo = "#43a2ca",
                            colourThree = "#8856a7")

traits<- colnames(pyn18[,4:ncol(pyn18)])
pyn18<- pyn18 %>% 
  filter(chrom != "UN")
manhattanp18<- myManhattans(dat = pyn18, traits = traits, 
                            saveFileFigures = "./Figures/GWAS/",
                            identifier = "pyn18",
                            colourOne = "#2ca25f", colourTwo = "#43a2ca",
                            colourThree = "#8856a7")

traits<- colnames(ayn19[,4:ncol(ayn19)])
ayn19<- ayn19 %>% 
  filter(chrom != "UN")
manhattana19<- myManhattans(dat = ayn19, traits = traits, 
                            saveFileFigures = "./Figures/GWAS/",
                            identifier = "ayn19",
                            colourOne = "#2ca25f", colourTwo = "#43a2ca",
                            colourThree = "#8856a7")

traits<- colnames(pyn19[,4:ncol(pyn19)])
pyn19<- pyn19 %>% 
  filter(chrom != "UN")
manhattanp19<- myManhattans(dat = pyn19, traits = traits, 
                            saveFileFigures = "./Figures/GWAS/",
                            identifier = "pyn19",
                            colourOne = "#2ca25f", colourTwo = "#43a2ca",
                            colourThree = "#8856a7")
###############################################################################
#### Additional diagnostics

geno<- t(geno)
geno<- setDT(as.data.frame(geno), keep.rownames = TRUE) %>% 
  row_to_names(1)

gndvi20180517RN_snp<- ayn18 %>% 
  select(snp, chrom, pos, GNDVI_2018.05.17_RN) %>% 
  filter(GNDVI_2018.05.17_RN == max(GNDVI_2018.05.17_RN))
gndvi20180517RN_snp<- geno %>% 
  select(snp, gndvi20180517RN_snp[1,1])

gndvi20180517RN<- blups18_ayn %>% 
  select(Variety,`GNDVI_2018-05-17_RN`, `GRYLD_2018-06-28_RN`) %>%
  left_join(gndvi20180517RN_snp, by = c("Variety" = "snp"))

p1<- gndvi20180517RN %>% 
  ggplot(aes(x = S7A_111223394, y = `GNDVI_2018-05-17_RN`)) +
  geom_boxplot() 
p1
p2<- gndvi20180517RN %>% 
  ggplot(aes(x = S7A_111223394, y =  `GRYLD_2018-06-28_RN`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p1,p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 22)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_gndvi20180517RN.png",
                 width = 1000,
                 height = 750)

ndre20180517RN_snp<- ayn18 %>% 
  select(snp, chrom, pos, NDRE_2018.05.17_RN) %>% 
  filter(NDRE_2018.05.17_RN == max(NDRE_2018.05.17_RN))
ndre20180517RN_snp<- geno %>% 
  select(snp, ndre20180517RN_snp[1,1])

ndre20180517RN<- blups18_ayn %>% 
  select(Variety,`NDRE_2018-05-17_RN`, `GRYLD_2018-06-28_RN`) %>%
  left_join(ndre20180517RN_snp, by = c("Variety" = "snp"))

p1<- ndre20180517RN %>% 
  ggplot(aes(x = S7A_111223394, y = `NDRE_2018-05-17_RN`)) +
  geom_boxplot() 
p1
p2<- ndre20180517RN %>% 
  ggplot(aes(x = S7A_111223394, y =  `GRYLD_2018-06-28_RN`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p1,p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 22)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_ndre20180517RN.png",
                 width = 1000,
                 height = 750)

ndvi20180517RN_snp<- ayn18 %>% 
  select(snp, chrom, pos, NDVI_2018.05.17_RN) %>% 
  filter(NDVI_2018.05.17_RN == max(NDVI_2018.05.17_RN))
ndvi20180517RN_snp<- geno %>% 
  select(snp, ndvi20180517RN_snp[1,1])

ndvi20180517RN<- blups18_ayn %>% 
  select(Variety,`NDVI_2018-05-17_RN`, `GRYLD_2018-06-28_RN`) %>%
  left_join(ndvi20180517RN_snp, by = c("Variety" = "snp"))

p1<- ndvi20180517RN %>% 
  ggplot(aes(x = S7A_111223394, y = `NDVI_2018-05-17_RN`)) +
  geom_boxplot() 
p1
p2<- ndvi20180517RN %>% 
  ggplot(aes(x = S7A_111223394, y =  `GRYLD_2018-06-28_RN`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p1,p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 22)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_ndvi20180517RN.png",
                 width = 1000,
                 height = 750)

ndvi20180606MP_snp<- ayn18 %>% 
  select(snp, chrom, pos, NDVI_2018.06.06_MP) %>% 
  filter(NDVI_2018.06.06_MP == max(NDVI_2018.06.06_MP))
ndvi20180606MP_snp<- geno %>% 
  select(snp, ndvi20180606MP_snp[1,1])

ndvi20180606MP<- blups18_ayn %>% 
  select(Variety,`NDVI_2018-06-06_MP`, `GRYLD_2018-06-29_MP`) %>%
  left_join(ndvi20180606MP_snp, by = c("Variety" = "snp"))

p1<- ndvi20180606MP %>% 
  ggplot(aes(x = S6B_719425028, y = `NDVI_2018-06-06_MP`)) +
  geom_boxplot() 
p1
p2<- ndvi20180606MP %>% 
  ggplot(aes(x = S6B_719425028, y =  `GRYLD_2018-06-29_MP`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p1,p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 22)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_ndvi20180606MP.png",
                 width = 1000,
                 height = 750)


ndvi20190530RN_snp<- ayn19 %>% 
  select(snp, chrom, pos, NDVI_2019.05.30_RN) %>% 
  filter(NDVI_2019.05.30_RN == max(NDVI_2019.05.30_RN))
ndvi20190530RN_snp<- geno %>% 
  select(snp, ndvi20190530RN_snp[1,1])

ndvi20190530RN<- blups19_ayn %>% 
  select(Variety,`NDVI_2019-05-30_RN`, `GRYLD_2019-07-01_RN`) %>%
  left_join(ndvi20190530RN_snp, by = c("Variety" = "snp"))

p1<- ndvi20190530RN %>% 
  ggplot(aes(x = S3D_27244923, y = `NDVI_2019-05-30_RN`)) +
  geom_boxplot() 
p1
p2<- ndvi20190530RN %>% 
  ggplot(aes(x = S3D_27244923, y =  `GRYLD_2019-07-01_RN`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p1,p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 22)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_ndvi20190530RN.png",
                 width = 1000,
                 height = 750)

ndvi20170506MP_snp<- pyn17 %>% 
  select(snp, chrom, pos, NDVI_2017.05.06_MP) %>% 
  filter(NDVI_2017.05.06_MP == max(NDVI_2017.05.06_MP))
ndvi20170506MP_snp<- geno %>% 
  select(snp, ndvi20170506MP_snp[1,1])

ndvi20170506MP<- blups17_pyn %>% 
  select(Variety,`NDVI_2017-05-06_MP`, `GRYLD_2017-06-20_MP`) %>%
  left_join(ndvi20170506MP_snp, by = c("Variety" = "snp"))

p1<- ndvi20170506MP %>% 
  ggplot(aes(x = S5B_559689190, y = `NDVI_2017-05-06_MP`)) +
  geom_boxplot() 
p1
p2<- ndvi20170506MP %>% 
  ggplot(aes(x = S5B_559689190, y =  `GRYLD_2017-06-20_MP`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p1,p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 22)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_ndvi20170506MP.png",
                 width = 1000,
                 height = 750)


gndvi20190522RN_snp<- pyn19 %>% 
  select(snp, chrom, pos, GNDVI_2019.05.22_RN) %>% 
  filter(GNDVI_2019.05.22_RN == max(GNDVI_2019.05.22_RN))
gndvi20190522RN_snp<- geno %>% 
  select(snp, gndvi20190522RN_snp[1,1])

gndvi20190522RN<- blups19_pyn %>% 
  select(Variety,`GNDVI_2019-05-22_RN`, `GRYLD_2019-07-01_RN`) %>%
  left_join(gndvi20190522RN_snp, by = c("Variety" = "snp"))

p1<- gndvi20190522RN %>% 
  ggplot(aes(x = S2A_17830617, y = `GNDVI_2019-05-22_RN`)) +
  geom_boxplot() 
p1
p2<- gndvi20190522RN %>% 
  ggplot(aes(x = S2A_17830617, y = `GRYLD_2019-07-01_RN`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p1,p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 22)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_gndvi20190522RN.png",
                 width = 1000,
                 height = 750)

gndvi20190530RN_snp<- pyn19 %>% 
  select(snp, chrom, pos, GNDVI_2019.05.30_RN) %>% 
  filter(GNDVI_2019.05.30_RN == max(GNDVI_2019.05.30_RN))
gndvi20190530RN_snp<- geno %>% 
  select(snp, gndvi20190530RN_snp[1,1])

gndvi20190530RN<- blups19_pyn %>% 
  select(Variety,`GNDVI_2019-05-30_RN`, `GRYLD_2019-07-01_RN`) %>%
  left_join(gndvi20190530RN_snp, by = c("Variety" = "snp"))

p1<- gndvi20190530RN %>% 
  ggplot(aes(x = S2A_19902461, y = `GNDVI_2019-05-30_RN`)) +
  geom_boxplot() 
p1
p2<- gndvi20190530RN %>% 
  ggplot(aes(x = S2A_19902461, y = `GRYLD_2019-07-01_RN`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p1,p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 30)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_gndvi20190530RN.png",
                 width = 1000,
                 height = 750)

gndvi20190604RN_snp<- pyn19 %>% 
  select(snp, chrom, pos, GNDVI_2019.06.04_RN) %>% 
  filter(GNDVI_2019.06.04_RN == max(GNDVI_2019.06.04_RN))
gndvi20190604RN_snp<- geno %>% 
  select(snp, gndvi20190604RN_snp[1,1])

gndvi20190604RN<- blups19_pyn %>% 
  select(Variety,`GNDVI_2019-06-04_RN`, `GRYLD_2019-07-01_RN`) %>%
  left_join(gndvi20190604RN_snp, by = c("Variety" = "snp"))

p1<- gndvi20190604RN %>% 
  ggplot(aes(x = S2A_19902461, y = `GNDVI_2019-06-04_RN`)) +
  geom_boxplot() 
p1
p2<- gndvi20190604RN %>% 
  ggplot(aes(x = S2A_19902461, y = `GRYLD_2019-07-01_RN`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p1,p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 04)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_gndvi20190604RN.png",
                 width = 1000,
                 height = 750)

gryld20190701RN_snp<- pyn19 %>% 
  select(snp, chrom, pos, GRYLD_2019.07.01_RN) %>% 
  filter(GRYLD_2019.07.01_RN == max(GRYLD_2019.07.01_RN))
gryld20190701RN_snp<- geno %>% 
  select(snp, gryld20190701RN_snp[1,1])

gryld20190701RN<- blups19_pyn %>% 
  select(Variety,`GRYLD_2019-07-01_RN`) %>%
  left_join(gryld20190701RN_snp, by = c("Variety" = "snp"))

p2<- gryld20190701RN %>% 
  ggplot(aes(x = S1A_568381115, y = `GRYLD_2019-07-01_RN`)) +
  geom_boxplot() 
p2
fig<- ggpubr::ggarrange(p2)
fig
fig<- ggpubr::annotate_figure(fig, 
                              fig.lab = "Differences in allele phenotype for significant SNP",
                              fig.lab.size = 04)
ggpubr::ggexport(fig, filename = "./Figures/GWAS/Allele_effects_gryld20190701RN.png",
                 width = 1000,
                 height = 750)
