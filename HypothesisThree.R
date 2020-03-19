rm(list = objects())
ls()

library(MatrixModels)
library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(beepr)
library(rrBLUP)
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
setwd("/Users/megzcalvert/Dropbox/Research_Poland_Lab/BreedingProgram")

set.seed(1642)

#### Load data ####

# geno <- fread("./GenoDatabase/SelectedGeno_Numeric.txt")
# geno[1:10, 1:10]
# hapgeno <- geno[, 4:ncol(geno)]
# hapgeno[hapgeno == 0] <- -1
# hapgeno[hapgeno == 0.5] <- 0
# 
# geno <- bind_cols(geno[, 1:3], hapgeno)
# geno[1:10, 1:10]
# 
# pheno <- fread("./PhenoDatabase/Blups_allYrs_allVarieties.txt")
# pheno[1:10, ]
# 
# genoVarieties <- as.data.frame(colnames(geno[, 4:ncol(geno)]))
# names(genoVarieties) <- "Variety"
# pheno <- pheno %>%
#   semi_join(genoVarieties)
# 
# # transpose and remove positions
# snpMatrix <- t(geno[, c(-1, -2, -3)])
# snpMatrix %>% glimpse()
# snpMatrix[1:10, 1:10]
# 
# snpMatrix <- setDT(as.data.frame(snpMatrix), keep.rownames = T)
# snpMatrix[1:10, 1:10]
# 
# snpMatrix <- snpMatrix %>%
#   semi_join(pheno, by = c("rn" = "Variety"))
# rownames(snpMatrix) <- snpMatrix[, 1]
# snpMatrix[, 1] <- NULL
# 
# snpMatrix <- as.matrix(snpMatrix)
# 
# ##### Defining the training and test populations ####
# 
# # define the training and test populations
# # training-80% validation-20%
# Pheno_train <- pheno %>%
#   dplyr::sample_frac(0.8)
# Pheno_valid <- pheno %>%
#   anti_join(Pheno_train, by = "Variety")
# 
# m_train <- snpMatrix[Pheno_train$Variety, ]
# m_valid <- snpMatrix[Pheno_valid$Variety, ]
# 
# ##### Predicting Phenotypes ####
# ## GRYLD
# yield <- (Pheno_train[, "GRYLD"])
# 
# yield_answer <- mixed.solve(yield,
#                             Z = m_train, 
#                             SE = TRUE
# )
# YLD <- yield_answer$u
# e <- as.matrix(YLD)
# pred_yield_valid <- m_valid %*% e
# pred_yield <- (pred_yield_valid[, 1]) + yield_answer$beta
# pred_yield
# yield_valid <- Pheno_valid[, "GRYLD"]
# YLD_accuracy <- cor(pred_yield_valid, yield_valid)
# YLD_accuracy
# 
# #### Cross-validation ####
# 
# traits=1
# cycles=100
# accuracy = matrix(nrow=cycles, ncol=traits)
# for(r in 1:cycles)
# {
# train= as.matrix(sample(1:96, 29))
# test<-setdiff(1:96,train)
# Pheno_train=Pheno[train,]
# m_train=Markers_impute2[train,]
# Pheno_valid=Pheno[test,]
# m_valid=Markers_impute2[test,]
# 
# yield=(Pheno_train[,1])
# yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
# YLD = yield_answer$u
# e = as.matrix(YLD)
# pred_yield_valid =  m_valid %*% e
# pred_yield=(pred_yield_valid[,1])+yield_answer$beta
# pred_yield
# yield_valid = Pheno_valid[,1]
# accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
# }
# mean(accuracy)
# 
# 
# write.csv(accuracy17,
#           "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy17.txt",
#           quote = F, row.names = F
# )
# colMeans(accuracy17)

accuracy<- fread("./BeocatScripts/GenomicSelection_80_100_accuracy_allyrs.txt")

accuracy %>% 
  ggplot(aes(V1)) + 
  geom_histogram(colour = "black", fill = NA) +
  geom_vline(xintercept = mean(accuracy$V1), colour = "blue") +
  labs(title = "Distribution of accuracy of GRYLD genomic prediction",
       subtitle = "rrBLUP",
       x = "Accuracy")
  
