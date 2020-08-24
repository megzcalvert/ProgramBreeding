rm(list = objects())
ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(broom)
library(reshape2)
library(corrplot)
library(stringr)
library(vcfR)

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

######### Reading in files as a list of data frames
# 
varietyNames <- fread("./PhenoDatabase/VarietyNames_Corrected.txt")
varietyNames <- varietyNames %>%
  select(Variety)
# 
write.table(varietyNames, "./GenoDatabase/MasterBreeding/varietiesSelected.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

## Key file
keyfile<- fread("./GenoDatabase/MasterBreedingProgram.txt")

keyfile<- keyfile %>% 
  mutate(Variety = str_replace(Variety, "-M-", "M-")) %>%
  mutate(Variety = str_replace(Variety, "-K-", "K-")) %>%
  mutate(Variety = tolower(Variety)) %>%
  mutate(Variety = str_squish(Variety)) %>%
  mutate(Variety = str_replace(Variety, "-", "_")) %>%
  mutate(Variety = str_replace(Variety, "wb_4458", "wb4458")) %>%
  mutate(Variety = str_replace(Variety, "symonument", "monument")) %>%
  mutate(Variety = str_replace(Variety, "smith'sgold", "smithsgold")) %>% 
  semi_join(varietyNames)

write.table(keyfile, "./GenoDatabase/HypothesisTwo/KeyFile_BP.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#### Numeric output ####
geno <- fread("./GenoDatabase/HypothesisTwo/numericImputed.txt")
# 
geno[1:10, 1:10]
# 
geno <- geno %>%
  rename(Variety = `<Marker>`) %>%
  mutate(Variety = tolower(Variety)) %>%
  mutate(Variety = str_remove_all(Variety, " ")) %>% 
  dplyr::arrange(Variety)
 
geno[1:10, 1:10]
 
genoSelected <- geno %>%
  semi_join(varietyNames, by = "Variety") %>% 
  arrange(Variety)

genoMissing <- varietyNames %>% 
  anti_join(geno, by = "Variety")

phenoMissing<- geno %>% 
  anti_join(varietyNames, by = "Variety")
  
genoSelected[1:10, 1:10]

rm(geno)

positions <- fread("./GenoDatabase/HypothesisTwo/HapMapFormat.hmp.txt",
  select = c(1:4)
)
# 
positions <- positions %>%
  rename(snp = `rs#`) %>%
  select(-alleles)
# 
genoSelected <- genoSelected %>%
  distinct(Variety, .keep_all = TRUE) %>% 
  arrange(Variety)
# 
genoSelected <- t(genoSelected)
genoSelected <- setDT(as.data.frame(genoSelected), keep.rownames = TRUE)
genoSelected[1:10,1:10]

write.table(genoSelected, "./GenoDatabase/HypothesisTwo/SelectedGeno_Numeric2.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

genoSelected <- fread("./GenoDatabase/HypothesisTwo/SelectedGeno_Numeric2.txt")

colnames(genoSelected)[1] <- "snp"

genoSelected[genoSelected == 0] <- -1
genoSelected[genoSelected == 0.5] <- 0

genoSelected <- genoSelected %>%
  inner_join(positions, by = "snp") %>%
  select(snp, chrom, pos, everything())

rm(positions)

genoSelected[1:10, 1:10]
# 
snpNames<- genoSelected[,1]
# 
write.table(genoSelected,
  "./GenoDatabase/SelectedGeno_Numeric.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)
# 
genoSelected <- fread("./GenoDatabase/SelectedGeno_Numeric.txt")
# 
str(genoSelected)
# 
positions<- genoSelected[,c(1,2,3)]
# 
snpMatrix <- t(genoSelected[, c(-1, -2, -3)])
# 
pcaMethods::checkData(snpMatrix) # Check PCA assumptions
#
pcaAM <- pcaMethods::pca(snpMatrix, nPcs = 5) # SVD PCA
# 
sumPCA <- as.data.frame(summary(pcaAM))
# #
Scores <- as.data.frame(pcaMethods::scores(pcaAM))
Scores <- setDT(Scores, keep.rownames = TRUE) %>%
  mutate(founders = if_else(
    rn %in% c(
      "everest", "kanmark", "joe", "wb4458", "zenda", "danby",
      "bobdole", "gallagher", "monument", "wb_cedar"
    ), "yes", "no"
  ))
# 
p <- Scores %>%
  ggplot(aes(x = PC1, y = PC3)) +
  geom_point() +
  gghighlight::gghighlight(founders == "yes", label_key = rn) +
  theme(aspect.ratio = 1) +
  labs(x = paste0("PC1 = ", sumPCA["R2","PC1"] * 100, "%"),
       y = paste0("PC3 = ", sumPCA["R2","PC3"] * 100, "%"))
p
ggsave("./Figures/PCA_2.png",
       height = 25,
       width = 25, units = "cm", dpi = 350
)
# 
# snpMatrix[1:10,1:10]
# colnames(snpMatrix)<- snpNames
# snpMatrix[1:10,1:10]
# # 
# Linkage<- cor(snpMatrix)
# 
# Linkage[1:5,1:5]
# 
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
# flattenCorrMatrix <- function(cormat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut]
#   )
# }
# 
# Linkage<- flattenCorrMatrix(cormat = Linkage)
# 
# Linkage[1:5,]
# 
# Linkage<- Linkage %>%
#  inner_join(positions, by = c("row" = "snp")) %>%
#   rename(chromRow = chrom,
#          posRow = pos) %>%
#   inner_join(positions, by = c("column" = "snp")) %>%
#   rename(chromCol = chrom,
#          posCol = pos) %>%
#   filter(chromRow == chromCol)
# 
# Linkage[1:5,]
# 
# write.table(Linkage, file = "./Results/snpCorrelations.txt", quote = FALSE,
#             sep = "\t", row.names = FALSE, col.names = TRUE)

snpCor<- fread("./Results/snpCorrelations.txt")

snpCor<- snpCor %>% 
  mutate(distance = abs(posRow - posCol) / 1000000,
         absoluteCor = abs(cor))

pdf("./Figures/absCorAbsDistance.pdf",
    width = 15, 
    height = 10)

snpCor %>% 
  ggplot(aes(x = distance, y = absoluteCor)) +
  geom_smooth() +
  labs(title = "Linkage disequilibrium for KSU Breeding Program",
       x = "Distance (Mb)")


