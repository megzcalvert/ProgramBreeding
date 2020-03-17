rm(list = objects())
ls()

library(tidyverse)
library(data.table)
library(tidylog)
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

######### Reading in files as a list of data frames

varietyNames <- fread("./PhenoDatabase/VarietyNames_Corrected.txt")
varietyNames<- varietyNames %>% 
  select(Variety)

write.table(varietyNames, "./GenoDatabase/MasterBreeding/varietiesSelected.txt",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

geno <- fread("./GenoDatabase/MasterBreeding/Numeric_output.txt")

geno[1:10,1:10]

geno<- geno %>% 
  rename(Variety = `<Marker>`) %>% 
  mutate(Variety = tolower(Variety)) %>% 
  mutate(Variety = str_remove_all(Variety," ")) 

geno[1:10,1:10]

genoSelected<- geno %>% 
  semi_join(varietyNames, by = "Variety")

genoSelected[1:10,1:10]

rm(geno)

positions<- fread("./GenoDatabase/MasterBreeding/Merged_Genotype_Table.hmp.txt",
                  select = c(1:4))

positions<- positions %>% 
  rename(snp = `rs#`) %>% 
  select(-alleles)

genoSelected<- genoSelected %>% 
  distinct(Variety, .keep_all = TRUE) 

genoSelected<- t(genoSelected)  
genoSelected<- setDT(as.data.frame(genoSelected), keep.rownames = TRUE) 

write.table(genoSelected, "./GenoDatabase/MasterBreeding/SelectedGeno_Numeric2.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

genoSelected<- fread("./GenoDatabase/MasterBreeding/SelectedGeno_Numeric2.txt")

colnames(genoSelected)[1] <- "snp"

genoSelected<- genoSelected %>% 
  inner_join(positions, by = "snp") %>% 
  select(snp, chrom, pos, everything())

heterozygosityTest <- genoSelected %>%
  select(-chrom, -pos) %>%
  pivot_longer(cols = -snp, names_to = "Variety") %>%
  group_by(snp, value) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = value, values_from = count) %>%
  mutate(Total = `0` + `0.5` + `1`,
         homoRecessive = `0` / Total,
         heterozygous = `0.5` / Total,
         homoDominant = `1` / Total,
         Check = homoRecessive + heterozygous + homoDominant) %>% 
  filter(homoRecessive > 0.05) %>% 
  filter(heterozygous > homoDominant)

genoSelected<- genoSelected %>% 
  semi_join(heterozygosityTest, by = "snp")

heterozygosityTest %>% 
  ggplot(aes(x = heterozygous)) +
  geom_histogram(colour = "black", fill = NA, binwidth = 0.001)

genoSelected[1:10,1:10]

write.table(genoSelected, 
            "./GenoDatabase/SelectedGeno_Numeric.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

genoSelected<- fread("./GenoDatabase/SelectedGeno_Numeric.txt")

str(genoSelected)

snpMatrix <- t(genoSelected[, c(-1, -2, -3)])

pcaMethods::checkData(snpMatrix) # Check PCA assumptions

pcaAM <- pcaMethods::pca(snpMatrix, nPcs = 5) # SVD PCA

sumPCA <- as.data.frame(summary(pcaAM))

Scores <- as.data.frame(pcaMethods::scores(pcaAM))
Scores <- setDT(Scores, keep.rownames = TRUE) %>%
  mutate(founders = if_else(
    rn %in% c("everest", "kanmark", "joe", "wb4458", "zenda", "danby",
              "bobdole", "gallagher", "monument"), "yes", "no"
  ))

p <- Scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  gghighlight::gghighlight(founders == "yes", label_key = rn)
p
