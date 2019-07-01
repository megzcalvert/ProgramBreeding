rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)

pheno<- fread("./PhenoDatabase/PhenoLong.txt")
phenoLines<- as.data.frame(unique(pheno$Variety))
colnames(phenoLines)<- "Variety"
phenoLines$Variety<- as.character(tolower(phenoLines$Variety))


genoMaster<- fread("./GenoDatabase/Master_2018_Marshall.txt")
genoMaster<- genoMaster %>% 
  mutate(Variety = tolower(FullSampleName)) %>% 
  glimpse()

genoFile<- genoMaster %>% 
  inner_join(phenoLines, by = c("Variety" = "Variety"))

genoLines<- as.data.frame(unique(genoFile$Variety))
colnames(genoLines)<- "Variety"

missinglines<- phenoLines %>% 
  anti_join(genoLines)

write.table(genoFile, "./GenoDatabase/BreedingProgramLines.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)

#############################################################################

snpChip <- read_delim(
  "./GenoDatabase/BreedingProgramSelected_Imputed.hmp.txt", 
  "\t", escape_double = FALSE, trim_ws = TRUE)
snpChip<- snpChip %>% 
  clean_names()
snpChip<- separate(snpChip,alleles,c("allele_a","allele_b"), sep = "/")

missAmbiguous = c('0', '+', '-')
hetCodes = c('R','Y','S','W','K','M','B','D','H','V')
hapgeno=as.matrix(snpChip[,13:ncol(snpChip)])
hapgeno[hapgeno %in% missAmbiguous]=NA
hapgeno[hapgeno=='N']=NA
hapgeno[hapgeno %in% hetCodes]='H'
snpChip=cbind(snpChip[,1:12], hapgeno)
rm(hapgeno)

write.table(snpChip, file="./GenoDatabase/SelectedImputedBeagle.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./GenoDatabase/SelectedImputedBeagle.txt", 
                header=TRUE, check.names=F, sep = "\t")

snpChip[snpChip == snpChip$allele_a] = -1
snpChip[snpChip == snpChip$allele_b] = 1
snpChip[snpChip == "H"] = 0
snpChip[snpChip == "C"] = NA
snpChip[snpChip == "A"] = NA
snpChip[snpChip == "T"] = NA
snpChip[snpChip == "G"] = NA
snpChip[snpChip == "-"] = NA
snpChip[snpChip == "."] = NA

snpChip<- snpChip[ ,c(1,4,5,13:826)]

write.table(snpChip, file="./GenoDatabase/SelectedImputedBeagleNumeric.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./GenoDatabase/SelectedImputedBeagleNumeric.txt", 
                header=TRUE, check.names=F, sep = "\t")

chrSum<- plyr::count(snpChip, vars = "chrom")
snpChip<- snpChip %>% 
  filter(chrom != "UN")
snpMatrix<- t(snpChip[ , c(-1, -2, -3)])

pcaMethods::checkData(snpMatrix)  #Check PCA assumptions

pcaAM<- pcaMethods::pca(snpMatrix, nPcs = 15) #SVD PCA

sumPCA<- as.data.frame(summary(pcaAM))

Scores<- as.data.frame(pcaMethods::scores(pcaAM)) 
Scores<- setDT(Scores, keep.rownames = TRUE)

Scores<- Scores %>% 
  mutate(Year = str_sub(rn,3,4),
         Year = as.numeric(Year),
         Year = replace_na(Year,"Unknown"),
         Year = str_replace(Year,"44","Unknown"),
         dh = str_detect(rn,"dh"))

Scores %>% 
  ggplot(aes(x = PC2, y = PC3, colour = factor(Year), shape = dh)) +
  geom_point(alpha = 1) +
  scale_color_manual(name = "Year", 
                     values = c('#1b9e77','#d95f02','#7570b3',
                                '#e7298a','#66a61e','#e6ab02',
                                '#a6761d','#666666','#000000')) +
  scale_shape_manual(name = "DoubleHaploid", 
                     values = c(17,19)) +
  theme_bw() +
  theme(axis.title = element_text(colour = "black", size = 16),
        axis.text = element_text(colour = "black", size = 14),
        aspect.ratio = 1:1,
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid.major = element_line(colour = "#b7b7b7",linetype = 2),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(colour = "black", size = 14)) +
  labs(title = "PCA of Breeding Program Markers",
       subtitle = "AYN and PYN in 2016-2019",
       x = expression(paste("PC2 ",R^2, " = 3.9%")),
       y =  expression(paste("PC3 ",R^2, " = 2.6%")))


