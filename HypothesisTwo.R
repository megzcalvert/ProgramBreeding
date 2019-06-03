rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)

pheno<- fread("./PhenoDatabase/PhenoLong.txt")
phenoLines<- as.data.frame(unique(pheno$Variety))
colnames(phenoLines)<- "Variety"

genoMaster<- fread("./GenoDatabase/Master_2018_Marshall.txt")

genoFile<- genoMaster %>% 
  inner_join(phenoLines, by = c("FullSampleName" = "Variety"))

genoLines<- as.data.frame(unique(genoFile$FullSampleName))
colnames(genoLines)<- "Variety"

missinglines<- phenoLines %>% 
  anti_join(genoLines)

write.table(genoFile, "./GenoDatabase/BreedingProgramLines.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)

#############################################################################

snpChip <- read_delim(
  "./GenoDatabase/BreedingProgramSelected_imputed.hmp.txt", 
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

snpChip<- snpChip[ ,c(1,4,5,13:311)]

write.table(snpChip, file="./GenoDatabase/SelectedImputedBeagleNumeric.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./GenoDatabase/SelectedImputedBeagleNumeric.txt", 
                header=TRUE, check.names=F, sep = "\t")

chrSum<- plyr::count(snpChip, vars = "chrom")
snpMatrix<- t(snpChip[ , c(-1, -2, -3)])

pcaMethods::checkData(snpMatrix)  #Check PCA assumptions

pcaAM<- pcaMethods::pca(snpMatrix, nPcs = 15) #SVD PCA

sumPCA<- as.data.frame(summary(pcaAM))

Scores<- as.data.frame(pcaMethods::scores(pcaAM)) 
Scores<- setDT(Scores, keep.rownames = TRUE)

Scores %>% 
  ggplot(aes(x = PC1, y = PC3)) +
  geom_point()


