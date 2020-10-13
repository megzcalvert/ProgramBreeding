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

##### Read in geno data ####
snpInfo<- fread("./GenoDatabase/SelectedGeno_Numeric.txt")

positions<- snpInfo[,1:3]
snpMatrix <- setDT(as.data.frame(t(snpInfo[,4:ncol(snpInfo)])),
                   keep.rownames = TRUE) %>% 
  dplyr::arrange(rn)

#### Pheno All ####
pheno<- fread("./Results/BlupsGRYLD_allYrs_allVarieties.txt") 

snpMatrix_all<- snpMatrix %>% 
  semi_join(pheno, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno<- pheno %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_all) == pheno$Variety), sep = ' '))

write.table(snpMatrix_all,
            "./GenoDatabase/snpMatrix_all.txt",
            quote = FALSE, col.names = TRUE, row.names = TRUE,
            sep = "\t")
geno_all<- fread("./GenoDatabase/snpMatrix_all.txt")

geno_all<- as.data.frame(t(geno_all)) %>% 
  row_to_names(row_number = 1) %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

write.table(geno_all,
            "./GenoDatabase/geno_all.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_all<- fread("./GenoDatabase/geno_all.txt")

rel_mat<- A.mat(snpMatrix_all)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(rel_mat) == pheno$Variety), sep = ' '))

gwas_gryld_all<- GWAS(pheno = pheno,
                      geno = geno_all,
                      K = rel_mat,
                      n.PC = 4,
                      P3D = FALSE,
                      plot = TRUE)

#### Pheno 16 ####
pheno_16<- fread("./Results/BlupsGRYLD_16.txt") 

snpMatrix_16<- snpMatrix %>% 
  semi_join(pheno_16, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_16<- pheno_16 %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_16) == pheno_16$Variety), sep = ' '))

write.table(snpMatrix_16,
            "./GenoDatabase/snpMatrix_16.txt",
            quote = FALSE, col.names = TRUE, row.names = TRUE,
            sep = "\t")
geno_16<- fread("./GenoDatabase/snpMatrix_16.txt")

geno_16<- as.data.frame(t(geno_16)) %>% 
  row_to_names(row_number = 1) %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

write.table(geno_16,
            "./GenoDatabase/geno_16.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_16<- fread("./GenoDatabase/geno_16.txt")

rel_mat<- A.mat(snpMatrix_16)

gwas_gryld_16<- GWAS(pheno = pheno_16,
                      geno = geno_16,
                      K = rel_mat,
                      n.PC = 3,
                      P3D = TRUE,
                      plot = TRUE)

#### Pheno 16 ayn ####
pheno_16_ayn<- fread("./Results/BlupsGRYLD_16_ayn.txt") 

snpMatrix_16_ayn<- snpMatrix %>% 
  semi_join(pheno_16_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_16_ayn<- pheno_16_ayn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_16_ayn) == pheno_16_ayn$Variety), sep = ' '))

write.table(snpMatrix_16_ayn,
            "./GenoDatabase/snpMatrix_16_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = TRUE,
            sep = "\t")
geno_16_ayn<- fread("./GenoDatabase/snpMatrix_16_ayn.txt")

geno_16_ayn<- as.data.frame(t(geno_16_ayn)) %>% 
  row_to_names(row_number = 1) %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

write.table(geno_16_ayn,
            "./GenoDatabase/geno_16_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_16_ayn<- fread("./GenoDatabase/geno_16_ayn.txt")

rel_mat<- A.mat(snpMatrix_16_ayn)

gwas_gryld_16_ayn<- GWAS(pheno = pheno_16_ayn,
                     geno = geno_16_ayn,
                     K = rel_mat,
                     n.PC = 3,
                     P3D = TRUE,
                     plot = TRUE)

#### Pheno 16 pyn ####
pheno_16_pyn<- fread("./Results/BlupsGRYLD_16_pyn.txt") 

snpMatrix_16_pyn<- snpMatrix %>% 
  semi_join(pheno_16_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_16_pyn<- pheno_16_pyn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_16_pyn) == pheno_16_pyn$Variety), sep = ' '))

write.table(snpMatrix_16_pyn,
            "./GenoDatabase/snpMatrix_16_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = TRUE,
            sep = "\t")
geno_16_pyn<- fread("./GenoDatabase/snpMatrix_16_pyn.txt")

geno_16_pyn<- as.data.frame(t(geno_16_pyn)) %>% 
  row_to_names(row_number = 1) %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

write.table(geno_16_pyn,
            "./GenoDatabase/geno_16_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_16_pyn<- fread("./GenoDatabase/geno_16_pyn.txt")

rel_mat<- A.mat(snpMatrix_16_pyn)

gwas_gryld_16_pyn<- GWAS(pheno = pheno_16_pyn,
                         geno = geno_16_pyn,
                         K = rel_mat,
                         n.PC = 3,
                         P3D = TRUE,
                         plot = TRUE)

#### Pheno ex16 ####
pheno_ex16<- fread("./Results/BlupsGRYLD_ex16.txt") 

snpMatrix_ex16<- snpMatrix %>% 
  semi_join(pheno_ex16, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex16<- pheno_ex16 %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex16) == pheno_ex16$Variety), sep = ' '))

write.table(snpMatrix_ex16,
            "./GenoDatabase/snpMatrix_ex16.txt",
            quote = FALSE, col.names = TRUE, row.names = TRUE,
            sep = "\t")
geno_ex16<- fread("./GenoDatabase/snpMatrix_ex16.txt")

geno_ex16<- as.data.frame(t(geno_ex16)) %>% 
  row_to_names(row_number = 1) %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

write.table(geno_ex16,
            "./GenoDatabase/geno_ex16.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex16<- fread("./GenoDatabase/geno_ex16.txt")

rel_mat<- A.mat(snpMatrix_ex16)

gwas_gryld_ex16<- GWAS(pheno = pheno_ex16,
                     geno = geno_ex16,
                     K = rel_mat,
                     n.PC = 3,
                     P3D = TRUE,
                     plot = TRUE)

#### Pheno ex16 ayn ####
pheno_ex16_ayn<- fread("./Results/BlupsGRYLD_ex16_ayn.txt") 

snpMatrix_ex16_ayn<- snpMatrix %>% 
  semi_join(pheno_ex16_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex16_ayn<- pheno_ex16_ayn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex16_ayn) == pheno_ex16_ayn$Variety), sep = ' '))

write.table(snpMatrix_ex16_ayn,
            "./GenoDatabase/snpMatrix_ex16_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = TRUE,
            sep = "\t")
geno_ex16_ayn<- fread("./GenoDatabase/snpMatrix_ex16_ayn.txt")

geno_ex16_ayn<- as.data.frame(t(geno_ex16_ayn)) %>% 
  row_to_names(row_number = 1) %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

write.table(geno_ex16_ayn,
            "./GenoDatabase/geno_ex16_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex16_ayn<- fread("./GenoDatabase/geno_ex16_ayn.txt")

rel_mat<- A.mat(snpMatrix_ex16_ayn)

gwas_gryld_ex16_ayn<- GWAS(pheno = pheno_ex16_ayn,
                         geno = geno_ex16_ayn,
                         K = rel_mat,
                         n.PC = 3,
                         P3D = TRUE,
                         plot = TRUE)

#### Pheno ex16 pyn ####
pheno_ex16_pyn<- fread("./Results/BlupsGRYLD_ex16_pyn.txt") 

snpMatrix_ex16_pyn<- snpMatrix %>% 
  semi_join(pheno_ex16_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex16_pyn<- pheno_ex16_pyn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex16_pyn) == pheno_ex16_pyn$Variety), sep = ' '))

write.table(snpMatrix_ex16_pyn,
            "./GenoDatabase/snpMatrix_ex16_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = TRUE,
            sep = "\t")
geno_ex16_pyn<- fread("./GenoDatabase/snpMatrix_ex16_pyn.txt")

geno_ex16_pyn<- as.data.frame(t(geno_ex16_pyn)) %>% 
  row_to_names(row_number = 1) %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

write.table(geno_ex16_pyn,
            "./GenoDatabase/geno_ex16_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex16_pyn<- fread("./GenoDatabase/geno_ex16_pyn.txt")

rel_mat<- A.mat(snpMatrix_ex16_pyn)

gwas_gryld_ex16_pyn<- GWAS(pheno = pheno_ex16_pyn,
                         geno = geno_ex16_pyn,
                         K = rel_mat,
                         n.PC = 3,
                         P3D = TRUE,
                         plot = TRUE)

#### Pheno 17 ####
pheno_17<- fread("./Results/BlupsGRYLD_17.txt") 

snpMatrix_17<- snpMatrix %>% 
  semi_join(pheno_17, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_17<- pheno_17 %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_17) == pheno_17$Variety), sep = ' '))

geno_17<-snpMatrix %>% 
  semi_join(pheno_17, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_17<- as.data.frame(t(geno_17)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_17,
            "./GenoDatabase/snpMatrix_17.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_17<- fread("./GenoDatabase/snpMatrix_17.txt")

geno_17<- geno_17 %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_17)

gwas_gryld_17<- GWAS(pheno = pheno_17,
                     geno = geno_17,
                     K = rel_mat,
                     n.PC = 3,
                     P3D = TRUE,
                     plot = TRUE)

#### Pheno 17 ayn ####
pheno_17_ayn<- fread("./Results/BlupsGRYLD_17_ayn.txt") 

snpMatrix_17_ayn<- snpMatrix %>% 
  semi_join(pheno_17_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_17_ayn<- pheno_17_ayn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_17_ayn) == pheno_17_ayn$Variety), sep = ' '))

geno_17_ayn<-snpMatrix %>% 
  semi_join(pheno_17_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_17_ayn<- as.data.frame(t(geno_17_ayn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_17_ayn,
            "./GenoDatabase/snpMatrix_17_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_17_ayn<- fread("./GenoDatabase/snpMatrix_17_ayn.txt")

geno_17_ayn<- geno_17_ayn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_17_ayn)

gwas_gryld_17_ayn<- GWAS(pheno = pheno_17_ayn,
                         geno = geno_17_ayn,
                         K = rel_mat,
                         n.PC = 3,
                         P3D = TRUE,
                         plot = TRUE)

#### Pheno 17 pyn ####
pheno_17_pyn<- fread("./Results/BlupsGRYLD_17_pyn.txt") 

snpMatrix_17_pyn<- snpMatrix %>% 
  semi_join(pheno_17_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_17_pyn<- pheno_17_pyn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_17_pyn) == pheno_17_pyn$Variety), sep = ' '))

geno_17_pyn<-snpMatrix %>% 
  semi_join(pheno_17_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_17_pyn<- as.data.frame(t(geno_17_pyn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_17_pyn,
            "./GenoDatabase/snpMatrix_17_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_17_pyn<- fread("./GenoDatabase/snpMatrix_17_pyn.txt")

geno_17_pyn<- geno_17_pyn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_17_pyn)

gwas_gryld_17_pyn<- GWAS(pheno = pheno_17_pyn,
                         geno = geno_17_pyn,
                         K = rel_mat,
                         n.PC = 3,
                         P3D = TRUE,
                         plot = TRUE)

#### Pheno ex17 ####
pheno_ex17<- fread("./Results/BlupsGRYLD_ex17.txt") 

snpMatrix_ex17<- snpMatrix %>% 
  semi_join(pheno_ex17, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex17<- pheno_ex17 %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex17) == pheno_ex17$Variety), sep = ' '))

geno_ex17<-snpMatrix %>% 
  semi_join(pheno_ex17, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_ex17<- as.data.frame(t(geno_ex17)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_ex17,
            "./GenoDatabase/snpMatrix_ex17.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex17<- fread("./GenoDatabase/snpMatrix_ex17.txt")

geno_ex17<- geno_ex17 %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_ex17)

gwas_gryld_ex17<- GWAS(pheno = pheno_ex17,
                       geno = geno_ex17,
                       K = rel_mat,
                       n.PC = 3,
                       P3D = TRUE,
                       plot = TRUE)

#### Pheno ex17 ayn ####
pheno_ex17_ayn<- fread("./Results/BlupsGRYLD_ex17_ayn.txt") 

snpMatrix_ex17_ayn<- snpMatrix %>% 
  semi_join(pheno_ex17_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex17_ayn<- pheno_ex17_ayn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex17_ayn) == pheno_ex17_ayn$Variety), sep = ' '))

geno_ex17_ayn<-snpMatrix %>% 
  semi_join(pheno_ex17_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_ex17_ayn<- as.data.frame(t(geno_ex17_ayn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_ex17_ayn,
            "./GenoDatabase/snpMatrix_ex17_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex17_ayn<- fread("./GenoDatabase/snpMatrix_ex17_ayn.txt")

geno_ex17_ayn<- geno_ex17_ayn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_ex17_ayn)

gwas_gryld_ex17_ayn<- GWAS(pheno = pheno_ex17_ayn,
                           geno = geno_ex17_ayn,
                           K = rel_mat,
                           n.PC = 3,
                           P3D = TRUE,
                           plot = TRUE)

#### Pheno ex17 pyn ####
pheno_ex17_pyn<- fread("./Results/BlupsGRYLD_ex17_pyn.txt") 

snpMatrix_ex17_pyn<- snpMatrix %>% 
  semi_join(pheno_ex17_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex17_pyn<- pheno_ex17_pyn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex17_pyn) == pheno_ex17_pyn$Variety), sep = ' '))

geno_ex17_pyn<-snpMatrix %>% 
  semi_join(pheno_ex17_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_ex17_pyn<- as.data.frame(t(geno_ex17_pyn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_ex17_pyn,
            "./GenoDatabase/snpMatrix_ex17_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex17_pyn<- fread("./GenoDatabase/snpMatrix_ex17_pyn.txt")

geno_ex17_pyn<- geno_ex17_pyn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_ex17_pyn)

gwas_gryld_ex17_pyn<- GWAS(pheno = pheno_ex17_pyn,
                           geno = geno_ex17_pyn,
                           K = rel_mat,
                           n.PC = 3,
                           P3D = TRUE,
                           plot = TRUE)

#### Pheno 18 ####
pheno_18<- fread("./Results/BlupsGRYLD_18.txt") 

snpMatrix_18<- snpMatrix %>% 
  semi_join(pheno_18, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_18<- pheno_18 %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_18) == pheno_18$Variety), sep = ' '))

geno_18<-snpMatrix %>% 
  semi_join(pheno_18, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_18<- as.data.frame(t(geno_18)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_18,
            "./GenoDatabase/snpMatrix_18.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_18<- fread("./GenoDatabase/snpMatrix_18.txt")

geno_18<- geno_18 %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_18)

gwas_gryld_18<- GWAS(pheno = pheno_18,
                     geno = geno_18,
                     K = rel_mat,
                     n.PC = 3,
                     P3D = TRUE,
                     plot = TRUE)

#### Pheno 18 ayn ####
pheno_18_ayn<- fread("./Results/BlupsGRYLD_18_ayn.txt") 

snpMatrix_18_ayn<- snpMatrix %>% 
  semi_join(pheno_18_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_18_ayn<- pheno_18_ayn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_18_ayn) == pheno_18_ayn$Variety), sep = ' '))

geno_18_ayn<-snpMatrix %>% 
  semi_join(pheno_18_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_18_ayn<- as.data.frame(t(geno_18_ayn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_18_ayn,
            "./GenoDatabase/snpMatrix_18_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_18_ayn<- fread("./GenoDatabase/snpMatrix_18_ayn.txt")

geno_18_ayn<- geno_18_ayn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_18_ayn)

gwas_gryld_18_ayn<- GWAS(pheno = pheno_18_ayn,
                         geno = geno_18_ayn,
                         K = rel_mat,
                         n.PC = 3,
                         P3D = TRUE,
                         plot = TRUE)

#### Pheno 18 pyn ####
pheno_18_pyn<- fread("./Results/BlupsGRYLD_18_pyn.txt") 

snpMatrix_18_pyn<- snpMatrix %>% 
  semi_join(pheno_18_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_18_pyn<- pheno_18_pyn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_18_pyn) == pheno_18_pyn$Variety), sep = ' '))

geno_18_pyn<-snpMatrix %>% 
  semi_join(pheno_18_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_18_pyn<- as.data.frame(t(geno_18_pyn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_18_pyn,
            "./GenoDatabase/snpMatrix_18_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_18_pyn<- fread("./GenoDatabase/snpMatrix_18_pyn.txt")

geno_18_pyn<- geno_18_pyn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_18_pyn)

gwas_gryld_18_pyn<- GWAS(pheno = pheno_18_pyn,
                         geno = geno_18_pyn,
                         K = rel_mat,
                         n.PC = 3,
                         P3D = TRUE,
                         plot = TRUE)

#### Pheno ex18 ####
pheno_ex18<- fread("./Results/BlupsGRYLD_ex18.txt") 

snpMatrix_ex18<- snpMatrix %>% 
  semi_join(pheno_ex18, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex18<- pheno_ex18 %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex18) == pheno_ex18$Variety), sep = ' '))

geno_ex18<-snpMatrix %>% 
  semi_join(pheno_ex18, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_ex18<- as.data.frame(t(geno_ex18)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_ex18,
            "./GenoDatabase/snpMatrix_ex18.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex18<- fread("./GenoDatabase/snpMatrix_ex18.txt")

geno_ex18<- geno_ex18 %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_ex18)

gwas_gryld_ex18<- GWAS(pheno = pheno_ex18,
                       geno = geno_ex18,
                       K = rel_mat,
                       n.PC = 3,
                       P3D = TRUE,
                       plot = TRUE)

#### Pheno ex18 ayn ####
pheno_ex18_ayn<- fread("./Results/BlupsGRYLD_ex18_ayn.txt") 

snpMatrix_ex18_ayn<- snpMatrix %>% 
  semi_join(pheno_ex18_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex18_ayn<- pheno_ex18_ayn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex18_ayn) == pheno_ex18_ayn$Variety), sep = ' '))

geno_ex18_ayn<-snpMatrix %>% 
  semi_join(pheno_ex18_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_ex18_ayn<- as.data.frame(t(geno_ex18_ayn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_ex18_ayn,
            "./GenoDatabase/snpMatrix_ex18_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex18_ayn<- fread("./GenoDatabase/snpMatrix_ex18_ayn.txt")

geno_ex18_ayn<- geno_ex18_ayn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_ex18_ayn)

gwas_gryld_ex18_ayn<- GWAS(pheno = pheno_ex18_ayn,
                           geno = geno_ex18_ayn,
                           K = rel_mat,
                           n.PC = 3,
                           P3D = TRUE,
                           plot = TRUE)

#### Pheno ex18 pyn ####
pheno_ex18_pyn<- fread("./Results/BlupsGRYLD_ex18_pyn.txt") 

snpMatrix_ex18_pyn<- snpMatrix %>% 
  semi_join(pheno_ex18_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex18_pyn<- pheno_ex18_pyn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex18_pyn) == pheno_ex18_pyn$Variety), sep = ' '))

geno_ex18_pyn<-snpMatrix %>% 
  semi_join(pheno_ex18_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_ex18_pyn<- as.data.frame(t(geno_ex18_pyn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_ex18_pyn,
            "./GenoDatabase/snpMatrix_ex18_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex18_pyn<- fread("./GenoDatabase/snpMatrix_ex18_pyn.txt")

geno_ex18_pyn<- geno_ex18_pyn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_ex18_pyn)

gwas_gryld_ex18_pyn<- GWAS(pheno = pheno_ex18_pyn,
                           geno = geno_ex18_pyn,
                           K = rel_mat,
                           n.PC = 3,
                           P3D = TRUE,
                           plot = TRUE)

#### Pheno 19 ####
pheno_19<- fread("./Results/BlupsGRYLD_19.txt") 

snpMatrix_19<- snpMatrix %>% 
  semi_join(pheno_19, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_19<- pheno_19 %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_19) == pheno_19$Variety), sep = ' '))

geno_19<-snpMatrix %>% 
  semi_join(pheno_19, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_19<- as.data.frame(t(geno_19)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_19,
            "./GenoDatabase/snpMatrix_19.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_19<- fread("./GenoDatabase/snpMatrix_19.txt")

geno_19<- geno_19 %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_19)

gwas_gryld_19<- GWAS(pheno = pheno_19,
                     geno = geno_19,
                     K = rel_mat,
                     n.PC = 3,
                     P3D = TRUE,
                     plot = TRUE)

#### Pheno 19 ayn ####
pheno_19_ayn<- fread("./Results/BlupsGRYLD_19_ayn.txt") 

snpMatrix_19_ayn<- snpMatrix %>% 
  semi_join(pheno_19_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_19_ayn<- pheno_19_ayn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_19_ayn) == pheno_19_ayn$Variety), sep = ' '))

geno_19_ayn<-snpMatrix %>% 
  semi_join(pheno_19_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_19_ayn<- as.data.frame(t(geno_19_ayn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_19_ayn,
            "./GenoDatabase/snpMatrix_19_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_19_ayn<- fread("./GenoDatabase/snpMatrix_19_ayn.txt")

geno_19_ayn<- geno_19_ayn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_19_ayn)

gwas_gryld_19_ayn<- GWAS(pheno = pheno_19_ayn,
                         geno = geno_19_ayn,
                         K = rel_mat,
                         n.PC = 3,
                         P3D = TRUE,
                         plot = TRUE)

#### Pheno 19 pyn ####
pheno_19_pyn<- fread("./Results/BlupsGRYLD_19_pyn.txt") 

snpMatrix_19_pyn<- snpMatrix %>% 
  semi_join(pheno_19_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_19_pyn<- pheno_19_pyn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_19_pyn) == pheno_19_pyn$Variety), sep = ' '))

geno_19_pyn<-snpMatrix %>% 
  semi_join(pheno_19_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_19_pyn<- as.data.frame(t(geno_19_pyn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_19_pyn,
            "./GenoDatabase/snpMatrix_19_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_19_pyn<- fread("./GenoDatabase/snpMatrix_19_pyn.txt")

geno_19_pyn<- geno_19_pyn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_19_pyn)

gwas_gryld_19_pyn<- GWAS(pheno = pheno_19_pyn,
                         geno = geno_19_pyn,
                         K = rel_mat,
                         n.PC = 3,
                         P3D = TRUE,
                         plot = TRUE)

#### Pheno ex19 ####
pheno_ex19<- fread("./Results/BlupsGRYLD_ex19.txt") 

snpMatrix_ex19<- snpMatrix %>% 
  semi_join(pheno_ex19, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex19<- pheno_ex19 %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex19) == pheno_ex19$Variety), sep = ' '))

geno_ex19<-snpMatrix %>% 
  semi_join(pheno_ex19, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_ex19<- as.data.frame(t(geno_ex19)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_ex19,
            "./GenoDatabase/snpMatrix_ex19.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex19<- fread("./GenoDatabase/snpMatrix_ex19.txt")

geno_ex19<- geno_ex19 %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_ex19)

gwas_gryld_ex19<- GWAS(pheno = pheno_ex19,
                       geno = geno_ex19,
                       K = rel_mat,
                       n.PC = 3,
                       P3D = TRUE,
                       plot = TRUE)

#### Pheno ex19 ayn ####
pheno_ex19_ayn<- fread("./Results/BlupsGRYLD_ex19_ayn.txt") 

snpMatrix_ex19_ayn<- snpMatrix %>% 
  semi_join(pheno_ex19_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex19_ayn<- pheno_ex19_ayn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex19_ayn) == pheno_ex19_ayn$Variety), sep = ' '))

geno_ex19_ayn<-snpMatrix %>% 
  semi_join(pheno_ex19_ayn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_ex19_ayn<- as.data.frame(t(geno_ex19_ayn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_ex19_ayn,
            "./GenoDatabase/snpMatrix_ex19_ayn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex19_ayn<- fread("./GenoDatabase/snpMatrix_ex19_ayn.txt")

geno_ex19_ayn<- geno_ex19_ayn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_ex19_ayn)

gwas_gryld_ex19_ayn<- GWAS(pheno = pheno_ex19_ayn,
                           geno = geno_ex19_ayn,
                           K = rel_mat,
                           n.PC = 3,
                           P3D = TRUE,
                           plot = TRUE)

#### Pheno ex19 pyn ####
pheno_ex19_pyn<- fread("./Results/BlupsGRYLD_ex19_pyn.txt") 

snpMatrix_ex19_pyn<- snpMatrix %>% 
  semi_join(pheno_ex19_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn) %>% 
  column_to_rownames(var = "rn") %>% 
  as.matrix()

pheno_ex19_pyn<- pheno_ex19_pyn %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

print(paste('Check that marker matrix and phenotypes align',  
            all(rownames(snpMatrix_ex19_pyn) == pheno_ex19_pyn$Variety), sep = ' '))

geno_ex19_pyn<-snpMatrix %>% 
  semi_join(pheno_ex19_pyn, by = c("rn" = "Variety")) %>% 
  arrange(rn)
geno_ex19_pyn<- as.data.frame(t(geno_ex19_pyn)) %>% 
  row_to_names(row_number = 1) 

write.table(geno_ex19_pyn,
            "./GenoDatabase/snpMatrix_ex19_pyn.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
geno_ex19_pyn<- fread("./GenoDatabase/snpMatrix_ex19_pyn.txt")

geno_ex19_pyn<- geno_ex19_pyn %>% 
  bind_cols(positions) %>% 
  select(snp, chrom, pos, everything())

rel_mat<- A.mat(snpMatrix_ex19_pyn)

gwas_gryld_ex19_pyn<- GWAS(pheno = pheno_ex19_pyn,
                           geno = geno_ex19_pyn,
                           K = rel_mat,
                           n.PC = 3,
                           P3D = TRUE,
                           plot = TRUE)

write.table(gwas_gryld_16,
            "./Results/gwas_GRYLD_16.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_16_ayn,
            "./Results/gwas_GRYLD_16_ayn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_16_pyn,
            "./Results/gwas_GRYLD_16_pyn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex16,
            "./Results/gwas_GRYLD_ex16.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex16_ayn,
            "./Results/gwas_GRYLD_ex16_ayn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex16_pyn,
            "./Results/gwas_GRYLD_ex16_pyn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_17,
            "./Results/gwas_GRYLD_17.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_17_ayn,
            "./Results/gwas_GRYLD_17_ayn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_17_pyn,
            "./Results/gwas_GRYLD_17_pyn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex17,
            "./Results/gwas_GRYLD_ex17.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex17_ayn,
            "./Results/gwas_GRYLD_ex17_ayn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex17_pyn,
            "./Results/gwas_GRYLD_ex17_pyn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_18,
            "./Results/gwas_GRYLD_18.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_18_ayn,
            "./Results/gwas_GRYLD_18_ayn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_18_pyn,
            "./Results/gwas_GRYLD_18_pyn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex18,
            "./Results/gwas_GRYLD_ex18.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex18_ayn,
            "./Results/gwas_GRYLD_ex18_ayn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex18_pyn,
            "./Results/gwas_GRYLD_ex18_pyn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_19,
            "./Results/gwas_GRYLD_19.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_19_ayn,
            "./Results/gwas_GRYLD_19_ayn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_19_pyn,
            "./Results/gwas_GRYLD_19_pyn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex19,
            "./Results/gwas_GRYLD_ex19.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex19_ayn,
            "./Results/gwas_GRYLD_ex19_ayn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(gwas_gryld_ex19_pyn,
            "./Results/gwas_GRYLD_ex19_pyn.txt",
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

#### Figures ####
don_rrB4_16 <- gwas_gryld_16 %>%
  # Compute chromosome size
  group_by(chrom) %>%
  summarise(chr_len = max(pos),
            .groups = "keep") %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwas_gryld_16, ., by = c("chrom" = "chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum = pos + tot)

axisdf <- don_rrB4_16 %>%
  group_by(chrom) %>%
  dplyr::summarize(
    center = (max(BPcum) + min(BPcum)) / 2,
    maximum = max(BPcum),
    .groups = "keep"
  )

ggplot(don_rrB4_16, aes(x = BPcum, 
                        y = GRYLD)) +
geom_point() +
  # Significance Threshold
  geom_hline(yintercept = -log10(0.05 / nrow(don_rrB4_16)), linetype = 2) +
  scale_x_continuous(label = axisdf$chrom, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0.05)) + 
  labs(
    title = "GWAS results GRYLD",
    subtitle = "Season 15/16, Bonferroni corection alpha = 0.05",
    x = "Chromosome",
    y = "-log10(P)"
  )


