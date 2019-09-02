rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(ggpubr)
library(broom)

#### Theme set
custom_theme<- theme_minimal() %+replace%
  theme(axis.title = element_text(colour = "black",
                                  size = rel(2)),
        axis.title.x = element_text(vjust = 0,
                                    margin = margin(t = 0, r = 0.25,
                                                    b = 0, l = 0,
                                                    unit = "cm")),
        axis.title.y = element_text(vjust = 1,
                                    angle = 90,
                                    margin = margin(t = 0, r = 0.25,
                                                    b = 0, l = 0.1,
                                                    unit = "cm")),
        axis.text = element_text(colour = "black",
                                 size = rel(1.5)),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(3,"pt"),
        axis.line = element_line(color = "black",
                                 size = 0.5),
        legend.key.size = unit(4,"lines"),
        # legend.background = element_rect(fill = NULL, colour = NULL),
        # legend.box = NULL,
        legend.margin = margin(t = 0, r = 0.75,
                               b = 0, l = 0.75,
                               unit = "cm"),
        legend.text = element_text(size = rel(2)),
        legend.title = element_text(size = rel(1.5)),
        panel.grid.major = element_line(colour = "#969696",
                                        linetype = 3),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = rel(2),
                                margin = margin(t = 0.1, r = 0.1, 
                                                b = 0.1, l = 0.1,
                                                unit = "cm")),
        plot.margin = margin(t = 0.5, r = 0.5, 
                             b = 0.5, l = 0,
                             unit = "cm"),
        plot.title = element_text(colour = "black",
                                  size = rel(3),
                                  vjust = 0,
                                  hjust = 0,
                                  margin = margin(t = 0.25, r = 0.25, 
                                                  b = 0.5, l = 0.25,
                                                  unit = "cm")),
        strip.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1),
        strip.text = element_text(colour = "black",
                                  size = rel(1)),
        complete = F)

theme_set(custom_theme)

#### Load pheno and blups data

blups17<- fread("./PhenoDatabase/BLUPs_2017.txt")
blups17<- blups17 %>% 
  clean_names()
blups18<- fread("./PhenoDatabase/BLUPs_2018.txt")
blupGryld<- fread("./PhenoDatabase/GRYLD_Blups_allyrs.txt")

pheno17<- fread("./PhenoDatabase/PhenoLong_vi17.txt")
pheno18<- fread("./PhenoDatabase/PhenoVI_18.txt")

traits17<- colnames(blups17[,2:ncol(blups17)])
traits18<- colnames(blups18[,2:ncol(blups18)])

pheno17<- pheno17 %>% 
  spread(key = trait_id, value = phenotypic_value) %>% 
  clean_names() %>% 
  glimpse()

pheno18<- pheno18 %>% 
  filter(trait_id != "PTHT") %>% 
  unite("trait_id", trait_id ,phenotype_date, sep = "_") %>% 
  spread(key = trait_id, value = phenotype_value) %>% 
  clean_names() %>% 
  glimpse()

#### Distributions

hist.plot <- function(traits, data, type, file, ...) {
  
  plotList = list()
  for (i in traits) {
    thisPlot <- ggplot(data = data, aes_string(x = i)) + 
      geom_histogram(colour="black", fill="white",
                     bins = 100) + 
      labs(title = paste0(type),
           x = i,
           y = "Frequency")
    
    plotList[[i]] = thisPlot
    ggsave(filename = paste0(i,"_",type,".pdf"),path = file)
    print(i)
    
  }
  return(plotList)
}

hist.plot(traits = traits17, data = pheno17, type = "Phenotypes2017",
          file = "./Figures/Distributions/")
hist.plot(traits = traits17, data = blups17, type = "BLUPs2017",
          file = "./Figures/Distributions/")
hist.plot(traits = traits18, data = pheno18, type = "Phenotypes2018",
          file = "./Figures/Distributions/")
hist.plot(traits = traits18, data = blups18, type = "BLUPs2018",
          file = "./Figures/Distributions/")

#### Selecting top 5% and drawing histogram

selectTop<- function(dat, traits, file, identifier, selecTrait, ...) {
  df <- data.frame(matrix(ncol = length(traits), 
                          nrow = nrow(dat) * 0.05 ))
  colnames(df) <- traits
  for (i in traits) {
    print(i)
    mT<- dat %>% 
      select(rn,paste(i),paste(selecTrait)) %>% 
      dplyr::rename(Trait = paste(i),
                    selecTrait = paste(selecTrait))
    topSelec<- dat %>% 
      select(rn,paste(i),paste(selecTrait)) %>% 
      dplyr::rename(Trait = paste(i),
                    selecTrait = paste(selecTrait)) %>% 
      top_frac(0.05, Trait)
    
    tidyT<-tidy(t.test(mT$selecTrait,topSelec$selecTrait))
    tidyP<-tidy(t.test(mT$Trait,topSelec$Trait))
    
    plot1<- ggplot(data = dat,
                   mapping = aes_string(i)) +
      geom_histogram(colour = "black",
                     fill = "white",
                     bins = 100) +
      geom_vline(xintercept = mean(mT$Trait), linetype = 2) +
      geom_histogram(data = topSelec, 
                     mapping = aes(Trait),
                     colour = "red",
                     bins = 100) +
      geom_vline(xintercept = mean(topSelec$Trait),
                 colour = "red", linetype = 2) +
      labs(subtitle = paste0("T-test for difference in mean p-value = ",
                             round(tidyP$p.value, 
                                   digits = 3)))
    
    plot2<- ggplot(data = dat,
                   mapping = aes_string(selecTrait)) +
      geom_histogram(colour = "black",
                     fill = "white",
                     bins = 100) + 
      geom_vline(xintercept = mean(mT$selecTrait), linetype = 2) +
      geom_histogram(data = topSelec, 
                     mapping = aes(selecTrait),
                     colour = "red",
                     bins = 100) +
      geom_vline(xintercept = mean(topSelec$selecTrait),
                 colour = "red", linetype = 2) +
      labs(subtitle = paste0("T-test for difference in mean p-value = ",
                             round(tidyT$p.value, 
                                   digits = 3)))
    
    plots<- ggarrange(plot1,plot2)
    ggexport(plots, filename = paste0(file,identifier,"_",i,".png"),
             width = 1000, height = 550)
    df[,paste(i)]<- topSelec$rn
  }
  return(df)
}

traits17<- colnames(blups17[,2:61])
top17<- selectTop(dat = blups17,traits = traits17, 
                  file = "./Figures/Selection/",
                  identifier = "Blups17",
                  selecTrait = "gryld")

blups18<- blups18 %>% 
  select(-ndvi_2017_12_01)

traits18<- colnames(blups18[,2:89])

top18<- selectTop(dat = blups18,traits = traits18, 
                  file = "./Figures/Selection/",
                  identifier = "Blups18",
                  selecTrait = "gryld_2018_06_13")

pheno17<- pheno17 %>% 
  rename(rn = variety) %>% 
  select(-gryld)
traits17<- colnames(pheno17[,11:ncol(pheno17)])

topPheno17<- selectTop(dat = pheno17,traits = traits17, 
                       file = "./Figures/Selection/",
                       identifier = "Pheno17")
pheno18<- pheno18 %>% 
  rename(rn = variety) 

topPheno18<- selectTop(dat = pheno18,traits = traits18, 
                       file = "./Figures/Selection/",
                       identifier = "Pheno18")


