rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(broom)
library(asreml)

set.seed(1990)
asreml.license.status()

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
        plot.title = element_text(colour = "#512888",
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

#### Read in data
pheno17<- fread("./PhenoDatabase/PhenoLong_vi17.txt")

pheno18<- fread("./PhenoDatabase/PhenoVI_18.txt")

#### Move from long to wide 
pheno17<- pheno17 %>% 
  spread(key = trait_id, value = phenotypic_value) %>% 
  mutate(Location = as.factor(Location),
         Variety = as.factor(Variety),
         Trial = as.factor(Trial),
         rep = as.factor(rep),
         range = as.factor(range),
         column = as.factor(column)) %>% 
  glimpse()

pheno18<- pheno18 %>% 
  filter(trait_id != "PTHT") %>% 
  unite("trait_id", trait_id, phenotype_date) %>% 
  spread(key = trait_id, value = phenotype_value) %>% 
  mutate(Location = as.factor(Location),
         Variety = as.factor(Variety),
         Trial = as.factor(Trial),
         rep = as.factor(rep),
         range = as.factor(range),
         column = as.factor(column)) %>% 
  glimpse()

#### Asreml trial

t17<- asreml(fixed = GRYLD ~ 1,
             random = ~ Variety 
             + Location 
             + Location:Trial 
             + Location:Trial:range
             + Location:Trial:column
             ,
             data = pheno17) 

plot(t17)
summary.asreml(t17)

blues<- setDT(as.data.frame(coef(t17)$random), keep.rownames = T)
blues$rn<- str_remove(blues$rn,"Variety_")
dat17<- blues %>% 
  rename(GRYLD = effect) %>% 
  glimpse() 

colnames(pheno17)
# pheno17<- pheno17 %>% 
#   select(-`Nir_2017-06-16`,-`RedEdge_2017-05-24`,-`RedEdge_2017-06-16`)

##### Generating all 2017 VI BLUEs

effectvars <- names(pheno17) %in% c("entity_id","Variety","Year","Trial",
                                    "Location","Treated","rep","range","column",
                                    "OverallTrial","GRYLD")
traits <- colnames(pheno17[ , !effectvars])
traits
fieldInfo<- pheno17 %>% 
  tidylog::select(Variety, Location, Trial, range, column)



for (i in traits) {
  print(paste("Working on trait", i))
  
  data<- cbind(fieldInfo, pheno17[,paste(i)])
  names(data)<- c("Variety","Location","Trial","range","column","Trait")
  print(colnames(data))
  
  t17<- asreml(fixed = Trait ~ 1,
               random = ~ Variety + Trial + column,
               data = data) 
  
  pdf(paste0("./Figures/AsremlPlots/ASREML_repBlock17_",i,".pdf"))
  plot(t17)
  dev.off()
  print(summary(t17))
  blues<- setDT(as.data.frame(coef(t17)$random), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat17<- blues %>% 
    inner_join(dat17)
  
}

warnings()

#### Asreml trial
pheno18<- pheno18 %>% 
  clean_names()

t18<- asreml(fixed = gryld_2018_06_13 ~ 1,
             random = ~ variety 
             + location 
             + location:trial 
             + location:trial:range
             + location:trial:column
             ,
             data = pheno18) 

plot(t18)
summary.asreml(t18)

blues<- setDT(as.data.frame(coef(t18)$random), keep.rownames = T)
blues$rn<- str_remove(blues$rn,"variety_")
dat18<- blues %>% 
  rename(gryld_2018_06_13 = effect) %>% 
  glimpse() 

colnames(pheno18)

##### Generating all 2017 VI BLUEs

effectvars <- names(pheno18) %in% c("entity_id","variety","year","trial",
                                    "location","treated","rep","range","column",
                                    "overall_trial","gryld_2018_06_13")
traits <- colnames(pheno18[ , !effectvars])
traits
fieldInfo<- pheno18 %>% 
  tidylog::select(variety, location, trial, range, column)



for (i in traits) {
  print(paste("Working on trait", i))
  
  data<- cbind(fieldInfo, pheno18[,paste(i)])
  names(data)<- c("Variety","Location","Trial","range","column","Trait")
  print(colnames(data))
  
  t17<- asreml(fixed = Trait ~ 1,
               random = ~ Variety + Trial + column,
               data = data) 
  
  pdf(paste0("./Figures/AsremlPlots/ASREML_repBlock18_",i,".pdf"))
  plot(t17)
  dev.off()
  print(summary(t17))
  blues<- setDT(as.data.frame(coef(t17)$random), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat18<- blues %>% 
    inner_join(dat18)
  
}

warnings()

write.table(dat17,"./PhenoDatabase/BLUPs_2017.txt", quote = F, sep = "\t",
            row.names = F, col.names = T)
write.table(dat18,"./PhenoDatabase/BLUPs_2018.txt", quote = F, sep = "\t",
            row.names = F, col.names = T)

#### PCA of phenotypes
dat17pca<- as.matrix(dat17[,-1])

pcaMethods::checkData(data = dat17pca)

pcaDat17<- prcomp(dat17pca)
summary(pcaDat17)

scores<- as.data.frame(pcaDat17$x)
ggplot(data = scores, aes(x = PC1,y = PC2)) +
  geom_point()
biplot(pcaDat17)

ggbiplot2(pcaDat17)

dat18pca<- as.matrix(dat18[,-1])

pcaMethods::checkData(data = dat18pca)

pcaDat18<- prcomp(dat18pca)
summary(pcaDat18)

scores<- as.data.frame(pcaDat18$x)
ggplot(data = scores, aes(x = PC1,y = PC2)) +
  geom_point()
biplot(pcaDat18)

ggbiplot2(pcaDat18)

##### GRYLD all years 

gryld17<- pheno17 %>% 
  select(entity_id:GRYLD) %>% 
  clean_names()
gryld18<- pheno18 %>% 
  select(location:overall_trial,gryld_2018_06_13) %>% 
  dplyr::rename(gryld = gryld_2018_06_13)

gryld<- gryld17 %>% 
  bind_rows(gryld18) %>% 
  mutate(variety = as.factor(variety),
         year = as.factor(year),
         location = as.factor(location),
         trial = as.factor(trial),
         range = as.factor(range),
         column = as.factor(column))

gryldBlup<- asreml(fixed = gryld ~ 1,
                   random = ~ variety 
                   + year 
                   + year:location
                   + year:location:trial
                   ,
                   data = gryld)
plot(gryldBlup)
summary(gryldBlup)

gryldBLUP<- setDT(as.data.frame(coef(gryldBlup)$random), keep.rownames = T)
gryldBLUP<- gryldBLUP %>% 
  filter(str_detect(rn,"variety_")) %>% 
  tidylog::mutate(rn = str_remove(rn,"variety_"))

write.table(gryldBLUP,"./PhenoDatabase/GRYLD_Blups_allyrs.txt", quote = F,
            sep = "\t", row.names = F,col.names = T)

