rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(readxl)
library(psych)
library(broom)

## End of season phenotypes
set.seed(1973)

pheno<- fread("./PhenoDatabase/PhenoLongRep.txt")

## Vegetation indices
# 2016-2017 season
fileNames<- list.files(path = "./PhenoDatabase/2016_2017",
                       full.names = T)
traitNames<- basename(fileNames) %>%
  str_remove_all(c("traits_2016-2017_|_no_fills|.xlsx"))

read_excel_sheets<- function(fileNames, traitNames, ...) {
  traitNames<- fileNames %>% 
    excel_sheets() %>% 
    set_names() %>% 
    map_df(~ read_excel(path = fileNames, sheet = .x), .id = "trait_id") 
}

data17<- lapply(fileNames, read_excel_sheets)
names(data17)<- traitNames

data17<- plyr::ldply(data17, data.frame, .id = "Field") %>% 
  gather(key = "phenotype_date", value = "phenotype_value",X42880:X42894) %>%
  rename(entity_id = Plot_ID) %>% 
  drop_na(phenotype_value) %>% 
  tidylog::filter(entity_id != "Fill") %>% 
  tidylog::filter( str_detect(entity_id, "AYN|PYN")) %>% 
  glimpse()

# Because excel dates are stupid
data17$phenotype_date<- as.numeric(sub('.', '', data17$phenotype_date))
data17$phenotype_date<- as.Date(data17$phenotype_date,origin = "1899-12-30")

data17<- data17 %>% 
  select(-Field) %>% 
  mutate(Sep = entity_id) %>% 
  separate(Sep, c("Year","Trial","Location","Treated","Plot"), sep = "-") 

# Checking data is there and makes sense

data17 %>% 
  ggplot(aes(x = factor(phenotype_date), y = phenotype_value, colour = Trial)) +
  geom_boxplot() +
  facet_wrap(trait_id ~ Location, scales = "free", ncol = 5) +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) + 
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a',
                                '#984ea3','#ff7f00','#ffff33',
                                '#a65628','#f781bf','#999999'))

data17<- data17 %>% 
  unite("trait_id", trait_id,phenotype_date, sep = "_") %>% 
  spread(key = trait_id, value = phenotype_value) 

# 2017-2018 Season
data18<- tibble()
htpPheno<-c("GNDVI","NDRE","NDVI","Nir","RE","CanopyHeight")

load.file<- function (filename) {
  d<- fread(file = filename,header = TRUE,check.names = F,data.table = F)
  d
}

for(i in htpPheno) {
  fileNames<- list.files(path = "./PhenoDatabase/2017_2018",
                         full.names = T,
                         pattern = paste0(i,".csv"))
  
  traitNames<- basename(fileNames) %>%
    str_remove_all(c(".csv"))
  
  data<- lapply(fileNames,load.file)
  names(data)<- traitNames 
  
  data<- plyr::ldply(data, data.frame, .id = "Phenotype")
  print(colnames(data))
  data18<- bind_rows(data18,data)
  rm(data)
}

data18L<- data18 %>% 
  separate(Phenotype,c("Date","Loc1","Trait"), sep = "_") %>% 
  unite("Trait", Trait, Date) %>% 
  select(-CanopyHeight_top3,-SoilHeight,-X96p,-X97p,-X98p,-X99p,-X100p) %>% 
  gather(key = "VI", value = "phenotypic_value",GNDVI:CanopyHeight_top5) %>%
  distinct() %>% 
  select(-VI) %>% 
  drop_na(phenotypic_value) %>% 
  glimpse()

data18L$Loc1<- str_remove(data18L$Loc1,"HE")
data18L$Loc1<- str_remove(data18L$Loc1,"HW")

data18L<- data18L %>% 
  separate(Plot_ID, c("Plot","range","column"), sep = ":") %>% 
  glimpse() %>% 
  filter(!str_detect(Plot,"18GYP_SA")) %>% 
  filter(!str_detect(Plot,"18BEL_RP")) %>% 
  filter(!str_detect(Plot,"18MP_")) %>% 
  filter(!str_detect(Plot,"Fill")) %>% 
  filter(!str_detect(Plot,"SRPN")) %>% 
  filter(!str_detect(Plot,"CKE")) %>% 
  filter(!str_detect(Plot,"ENEV")) %>% 
  filter(!str_detect(Plot,"KIN"))

data18gg<- data18L %>% 
  separate(Trait, c("Trait","Date"), sep = "_") %>% 
  dplyr::mutate(Sep = Plot) %>% 
  separate(Sep, c("Year","Trial","Location","Treated","Plot"))

data18gg$Date<- as.Date(data18gg$Date, format = "%Y%m%d")
data18gg$phenotypic_value<- as.numeric(data18gg$phenotypic_value)

data18gg %>% 
  ggplot(aes(x = factor(Date), y = phenotypic_value, colour = Trial)) +
  geom_boxplot() +
  facet_wrap(Loc1~Trait, scales = "free")

## Joining to pheno column
pheno17<- pheno %>% 
  filter(Year == "17" & trait_id == "GRYLD") %>% 
  select(-ID,-phenotype_date,-phenotype_person) %>% 
  spread(key = trait_id, value = phenotype_value) %>% 
  drop_na(GRYLD) 

plots_missingVI<- pheno17 %>% 
  anti_join(data17, by = "entity_id")

pheno17W<- data17 %>% 
  select(-Year,-Trial,-Location,-Treated,-Plot) %>% 
  left_join(x = pheno17, by = "entity_id") 

### Correlation 
pheno17Matrix<- pheno17W %>% 
  select(11:75)

corMat17<- corr.test(as.matrix(pheno17Matrix), method = "pearson",
                     adjust = "holm")

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

flatCor17<- flattenCorrMatrix(corMat17$r,corMat17$p)

flatCor17<- flatCor17 %>% 
  drop_na() %>% 
  tidylog::filter(row == "GRYLD") %>% 
  separate(column, c("trait_id","Date"), sep = "_")

flatCor17$Date <- as.Date(flatCor17$Date)

flatCor17 %>% 
  ggplot(aes(x = Date, y = cor)) +
  geom_point() +
  facet_wrap(~trait_id, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  coord_cartesian(ylim = c(-1,1)) +
  labs(title = "Correlation between GRYLD and VI")


colnames(pheno17W)

#### Long format

pheno17L<- pheno17W %>% 
  gather(key = "trait_id",value = "phenotypic_value",12:76) %>% 
  drop_na(phenotypic_value)
colnames(pheno17L)

#write_delim(pheno17L,"./PhenoDatabase/PhenoLong_vi17.txt", delim = "\t")

nested17<- pheno17L %>% 
  filter(Location != "RP" 
         #& #Trial!= "GSPYN" 
         #& trait_id != "GNDVI_2017-05-24"
  ) %>% 
  tidylog::select(-entity_id, -Variety, -range, 
                  -column, -Year, -Treated,-OverallTrial) %>% 
  group_by(Location, Trial, trait_id) %>% 
  nest() %>% 
  mutate(correlation = map(data, ~ cor.test(.x$GRYLD,.x$phenotypic_value)),
         tidyCor = map(correlation,glance)) %>% 
  unnest(tidyCor) %>% 
  tidylog::select(-data,-correlation) %>% 
  separate(trait_id,c("trait_id","Date"), sep = "_")

nested17$Date <- as.Date(nested17$Date)

nested17 %>% 
  ggplot(aes(x = Date, y = estimate, colour = Trial)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(trait_id ~ Location, scales = "free", ncol = 3) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  coord_cartesian(ylim = c(-1,1)) +
  labs(title = "Correlation between GRYLD and VI 2017")





