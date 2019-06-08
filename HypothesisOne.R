rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(RMySQL)
library(janitor)
library(broom)
library(reshape2)

getwd()
set.seed(1964)

wheatgenetics = dbConnect(MySQL( ),
                          user=rstudioapi::askForPassword("Database user"),
                          dbname='wheatgenetics', host='beocat.cis.ksu.edu',
                          password = 
                            rstudioapi::askForPassword("Database password"),
                          port = 6306) 

#SQL Query to get all AM Panel phenotype data

pheno_query <- "SELECT phenotype.entity_id,
phenotype.trait_id,
phenotype.phenotype_value,
phenotype.phenotype_date,
phenotype.phenotype_person,
plot.plot_name AS 'Variety',
plot.location,
plot.range,
plot.column,
plot.rep
FROM wheatgenetics.phenotype LEFT JOIN 
wheatgenetics.plot ON plot.plot_id = phenotype.entity_id 
WHERE wheatgenetics.phenotype.entity_id LIKE '%YN-MP-%' OR
wheatgenetics.phenotype.entity_id LIKE '%YN%-RN-%' OR 
wheatgenetics.phenotype.entity_id LIKE '%YN%-RP-%' OR
wheatgenetics.phenotype.entity_id LIKE '%YN%-SA-%' OR
wheatgenetics.phenotype.entity_id LIKE '%YN%-TH-%' OR
wheatgenetics.phenotype.entity_id LIKE '%YN%-RL-%' ;"

#run the query to get plot information

pheno <- dbGetQuery(wheatgenetics, pheno_query)
str(pheno)

#save original data
getwd( ) #get working directory set if needed

saveRDS(pheno, "./PhenoDatabase/Pheno.RDS") 

dbDisconnect(wheatgenetics) #disconnect from database

sessionInfo() # useful infos **reproducible research**
rm(wheatgenetics, pheno_query, pheno)

#Read in original data
pheno_long<- readRDS("./PhenoDatabase/Pheno.RDS")

glimpse(pheno_long)

## Divide out all important info from the entity_id
# remove awns and pcthead as they were taken sporadically

pheno_long<- pheno_long %>% 
  mutate(Sep = entity_id) %>% 
  separate(Sep, c("Year","Trial","Location","Treated","Plot"), sep = "-") %>% 
  tidylog::filter(trait_id != "AWNS") %>% 
  tidylog::filter(trait_id != "PCTHEAD") %>% 
  tidylog::filter(phenotype_value >= 0) %>% 
  tidylog::filter(phenotype_value < 6500) %>% 
  tidylog::select(-location,-rep,-Plot) %>% 
  mutate(ID = row_number()) %>% 
  glimpse()

unique(pheno_long$Trial)

#Function to add a column based on a portion of text in another column
ff = function(x, patterns, replacements = patterns, fill = NA, ...)
{
  stopifnot(length(patterns) == length(replacements))
  
  ans = rep_len(as.character(fill), length(x))    
  empty = seq_along(x)
  
  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]  
    empty = empty[!greps]
  }
  
  return(ans)
}

#Adding an overall trial column
pheno_long$OverallTrial<- ff(pheno_long$Trial, 
                             c("AYN1","AYN2","AYN3","AYN4","DHAYN1","DHAYN2",
                               "PYN","PYN1","PYN2","PYNA","PYNB",
                               "GSPYN","DHPYN","DHPYN1","DHPYN2"), 
                             c("AYN","AYN","AYN","AYN","DHAYN","DHAYN",
                               "PYN","PYN","PYN","PYN","PYN",
                               "GSPYN","DHPYN","DHPYN","DHPYN"),
                             "NA", ignore.case = TRUE)

pheno_long$phenotype_value<- as.numeric(pheno_long$phenotype_value)
pheno_long$Location<- as.factor(pheno_long$Location)
pheno_long$Trial<- as.factor(pheno_long$Trial)
pheno_long$Year<- as.factor(pheno_long$Year)
traits<- unique(pheno_long$trait_id)

for (i in traits) {
  p<- pheno_long %>% 
    tidylog::filter(trait_id == paste(i)) %>% 
    ggplot(aes(x = phenotype_value, colour = Location)) +
    geom_density() +
    facet_wrap(Year ~ Trial, scales = "free") +
    theme_bw() +
    labs(title = paste(i))
  print(p)
}


write.table(pheno_long, "./PhenoDatabase/PhenoLong.txt", quote = F, sep = "\t",
            col.names = T, row.names = F)

# separating into traits
unique(pheno_long$trait_id)
phenoGryld<- pheno_long %>%    
  tidylog::filter(trait_id == "GRYLD")
phenoPtht<- pheno_long %>%    
  tidylog::filter(trait_id == "PTHT")
phenoMoist<- pheno_long %>%    
  tidylog::filter(trait_id == "MOIST")
phenoTestwt<- pheno_long %>%    
  tidylog::filter(trait_id == "TESTWT")

gryld<- phenoGryld %>%
  group_by(Year) %>%
  do(tidy(anova(lm(phenotype_value ~ Trial * Location, data = .))))

ptht<- phenoPtht %>%
  group_by(Year) %>%
  do(tidy(anova(lm(phenotype_value ~ Trial * Location, data = .))))

moist<- phenoMoist %>%
  group_by(Year) %>%
  do(tidy(anova(lm(phenotype_value ~ Trial * Location, data = .))))

testwt<- phenoTestwt %>%
  group_by(Year) %>%
  do(tidy(anova(lm(phenotype_value ~ Trial * Location, data = .))))

pheno_longRep<- pheno_long %>% 
  tidylog::select(-trait_id,-phenotype_value,-phenotype_date,
                  -phenotype_person, -ID,-range,-column,-OverallTrial) %>% 
  group_by(Year, Location, Trial, Treated, Variety) %>% 
  distinct() %>% 
  mutate(rep = row_number()) %>% 
  inner_join(pheno_long, by = c("entity_id","Variety","Year","Trial",
                                "Location","Treated")) %>%
  select(ID,entity_id,Variety,Year,Trial,
         Location,Treated,rep,range,column,
         OverallTrial,trait_id,phenotype_value,
         phenotype_date,phenotype_person) %>% 
  glimpse 

write.table(pheno_longRep, "./PhenoDatabase/PhenoLongRep.txt", quote = F, sep = "\t",
            col.names = T, row.names = F)
