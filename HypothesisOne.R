rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(RMySQL)
library(janitor)

getwd()

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
plot.plot_name AS 'Variety',
plot.location,
plot.range,
plot.column,
plot.rep
FROM wheatgenetics.phenotype LEFT JOIN wheatgenetics.plot ON plot.plot_id = phenotype.entity_id 
WHERE wheatgenetics.phenotype.entity_id LIKE '%PYN-MP-%' OR
wheatgenetics.phenotype.entity_id LIKE '%PYN-RN-%' OR 
wheatgenetics.phenotype.entity_id LIKE '%PYN-RP-%' OR
wheatgenetics.phenotype.entity_id LIKE '%PYN-SA-%' OR
wheatgenetics.phenotype.entity_id LIKE '%PYN-TH-%' OR
wheatgenetics.phenotype.entity_id LIKE '%PYN-RL-%' OR 
wheatgenetics.phenotype.entity_id LIKE '%AYN%-RN-%' OR 
wheatgenetics.phenotype.entity_id LIKE '%AYN%-RL-%' OR 
wheatgenetics.phenotype.entity_id LIKE '%AYN%-RP-%' OR 
wheatgenetics.phenotype.entity_id LIKE '%AYN%-SA-%' OR 
wheatgenetics.phenotype.entity_id LIKE '%AYN%-MP-%' OR 
wheatgenetics.phenotype.entity_id LIKE '%AYN%-TH-%' ;"

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
  tidylog::filter(phenotype_value < 7000) %>% 
  group_by(Year, Trial, Location) %>% 
  glimpse()

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
                             c("AYN3","AYN2","AYN1","AYN4","DHAYN2",
                               "DHAYN1","PYN","GSPYN","DHPYN"), 
                             c("AYN","AYN","AYN","AYN","DHAYN",
                               "DHAYN","PYN","GSPYN","DHPYN"),
                     "NA", ignore.case = TRUE)
  
pheno_long$phenotype_value<- as.numeric(pheno_long$phenotype_value)
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


