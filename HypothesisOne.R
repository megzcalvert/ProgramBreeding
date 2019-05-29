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

#save original data
getwd( ) #get working directory set if needed
setwd("~/Dropbox/Research_Poland_Lab/AM Panel")

saveRDS(pheno, "./Phenotype_Database/Pheno1718.RDS") 

dbDisconnect(wheatgenetics) #disconnect from database