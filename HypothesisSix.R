rm(list = objects()); ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(broom)

##### Set up work space ####
#### Theme set
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
    legend.key.size = unit(1.75, "lines"),
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
    plot.subtitle = element_text(
      colour = "black",
      size = rel(2)
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

set.seed(1935)

##Load data 

pheno17<- fread("./PhenoDatabase/PhenoLong_vi17.txt")
pheno18<- fread("./PhenoDatabase/PhenoVI_18.txt")

#### linear regression with purrr ####

nested17<- pheno17 %>% 
  tidylog::select(-entity_id, -Variety, -range, 
                  -column, -Year, -Treated,-OverallTrial) %>% 
  group_by(Location, Trial, trait_id) %>% 
  nest() %>% 
  mutate(correlation = map(data, ~ lm(GRYLD ~ phenotypic_value, data = .x)),
         tidyCor = map(correlation,glance)) %>% 
  unnest(tidyCor) %>% 
  tidylog::select(-data,-correlation) %>% 
  separate(trait_id,c("trait_id","Date"), sep = "_") %>% 
  mutate(Date = as.Date(Date)) %>% 
  glimpse()

nested18<- pheno18 %>% 
  filter(trait_id != "PTHT") %>% 
  unite("trait_id",trait_id,phenotype_date, sep = "_") %>% 
  spread(key = trait_id, value = phenotype_value) %>% 
  gather(key = "trait_id",value = "phenotype_value", 
         `CanopyHeight_2017-12-01`:`GNDVI_2018-06-11`,
         `NDRE_2017-12-01`:`RE_2018-06-11`) %>% 
  drop_na(phenotype_value) %>% 
  group_by(Location, Trial, trait_id) %>% 
  nest() %>% 
  mutate(correlation = map(data, ~ lm(`GRYLD_2018-06-13` ~ phenotype_value,
                                      data = .x)),
         tidyCor = map(correlation,glance)) %>% 
  unnest(tidyCor) %>% 
  tidylog::select(-data,-correlation) %>% 
  separate(trait_id,c("trait_id","Date"), sep = "_") %>% 
  mutate(Date = as.Date(Date)) %>% 
  glimpse()

nested17 %>% 
  ggplot(aes(x = Date,
             y = adj.r.squared, 
             colour = Trial)) +
  geom_point() +
  facet_wrap(trait_id ~ Location, ncol = 4, scales = "free") +
  scale_x_date(breaks = "1 week", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(-1,1)) +
  theme(strip.text = element_text(colour = "black",
                                   size = rel(2))) +
  labs(title = "Linear regression for GRYLD 2016-2017",
       y = expression(paste(" adjusted for multiple tests r"^{2})))

nested18 %>% 
  ggplot(aes(x = Date,
             y = adj.r.squared, 
             colour = Trial)) +
  geom_point(size = 1.5) +
  facet_wrap(trait_id ~ Location, ncol = 4, scales = "free") +
  scale_x_date(breaks = "20 days", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(-1,1)) +
  theme(strip.text = element_text(colour = "black",
                                  size = rel(2))) +
  labs(title = "Linear regression for GRYLD 2017-2018",
       y = expression(paste(" adjusted for multiple tests r"^{2})))
