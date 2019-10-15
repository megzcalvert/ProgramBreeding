rm(list = objects())
ls()

library(tidyverse)
library(data.table)
library(tidylog)
library(janitor)
library(ggpubr)
library(broom)

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
    legend.key.size = unit(4, "lines"),
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

#### Load pheno and blups data

blups17 <- fread("./PhenoDatabase/BLUPs_TrialColumnRange_2017.txt")
blups17 <- blups17 %>%
  unite("trait_id", location:phenotype_date, sep = "_") %>%
  tidylog::mutate(phenotype_value = as.numeric(phenotype_value)) %>% 
  pivot_wider(names_from = trait_id, values_from = phenotype_value) %>%
  clean_names() 

blups18 <- fread("./PhenoDatabase/BLUPs_TrialColumnRange_2018.txt")
blups18 <- blups18 %>%
  unite("trait_id", location:phenotype_date, sep = "_") %>%
  tidylog::mutate(phenotype_value = as.numeric(phenotype_value)) %>% 
  pivot_wider(names_from = trait_id, values_from = phenotype_value) %>%
  clean_names() 

blupGryld <- fread("./PhenoDatabase/GRYLD_Blups_allyrs.txt")

### raw phenotype data
pheno17 <- fread("./PhenoDatabase/PhenoLong_vi17.txt")
pheno18 <- fread("./PhenoDatabase/PhenoVI_18.txt")

traits17 <- colnames(blups17[, 2:ncol(blups17)])
traits18 <- colnames(blups18[, 2:ncol(blups18)])

pheno17 <- pheno17 %>%
  separate(trait_id, c("trait_id", "phenotype_date"), sep = "_") %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = format(phenotype_date, "%Y%m%d"),
    trait_id = str_replace(trait_id, "RedEdge", "rededge")
  ) %>%
  unite("trait_id", trait_id, phenotype_date, sep = "_") %>%
  unite("trait_id", Location, trait_id, sep = "_") %>%
  spread(key = trait_id, value = phenotypic_value) %>%
  clean_names()

pheno18 <- pheno18 %>%
  filter(trait_id != "PTHT") %>%
  mutate(
    phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"),
    phenotype_date = format(phenotype_date, "%Y%m%d")
  ) %>%
  unite("trait_id", trait_id, phenotype_date, sep = "_") %>%
  unite("trait_id", Location, trait_id, sep = "_") %>%
  spread(key = trait_id, value = phenotype_value) %>%
  clean_names()

#### Distributions

hist.plot <- function(traits, data, type, file, ...) {
  plotList <- list()
  for (i in traits) {
    thisPlot <- ggplot(data = data, aes_string(
      x = i
    )) +
      geom_histogram(
        colour = "black",
        fill = NA,
        bins = 100
      ) +
      labs(
        title = paste0(type),
        x = i,
        y = "Frequency"
      )
    
    plotList[[i]] <- thisPlot
    ggsave(filename = paste0(i, "_", type, ".pdf"), path = file)
    print(i)
  }
  return(plotList)
}


hist.plot(
  traits = traits17[1:63], data = pheno17, type = "Phenotypes2017",
  file = "./Figures/Distributions/"
)

hist.plot(
  traits = traits17, data = blups17, type = "BLUPs2017",
  file = "./Figures/Distributions/"
)

phenotTraits18 <- colnames(pheno18[, 10:114])

hist.plot(
  traits = phenotTraits18, data = pheno18, type = "Phenotypes2018",
  file = "./Figures/Distributions/"
)

hist.plot(
  traits = traits18, data = blups18, type = "BLUPs2018",
  file = "./Figures/Distributions/"
)

#### Selecting top #% and drawing histogram

selectTop <- function(dat, traits, file, identifier,
                      selecTrait, selectionIndex, ...) {
  df <- data.frame(matrix(
    ncol = length(traits),
    nrow = nrow(dat) * selectionIndex
  ))
  
  colnames(df) <- traits
  
  for (i in traits) {
    print(i)
    
    mT <- dat %>%
      select(variety, paste(i), paste(selecTrait)) %>%
      dplyr::rename(
        Trait = paste(i),
        selecTrait = paste(selecTrait)
      ) %>%
      drop_na(Trait) %>% 
      tidylog::mutate(Trait = as.numeric(Trait))
    
    topSelec <- dat %>%
      select(variety, paste(selecTrait), paste(i)) %>%
      dplyr::rename(
        Trait = paste(i),
        selecTrait = paste(selecTrait)
      ) %>%
      drop_na(Trait) %>% 
      arrange(desc(Trait))
    
    if (nrow(mT) > nrow(df)) {
      topSelec <- topSelec[1:nrow(df), ]
      
      mTafter <- mT %>%
        anti_join(topSelec, by = "variety")
      
      tidyT <- tidy(t.test(mTafter$selecTrait, topSelec$selecTrait))
      tidyP <- tidy(t.test(mT$Trait, topSelec$Trait))
      
      plot1 <- ggplot(
        data = mT,
        mapping = aes(Trait)
      ) +
        geom_histogram(
          colour = "black",
          fill = "white",
          bins = 100
        ) +
        geom_vline(
          xintercept = tidyP$estimate1,
          linetype = "solid"
        ) +
        geom_histogram(
          data = topSelec,
          mapping = aes(Trait),
          colour = "red",
          fill = "white",
          bins = 100
        ) +
        geom_vline(
          xintercept = tidyP$estimate2,
          colour = "red", linetype = "solid"
        ) +
        labs(subtitle = paste0(
          "T-test for difference in mean p-value = ",
          round(tidyP$p.value,
                digits = 3
          )
        ),
        x = paste0(i))
      
      plot2 <- ggplot(
        data = mTafter,
        mapping = aes(selecTrait)
      ) +
        geom_histogram(
          colour = "black",
          fill = "white",
          bins = 100
        ) +
        geom_vline(
          xintercept = mean(mTafter$selecTrait),
          linetype = "solid"
        ) +
        geom_histogram(
          data = topSelec,
          mapping = aes(selecTrait),
          colour = "red",
          fill = "white",
          bins = 100
        ) +
        geom_vline(
          xintercept = mean(topSelec$selecTrait),
          colour = "red", linetype = "solid"
        ) +
        labs(
          x = paste(selecTrait),
          subtitle = paste0(
            "T-test for difference in mean p-value = ",
            round(tidyT$p.value,
                  digits = 3
            )
          )
        )
      
      plots <- ggarrange(plot1, plot2)
      plots <- annotate_figure(plots,
                               top = text_grob(
                                 paste0(
                                   "Indirect selection on ", i,
                                   " with a selection index of ",
                                   selectionIndex
                                 )
                               )
      )
      print(plots)
      ggexport(plots,
               filename = paste0(file, identifier, "_", i,"_", 
                                 selectionIndex, ".png"),
               width = 1000, height = 550
      )
      df[, paste(i)] <- topSelec$variety
      dev.off()
    }
    else{
      print("Not enough observartions in field")
    }
    
  }
  return(df)
}

traits17 <- colnames(blups17[, c(2:64)])
top17_60 <- selectTop(
  dat = blups17, traits = traits17,
  file = "./Figures/Selection/",
  identifier = "Blups17_",
  selecTrait = "all_gryld_2017",
  selectionIndex = 0.6
)

top17_40 <- selectTop(
  dat = blups17, traits = traits17,
  file = "./Figures/Selection/",
  identifier = "Blups17_",
  selecTrait = "all_gryld_2017",
  selectionIndex = 0.4
)

traits18 <- colnames(blups18[, 2:102])

top18_60 <- selectTop(
  dat = blups18, traits = traits18,
  file = "./Figures/Selection/",
  identifier = "Blups18",
  selecTrait = "all_gryld_2018",
  selectionIndex = 0.6
)

top18_40 <- selectTop(
  dat = blups18, traits = traits18,
  file = "./Figures/Selection/",
  identifier = "Blups18",
  selecTrait = "all_gryld_2018",
  selectionIndex = 0.4
)

top17d_60 <- setDT(top17_60, keep.rownames = T)
top17d_60 <- top17d_60 %>%
  gather(
    key = "trait_id", value = "Variety",
    rl_rededge_20170605:mp_gndvi_20170506
  ) %>%
  dplyr::group_by(Variety) %>%
  summarise(count = n()) %>% 
  drop_na()

write.table(top17d_60, "./PhenoDatabase/Count_in_Selection_17_60.txt",
            quote = F,
            sep = "\t", row.names = F, col.names = T
)

top17d_40 <- setDT(top17_40, keep.rownames = T)
top17d_40 <- top17d_40 %>%
  gather(
    key = "trait_id", value = "Variety",
    rl_rededge_20170605:mp_gndvi_20170506
  ) %>%
  dplyr::group_by(Variety) %>%
  summarise(count = n()) %>% 
  drop_na()

write.table(top17d_40, "./PhenoDatabase/Count_in_Selection_17_40.txt",
            quote = F,
            sep = "\t", row.names = F, col.names = T
)

top18d_60 <- setDT(top18_60, keep.rownames = T)
top18d_60 <- top18d_60 %>%
  gather(
    key = "trait_id", value = "Variety",
    sa_re_20180522:mp_canopyheight_20180606
  ) %>%
  dplyr::group_by(Variety) %>%
  summarise(count = n()) %>% 
  drop_na()

write.table(top18d_60, "./PhenoDatabase/Count_in_Selection_18_60.txt",
            quote = F,
            sep = "\t", row.names = F, col.names = T
)

top18d_40 <- setDT(top18_40, keep.rownames = T)
top18d_40 <- top18d_40 %>%
  gather(
    key = "trait_id", value = "Variety",
    sa_re_20180522:mp_canopyheight_20180606
  ) %>%
  dplyr::group_by(Variety) %>%
  summarise(count = n()) %>% 
  drop_na()

write.table(top18d_40, "./PhenoDatabase/Count_in_Selection_18_40.txt",
            quote = F,
            sep = "\t", row.names = F, col.names = T
)

##### Comparing to Alan's Selection ####

lineSummary <- fread("./PhenoDatabase/VarietySummary.txt")

# Function to add a column based on a portion of text in another column
ff <- function(x, patterns, replacements = patterns, fill = NA, ...) {
  stopifnot(length(patterns) == length(replacements))
  
  ans <- rep_len(as.character(fill), length(x))
  empty <- seq_along(x)
  
  for (i in seq_along(patterns)) {
    greps <- grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] <- replacements[[i]]
    empty <- empty[!greps]
  }
  
  return(ans)
}

# Adding an overall trial column
lineSummary$OverallTrial <- ff(lineSummary$Trial,
                               c(
                                 "AYN1", "AYN2", "AYN3", "AYN4", "DHAYN1", "DHAYN2",
                                 "PYN", "PYN1", "PYN2", "PYNA", "PYNB",
                                 "GSPYN", "DHPYN", "DHPYN1", "DHPYN2"
                               ),
                               c(
                                 "AYN", "AYN", "AYN", "AYN", "DHAYN", "DHAYN",
                                 "PYN", "PYN", "PYN", "PYN", "PYN",
                                 "GSPYN", "DHPYN", "DHPYN", "DHPYN"
                               ),
                               "NA",
                               ignore.case = TRUE
)

lineSummary <- lineSummary %>% 
  filter(Variety != "CentralBlend") %>% 
  filter(Variety != "Zenda") %>% 
  filter(Variety != "Bob Dole") %>% 
  filter(Variety != "Gallagher") %>% 
  filter(Variety != "Joe") %>% 
  filter(Variety != "SYMonument") %>% 
  filter(Variety != "LCSChrome") %>% 
  filter(Variety != "Larry") %>% 
  filter(Variety != "Everest") %>% 
  filter(Variety != "FILL") %>% 
  filter(Variety != "KS031009K-6")

top17d_60 <- top17d_60 %>%
  left_join(lineSummary, by = "Variety") %>% 
  select(-n,-Trial, -Location) %>% 
  distinct() %>% 
  pivot_wider(names_from = Year, values_from = OverallTrial) %>% 
  select(-`NA`)

top17d_40 <- top17d_40 %>%
  left_join(lineSummary, by = "Variety") %>% 
  select(-n,-Trial, -Location) %>% 
  distinct() %>% 
  pivot_wider(names_from = Year, values_from = OverallTrial) %>% 
  select(-`NA`)

top18d_60 <- top18d_60 %>%
  left_join(lineSummary, by = "Variety")  %>% 
  select(-n,-Trial,-Location) %>% 
  distinct() %>% 
  pivot_wider(names_from = Year, values_from = OverallTrial) %>% 
  select(-`NA`)

top18d_40 <- top18d_40 %>%
  left_join(lineSummary, by = "Variety")  %>% 
  select(-n,-Trial,-Location) %>% 
  distinct() %>% 
  pivot_wider(names_from = Year, values_from = OverallTrial) %>% 
  select(-`NA`)

write.table(top17d_60, "./PhenoDatabase/Count_in_Selection_Alan_17_60.txt",
            quote = F,
            sep = "\t", row.names = F, col.names = T
)
write.table(top17d_40, "./PhenoDatabase/Count_in_Selection_Alan_17_40.txt",
            quote = F,
            sep = "\t", row.names = F, col.names = T
)
write.table(top18d_60, "./PhenoDatabase/Count_in_Selection_Alan_18_60.txt",
            quote = F,
            sep = "\t", row.names = F, col.names = T
)
write.table(top18d_40, "./PhenoDatabase/Count_in_Selection_Alan_18_40.txt",
            quote = F,
            sep = "\t", row.names = F, col.names = T
)

