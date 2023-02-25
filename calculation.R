library("tidyverse")
library("dplyr")
library("ggplot2")
library("stringr")

calculateDDT <- function(file, housekeeping, goi) {
  
  file <- read.csv("data/bna5 iron 6hours_data_rox.csv", skip = 7)
  housekeeping <- "gapdh"
  goi <- "cft1"
  
  data <- file %>%
    select(Sample.Name, Target.Name, Cт) %>%
    na_if("") %>%
    na.omit() %>%
    mutate(Cт = na_if(Cт, "Undetermined")) %>%
    mutate(Cт = as.numeric(Cт))
  
  data$replicate <- str_extract(data$Sample.Name, "\\d{1,3}$")
  data$Sample.Name <- str_replace(data$Sample.Name, "\\d{1,3}$", "")
  
  ##dataset specific
  
  analysis <- data %>%
    filter(replicate == c(1,2,3)) %>%
    filter(Target.Name == housekeeping | Target.Name == goi)
  analysis$Sample.Name <- str_c(analysis$Sample.Name, '', analysis$replicate)
  
  analysis <- analysis %>%
    select(!replicate) %>%
    group_by(Sample.Name) %>%
    spread(Target.Name, Cт)
  colnames(analysis)[2:3] <- c("goiCT", "housekeepingCT")
  analysis$deltaCT <- analysis$goiCT - analysis$housekeepingCT
  control_average = mean(analysis$housekeepingCT)
  analysis$deltadeltaCT = analysis$deltaCT - control_average
  analysis$foldGeneExp = 2^-(analysis$deltadeltaCT)
  analysis$logFC = log(analysis$foldGeneExp)
  
  analysis
  
}


calculatesample <- function(file, housekeeping, goi, vtreated, vuntreated) {
  
  # file <- read.csv("data/bna5 iron 6hours_data_rox.csv", skip = 7)
  # housekeeping <- "gapdh"
  # goi <- c("cft1 cfo1")
  # vtreated <- "6h"
  # vuntreated <- "6b"
  
  genes <- unlist(strsplit(goi, split = "\\s+"))
  combined<- data.frame("x")
  
  for (i in genes){
    
    data <- file %>%
      select(Sample.Name, Target.Name, Cт) %>%
      na_if("") %>%
      na.omit() %>%
      mutate(Cт = na_if(Cт, "Undetermined")) %>%
      mutate(Cт = as.numeric(Cт))
    
    data$replicate <- str_extract(data$Sample.Name, "\\d{1,3}$")
    data$Sample.Name <- str_replace(data$Sample.Name, "\\d{1,3}$", "")
    
    analysis <- data %>%
      filter(replicate == c(1,2,3)) %>%
      filter(Target.Name == housekeeping | Target.Name == i) 
    
    treated <- filter(analysis, Sample.Name == vtreated)
    untreated <- filter(analysis, Sample.Name == vuntreated)
    
    treated$Sample.Name <- str_c(treated$Sample.Name, '', treated$replicate)
    untreated$Sample.Name <- str_c(untreated$Sample.Name, '', untreated$replicate)
    
    treated <- treated %>%
      select(!replicate) %>%
      group_by(Sample.Name) %>%
      spread(Target.Name, Cт)
    
    untreated <- untreated %>%
      select(!replicate) %>%
      group_by(Sample.Name) %>%
      spread(Target.Name, Cт)
    
    colnames(treated)[2:3] <- c("goiCT", "housekeepingCT")
    colnames(untreated)[2:3] <- c("goiCT", "housekeepingCT")
    
    treated$deltaCT <- treated$goiCT - treated$housekeepingCT
    untreated$deltaCT <- untreated$goiCT - untreated$housekeepingCT
    
    control_average = mean(untreated$deltaCT)
    treated$deltadeltaCT = treated$deltaCT - control_average
    treated$foldGeneExp = 2^-(treated$deltadeltaCT)
    treated$logFC = log(treated$foldGeneExp)
    
    untreated$deltadeltaCT = untreated$deltaCT - control_average
    untreated$foldGeneExp = 2^-(untreated$deltadeltaCT)
    untreated$logFC = log(untreated$foldGeneExp)
    
    final <- rbind(untreated, treated)
    colnames(final)[6] <- paste(i, "fold gene expression")
    
    final <- final[,c(1,6)]
    combined <- cbind(combined, final)
      
    
  }
  
  colnames(combined)[2] <- "Sample"
  combined <- combined %>%
    select(-c(X.x., Sample.Name))
  
  combined
  
}

