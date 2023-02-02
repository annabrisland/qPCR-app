library("tidyverse")
library("dplyr")
library("ggplot2")
library("stringr")

calculateDDT <- function(housekeeping, goi) {
  
  data <- read.csv("data/bna5 iron 6hours_data_rox.csv", skip = 7) %>%
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
    filter(Target.Name == housekeeping | Target.Name == goi) #%>%
  #filter(Sample.Name == "6h")
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