library("tidyverse")
library("dplyr")
library("ggplot2")
library("stringr")


calculatesample <- function(file, housekeeping, goi, vtreated, vuntreated) {
  
  # file <- read.csv("data/bna5 iron 6hours_data_rox.csv", skip = 7)
  # housekeeping <- "gapdh"
  # goi <- c("cft1 cfo1")
  # #i = "cft1"
  #vtreated <- "6b, 6f"
  # vuntreated <- "6b"

  genes <- unlist(strsplit(goi, split = "\\s+"))
  vtreated <- unlist(strsplit(vtreated, split = ",\\s+"))
  combined <- data.frame("x")
  treatments <- data.frame()
  
    
  for (i in genes){
      
      data <- file %>%
        select(Sample.Name, Target.Name, Cт) %>%
        mutate(across(where(is.character), ~na_if(., ""))) %>%
        mutate(across(where(is.numeric), ~na_if(., ""))) %>%
        na.omit() %>%
        mutate(Cт = na_if(Cт, "Undetermined")) %>%
        mutate(Cт = as.numeric(Cт))
      
      data$replicate <- str_extract(data$Sample.Name, "\\d{1,3}$")
      data$Sample.Name <- str_replace(data$Sample.Name, "\\d{1,3}$", "")
      
      analysis <- data %>%
        filter(replicate == c(1,2,3)) %>%
        filter(Target.Name == housekeeping | Target.Name == i)
      
      for(j in vtreated){
        
        treated <- filter(analysis, Sample.Name == j)
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
        
        treatments <- rbind(treatments, treated)
        
      }
      
      
      final <- rbind(untreated, treatments)
      colnames(final)[6] <- paste(i)
      
      final <- final[,c(1,6)]
      combined <- merge(combined, final)
      
      
    }
    
  combined <- combined %>%
    select(-X.x.)
  colnames(combined)[1] <- "Sample"

  combined
  
}


qPCRplot <- function(file, housekeeping, goi, vtreated, vuntreated){
  
  # file <- read.csv("data/bna5 iron 6hours_data_rox.csv", skip = 7)
  # housekeeping <- "gapdh"
  # goi <- c("cft1")
  # #i = "cft1"
  # vtreated <- "6h"
  # vuntreated <- "6b"
  
  final <- calculatesample(file, housekeeping, goi, vtreated, vuntreated) %>%
    gather(gene, foldGeneExp, -1)
  
  colnames(final)[1] <- "Sample"

  final$Sample <- str_replace(final$Sample, "\\d{1,3}$", "")
  
  summary <- final %>%
    group_by(Sample, gene) %>%
    summarise(mean = mean(foldGeneExp), sd = sd(foldGeneExp))
  
  
  #### Plotting
  plot1 = ggplot(summary,aes(x=gene,y=mean, fill = Sample))+
    geom_bar(stat= "identity", position = position_dodge(0.5), color = "black", width = 0.5, linewidth = 1.5)+
    geom_errorbar(aes(ymin = mean-sd,ymax = mean+sd),width =0.3, position = position_dodge(0.5), linewidth = 1.5)+
    theme_classic()+
    ggtitle("")+
    theme(axis.text.x = element_text(face = "bold", size = 20,angle = -45),
          axis.text.y = element_text(face = "bold", size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold", size = 20),
          legend.text = element_text(face = "bold",size = 12),
          title = element_text(face = "bold",size = 14))+
    labs(y = "Fold gene expression")+
    scale_fill_manual(values =c("#E1999E","#7D0008","#4D0007","#A1CF9C","#117A07","#004D0E"))
  
  plot1
  
}
