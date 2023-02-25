library("tidyverse")
library("ggplot2")

qPCRplot <- function(file, housekeeping, goi, treated, untreated){
  
  # housekeeping <- "gapdh"
  # goi <- "cft1"
  
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
    filter(Target.Name == housekeeping | Target.Name == goi) 
  
  treated <- filter(analysis, Sample.Name == treated)
  untreated <- filter(analysis, Sample.Name == untreated)
  
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
  
  
  final$replicate <- str_extract(final$Sample.Name, "\\d{1,3}$")
  final$Sample.Name <- str_replace(final$Sample.Name, "\\d{1,3}$", "")
  
  summary <- final %>%
    group_by(Sample.Name) %>%
    summarise(mean = mean(foldGeneExp), sd = sd(foldGeneExp))
  
  
  #### Plotting
  plot1 = ggplot(summary,aes(x=Sample.Name,y=mean, fill = Sample.Name))+
    geom_bar(stat= "identity", position = position_dodge(0.9), color = "black",size = 1.5)+
    geom_errorbar(aes(ymin = mean-sd,ymax = mean+sd),size =1.5,width =0.5,position = position_dodge(0.9))+
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