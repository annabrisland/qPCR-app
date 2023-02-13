library("tidyverse")
library("ggplot2")

qPCRplot <- function(file, housekeeping, goi){
  
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
  
  analysis$replicate <- str_extract(analysis$Sample.Name, "\\d{1,3}$")
  analysis$Sample.Name <- str_replace(analysis$Sample.Name, "\\d{1,3}$", "")
  
  summary <- analysis %>%
    group_by(Sample.Name) %>%
    summarise(mean = mean(logFC), sd = sd(logFC))
  
  
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
    labs(y = "logFC")+
    scale_fill_manual(values =c("#E1999E","#7D0008","#4D0007","#A1CF9C","#117A07","#004D0E"))
  
  plot1
  
}