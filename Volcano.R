library(readxl)
library(tidyverse)
library(ggrepel)
library(ggpubr)

##### Configuration ####

# this is the path to your data
dir <- "path/to/your/data"

# adjust FC and pval threshold if needed
settings <- list("FC" = c(0.3),
                 "pval" = c(0.01))


#### Main script #####
files <- list.files(dir)
for(i in 1:length(files)) {
  
  data <- read_excel(file.path(dir, files[i]), 
                     col_types = c("numeric"))
  
  for(j in 1:length(settings$FC)) {
    FC_ <- settings$FC[j]
    pval_ <- settings$pval[j]
    data_ready <- data %>%
      gather("key", "val", -mz) %>%
      separate(key, into = c("class", "rep")) %>%
      select(-rep) %>%
      group_by(mz, class) %>%
      summarise(int = list(val)) %>%
      group_by(mz) %>%
      spread(class, int) %>%
      mutate(pval = t.test(x = unlist(Treated),
                           y = unlist(Untreated))$p.value,
             meanTreat = mean(unlist(Treated)),
             meanUntreat = mean(unlist(Untreated)),
             logFC = log(meanTreat/meanUntreat, base = 2)) %>%
      ungroup() %>%
      mutate(padj = p.adjust(pval, method = "BH", n = length(mz)),
             sig = ifelse(padj < pval_, TRUE, FALSE),
             lab = as.numeric(ifelse(sig & ( logFC > FC_ | logFC < -FC_), mz, ""))) %>%
      select(-Treated, -Untreated)
    
    
    ggplot(data_ready, 
           aes(x = logFC, y = -log10(padj), 
               label = round(lab, digits = 1))) +
      geom_point(aes(col = sig), show.legend = FALSE) +
      geom_text_repel(min.segment.length = 0) +
      theme_classic2() +
      theme(text = element_text(size = 16)) +
      geom_vline(aes(xintercept = -FC_), linetype = "dotted") +
      geom_vline(aes(xintercept = FC_), linetype = "dotted") +
      geom_hline(aes(yintercept = -log10(pval_)), linetype = "dotted") +
      scale_x_continuous(limits =  c(-3,3), breaks = seq(-3, 3, 1)) +
      labs(x = "Log2(Fold-Change)",
           y = "-log10(adjusted p-value)",
           col = "p < 0.05") -> p
    
    filename <- file.path(getwd(), paste0("volc_", tools::file_path_sans_ext( files[i]), "_p", settings$pval[j], "_FC", settings$FC[j]))
    ggsave(plot = p, filename = paste0(filename, ".png"), dpi = 600)
    write_excel_csv(x = data_ready, path = paste0(filename, ".csv"))  
  }
}

