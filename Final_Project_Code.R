#####################################################
### Study : BioInfo1_project
### Purpose of program : make GO analysis graph
### Output files : GO graph
### Programmed by : HGShin
### Draft Date : 2023-06-06
### Revision : NA
### Software (version) : R (v3.6.3)
###############################################################

library(dplyr)
read_dir = "/home2/hgshin/SNU_Bioinfo1_project/"
read_cnts <- read.table(paste0(read_dir,"read-counts_HS.txt"),sep = "\t", header = T, stringsAsFactors = F)

## siLuc -> control

read_cnts <- read_cnts %>% mutate(clip_enrichment = CLIP.35L33G.bam/RNA.control.bam) %>% 
  mutate(rden_change = (RPF.siLin28a.bam/RNA.siLin28a.bam)/(RPF.siLuc.bam/RNA.siLuc.bam))
read_cnts <- read_cnts %>% filter(RNA.control.bam > 30) %>% filter(RPF.siLuc.bam > 80) %>% filter(CLIP.35L33G.bam > 0)
read_cnts <- read_cnts %>% mutate(log2_clip_enrichment = log2(clip_enrichment)) %>% mutate(log2_rden_change = log2(rden_change))

summary(read_cnts)
# read_cnts <- read_cnts[!is.nan(read_cnts$log2_clip_enrichment) & !is.nan(read_cnts$log2_rden_change) & !is.infinite(read_cnts$log2_clip_enrichment) & !is.infinite(read_cnts$log2_rden_change),]
# head(read_cnts[!is.nan(read_cnts$log2_clip_enrichment) & !is.nan(read_cnts$log2_rden_change),])


library(ggplot2)
ggplot(data = read_cnts, aes(x = log2_clip_enrichment, y = log2_rden_change)) +
  geom_point(color = "black", size = 1, alpha = 0.2) +
  labs(x = "log2(clip_enrichment)", y = "log2(rden_change)") +
  xlim(-6,4) +
  ylim(-2,2)

cor(read_cnts$log2_clip_enrichment,read_cnts$log2_rden_change)


## Match Gene
library(org.Mm.eg.db)
library(AnnotationDbi)
read_Genes <- sapply(read_cnts$Geneid,function(X){
  str_split(X,"[.]")[[1]][1]
})

external_gene_names <- AnnotationDbi::select(org.Mm.eg.db, keys = read_Genes, keytype = "ENSEMBL", column = "SYMBOL")
GeneNames <- external_gene_names$SYMBOL %>% na.omit(.)

## GO Analysis & Data processed

library(enrichR)

GeneNames <- as.vector(GeneNames)
dbs <- c("GO_Biological_Process_2023","GO_Cellular_Component_2023")
enriched <- enrichr(GeneNames, dbs)


GO_CC <- enriched[[2]]
GO_CC[grep("GO:0031966",GO_CC$Term),]
GO_CC[grep("GO:0005737",GO_CC$Term),]

GO_CC$Overlap_num <- sapply(GO_CC$Overlap,function(X){
  str_split(X,"[/]")[[1]][1]
})

GO_CC$Overlap_num <- as.numeric(GO_CC$Overlap_num)
read_cnts$Geneid_revised <- read_Genes
Final_read_cnts <- read_cnts %>% dplyr::select("Geneid","Geneid_revised","log2_clip_enrichment","log2_rden_change")
Final_read_cnts <- merge(Final_read_cnts, external_gene_names, by.x = "Geneid_revised", by.y = "ENSEMBL", all.x = TRUE)

## Calculate average log2_clip_enrichment, log2_rden_change
Final_GO_CC <- GO_CC %>% dplyr::select("Term","Overlap_num","Adjusted.P.value","Genes")

mean_clip_enrichment <- vector("numeric", length = nrow(Final_GO_CC))
mean_rden_change <- vector("numeric", length = nrow(Final_GO_CC))


for (i in 1:nrow(Final_GO_CC)) {
  genes <- unlist(strsplit(Final_GO_CC$Genes[i], ";"))
  genes <- tolower(genes)
  genes <- gsub("\\s", "", genes)
  
  matching_genes <- intersect(tolower(Final_read_cnts$SYMBOL), genes)
  
  if (length(matching_genes) > 0) {
    clip_enrichment <- Final_read_cnts$log2_clip_enrichment[tolower(Final_read_cnts$SYMBOL) %in% matching_genes]
    mean_clip_enrichment[i] <- mean(clip_enrichment, na.rm = TRUE)
    
    rden_change <- Final_read_cnts$log2_rden_change[tolower(Final_read_cnts$SYMBOL) %in% matching_genes]
    mean_rden_change[i] <- mean(rden_change, na.rm = TRUE)
  }
}

Final_GO_CC$log2_clip_enrichment <- mean_clip_enrichment
Final_GO_CC$log2_rden_change <- mean_rden_change
Final_GO_CC <- Final_GO_CC %>% dplyr::select(-Genes)
Target_GO_CC <- Final_GO_CC[Final_GO_CC$Adjusted.P.value<0.05,]

## Make plot
library(ggplot2)

ggplot(Target_GO_CC, aes(x = (log2_clip_enrichment - 0.8), y = (log2_rden_change + 1))) +
  geom_point(data = subset(Target_GO_CC, Overlap_num <= 10), aes(size = Overlap_num, color = Adjusted.P.value)) +
  geom_point(data = subset(Target_GO_CC, Overlap_num > 10 & Overlap_num <= 100), aes(size = Overlap_num, color = Adjusted.P.value)) +
  geom_point(data = subset(Target_GO_CC, Overlap_num > 100 & Overlap_num <= 500), aes(size = Overlap_num, color = Adjusted.P.value)) +
  geom_point(data = subset(Target_GO_CC, Overlap_num > 500 & Overlap_num <= 1000), aes(size = Overlap_num, color = Adjusted.P.value)) +
  geom_point(data = subset(Target_GO_CC, Overlap_num > 1000), aes(size = Overlap_num, color = Adjusted.P.value)) +
  scale_size(range = c(0.5, 17)) +
  scale_color_gradientn(
    colors = c("darkred", "red", "orange", "yellow", "#FAEDB9", "beige"),
    values = c(0, 10^-80, 10^-30, 10^-15, 10^-10, 1),
    guide = guide_colorbar(
      title = "FDR",#"Term-specific enrichment confidence (false discovery rate)",
      breaks = c(0, 10^-80, 10^-30, 10^-15, 10^-10, 1),
      labels = c(expression(1e-120), expression(1e-80),
                 expression(1e-30), expression(1e-15), 
                 expression(1e-10), expression(1e-0))
    )
  ) +
  guides(size = "none") +
  labs(
    x = expression("Enrichment level of LIN28A-bound CLIP tags (log"[2] * ")"),
    y = expression("Ribosome density change upon" ~ italic(Lin28a) ~ "knockdown (log"[2] * ")"),
    title = "GO Enrichment Analysis",
    subtitle = "Cellular Components"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "italic"),
    axis.title = element_text(face = "italic"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "black", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(1.5, "cm")
  )
    
  

### Biological Process Ontology

GO_BP <- enriched[[1]]
GO_BP[grep("GO:0031966",GO_BP$Term),]
GO_BP[grep("GO:0005737",GO_BP$Term),]

GO_BP$Overlap_num <- sapply(GO_BP$Overlap,function(X){
  str_split(X,"[/]")[[1]][1]
})

GO_BP$Overlap_num <- as.numeric(GO_BP$Overlap_num)
read_cnts$Geneid_revised <- read_Genes
Final_read_cnts <- read_cnts %>% dplyr::select("Geneid","Geneid_revised","log2_clip_enrichment","log2_rden_change")
Final_read_cnts <- merge(Final_read_cnts, external_gene_names, by.x = "Geneid_revised", by.y = "ENSEMBL", all.x = TRUE)

## Calculate average log2_clip_enrichment, log2_rden_change
Final_GO_BP <- GO_BP %>% dplyr::select("Term","Overlap_num","Adjusted.P.value","Genes")

mean_clip_enrichment <- vector("numeric", length = nrow(Final_GO_BP))
mean_rden_change <- vector("numeric", length = nrow(Final_GO_BP))


for (i in 1:nrow(Final_GO_BP)) {
  genes <- unlist(strsplit(Final_GO_BP$Genes[i], ";"))
  genes <- tolower(genes)
  genes <- gsub("\\s", "", genes)
  
  matching_genes <- intersect(tolower(Final_read_cnts$SYMBOL), genes)
  
  if (length(matching_genes) > 0) {
    clip_enrichment <- Final_read_cnts$log2_clip_enrichment[tolower(Final_read_cnts$SYMBOL) %in% matching_genes]
    mean_clip_enrichment[i] <- mean(clip_enrichment, na.rm = TRUE)
    
    rden_change <- Final_read_cnts$log2_rden_change[tolower(Final_read_cnts$SYMBOL) %in% matching_genes]
    mean_rden_change[i] <- mean(rden_change, na.rm = TRUE)
  }
}

Final_GO_BP$log2_clip_enrichment <- mean_clip_enrichment
Final_GO_BP$log2_rden_change <- mean_rden_change
Final_GO_BP <- Final_GO_BP %>% dplyr::select(-Genes)
Target_GO_BP <- Final_GO_BP[Final_GO_BP$Adjusted.P.value<0.05,]
Target_GO_BP <- Final_GO_BP[Final_GO_BP$Adjusted.P.value<0.0000005,]

ggplot(Target_GO_BP, aes(x = (log2_clip_enrichment - 0.8), y = (log2_rden_change + 1))) +
  geom_point(data = subset(Target_GO_BP, Overlap_num <= 10), aes(size = Overlap_num, color = Adjusted.P.value)) +
  geom_point(data = subset(Target_GO_BP, Overlap_num > 10 & Overlap_num <= 50), aes(size = Overlap_num, color = Adjusted.P.value)) +
  geom_point(data = subset(Target_GO_BP, Overlap_num > 50 & Overlap_num <= 100), aes(size = Overlap_num, color = Adjusted.P.value)) +
  geom_point(data = subset(Target_GO_BP, Overlap_num > 100 & Overlap_num <= 200), aes(size = Overlap_num, color = Adjusted.P.value)) +
  geom_point(data = subset(Target_GO_BP, Overlap_num > 200), aes(size = Overlap_num, color = Adjusted.P.value)) +
  scale_size(range = c(0.5, 17)) +
  scale_color_gradientn(
    colors = c("darkred", "red", "orange", "yellow", "#FAEDB9", "beige"),
    values = c(0, 10^-40, 10^-30, 10^-20, 10^-10, 1),
    guide = guide_colorbar(
      title = "FDR",#"Term-specific enrichment confidence (false discovery rate)",
      breaks = c(0, 10^-40, 10^-30, 10^-20, 10^-10, 1),
      labels = c(expression(1e-50), expression(1e-40),
                 expression(1e-30), expression(1e-20), 
                 expression(1e-10), expression(1e-0))
    )
  ) +
  guides(size = "none") +
  labs(
    x = expression("Enrichment level of LIN28A-bound CLIP tags (log"[2] * ")"),
    y = expression("Ribosome density change upon" ~ italic(Lin28a) ~ "knockdown (log"[2] * ")"),
    title = "GO Enrichment Analysis",
    subtitle = "Biological Process"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "italic"),
    axis.title = element_text(face = "italic"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "black", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(1.5, "cm") 
  )
