## Import and load libraries needed
library(readr)
library(tibble)
library(tidyverse)
library(reshape2)
library(ggforce)
library("FactoMineR")
library("factoextra")
library("Hmisc")
library(corrplot)
library(gplots)
library(readxl)


setwd("/Volumes/GoogleDrive/My Drive/Anglian project/PCA and correlation/All cambi/")

cambi_rel <- read_csv("cambi_997_gtdb_relative_counts.csv")

## Transpose table
cambi_rel_t <- t(cambi_rel)

## Rename 
cambi_t2 <- cambi_rel_t

## Set colnames
colnames(cambi_t2) <- cambi_rel_t[1,]

## Remove first row now set to colnames
cambi_t2 <- cambi_t2[-1,]

## Convert matrix back into dataframe 
cambi_t2 = as.data.frame(cambi_t2)

## Add in Sample name
cambi_t2 <- add_column(cambi_t2, Sample_name = "NA", .before=1)

## Add in date
cambi_t2 <- add_column(cambi_t2, Date = "NA", .before=2)

## Add in date as a string to facet on plots
cambi_t2 <- add_column(cambi_t2, Date_2 = "NA", .before=3)

## Add in conversion rate 
cambi_t2 <- add_column(cambi_t2, "MWh/tDS" = "NA", .before=4)

## Add in site 
cambi_t2 <- add_column(cambi_t2, Site = "NA", .before=5)

## Add in method - all sites for this are cambi 
cambi_t2$Method = "CAMBI"

## Add in group values based on sample names - only whitlingham and cottonvalley in here 
cambi_t2 <- mutate(cambi_t2, Site = case_when(
  str_detect(rownames(cambi_t2), "Cotton_Valley") ~ "Cotton_Valley",
  str_detect(rownames(cambi_t2), "Whitlingham") ~ "Whitlingham"))

## Set date as string 
cambi_t2 <- mutate(cambi_t2, Date_2 = case_when(
  str_detect(rownames(cambi_t2), "20201014") ~ "20201014",
  str_detect(rownames(cambi_t2), "20201021") ~ "20201021",
  str_detect(rownames(cambi_t2), "20201028") ~ "20201028",
  str_detect(rownames(cambi_t2), "2020114") ~ "2020114",
  str_detect(rownames(cambi_t2), "20201111") ~ "20201111",
  str_detect(rownames(cambi_t2), "20201118") ~ "20201118",
  str_detect(rownames(cambi_t2), "20201125") ~ "20201125",
  str_detect(rownames(cambi_t2), "2020122") ~ "2020122",
  str_detect(rownames(cambi_t2), "20201112") ~ "20201111"))

## Set Date 
cambi_t2 <- mutate(cambi_t2, Date = case_when(
  str_detect(rownames(cambi_t2), "20201014") ~ "2020-10-14",
  str_detect(rownames(cambi_t2), "20201021") ~ "2020-10-21",
  str_detect(rownames(cambi_t2), "20201028") ~ "2020-10-28",
  str_detect(rownames(cambi_t2), "2020114") ~ "2020-11-04",
  str_detect(rownames(cambi_t2), "20201111") ~ "2020-11-11",
  str_detect(rownames(cambi_t2), "20201118") ~ "2020-11-18",
  str_detect(rownames(cambi_t2), "20201125") ~ "2020-11-25",
  str_detect(rownames(cambi_t2), "2020122") ~ "2020-12-02",
  str_detect(rownames(cambi_t2), "20201112") ~ "2020-11-11"))

## Set to date 
cambi_t2$Date = as.Date(cambi_t2$Date)

## Set sample name 
cambi_t2$Sample_name = paste(cambi_t2$Site, cambi_t2$Date, sep = '_')

## Shorten sample name 
cambi_t2$Sample_name <- str_replace(cambi_t2$Sample_name, "Whitlingham", "W")
cambi_t2$Sample_name <- str_replace(cambi_t2$Sample_name, "Cotton_Valley", "CV")

## Read in operational info 
op_data <- read_csv("../../Operational params/all_sites_averaged_columns_kings_lynn_both_corrected_cut.csv")

## Set id for merging tables 
op_data <- add_column(op_data, Sample_name = "NA", .before=1)

op_data$Sample_name = paste(op_data$Site, op_data$Date, sep = '_')

op_data$Sample_name <- str_replace(op_data$Sample_name, "Whitlingham", "W")
op_data$Sample_name <- str_replace(op_data$Sample_name, "Cotton_valley", "CV")
op_data$Sample_name <- str_replace(op_data$Sample_name, "Basildon", "B")
op_data$Sample_name <- str_replace(op_data$Sample_name, "Cambridge", "Ca")
op_data$Sample_name <- str_replace(op_data$Sample_name, "Cliff_quay", "CQ")
op_data$Sample_name <- str_replace(op_data$Sample_name, "Great_billing", "GB")
op_data$Sample_name <- str_replace(op_data$Sample_name, "Colchester", "Co")
op_data$Sample_name <- str_replace(op_data$Sample_name, "Kings_lynn", "KL")
op_data$Sample_name <- str_replace(op_data$Sample_name, "Pyewipe", "P")

cambi_t2 <- add_column(cambi_t2, ID = "NA", .before=1)

cambi_t2$ID = row.names(cambi_t2)

cambi_t2$ID <- str_replace(cambi_t2$ID, "Whitlingham", "W")
cambi_t2$ID <- str_replace(cambi_t2$ID, "Cotton_Valley", "CV")
cambi_t2$ID <- str_replace(cambi_t2$ID, "Basildon", "B")
cambi_t2$ID <- str_replace(cambi_t2$ID, "Cambridge", "Ca")
cambi_t2$ID <- str_replace(cambi_t2$ID, "Cliff_quay", "CQ")
cambi_t2$ID <- str_replace(cambi_t2$ID, "Great_billing", "GB")
cambi_t2$ID <- str_replace(cambi_t2$ID, "Colchester", "Co")
cambi_t2$ID <- str_replace(cambi_t2$ID, "Kings_lynn", "KL")
cambi_t2$ID <- str_replace(cambi_t2$ID, "Pyewipe", "P")

cambi_t2$`MWh/tDS` = NULL

cambi_t2 <- add_column(cambi_t2, Stage = "NA", .before=2)

cambi_t2 <- mutate(cambi_t2, Stage = case_when(
  str_detect(rownames(cambi_t2), "pre") ~ "Pre",
  str_detect(rownames(cambi_t2), "hyd") ~ "Hyd",
  str_detect(rownames(cambi_t2), "dig") ~ "Dig",
  str_detect(rownames(cambi_t2), "post") ~ "Post",
  str_detect(rownames(cambi_t2), "_x_") ~ "X"))



## Merge, to have conversion rate, drop redundant columns 
merge_op_data <- merge(op_data, cambi_t2, by ="Sample_name", suffixes = c("",".y"))

## NA values for PCA so remove column 
merge_op_data$Average_alkalinity = NULL 

merge_op_data <- merge_op_data %>% remove_rownames %>% column_to_rownames(var="ID")

merge_pca <- merge_op_data

merge_pca <- merge_pca[-c(1:3,14:18)]

merge_pca = merge_pca %>% mutate_if(is.character,as.numeric)

res.pca <- PCA(merge_pca, scale.unit =TRUE, graph = TRUE, ncp=5)

eig.val <- get_eigenvalue(res.pca)

## Make screeplot 
pdf("Screeplot_cambi_only.pdf")
fviz_eig(res.pca, addlabels = TRUE)
dev.off()

### Get var pca values and save to object 
var <- get_pca_var(res.pca)

## Get individuals pca values and save to object 
ind <- get_pca_ind(res.pca)

## Default black pca plot 
pdf("pca_var.pdf")
fviz_pca_var(res.pca, col.var = "black")
dev.off()

pdf("pca_ind_all_variables_taxon_by_site_cambi_only.pdf")
fviz_pca_ind(res.pca, pointshape= 21, geom.ind ="point", fill=merge_op_data$Site, addEllipses=TRUE, col.ind= merge_op_data$Site,  legend.title = "Site", mean.point= FALSE, pointsize=3)
dev.off()

pdf("pca_ind_all_variables_taxon_by_date_cambi_only.pdf")
fviz_pca_ind(res.pca, pointshape= 21, geom.ind ="point", fill=merge_op_data$Date_2, addEllipses=TRUE, col.ind= merge_op_data$Date_2,  legend.title = "Date", mean.point= FALSE, pointsize=3)
dev.off()

pdf("pca_ind_all_variables_taxon_by_stage_cambi_only.pdf")
fviz_pca_ind(res.pca, pointshape= 21, geom.ind ="point", fill=merge_op_data$Stage, addEllipses=TRUE, col.ind= merge_op_data$Stage,  legend.title = "Date", mean.point= FALSE, pointsize=3)
dev.off()


pdf("pca_var_cos2_cambi_only.pdf")
fviz_pca_var(res.pca, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, col.var = "cos2")
dev.off()

pdf("pca_var_cos2_top50_cambi_only.pdf")
fviz_pca_var(res.pca, col.var = "cos2", 
             select.var = list(contrib = 50),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping 
)
dev.off()

pdf("pca_ind_cambi_only.pdf")
fviz_pca_ind(res.pca, pointshape= 21,  fill=merge_op_data$Stage, col.ind= merge_op_data$Stage,  legend.title = "Stage", mean.point= FALSE, pointsize=3, repel= TRUE)  # Avoid text overlapping 
dev.off()


pdf("pca_biplot_cambi_only.pdf")
fviz_pca_biplot(res.pca, label="all", col.ind = merge_op_data$Stage,legend.title = "Stage", mean.point= FALSE,
                repel = TRUE)
dev.off()


pdf("pca_var_contrib_top50_cambi_only.pdf")
fviz_pca_var(res.pca, col.var = "contrib", 
             select.var = list(contrib = 50),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping 
)
dev.off()

pdf("top20_var.contrib_cambi_only.pdf")
fviz_contrib(res.pca, choice = "var", axes=1, top = 20)
dev.off()

ind_contrib = res.pca$ind$contrib

ind_dim1_contribution <- cbind(merge_op_data[c(1:3,14)], ind_contrib)

var_contrib = res.pca$var$contrib

write.table(var_contrib, file="variable_dim1_dim5_contribs.txt", sep="\t")

long_ind_contrib <- gather(ind_dim1_contribution, contrib, value, Dim.1:Dim.5 )

long_ind_contrib$value <- factor(long_ind_contrib$value,                                    # Factor levels in decreasing order
                                 levels =long_ind_contrib$value[order(long_ind_contrib$Sample_name, decreasing = TRUE)])

Dim1 = subset(long_ind_contrib, contrib == "Dim.1")

Dim1$value= as.numeric(Dim1$value)

ggplot(Dim1, aes(Sample_name, value, col=Site, fill=Site)) + geom_bar(stat="identity")

pdf("individual_contribution_barplot_per_stage_cambi_only.pdf")
Dim1 %>% 
  ggplot(aes(fct_rev(fct_reorder(Sample_name,
                                 value)), value), fill=Stage) + geom_col(aes(fill=Stage))  + ylab("Contribution") + xlab("Sample") + facet_grid(Stage~Site, scales = "free_y") + theme(axis.text.x=element_text(angle=90))

dev.off()


