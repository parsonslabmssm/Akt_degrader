### Need to download data from CCLE for CCLE_RPPA_20181003.csv and  CCLE_RPPA_Ab_info_20181226.csv
library("DescTools")
library("dplyr")
library("tidyverse")
library("stringr")

DirJia="/AKTDegrader/CCLE/"

FileCCLE_RPPA=paste0(DirJia,"data/CCLE_RPPA_20181003.csv")
rppa_ccle <- read.csv(file=FileCCLE_RPPA,header=TRUE)
rownames(rppa_ccle) <- rppa_ccle[,1]

FileCCLE_RPPA_info=paste0(DirJia,"data/CCLE_RPPA_Ab_info_20181226.csv")
rppa_ccle_info <- read.csv(file=FileCCLE_RPPA_info,header=TRUE)
rppa_ccle_info_akt <- rppa_ccle_info[grep("^AKT|^EGFR|ERBB|^MAPK|^VHL", rppa_ccle_info[,2]),]

rppa_ccle_akt <- rppa_ccle[,which(colnames(rppa_ccle) %in% rppa_ccle_info_akt$Antibody_Name)]
dim(rppa_ccle_akt)
rppa_ccle_akt$tissue <- rownames(rppa_ccle_akt)
rppa_ccle_akt_ID<- rppa_ccle_akt %>% separate(tissue, "Cell.Line")
rppa_ccle_akt_ID$tissue <- rownames(rppa_ccle_akt)

num_ccle <- nrow(rppa_ccle_akt_ID)
for (i in 1:num_ccle){
  name_tmp <- str_split(rownames(rppa_ccle_akt_ID)[i], "_", simplify = TRUE)
  name_l <-length(name_tmp)
  new_name <-paste(name_tmp[2:name_l],collapse = ' ')
  rppa_ccle_akt_ID[i,"tissue"] <- new_name
}

rppa_ccle_akt_ID_order <- rppa_ccle_akt_ID

#Heatmap
library(ComplexHeatmap)
library(circlize)
library(data.table)

df_Site.Primary <-data.frame(rppa_ccle_akt_ID_order$tissue)
colnames(df_Site.Primary)[1]<- "Tissue"
ha = HeatmapAnnotation(df = df_Site.Primary)

rownames(rppa_ccle_akt_ID_order) <- rppa_ccle_akt_ID_order$ID
heatmap_short_order <- Heatmap(t(rppa_ccle_akt_ID_order[,c(1:19)]),
                               cluster_columns = FALSE,
                               cluster_rows = FALSE,
                               name = "Protein expression of AKT genes in CCLE",
                               column_names_side = "bottom", column_dend_side = "bottom",
                               bottom_annotation = ha)

filename_heatmap_short=paste0(DirJia,"/results/AKT_CCLE_RPPA_heatmap_largeProtein_20191119.svg")
svg(filename=filename_heatmap_short, 
    width=100, 
    height=9, 
    pointsize=12)
heatmap_short_order
dev.off()

### Mutations Cell lines
folder_visu<-"/AKTDegrader/Visualisation/"
filename <-"Table_celllines_nopAKT_20191115_order.csv"

data_ccle <- read.csv(file=paste0(folder_visu,"data/",filename),header=TRUE)

rppa_mutations_ccle_akt_all <- merge(data_ccle,rppa_ccle_akt_ID_order,by="Cell.Line",all=TRUE)
dim(data_ccle)
dim(rppa_ccle_akt_ID_order)
dim(rppa_mutations_ccle_akt_all)
write_csv(rppa_mutations_ccle_akt_all,path=paste0(DirJia,"data/RPPA_mutation_CCLE_4AKT_all_updatedNameCell_20191119.csv"))

rppa_mutations_ccle_akt <- merge(data_ccle,rppa_ccle_akt_ID_order,by="Cell.Line")
dim(rppa_mutations_ccle_akt)
write_csv(rppa_mutations_ccle_akt,path=paste0(DirJia,"data/RPPA_mutation_CCLE_4AKT_updatedNameCell_20191119.csv"))

##order to put resistance/weak sentitive/sentitive
rppa_mutations_ccle_akt_order <- read.csv(file=paste0(DirJia,"data/RPPA_mutation_CCLE_4AKT_updatedNameCell_20191119_order.csv"),header=TRUE)
rownames(rppa_mutations_ccle_akt_order) <- rppa_mutations_ccle_akt_order$Cell.Line_old_me
ccle_heatmap<-HeatmapAnnotation(KRAS = rppa_mutations_ccle_akt_order$KRAS,
                                BRAF = rppa_mutations_ccle_akt_order$BRAF,
                                NF1 = rppa_mutations_ccle_akt_order$NF1,
                                PIK3CA = rppa_mutations_ccle_akt_order$PIK3CA,
                                PTEN = rppa_mutations_ccle_akt_order$PTEN,
                                "HER2+" = rppa_mutations_ccle_akt_order$HER2.,
                                Tissue = rppa_mutations_ccle_akt_order$Tissue,
                                annotation_name_side = "left")
draw(ccle_heatmap)
col_fun = colorRamp2(c(-4, 0, 4), c("#5e3c99", "#f7f7f7", "#e66101"))
heatmap_rppa_mutation <- Heatmap(t(rppa_mutations_ccle_akt_order[,c(13:31)]),
                               cluster_columns = FALSE,
                               cluster_rows = FALSE,
                               name = "Protein expression\n of genes in CCLE",
                               column_names_side = "bottom", column_dend_side = "bottom",
                               bottom_annotation = ccle_heatmap, col = col_fun)
heatmap_rppa_mutation
filename_heatmap_short=paste0(DirJia,"/results/AKT_CCLE_RPPA_mutation_heatmap_order_updatedNameCell_2019119.svg")
svg(filename=filename_heatmap_short, 
    width=15, 
    height=7, 
    pointsize=12)
heatmap_rppa_mutation
dev.off()

filename=paste0(DirJia,"/results/AKT_CCLE_RPPA_mutation_heatmap_order_updatedNameCell_20191119.png")
png(filename=filename,width = 900, height = 550)
heatmap_rppa_mutation
dev.off()

