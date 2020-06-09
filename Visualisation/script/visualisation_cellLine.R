library("ComplexHeatmap")
folder_visu<-"/AKTDegrader/Visualisation/"
filename <-"Table_celllines_nopAKT_20191115_order.csv"

data_ccle <- read.csv(file=paste0(folder_visu,"data/",filename),header=TRUE)

mat_ccle <-t(as.matrix(rep(1,38)))
colnames(mat_ccle)<-data_ccle$Cell.Line_old_me
ccle_heatmap<-HeatmapAnnotation(KRAS = data_ccle$KRAS,
                                BRAF = data_ccle$BRAF,
                                NF1 = data_ccle$NF1,
                                PIK3CA = data_ccle$PIK3CA,
                                PTEN = data_ccle$PTEN,
                                "HER2+" = data_ccle$HER2.,
                                Tissue = data_ccle$Tissue,
                                annotation_name_side = "left")
dev.off()
draw(ccle_heatmap)

filename=paste0(folder_visu,"/results/Heatmap_cellLines_KRAS_vs_AKT_updateNF1_status1.svg")
svg(filename=filename,width=18,height=12)
Heatmap(mat_ccle, name = "mat", 
        bottom_annotation = ccle_heatmap,
        cluster_columns = FALSE)
dev.off()

filename=paste0(folder_visu,"/results/Heatmap_cellLines_KRAS_vs_AKT_updateNF1_status1.png")
png(filename=filename,width = 600, height = 480)
Heatmap(mat_ccle, name = "mat", 
        bottom_annotation = ccle_heatmap,
        cluster_columns = FALSE)
dev.off()

