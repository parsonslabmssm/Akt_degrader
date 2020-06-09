###Please download data from GEO #

library(ComplexHeatmap)
library(circlize)
library(data.table)

dirRPPA <- "/AKTDegrader/RPPA"
fileRppa <- paste0(dirRPPA,"RPPA_Akt_normlog2.csv")
data<- read.csv(file=fileRppa,header=TRUE)
rownames(data) <- data[,1]
data <- data[,-1]

col_fun = colorRamp2(c(-4, 0, 4), c("#5e3c99", "#f7f7f7", "#e66101"))

### Heatmap 1
listassay <- c("468.DMSO", "468.DMSO-R","468.AZD","468.AZD-R", "468.021", "468.021-R",
               "PC-3.DMSO", "PC-3.DMSO-R", "PC-3.AZD",  "PC-3.AZD-R", "PC-3.021", "PC-3.021-R")

#reorder according to listassay as Jia asks
mydf<-t(data)
mydf <-as.data.frame(mydf)
setDT(mydf)
mydf_order<-setcolorder(mydf, listassay)
data_order <- t(mydf_order)
colnames(data_order) <- colnames(data)

heatmap_all <- Heatmap(data_order,cluster_rows =FALSE,
        col =col_fun ,
        name = "Log2 median centered",
        column_names_side = "top", column_dend_side = "bottom")

filename_heatmapall=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_heatmap.svg")
svg(filename=filename_heatmapall, 
    width=70, 
    height=9, 
    pointsize=12)
heatmap_all
dev.off()

######### 468
## Wilcoxon
data_ctrl <- data[which(rownames(data) %in% c("468.DMSO", "468.DMSO-R")),]
data_azd <- data[which(rownames(data) %in% c("468.AZD","468.AZD-R")),]
data_021 <- data[which(rownames(data) %in% c("468.021", "468.021-R")),]

## Wilcoxon 468-DMSO vs 468-AZD
data_results_wilcoxon <- c("Protein",
                           "mean_468_DMSO ","sd_468_DMSO ",
                           "mean_468_AZD","sd_468_AZD", "Beta",
                           "Pvalue_Wilcoxon")
num_prot <- ncol(data)
for(j in 1:num_prot){
  
  name_prot <- colnames(data_ctrl)[j]
  prot_ctrl <-data_ctrl[,j]
  prot_case <-data_azd[,j]
  
  mean_control <- mean(prot_ctrl)
  sd_control <- sd(prot_ctrl)
  mean_case <- mean(prot_case)
  sd_case <- sd(prot_case)
  
  prot_test_mb <- wilcox.test(prot_ctrl, prot_case) 
  prot_pvalue <- prot_test_mb$p.value
  
  ### anova
  protein_value <- c(prot_ctrl,prot_case)
  sample_status <- as.factor(c("1DMSO","1DMSO","2AZD","2AZD"))
  #null model (only age)
  null_model <- lm(protein_value ~ 1 )
  
  ##full model (add status)
  full_model <- lm(protein_value ~ 1 + sample_status)
  
  compare <- anova(full_model,null_model)
  fitted <- summary(full_model)$coefficients
  beta_protein <- fitted[2,1]
  prot_pvalue<- compare[2,6]
  
  prot_data <- c(name_prot,
                 mean_control,sd_control,mean_case,sd_case,
                 beta_protein,prot_pvalue)
  data_results_wilcoxon <-rbind(data_results_wilcoxon,prot_data)
}
colnames(data_results_wilcoxon) <- data_results_wilcoxon[1,]
data_results_wilcoxon <- data_results_wilcoxon[-1,]
data_results_wilcoxon <- as.data.frame(data_results_wilcoxon)
pvalue_all <- as.numeric(as.character(data_results_wilcoxon$Pvalue_Wilcoxon))

prot_pvalue.adjust_Wilcoxon <- p.adjust(pvalue_all,method ="fdr")
data_results_wilcoxon_DMSO_AZD <- cbind(data_results_wilcoxon,prot_pvalue.adjust_Wilcoxon)
data_results_wilcoxon_DMSO_AZD$Pvalue_Wilcoxon <- as.numeric(as.character(data_results_wilcoxon_DMSO_AZD$Pvalue_Wilcoxon))

data_wilcoxon_DMSO_AZD_signFDR <- data_results_wilcoxon_DMSO_AZD[which(data_results_wilcoxon_DMSO_AZD$prot_pvalue.adjust_Wilcoxon <= 0.05),]
dim(data_wilcoxon_DMSO_AZD_signFDR)
data_wilcoxon_DMSO_AZD_Nomsign <- data_results_wilcoxon_DMSO_AZD[which(data_results_wilcoxon_DMSO_AZD$Pvalue_Wilcoxon <= 0.05),]
dim(data_wilcoxon_DMSO_AZD_Nomsign)
data_wilcoxon_DMSO_AZD_Nomsign$Beta <- as.numeric(as.character(data_wilcoxon_DMSO_AZD_Nomsign$Beta))
data_wilcoxon_DMSO_AZD_Nomsign_1 <- data_wilcoxon_DMSO_AZD_Nomsign[which(data_wilcoxon_DMSO_AZD_Nomsign$Beta >=1 | data_wilcoxon_DMSO_AZD_Nomsign$Beta <=-1),]
dim(data_wilcoxon_DMSO_AZD_Nomsign_1)

fileRppa <- paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_DMSO_AZD_20191119.csv")
write.csv(data_results_wilcoxon_DMSO_AZD,file=fileRppa,row.names = FALSE)


### Heatmap 2
listcellline2 <- c("468.DMSO", "468.DMSO-R", "468.AZD","468.AZD-R")
listantibody2<- data_wilcoxon_DMSO_AZD_Nomsign$Protein

data_heatmap2 <- data[which(rownames(data) %in% listcellline2),
                      which(colnames(data) %in% listantibody2)]
dim(data_heatmap2)

heatmapDMSOvsAZD<- Heatmap(t(data_heatmap2),col = col_fun,
        name = "Log2 median centered",
        column_names_side = "top", column_dend_side = "bottom")

filename_heatmapDMSOvsAZD=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_DMSOvsAZD_20191119.svg")
svg(filename=filename_heatmapDMSOvsAZD, 
    width=9, 
    height=9, 
    pointsize=12)
heatmapDMSOvsAZD
dev.off()

######### Wilcoxon 468-DMSO vs 468.021
data_results_wilcoxon <- c("Protein",
                           "mean_468_DMSO","sd_468_DMSO",
                           "mean_468_021","sd_468_021", "Beta",
                           "Pvalue_Wilcoxon")
num_prot <- ncol(data)
for(j in 1:num_prot){
  
  name_prot <- colnames(data_ctrl)[j]
  prot_ctrl <-data_ctrl[,j]
  prot_case <-data_021[,j]
  
  mean_control <- mean(prot_ctrl)
  sd_control <- sd(prot_ctrl)
  mean_case <- mean(prot_case)
  sd_case <- sd(prot_case)
  
  prot_test_mb <- wilcox.test(prot_ctrl, prot_case) 
  prot_pvalue <- prot_test_mb$p.value
  
  ### anova
  protein_value <- c(prot_ctrl,prot_case)
  sample_status <- as.factor(c("1DMSO","1DMSO","2021","2021"))
  #null model (only age)
  null_model <- lm(protein_value ~ 1 )
  
  ##full model (add status)
  full_model <- lm(protein_value ~ 1 + sample_status)
  
  compare <- anova(full_model,null_model)
  fitted <- summary(full_model)$coefficients
  beta_protein <- fitted[2,1]
  prot_pvalue<- compare[2,6]
  
  prot_data <- c(name_prot,
                 mean_control,sd_control,mean_case,sd_case,
                 beta_protein,prot_pvalue)
  data_results_wilcoxon <-rbind(data_results_wilcoxon,prot_data)
}
colnames(data_results_wilcoxon) <- data_results_wilcoxon[1,]
data_results_wilcoxon <- data_results_wilcoxon[-1,]
data_results_wilcoxon <- as.data.frame(data_results_wilcoxon)
pvalue_all <- as.numeric(as.character(data_results_wilcoxon$Pvalue_Wilcoxon))

prot_pvalue.adjust_Wilcoxon <- p.adjust(pvalue_all,method ="fdr")
data_results_wilcoxon_DMSO_021 <- cbind(data_results_wilcoxon,prot_pvalue.adjust_Wilcoxon)
data_results_wilcoxon_DMSO_021$Pvalue_Wilcoxon <- as.numeric(as.character(data_results_wilcoxon_DMSO_021$Pvalue_Wilcoxon))

data_wilcoxon_DMSO_021_signFDR <- data_results_wilcoxon_DMSO_021[which(data_results_wilcoxon_DMSO_021$prot_pvalue.adjust_Wilcoxon <= 0.05),]
dim(data_wilcoxon_DMSO_021_signFDR)
data_wilcoxon_DMSO_021_Nomsign <- data_results_wilcoxon_DMSO_021[which(data_results_wilcoxon_DMSO_021$Pvalue_Wilcoxon <= 0.05),]
dim(data_wilcoxon_DMSO_021_Nomsign)
data_wilcoxon_DMSO_021_Nomsign$Beta <- as.numeric(as.character(data_wilcoxon_DMSO_021_Nomsign$Beta))
data_wilcoxon_DMSO_021_Nomsign_1 <- data_wilcoxon_DMSO_021_Nomsign[which(data_wilcoxon_DMSO_021_Nomsign$Beta <= -1 | 
                                                                           data_wilcoxon_DMSO_021_Nomsign$Beta >= 1),]
dim(data_wilcoxon_DMSO_021_Nomsign_1)

fileRppa_DMSO_021 <- paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_DMSO_021_20191119.csv")
write.csv(data_results_wilcoxon_DMSO_021,file=fileRppa_DMSO_021,row.names = FALSE)


### Heatmap 2
listcellline2 <- c("468.DMSO", "468.DMSO-R", "468.021","468.021-R")
listantibody2<- data_wilcoxon_DMSO_021_Nomsign$Protein

data_heatmap2 <- data[which(rownames(data) %in% listcellline2),
                      which(colnames(data) %in% listantibody2)]
dim(data_heatmap2)

heatmapDMSOvs021<- Heatmap(t(data_heatmap2),col = col_fun,
                           name = "Log2 median centered",
                           column_names_side = "top", 
                           column_dend_side = "bottom",
                           cluster_columns = TRUE,
                           cluster_rows = FALSE)

filename_heatmapDMSOvs021=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_DMSOvs021_20191119.svg")
svg(filename=filename_heatmapDMSOvs021, 
    width=9, 
    height=9, 
    pointsize=12)
heatmapDMSOvs021
dev.off()

### Heatmap 4
listcellline2 <- c("468.DMSO", "468.DMSO-R", "468.021","468.021-R", "468.AZD","468.AZD-R")
listantibody2<- c(as.character(data_wilcoxon_DMSO_AZD_Nomsign$Protein),
                  as.character(data_wilcoxon_DMSO_021_Nomsign$Protein))

data_heatmap2 <- data_order[which(rownames(data_order) %in% listcellline2),
                      which(colnames(data_order) %in% listantibody2)]
dim(data_heatmap2)

heatmapDMSOvs021vsAZ<- Heatmap(t(data_heatmap2),col = col_fun,
                           name = "Log2 median centered",
                           column_names_side = "top", 
                           column_dend_side = "bottom",
                           cluster_columns = FALSE,
                           cluster_rows = FALSE)

filename_heatmapDMSOvs021vsAZ=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_DMSOvs021vsAZD_20191119.svg")
svg(filename=filename_heatmapDMSOvs021vsAZ, 
    width=9, 
    height=9, 
    pointsize=12)
heatmapDMSOvs021vsAZ
dev.off()

### Heatmap 5
listcellline2 <- c("468.DMSO", "468.DMSO-R", "468.021","468.021-R", "468.AZD","468.AZD-R")
listantibody3<- c("GSK.3a.b_pS21_S9","S6_pS235_S236" ,"S6_pS240_S244",
                  "Akt1_pS473","Akt_pT308","Akt_pS473",
                  "Akt","Akt1","Akt2","Akt2_pS474")

data_heatmap2_458_tmp <- data_order[which(rownames(data_order) %in% listcellline2),
                      which(colnames(data_order) %in% listantibody3)]
dim(data_heatmap2_458_tmp)
data_heatmap2_458_tmp_dt<-setDT(as.data.frame(data_heatmap2_458_tmp))
data_heatmap2_458 <- setcolorder(data_heatmap2_458_tmp_dt, as.character(listantibody3))
dim(data_heatmap2_458)
data_heatmap2_458_mat <- as.matrix(data_heatmap2_458)
rownames(data_heatmap2_458_mat) <- rownames(data_heatmap2_458_tmp)

heatmapDMSOvs021vsAZ_Beta1_458<- Heatmap(t(data_heatmap2_458_mat),col = col_fun,
                               name = "Log2 median centered",
                               column_names_side = "top", 
                               column_dend_side = "bottom",
                               cluster_columns = FALSE,
                               cluster_rows = FALSE)

filename_heatmapDMSOvs021vsAZ_Beta1_458=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_DMSOvs021vsAZD_Beta1_noorder_20191119.svg")
svg(filename=filename_heatmapDMSOvs021vsAZ_Beta1_458, 
    width=9, 
    height=9, 
    pointsize=12)
heatmapDMSOvs021vsAZ_Beta1_458
dev.off()

filename_heatmapDMSOvs021vsAZ_Beta1_458_png=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_458_DMSOvs021vsAZD_Beta1_20191119.png")
png(filename=filename_heatmapDMSOvs021vsAZ_Beta1_458_png, 
    width=500, 
    height=500, 
    pointsize=12)
heatmapDMSOvs021vsAZ_Beta1_458
dev.off()

######### 468 vs PC3
## Wilcoxon
data_468 <- data[which(rownames(data) %in% c("468.DMSO", "468.DMSO-R","468.AZD","468.AZD-R", "468.021", "468.021-R")),]
data_PC3 <- data[which(rownames(data) %in% c("PC-3.DMSO", "PC-3.DMSO-R", "PC-3.AZD",  "PC-3.AZD-R", "PC-3.021", "PC-3.021-R")),]

## Wilcoxon 468-DMSO vs 468-AZD
data_results_wilcoxon <- c("Protein",
                           "mean_468 ","sd_468 ",
                           "mean_PC-3","sd_PC-3", "Beta",
                           "Pvalue_Wilcoxon")
num_prot <- ncol(data)
for(j in 1:num_prot){
  
  name_prot <- colnames(data_468)[j]
  prot_ctrl <-data_468[,j]
  prot_case <-data_PC3[,j]
  
  mean_control <- mean(prot_ctrl)
  sd_control <- sd(prot_ctrl)
  mean_case <- mean(prot_case)
  sd_case <- sd(prot_case)
  
  prot_test_mb <- wilcox.test(prot_ctrl, prot_case) 
  prot_pvalue <- prot_test_mb$p.value
  
  ### anova
  protein_value <- c(prot_ctrl,prot_case)
  sample_status <- as.factor(c("A468","A468","A468","A468","A468","A468",
                               "PC3","PC3","PC3","PC3","PC3","PC3"))
  #null model (only age)
  null_model <- lm(protein_value ~ 1 )
  
  ##full model (add status)
  full_model <- lm(protein_value ~ 1 + sample_status)
  
  compare <- anova(full_model,null_model)
  fitted <- summary(full_model)$coefficients
  beta_protein <- fitted[2,1]
  prot_pvalue<- compare[2,6]
  
  prot_data <- c(name_prot,
                 mean_control,sd_control,mean_case,sd_case,
                 beta_protein,prot_pvalue)
  data_results_wilcoxon <-rbind(data_results_wilcoxon,prot_data)
}
colnames(data_results_wilcoxon) <- data_results_wilcoxon[1,]
data_results_wilcoxon <- data_results_wilcoxon[-1,]
data_results_wilcoxon <- as.data.frame(data_results_wilcoxon)
pvalue_all <- as.numeric(as.character(data_results_wilcoxon$Pvalue_Wilcoxon))

prot_pvalue.adjust_Wilcoxon <- p.adjust(pvalue_all,method ="fdr")
data_results_wilcoxon_468_PC3<- cbind(data_results_wilcoxon,prot_pvalue.adjust_Wilcoxon)
data_results_wilcoxon_468_PC3$Pvalue_Wilcoxon <- as.numeric(as.character(data_results_wilcoxon_468_PC3$Pvalue_Wilcoxon))

data_wilcoxon_468_PC3_signFDR <- data_results_wilcoxon_468_PC3[which(data_results_wilcoxon_468_PC3$prot_pvalue.adjust_Wilcoxon <= 0.05),]
dim(data_wilcoxon_468_PC3_signFDR)
data_wilcoxon_468_PC3_Nomsign <- data_results_wilcoxon_468_PC3[which(data_results_wilcoxon_468_PC3$Pvalue_Wilcoxon <= 0.05),]
dim(data_wilcoxon_468_PC3_Nomsign)
data_wilcoxon_468_PC3_signFDR$Beta <- as.numeric(as.character(data_wilcoxon_468_PC3_signFDR$Beta))
data_wilcoxon_468_PC3_signFDR_beta1 <- data_wilcoxon_468_PC3_signFDR[which(data_wilcoxon_468_PC3_signFDR$Beta <= -1 |
                                                                             data_wilcoxon_468_PC3_signFDR$Beta >= 1),]
dim(data_wilcoxon_468_PC3_signFDR_beta1)

fileRppa_468_PC3 <- paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_PC3_20191119.csv")
write.csv(data_results_wilcoxon_468_PC3,file=fileRppa_468_PC3,row.names = FALSE)

### Heatmap 10
listantibody_FDR<- data_wilcoxon_468_PC3_signFDR$Protein
listantibody_Nomsign<- data_wilcoxon_468_PC3_Nomsign$Protein

data_heatmap2 <- data_order[, which(colnames(data_order) %in% listantibody_FDR)]
dim(data_heatmap2)

heatmap_468_PC3<- Heatmap(data_heatmap2,col = col_fun,
                           name = "Log2 median centered",
                           column_names_side = "top", column_dend_side = "bottom",
                          cluster_columns = FALSE,
                          cluster_rows = FALSE)

filename_heatmap_468_PC3=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_PC3_FDR.svg")
svg(filename=filename_heatmap_468_PC3, 
    width=40, 
    height=9, 
    pointsize=12)
heatmap_468_PC3
dev.off()

heatmap_468_PC3_noorder<- Heatmap(data_heatmap2,col = col_fun,
                          name = "Log2 median centered",
                          column_names_side = "top", column_dend_side = "bottom",
                          cluster_columns = FALSE,
                          cluster_rows = FALSE)

filename_heatmap_468_PC3_noorder=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_PC3_noorder_FDR_20191119.svg")
svg(filename=filename_heatmap_468_PC3_noorder, 
    width=40, 
    height=9, 
    pointsize=12)
heatmap_468_PC3_noorder
dev.off()

### Heatmap 11
listantibody_FDR_beta1<- data_wilcoxon_468_PC3_signFDR_beta1$Protein


data_heatmap2_beta1 <- data_order[, which(colnames(data_order) %in% listantibody_FDR_beta1)]
dim(data_heatmap2_beta1)

heatmap_468_PC3_noorder_beta1<- Heatmap(data_heatmap2_beta1,col = col_fun,
                                  name = "Log2 median centered",
                                  column_names_side = "top", column_dend_side = "bottom",
                                  cluster_columns = FALSE,
                                  cluster_rows = FALSE)

filename_heatmap_468_PC3_noorder_beta1=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_PC3_noorder_FDR_beta1_20191119.svg")
svg(filename=filename_heatmap_468_PC3_noorder_beta1, 
    width=40, 
    height=9, 
    pointsize=12)
heatmap_468_PC3_noorder_beta1
dev.off()

select_antibody <-c("Connexin.43","Axl","MCT4","EphA2",
                    "Caveolin.1","MMP14","NDRG1_pT346","Hif.1.alpha","HMHA1",
                    "IDO","Paxillin","Vimentin","EphA2_pY588","SHP.2_pY542",
                    "Atg7","VHL")
order_samples<-c("PC-3.DMSO","PC-3.DMSO-R","PC-3.AZD","PC-3.AZD-R","PC-3.021","PC-3.021-R",
                 "468.DMSO","468.DMSO-R","468.AZD","468.AZD-R","468.021","468.021-R")
data_heatmap2_beta1_selected<-t(data_heatmap2_beta1[,which(colnames(data_heatmap2_beta1) %in% select_antibody)])
data_heatmap2_beta1_selectedtmp <- data_heatmap2_beta1_selected[match(select_antibody, rownames(data_heatmap2_beta1_selected)),]
data_heatmap2_beta1_selected_dt <- setDT(as.data.frame(data_heatmap2_beta1_selectedtmp))
setcolorder(data_heatmap2_beta1_selected_dt, as.character(order_samples))
df_data_heatmap2_beta1_selected <- as.matrix(data_heatmap2_beta1_selected_dt)
rownames(df_data_heatmap2_beta1_selected) <- rownames(data_heatmap2_beta1_selectedtmp)
heatmap_468_PC3_noorder_beta1_selected<- Heatmap(df_data_heatmap2_beta1_selected,col = col_fun,
                                        name = "Log2 median centered",
                                        column_names_side = "top", column_dend_side = "bottom",
                                        cluster_columns = FALSE,
                                        cluster_rows = FALSE)

filename_heatmap_468_PC3_noorder_beta1=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_PC3_noorder_FDR_selected_20191119.svg")
svg(filename=filename_heatmap_468_PC3_noorder_beta1, 
    width=9, 
    height=12, 
    pointsize=12)
heatmap_468_PC3_noorder_beta1_selected
dev.off()

filename_heatmap_468_PC3_noorder_selected_png=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_468_PC3_noorder_FDR_selected_20191119.png")
png(filename=filename_heatmap_468_PC3_noorder_selected_png, 
    width=500, 
    height=600, 
    pointsize=12)
heatmap_468_PC3_noorder_beta1_selected
dev.off()


######### PC-3
## Wilcoxon
data_ctrl <- data[which(rownames(data) %in% c("PC-3.DMSO", "PC-3.DMSO-R")),]
data_azd <- data[which(rownames(data) %in% c("PC-3.AZD","PC-3.AZD-R")),]
data_021 <- data[which(rownames(data) %in% c("PC-3.021", "PC-3.021-R")),]

## Wilcoxon PC-3-DMSO vs PC-3-AZD
data_results_wilcoxon <- c("Protein",
                           "mean_PC-3_DMSO ","sd_PC-3_DMSO ",
                           "mean_PC-3_AZD","sd_PC-3_AZD", "Beta",
                           "Pvalue_Wilcoxon")
num_prot <- ncol(data)
for(j in 1:num_prot){
  
  name_prot <- colnames(data_ctrl)[j]
  prot_ctrl <-data_ctrl[,j]
  prot_case <-data_azd[,j]
  
  mean_control <- mean(prot_ctrl)
  sd_control <- sd(prot_ctrl)
  mean_case <- mean(prot_case)
  sd_case <- sd(prot_case)
  
  prot_test_mb <- wilcox.test(prot_ctrl, prot_case) 
  prot_pvalue <- prot_test_mb$p.value
  
  ### anova
  protein_value <- c(prot_ctrl,prot_case)
  sample_status <- as.factor(c("1DMSO","1DMSO","2AZD","2AZD"))
  #null model (only age)
  null_model <- lm(protein_value ~ 1 )
  
  ##full model (add status)
  full_model <- lm(protein_value ~ 1 + sample_status)
  
  compare <- anova(full_model,null_model)
  fitted <- summary(full_model)$coefficients
  beta_protein <- fitted[2,1]
  prot_pvalue<- compare[2,6]
  
  prot_data <- c(name_prot,
                 mean_control,sd_control,mean_case,sd_case,
                 beta_protein,prot_pvalue)
  data_results_wilcoxon <-rbind(data_results_wilcoxon,prot_data)
}
colnames(data_results_wilcoxon) <- data_results_wilcoxon[1,]
data_results_wilcoxon <- data_results_wilcoxon[-1,]
data_results_wilcoxon <- as.data.frame(data_results_wilcoxon)
pvalue_all <- as.numeric(as.character(data_results_wilcoxon$Pvalue_Wilcoxon))

prot_pvalue.adjust_Wilcoxon <- p.adjust(pvalue_all,method ="fdr")
data_results_wilcoxon_DMSO_AZD_PC3 <- cbind(data_results_wilcoxon,prot_pvalue.adjust_Wilcoxon)
data_results_wilcoxon_DMSO_AZD_PC3$Pvalue_Wilcoxon <- as.numeric(as.character(data_results_wilcoxon_DMSO_AZD_PC3$Pvalue_Wilcoxon))

data_wilcoxon_DMSO_AZD_signFDR_PC3 <- data_results_wilcoxon_DMSO_AZD[which(data_results_wilcoxon_DMSO_AZD_PC3$prot_pvalue.adjust_Wilcoxon <= 0.05),]
dim(data_wilcoxon_DMSO_AZD_signFDR_PC3)
data_wilcoxon_DMSO_AZD_Nomsign_PC3 <- data_results_wilcoxon_DMSO_AZD[which(data_results_wilcoxon_DMSO_AZD_PC3$Pvalue_Wilcoxon <= 0.05),]
dim(data_wilcoxon_DMSO_AZD_Nomsign_PC3)
data_wilcoxon_DMSO_AZD_Nomsign_PC3$Beta <- as.numeric(as.character(data_wilcoxon_DMSO_AZD_Nomsign_PC3$Beta))
data_wilcoxon_DMSO_AZD_Nomsign_1_PC3 <- data_wilcoxon_DMSO_AZD_Nomsign_PC3[which(data_wilcoxon_DMSO_AZD_Nomsign_PC3$Beta >=1 | data_wilcoxon_DMSO_AZD_Nomsign_PC3$Beta <=-1),]
dim(data_wilcoxon_DMSO_AZD_Nomsign_1_PC3)

fileRppa <- paste0(dirRPPA,"/results/RPPA_Akt_normlog2_PC-3_DMSO_AZD_20191119.csv")
write.csv(data_results_wilcoxon_DMSO_AZD_PC3,file=fileRppa,row.names = FALSE)


### Heatmap 2
listcellline2 <- c("PC-3.DMSO", "PC-3.DMSO-R", "PC-3.AZD","PC-3.AZD-R")
listantibody2<- data_wilcoxon_DMSO_AZD_Nomsign$Protein

data_heatmap2 <- data[which(rownames(data) %in% listcellline2),
                      which(colnames(data) %in% listantibody2)]
dim(data_heatmap2)

heatmapDMSOvsAZD_PC3<- Heatmap(t(data_heatmap2),col = col_fun,
                           name = "Log2 median centered",
                           column_names_side = "top", column_dend_side = "bottom",
                           cluster_columns = FALSE,
                           cluster_rows = FALSE)

filename_heatmapDMSOvsAZD_PC3=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_PC-3_DMSOvsAZD_20191119.svg")
svg(filename=filename_heatmapDMSOvsAZD_PC3, 
    width=9, 
    height=9, 
    pointsize=12)
heatmapDMSOvsAZD_PC3
dev.off()

######### Wilcoxon PC-3-DMSO vs PC-3.021
data_results_wilcoxon <- c("Protein",
                           "mean_PC-3_DMSO","sd_PC-3_DMSO",
                           "mean_PC-3_021","sd_PC-3_021", "Beta",
                           "Pvalue_Wilcoxon")
num_prot <- ncol(data)
for(j in 1:num_prot){
  
  name_prot <- colnames(data_ctrl)[j]
  prot_ctrl <-data_ctrl[,j]
  prot_case <-data_021[,j]
  
  mean_control <- mean(prot_ctrl)
  sd_control <- sd(prot_ctrl)
  mean_case <- mean(prot_case)
  sd_case <- sd(prot_case)
  
  prot_test_mb <- wilcox.test(prot_ctrl, prot_case) 
  prot_pvalue <- prot_test_mb$p.value
  
  ### anova
  protein_value <- c(prot_ctrl,prot_case)
  sample_status <- as.factor(c("1DMSO","1DMSO","2021","2021"))
  #null model (only age)
  null_model <- lm(protein_value ~ 1 )
  
  ##full model (add status)
  full_model <- lm(protein_value ~ 1 + sample_status)
  
  compare <- anova(full_model,null_model)
  fitted <- summary(full_model)$coefficients
  beta_protein <- fitted[2,1]
  prot_pvalue<- compare[2,6]
  
  prot_data <- c(name_prot,
                 mean_control,sd_control,mean_case,sd_case,
                 beta_protein,prot_pvalue)
  data_results_wilcoxon <-rbind(data_results_wilcoxon,prot_data)
}
colnames(data_results_wilcoxon) <- data_results_wilcoxon[1,]
data_results_wilcoxon <- data_results_wilcoxon[-1,]
data_results_wilcoxon <- as.data.frame(data_results_wilcoxon)
pvalue_all <- as.numeric(as.character(data_results_wilcoxon$Pvalue_Wilcoxon))

prot_pvalue.adjust_Wilcoxon <- p.adjust(pvalue_all,method ="fdr")
data_results_wilcoxon_DMSO_021_PC3 <- cbind(data_results_wilcoxon,prot_pvalue.adjust_Wilcoxon)
data_results_wilcoxon_DMSO_021_PC3$Pvalue_Wilcoxon <- as.numeric(as.character(data_results_wilcoxon_DMSO_021_PC3$Pvalue_Wilcoxon))

data_wilcoxon_DMSO_021_signFDR_PC3 <- data_results_wilcoxon_DMSO_021_PC3[which(data_results_wilcoxon_DMSO_021_PC3$prot_pvalue.adjust_Wilcoxon <= 0.05),]
dim(data_wilcoxon_DMSO_021_signFDR_PC3)
data_wilcoxon_DMSO_021_Nomsign_PC3 <- data_results_wilcoxon_DMSO_021_PC3[which(data_results_wilcoxon_DMSO_021_PC3$Pvalue_Wilcoxon <= 0.05),]
dim(data_wilcoxon_DMSO_021_Nomsign_PC3)
data_wilcoxon_DMSO_021_Nomsign_PC3$Beta <- as.numeric(as.character(data_wilcoxon_DMSO_021_Nomsign_PC3$Beta))
data_wilcoxon_DMSO_021_Nomsign_1_PC3 <- data_wilcoxon_DMSO_021_Nomsign_PC3[which(data_wilcoxon_DMSO_021_Nomsign_PC3$Beta <= -1 | 
                                                                           data_wilcoxon_DMSO_021_Nomsign_PC3$Beta >= 1),]
dim(data_wilcoxon_DMSO_021_Nomsign_1_PC3)

fileRppa_DMSO_021 <- paste0(dirRPPA,"/results/RPPA_Akt_normlog2_PC-3_DMSO_021_20191119.csv")
write.csv(data_results_wilcoxon_DMSO_021_PC3,file=fileRppa_DMSO_021,row.names = FALSE)


### Heatmap 2
listcellline2 <- c("PC-3.DMSO", "PC-3.DMSO-R", "PC-3.021","PC-3.021-R")
listantibody2<- data_wilcoxon_DMSO_021_Nomsign_PC3$Protein

data_heatmap2 <- data[which(rownames(data) %in% listcellline2),
                      which(colnames(data) %in% listantibody2)]
dim(data_heatmap2)

heatmapDMSOvs021_PC3<- Heatmap(t(data_heatmap2),col = col_fun,
                           name = "Log2 median centered",
                           column_names_side = "top", column_dend_side = "bottom",
                           cluster_columns = TRUE,
                           cluster_rows = FALSE)

filename_heatmapDMSOvs021_PC3=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_PC-3_DMSOvs021.svg")
svg(filename=filename_heatmapDMSOvs021_PC3, 
    width=9, 
    height=9, 
    pointsize=12)
heatmapDMSOvs021_PC3
dev.off()

### Heatmap 4
listcellline2 <- c("PC-3.DMSO", "PC-3.DMSO-R", "PC-3.021","PC-3.021-R", "PC-3.AZD","PC-3.AZD-R")
listantibody2<- c(as.character(data_wilcoxon_DMSO_AZD_Nomsign_PC3$Protein),
                  as.character(data_wilcoxon_DMSO_021_Nomsign_PC3$Protein))

data_heatmap2 <- data[which(rownames(data) %in% listcellline2),
                      which(colnames(data) %in% listantibody2)]
dim(data_heatmap2)

heatmapDMSOvs021vsAZ_PC3<- Heatmap(t(data_heatmap2),col = col_fun,
                               name = "Log2 median centered",
                               column_names_side = "top", column_dend_side = "bottom",
                               cluster_columns = TRUE,
                               cluster_rows = TRUE)

filename_heatmapDMSOvs021vsAZ=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_PC-3_DMSOvs021vsAZD_20191119.svg")
svg(filename=filename_heatmapDMSOvs021vsAZ, 
    width=9, 
    height=9, 
    pointsize=12)
heatmapDMSOvs021vsAZ_PC3
dev.off()

### Heatmap 5
listcellline2 <- c("PC-3.DMSO", "PC-3.DMSO-R", "PC-3.021","PC-3.021-R", "PC-3.AZD","PC-3.AZD-R")
listantibody3<- c("GSK.3a.b_pS21_S9","S6_pS235_S236" ,"S6_pS240_S244",
                  "Akt1_pS473","Akt_pT308","Akt_pS473",
                  "Akt","Akt1","Akt2","Akt2_pS474")

data_heatmap2_PC3_tmp <- data[which(rownames(data) %in% listcellline2),
                      which(colnames(data) %in% listantibody3)]
dim(data_heatmap2_PC3_tmp)
name_row_PC3 <-  rownames(data_heatmap2_PC3_tmp)
data_heatmap2_PC3_tmp_dt<-setDT(data_heatmap2_PC3_tmp)
data_heatmap2_PC3 <- setcolorder(data_heatmap2_PC3_tmp_dt, as.character(listantibody3))
dim(data_heatmap2_PC3)
data_heatmap2_PC3_mat_tmp <- as.matrix(data_heatmap2_PC3)
rownames(data_heatmap2_PC3_mat_tmp) <- name_row_PC3
order_samples_PC3 <- c("PC-3.DMSO","PC-3.DMSO-R",
                       "PC-3.AZD","PC-3.AZD-R",
                       "PC-3.021","PC-3.021-R")
data_heatmap2_PC3_mat <- data_heatmap2_PC3_mat_tmp[match(order_samples_PC3,rownames(data_heatmap2_PC3_mat_tmp)),]
#
heatmapDMSOvs021vsAZ_Beta1_PC3<- Heatmap(t(data_heatmap2_PC3_mat),col = col_fun,
                                     name = "Log2 median centered",
                                     column_names_side = "top", 
                                     column_dend_side = "bottom",
                                     cluster_columns = FALSE,
                                     cluster_rows = FALSE)

filename_heatmapDMSOvs021vsAZ_Beta1_PC3=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_PC-3_DMSOvs021vsAZD_Beta1_20191119.svg")
svg(filename=filename_heatmapDMSOvs021vsAZ_Beta1_PC3, 
    width=9, 
    height=9, 
    pointsize=12)
heatmapDMSOvs021vsAZ_Beta1_PC3
dev.off()

filename_heatmapDMSOvs021vsAZ_Beta1_PC3_png=paste0(dirRPPA,"/results/RPPA_Akt_normlog2_PC-3_DMSOvs021vsAZD_Beta1_20191119.png")
png(filename=filename_heatmapDMSOvs021vsAZ_Beta1_PC3_png, 
    width=500, 
    height=500, 
    pointsize=12)
heatmapDMSOvs021vsAZ_Beta1_PC3
dev.off()

