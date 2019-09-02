#GTEX COMP PROJECT


#GTEX analysis process with tpm

#load("~/GRIPS/gtex_tpm.txt")
#gtex_tpm <- gtex_tpm[-24493,]
#rownames(gtex_tpm) <- gtex_tpm$gene
#gtex_tpm <- gtex_tpm[,-1]

gtex_an<- gtex_tpm
colnames(gtex_an)<- gtexat_pt$tissue
#medians for each tissue type
gtex_anmeans<- data.frame(sapply(split.default(gtex_an, names(gtex_an)), rowMedian))
#changing zero values to NA
gtex_anmeans[gtex_anmeans == 0]<- NA
#log2 transforming df
gtex_anmeans<- log2(gtex_anmeans)
#row median centering by subtracting rowmedians
gtex_anmeans <- gtex_anmeans - rowMedian(gtex_anmeans, na.rm = TRUE)

#pearson correlation matric (pairwise deletion)
gtex_corr <- cor(gtex_anmeans, method = "pearson", use = "complete.obs")

#compute euclidean distances in matric and then create cluster object with Ward's hierarchical clustering
gtex_euc<- dist(gtex_corr)
gtex_clust<- hclust(gtex_euc, method = "ward.D2")

#pca plot for gtex norm + gtex hashimotos_________USING DESEQ + NON TPM_______________________

thyroid_deseq<- DESeqDataSetFromMatrix(gtexdf_ceil , gtexdf_pt, design= ~condition)
thyroid_vst <- vst(thyroid_deseq, blind = FALSE)
plotPCA(thyroid_vst, intgroup=c("condition"))

#PCA PLOT FOR gtex norm + thyroid hashimotos using tpm
thyroid_nids<- gsub('-','\\.',gtex_normids)
thyroid_aids<- gsub('-','\\.',gtex_autoids)
thyroid_p1 <-gtex_tpm[, colnames(gtex_tpm) %in% thyroid_nids]
thyroid_p2 <-gtex_tpm[, colnames(gtex_tpm) %in% thyroid_aids]


thyroid_dcols <- (grep(paste0( "^",thyroid_dids,"$",  collapse = "|"), colnames(gtex_tpm)))
thyroid_p<- cbind(thyroid_p1, thyroid_p2)

thyroid_tpca<- gtex_tpm[,colnames(gtex_tpm) %in% thyroid_dids]

thyroid_pca<- prcomp(thyroid_p)


plot(thyroid_pca)
ggbiplot(thyroid_pca, varname.size = 0, obs.scale = 1, var.scale = 1)
ggbiplot(thyroid_pca, obs.scale = 1, var.scale = 1)


thyroid_pc<- data.frame(t(thyroid_tpca))
b <- thyroid_pc[,apply(thyroid_pc, 2, var, na.rm=TRUE) != 0]
c<- prcomp(b)
ggbiplot(c)

GTEX_Thyroid_AS <-read.xlsx("/Users/shreyashekhar/GRIPS/GTEX_Thyroid_Ann.xlsx")
GTEX_Liver_AS <-read.xlsx("/Users/shreyashekhar/GRIPS/GTEX_Liver_Ann.xlsx")
GTEX_Prostate_AS <-read.xlsx("/Users/shreyashekhar/GRIPS/GTEX_Prostate_Ann.xlsx")
GTEX_Stomach_AS <-read.xlsx("/Users/shreyashekhar/GRIPS/GTEX_Stomach_Ann.xlsx")

#Datasets for Proteomic Clinical Data
GTEX_Proteomic_Path <-read.xlsx("/Users/shreyashekhar/GRIPS/PhaseIAllSample_pathology notes.xlsx")
colnames(GTEX_Proteomic_Path)<- GTEX_Proteomic_Path[1,]
GTEX_Proteomic_Path<- as.data.frame(GTEX_Proteomic_Path[-1,])

GTEX_Proteomic_Ann <-read.xlsx("/Users/shreyashekhar/GRIPS/GTEx_v7_Annotations_SampleAttributesDS.xlsx")

ProPath_ids<- intersect(colnames(gtexat_nzd), GTEX_Proteomic_Ann$SAMPID)
GTEX_Proteomic_Sub <-  GTEX_Proteomic_Ann[GTEX_Proteomic_Ann$SAMPID %in% ProPath_ids,]


colnames(gtexat) = gsub('[\\.]','-',colnames(gtexat))


#Make df of clinical descriptions and conditions for gtex dataset
#using thyroid, liver, prostate, stomach-- 
#one row for patient ids, description, and condition
patientids<- append(append(append(GTEX_Thyroid_AS$SAMPID, GTEX_Liver_AS$SAMPID), GTEX_Prostate_AS$SAMPID), GTEX_Stomach_AS$SAMPID)
tisdescriptions<- append(append(append(GTEX_Thyroid_AS$SMPTHNTS, GTEX_Liver_AS$SMPTHNTS), GTEX_Prostate_AS$SMPTHNTS), GTEX_Stomach_AS$SMPTHNTS)
sampcondition<- append(append(append(GTEX_Thyroid_AS$SMUBRID, GTEX_Liver_AS$Status), GTEX_Prostate_AS$Status), GTEX_Stomach_AS$Status)
tistype<- append(append(append(rep("Thyroid", 564), rep("Liver", 188)), rep("Prostate", 159)), rep("Stomach", 272))

classifier_df <- data.frame(patientids,tistype, tisdescriptions, sampcondition, stringsAsFactors=FALSE)
classifier_df$sampcondition[grepl(c("N|H"), classifier_df$sampcondition)] = "1"
classifier_df$sampcondition[!grepl(c("1"), classifier_df$sampcondition)] = "0"

write.csv(classifier_df, file = "/Users/shreyashekhar/GRIPS/classifier_df.csv")


#makes & returns df for input to DESEQ object based on annotation + whole gtex df(gtexat)
#df is all against healthy/norm, where 1 added to all part
#will need to split df to enter into EBSEQ
#full function included in createdeseq function
makecompdf = function(GTEX_AS){
  healthy_ids <- GTEX_AS$SAMPID[grepl("H", GTEX_AS$Status)]
  allids<- GTEX_AS$SAMPID
  
  alldf<- gtexat[,colnames(gtexat) %in% allids]
  colnames(alldf) <- paste( colnames(alldf), "1", sep = "")
  
  normdf<- gtexat[,colnames(gtexat) %in% healthy_ids]
  
  compdf<- cbind(alldf, normdf)
  
  
}

#GTEX ANNOTATION files with only 
GTEX_Sann<- GTEX_Ann[,colnames(GTEX_Ann) %in% c("SAMPID","SMPTHNTS", "SMTS")]


#Colon
GTEX_Colon_Ann<- subset(GTEX_Sann, SMTS =="Colon")
write.xlsx(GTEX_Colon_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Colon_Ann.xlsx", sheetName="GTEX_Colon_Ann")

#Liver
GTEX_Liver_Ann<- subset(GTEX_Sann, SMTS =="Liver")
write.xlsx(GTEX_Liver_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Liver_Ann.xlsx", sheetName="GTEX_Liver_Ann")


#Prostate
GTEX_Prostate_Ann<- subset(GTEX_Sann, SMTS =="Prostate")
write.xlsx(GTEX_Prostate_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Prostate_Ann.xlsx", sheetName="GTEX_Prostate_Ann")


#Skin
GTEX_Skin_Ann<- subset(GTEX_Sann, SMTS =="Skin")
write.xlsx(GTEX_Skin_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Skin_Ann.xlsx", sheetName="GTEX_Skin_Ann")



#Stomach
GTEX_Stomach_Ann<- subset(GTEX_Sann, SMTS =="Stomach")
write.xlsx(GTEX_Stomach_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Stomach_Ann.xlsx", sheetName="GTEX_Stomach_Ann")

#____________FUNCTIONS____________________________________________________

#DESEQ PL

#CREATE DESEQ OBJ For Comp
createdeseq = function(GTEX_AS, tissue){
  healthy_ids <- GTEX_AS$SAMPID[grepl("H", GTEX_AS$Status)]
  allids<- GTEX_AS$SAMPID
  
  #print healthy + all
  print(paste0("Healthy ids ", length(healthy_ids)))
  print(paste0("All ids ", length(allids)))
  
  
  alldf<- gtexat[,colnames(gtexat) %in% allids]
  allsum<-sum(colnames(gtexat) %in% allids)
  colnames(alldf) <- paste( colnames(alldf), "1", sep = "")
  
  normdf<- gtexat[,colnames(gtexat) %in% healthy_ids]
  normsum<- sum(colnames(gtexat) %in% healthy_ids)
  
  
  #print healthy + all found
  print(paste0("Healthy ids found ", normsum))
  print(paste0("All ids found ", allsum))
  
  compdf<- cbind(alldf, normdf)
  
  compdf_tmm<- calcNormFactors(compdf , method = "TMM", refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,Acutoff=-1e10)
  
  compdf_ceil <- ceiling(compdf)
  
  compdf_pt <- data.frame(colnames(compdf))
  compdf_pt$condition <- c(rep(paste0(tissue, "_all"), allsum),rep(paste0(tissue , "_norm"), normsum))
  compdf_deseq<- DESeqDataSetFromMatrix(compdf_ceil , compdf_pt, design= ~condition)
  sizeFactors(compdf_deseq)<- compdf_tmm
  
  return(compdf_deseq)
  
  
}

#GET DEG- DEres is a deseq object here
getDEG = function(DEres, comp1, comp2){
  
  DEres<- data.frame(DEres[order(DEres$padj),])
  
  DEres<- na.omit(DEres)
  
  #pdf with histogram of LFC and padj values
  pdf(paste0("DESeq plots", comp1, " vs ", comp2, ".pdf" ), height = 15, width = 5)
  
  hist(DEres$log2FoldChange, breaks = 60, col = "red")
  hist(DEres$padj, breaks = 50, col = "red")
  
  dev.off()
  
  LFC_t<- t.test(DEres$log2FoldChange, conf.level = 0.95)
  
  u_thresh = LFC_t$estimate + prod(sd(DEres$log2FoldChange, na.rm=TRUE),2)
  l_thresh = LFC_t$estimate - prod(sd(DEres$log2FoldChange, na.rm=TRUE),2)
  
  DEres_de <- DEres[DEres$log2FoldChange >u_thresh,]
  DEres_de<- rbind(DEres_de, DEres[DEres$log2FoldChange < l_thresh,])
  
  DEres_de <- DEres_de[DEres_de$padj<0.05,]
  
  return(DEres_de)
}

#IDENTIFY IMMUNE GENES
isimmune = function(DE_genes){
  
  iclist <- list()
  i=1
  for(g in DE_genes){
    icstr = ""
    #1
    if(g %in% T_vec){
      icstr = paste(icstr, "Tcell")
    }
    #2
    if(g %in% TCD8_vec){
      icstr = paste(icstr, "T+CD8cells")
    }
    #3
    if(g %in% Treg_vec){
      icstr = paste(icstr, "Tregcells")
    }
    #4
    if(g %in% B_vec){
      icstr = paste(icstr, "Bcells")
    }
    #5
    if(g %in% NK_vec){
      icstr = paste(icstr, "NKcells")
    }
    #6
    if(g %in% ClasMon_vec){
      icstr = paste(icstr, "ClasMon")
    }
    #7
    if(g %in% IntMon_vec){
      icstr = paste(icstr, "IntMon")
    }
    #8
    if(g %in% NonclasMon_vec){
      icstr = paste(icstr, "NonclasMon")
    }
    #9
    if(g %in% MatDC_vec){
      icstr = paste(icstr, "MatDC")
    }
    #10
    if(g %in% MonDC_vec){
      icstr = paste(icstr, "MonDC")
    }
    #11
    if(g %in% MatNeut_vec){
      icstr = paste(icstr, "MatNeut")
    }
    #12
    if(g %in% ImNeut_vec){
      icstr = paste(icstr, "ImNeut")
    }
    
    iclist[[i]] = icstr
    i = i +1
    
    
  }
  
  DE_df <- data.frame(DE_genes, unlist(iclist))
  return(DE_df)
  
  
}

#heatmaps
#use gtexat_nzd as nzdf, gtexdf_generef as generef, gtexat_pt as tispt
heatcomp = function(nzdf, DEG,GTEX_AS, generef, tispt, tissue){
  #finding healthy ids and making df of only those, adding 1 to colnames(ids)
  healthy_ids <- GTEX_AS$SAMPID[grepl("H", GTEX_AS$Status)]
  normsum<- sum(colnames(nzdf) %in% healthy_ids)
  
  compdf<- nzdf[,colnames(nzdf) %in% healthy_ids]
  colnames(compdf) <- paste( colnames(compdf), "1", sep = "")
  
  nzdf<-cbind(nzdf, compdf)
  
  rownums = unlist(lapply(DEG$DE_genes,function(x) grep( paste0("^", x ,"$"), gtexdf_generef$gene)))
  
  #
  tisname<- append(gtexat_pt$tissue, rep(paste0(tissue, "_norm"), normsum))
  
  nzdf<- nzdf[as.numeric(rownums),]
  colnames(nzdf)<- tisname
  #taking mean of all columns with same tissue name based on gene, putting in df
  tis_means<- data.frame(sapply(split.default(nzdf, names(nzdf)), rowMeans))
  
  rownames(tis_means)<- paste(DEG$DE_genes,DEG$unlist.iclist.)
  
  meta_tm<- data.frame(colnames(tis_means))
  colnames(meta_tm)<- "tissue"
  rownames(meta_tm)<- meta_tm$tissue
  
  
  pdf(paste0("pheatmap_gtex", tissue, "comp_ztest.pdf" ), height = 50, width = 7, pointsize = 1)
  
    pheatmap(tis_means, annotation_col=meta_tm, cellheight = 7, cluster_cols = FALSE, scale = "row")
  
  dev.off()
  
}


#heatmaps-- compare with another category from dataset
#use gtexat_nzd as nzdf, gtexdf_generef as generef, gtexat_pt as tispt
heatcomp3 = function(nzdf, DEG,GTEX_AS, generef, tispt, tissue, cat2){
  #finding healthy ids and making df of only those, adding 1 to colnames(ids)
  healthy_ids <- GTEX_AS$SAMPID[grepl("H", GTEX_AS$Status)]
  normsum<- sum(colnames(nzdf) %in% healthy_ids)
  
  compdf<- nzdf[,colnames(nzdf) %in% healthy_ids]
  colnames(compdf) <- paste( colnames(compdf), "1", sep = "")
  
  nzdf<-cbind(nzdf, compdf)
  
  #finding cat2 ids and making df of only those, adding 1 to colnames(ids)
  cat2_ids <- GTEX_AS$SAMPID[grepl(cat2, GTEX_AS$Status)]
  cat2sum<- sum(colnames(nzdf) %in% cat2_ids)
  
  cat2df<- nzdf[,colnames(nzdf) %in% cat2_ids]
  colnames(cat2df) <- paste( colnames(cat2df), "2", sep = "")
  
  nzdf<-cbind(nzdf, cat2df)
  
  #finding DEG genes rownums to extract in nzdf
  rownums = unlist(lapply(DEG$DE_genes,function(x) grep( paste0("^", x ,"$"), gtexdf_generef$gene)))
  
  #appending the metadf with norm + cat2 categories
  tisname<- append(gtexat_pt$tissue, rep(paste0(tissue, "_norm"), normsum))
  tisname<- append(tisname, rep(paste0(tissue, "_cat2"), cat2sum))
  
  nzdf<- nzdf[as.numeric(rownums),]
  colnames(nzdf)<- tisname
  #taking mean of all columns with same tissue name based on gene, putting in df
  tis_means<- data.frame(sapply(split.default(nzdf, names(nzdf)), rowMeans))
  
  rownames(tis_means)<- paste(DEG$DE_genes,DEG$unlist.iclist.)
  
  meta_tm<- data.frame(colnames(tis_means))
  colnames(meta_tm)<- "tissue"
  rownames(meta_tm)<- meta_tm$tissue
  
  
  pdf(paste0("pheatmap_gtex", tissue, "compcat2", cat2, "_ztest.pdf" ), height = 50, width = 7, pointsize = 1)
  
  pheatmap(tis_means, annotation_col=meta_tm, cellheight = 7, cluster_cols = FALSE, scale = "row")
  
  dev.off()
  
}

#FULL GTEX COMP BY Z SCORES:
gtex_fc<- gtexat_nzd
colnames(gtex_fc)<- gtexat_pt$tissue
gtex_fcmeans<- data.frame(sapply(split.default(gtex_fc, names(gtex_fc)), rowMeans))
tnames<- colnames(gtex_fcmeans)
gtex_fcscaled<- apply(gtex_fcmeans, 1, scale)
gtex_fcscaled<-data.frame(t(gtex_fcscaled))
colnames(gtex_fcscaled) <- tnames

#ADD THOSE WITH Z SCORES BELOW -4.5 AS WELL!!
#top genes for EACH tissue
adipose_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Adipose.Tissue>=4.5))
adrenal_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Adrenal.Gland>=4.5))
bladder_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Bladder>=4.5))
blood_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Blood >=4.5))
bloodvessel_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Blood.Vessel>=4.5))
brain_genes<-rownames(subset(gtex_fcscaled, gtex_fcscaled$Brain >=4.5))
breast_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Breast>=4.5))
cervixuteri_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Cervix.Uteri>=4.5))
colon_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Colon>=4.5))
esophagus_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Esophagus>=4.5))
fallopian_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Fallopian.Tube>=4.5))
heart_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Heart>=4.5))
kidney_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Kidney>=4.5))
liver_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Liver>=4.5))
lung_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Lung>=4.5))
muscle_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Muscle >=4.5))
nerve_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Nerve>=4.5))
ovary_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Ovary>=4.5))
pancreas_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Pancreas>=4.5))
pituitary_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Pituitary>=4.5))
prostate_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Prostate>=4.5))
salivary_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Salivary.Gland>=4.5))
skin_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Skin>=4.5))
smallintestine_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Small.Intestine>=4.5))
spleen_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Spleen>=4.5))
stomach_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Stomach>=4.5))
testis_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Testis>=4.5))
thyroid_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Thyroid>=4.5))
uterus_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Uterus>=4.5))
vagina_genes<- rownames(subset(gtex_fcscaled, gtex_fcscaled$Vagina>=4.5))

#making txt file with thyroid + other tissue genes
lapply(thyroid_genes, write, "thyroid_genes.txt", append=TRUE, ncolumns=1000)
lapply(liver_genes, write, "liver_genes.txt", append=TRUE, ncolumns=1000)
lapply(prostate_genes, write, "prostate_genes.txt", append=TRUE, ncolumns=1000)
lapply(stomach_genes, write, "stomach_genes.txt", append=TRUE, ncolumns=1000)


gtex_tisgenes<- list(adipose_genes, adrenal_genes, bladder_genes, blood_genes, bloodvessel_genes, brain_genes, breast_genes, cervixuteri_genes, colon_genes, esophagus_genes, fallopian_genes, heart_genes, kidney_genes, liver_genes,lung_genes, muscle_genes, nerve_genes, ovary_genes, pancreas_genes, pituitary_genes, prostate_genes, salivary_genes, skin_genes, smallintestine_genes, spleen_genes, stomach_genes, testis_genes, thyroid_genes, uterus_genes, vagina_genes)
names(gtex_tisgenes)<- c("adipose_genes", "adrenal_genes", "bladder_genes", "blood_genes", "bloodvessel_genes", "brain_genes", "breast_genes", "cervixuteri_genes", "colon_genes", "esophagus_genes", "fallopian_genes", "heart_genes", "kidney_genes", "liver_genes","lung_genes", "muscle_genes", "nerve_genes", "ovary_genes", "pancreas_genes", "pituitary_genes", "prostate_genes", "salivary_genes", "skin_genes", "smallintestine_genes", "spleen_genes", "stomach_genes", "testis_genes", "thyroid_genes", "uterus_genes", "vagina_genes")

#find unique genes for each tissue by comparing to other tissue genes and eliminating overlaps
#tlist <- gtex_tisgenes, tgenes is the tissue gene chr vector, and tname is that chr vector name
ugenes = function(tlist, tgenes, tname){
  
  tlist[paste0(tname)]<- NULL
  
  alltgenes<- unlist(tlist)
  
  tgenes_unique<- setdiff(tgenes, alltgenes)
  
  return(tgenes_unique)
}

adipose_ugenes<- ugenes(gtex_tisgenes, adipose_genes, "adipose_genes")
adrenal_ugenes<- ugenes(gtex_tisgenes, adrenal_genes, "adrenal_genes")
muscle_ugenes<- ugenes(gtex_tisgenes, muscle_genes, "muscle_genes")
testis_ugenes<- ugenes(gtex_tisgenes, testis_genes, "testis_genes")
prostate_ugenes<- ugenes(gtex_tisgenes, prostate_genes, "prostate_genes")


#FULL GTEX comp by z scores with modified: thyroid, liver, prostate, stomach 
gtex_hc<- gtexat_nzd
gtex_hcpt <- gtexat_pt
gtex_hcpt$colnames.gtexat. = gsub('[\\.]','-',gtex_hcpt$colnames.gtexat.)

#removing non-norm thyroid genes from both gtex_hc & gtex_hcpt -- 69
#diffid is reused as the different non-norm ids removed in each tissue
diffid<- setdiff(colnames(GTEX_genes_data), colnames(GTEX_normag))
gtex_hc<- gtex_hc[,colnames(gtex_hc) %ni% diffid]
gtex_hcpt<- subset(gtex_hcpt, gtex_hcpt$colnames.gtexat. %ni% diffid)

#removing non-norm liver genes from both gtex_hc & gtex_hcpt
diffid<- GTEX_Liver_AS$SAMPID[grepl("H", GTEX_Liver_AS$Status)]
diffid<- setdiff(GTEX_Liver_AS$SAMPID, diffid)
gtex_hc<- gtex_hc[,colnames(gtex_hc) %ni% diffid]
gtex_hcpt<- subset(gtex_hcpt, gtex_hcpt$colnames.gtexat. %ni% diffid)

#removing non-norm prostate genes from both gtex_hc & gtex_hcpt
diffid<- GTEX_Prostate_AS$SAMPID[grepl("H", GTEX_Prostate_AS$Status)]
diffid<- setdiff(GTEX_Prostate_AS$SAMPID, diffid)
gtex_hc<- gtex_hc[,colnames(gtex_hc) %ni% diffid]
gtex_hcpt<- subset(gtex_hcpt, gtex_hcpt$colnames.gtexat. %ni% diffid)

#removing non-norm stomach genes from both gtex_hc & gtex_hcpt
diffid<- GTEX_Stomach_AS$SAMPID[grepl("H", GTEX_Stomach_AS$Status)]
diffid<- setdiff(GTEX_Stomach_AS$SAMPID, diffid)
gtex_hc<- gtex_hc[,colnames(gtex_hc) %ni% diffid]
gtex_hcpt<- subset(gtex_hcpt, gtex_hcpt$colnames.gtexat. %ni% diffid)

#Z-score analysis
colnames(gtex_hc)<- gtex_hcpt$tissue
gtex_hcmeans<- data.frame(sapply(split.default(gtex_hc, names(gtex_hc)), rowMeans))
hcnames<- colnames(gtex_hcmeans)
gtex_hcscaled<- apply(gtex_hcmeans, 1, scale)
gtex_hcscaled<-data.frame(t(gtex_hcscaled))
colnames(gtex_hcscaled) <- hcnames

#tissue-specific genes for those four:
thyroid_hcgenes<- rownames(subset(gtex_hcscaled, gtex_hcscaled$Thyroid>=4.5))
liver_hcgenes<- rownames(subset(gtex_hcscaled, gtex_hcscaled$Liver>=4.5))
prostate_hcgenes<- rownames(subset(gtex_hcscaled, gtex_hcscaled$Prostate>=4.5))
stomach_hcgenes<- rownames(subset(gtex_hcscaled, gtex_hcscaled$Stomach>=4.5))

#txt file with hc thyroid + other tissue genes
lapply(thyroid_hcgenes, write, "thyroid_hcgenes.txt", append=TRUE, ncolumns=1000)
lapply(liver_hcgenes, write, "liver_hcgenes.txt", append=TRUE, ncolumns=1000)
lapply(prostate_hcgenes, write, "prostate_hcgenes.txt", append=TRUE, ncolumns=1000)
lapply(stomach_hcgenes, write, "stomach_hcgenes.txt", append=TRUE, ncolumns=1000)


#______________________________________________

#DESeq for all GTEX with GTEX_norm~ Thyroid

#TMM for all GTEX comparison with normal gtexnc(normcomp)df
gtexncdf<-cbind(GTEX_allag, GTEX_normag, by = 'gene')
rownames(gtexncdf)<- gtexncdf$gene
gtexncdf <- data.frame(gtexncdf[grepl(paste(c("GTEX", "TCGA", "expected"),collapse="|"), colnames(gtexncdf))])
gtexnc_tmm <- calcNormFactors(gtexncdf , method = "TMM", refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,Acutoff=-1e10)

gtexncdf_ceil <- ceiling(gtexncdf)

gtexncdf_pt <- data.frame(colnames(gtexncdf))
gtexncdf_pt$condition <- c(rep("GTEX_all", 446),rep("GTEX_norm", 377))
gtexncdf_deseq<- DESeqDataSetFromMatrix(gtexncdf_ceil , gtexncdf_pt, design= ~condition)
sizeFactors(gtexncdf_deseq)<- gtexnc_tmm

gtexncdf_deseq <- DESeq(gtexncdf_deseq, parallel = TRUE, BPPARAM=MulticoreParam(6))
gall_gn <- results(gtexncdf_deseq, contrast=c("condition","GTEX_all","GTEX_norm"), parallel = TRUE)

gtexncdf_vst <- varianceStabilizingTransformation(gtexncdf_deseq, blind=FALSE)





#DESeq for all GTEX with GTEX_healthy~ Liver
liver_deseq<- createdeseq(GTEX_Liver_AS, "Liver")
save(liver_deseq, file = "/Users/shreyashekhar/Downloads/liver_deseq.csv")

load("/Users/shreyashekhar/Downloads/liver_deseq.csv")
load("/Users/shreyashekhar/Downloads/liver_comp.csv")

liver_de<- getDEG(liver_comp, "Liver_all", "Liver_norm")

heatcomp(gtexat_nzd, isimmune(rownames(liver_de)), GTEX_Liver_AS, gtexdf_generef, gtexat_pt, "Liver")

#DESeq for all GTEX with GTEX_healthy~ Prostate
prostate_deseq<- createdeseq(GTEX_Prostate_AS, "Prostate")
save(prostate_deseq, file = "/Users/shreyashekhar/Downloads/prostate_deseq.csv")

load("/Users/shreyashekhar/Downloads/prostate_deseq.csv")
load("/Users/shreyashekhar/Downloads/prostate_comp.csv")

prostate_de<- getDEG(prostate_comp, "Prostate_all", "Prostate_norm")
heatcomp(gtexat_nzd, isimmune(rownames(prostate_de)), GTEX_Prostate_AS, gtexdf_generef, gtexat_pt, "Prostate")
p_deseq<- createdeseq(GTEX_Prostate_AS, "Prostate")

heatcomp3(gtexat_nzd, isimmune(rownames(prostate_de)), GTEX_Prostate_AS, gtexdf_generef, gtexat_pt, "Prostate", "P")
heatcomp3(gtexat_nzd, isimmune(rownames(prostate_de)), GTEX_Prostate_AS, gtexdf_generef, gtexat_pt, "Prostate", "A")


#DESeq for all GTEX with GTEX_healthy~ Stomach
stomach_deseq<- createdeseq(GTEX_Stomach_AS, "Stomach")
save(stomach_deseq, file = "/Users/shreyashekhar/Downloads/stomach_deseq.csv")

load("/Users/shreyashekhar/Downloads/stomach_deseq.csv")
load("/Users/shreyashekhar/Downloads/stomach_comp.csv")

stomach_de<- getDEG(stomach_comp, "Stomach_all", "Stomach_norm")
heatcomp(gtexat_nzd, isimmune(rownames(stomach_de)), GTEX_Stomach_AS, gtexdf_generef, gtexat_pt, "Stomach")


#save them, send to cluster, do deseq + results, return data


