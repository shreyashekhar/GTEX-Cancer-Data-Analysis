echofn <-function(v) {
  deparse(substitute(v))
}

library(plotROC)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(glmnet)
library(stringr)
library(openxlsx)
library(plyr)
library(EBSeq)
library(grid)
library(edgeR)
library(dplyr)
library(DESeq2)
library(BiocParallel)
library(pheatmap)
library(GenomicRanges)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(DOSE)
library(biomaRt)
library(Hmisc)
library(devtools)
library(ggbiplot)
library(vqv)


THCA.ge <- read.table("/Users/shreyashekhar/Downloads/HiSeqV2_THCA.txt",row.names = 1,header=TRUE)
THCA.pt <- read.table("/Users/shreyashekhar/Downloads/THCA_clinicalMatrix.txt",row.names = 1,header=TRUE,sep='\t')
colnames(THCA.ge) = gsub('[\\.]','-',colnames(THCA.ge))
A = intersect(colnames(THCA.ge),rownames(THCA.pt))
CM = THCA.pt[A,]
THCA_Clin.ge <- data.frame(THCA.ge[,A])
colnames(THCA_Clin.ge) = gsub('[\\.]','-',colnames(THCA_Clin.ge))

path_m = CM$pathologic_M
path_n = CM$pathologic_N
path_df <- data.frame(path_m, path_n)
row.names(path_df)<-A
path_df$path_clas='0'
#based on M1 & N1b
path_df$path_clas[grepl("M1", path_df$path_m)] = '1'
path_df$path_clas[grepl("N1b", path_df$path_n)] = '1'

path_df$path_clas_m='0'
path_df$path_clas_n='0'
#based on only M1
path_df$path_clas_m[grepl("M1", path_df$path_m)] = '1'
#based on only N1b
path_df$path_clas_n[grepl("N1b", path_df$path_n)] = '1'

#based on N1b with M1 samples removed
pm = CM$pathologic_M
pn = CM$pathologic_N
npath_df <- data.frame(pm, pn)
row.names(npath_df) <- A
#find & store rownames of rows with pm = M1
npath_del = rownames(subset(npath_df, pm =="M1"))
#delete all rows in npath_del(have M1)
npath_df= npath_df[!rownames(npath_df) %in% npath_del,]
#Make .ge df without colnames in npath_del
`%ni%` <- Negate(`%in%`)
THCA_N.ge = subset(THCA_Clin.ge,select = names(THCA_Clin.ge) %ni% npath_del)

npath_df$nclas = '0'
npath_df$nclas[grepl("N1b", npath_df$pn)] = '1'

THCA_iso_norm.ge <- read.table("/Users/shreyashekhar/Downloads/newnormal_thyroid.isoforms.results.txt",row.names = 1,header=TRUE)
THCA_iso_tcga.ge <- read.table("/Users/shreyashekhar/Downloads/THCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt",row.names = 1,header=TRUE, sep='\t', stringsAsFactors = FALSE)
THCA_isoform_id <- read.table("/Users/shreyashekhar/Downloads/isoformids.txt",col.names = 1, header=TRUE, sep='\t')
THCA_transtogene <- read.csv("/Users/shreyashekhar/Downloads/transtogenexl.csv",header=TRUE)
HK_Genes_Data <-read.csv("/Users/shreyashekhar/Downloads/housekeeping_genes_data.csv",header=TRUE)
GTEX_EX_Count <- read.table("/Users/shreyashekhar/Downloads/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_expected_count.txt",row.names = 1,header=TRUE)
#THCA_norm <- read.csv("/Users/shreyashekhar/Downloads/newnormal_thyroid.genes.fpkm.csv",header=TRUE, stringsAsFactors = FALSE)
#THCA_tum <- read.csv("/Users/shreyashekhar/Downloads/newthyroid_tumor.genes.fpkm.csv",header=TRUE, stringsAsFactors = FALSE)
GTEX_Ann <- read.csv("/Users/shreyashekhar/Downloads/GTEx_v7_Annotations_SampleAttributesDS.csv",header=TRUE, stringsAsFactors = FALSE)
GTEX_thyroid_Grouped <- read.csv("/Users/shreyashekhar/Downloads/GTEX_thyroid_Grouped.csv",header=TRUE, stringsAsFactors = FALSE)
GTEX_genes_data<- read.csv("/Users/shreyashekhar/Downloads/gtex_genes_dataxl.csv",header=TRUE, stringsAsFactors = FALSE)
imgenes_stlength <-read.csv("/Users/shreyashekhar/Downloads/imgenes_efflength.csv",header=TRUE, stringsAsFactors = FALSE)
THCA_normex <- read.csv("/Users/shreyashekhar/Downloads/newnormal_thyroid.genes.exp.csv",header=TRUE, stringsAsFactors = FALSE)
THCA_tumex <- read.csv("/Users/shreyashekhar/Downloads/newthyroid_tumor.genes.exp.csv",header=TRUE, stringsAsFactors = FALSE)
THCA_tcga_raw<- read.csv("/Users/shreyashekhar/Downloads/THCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.csv",header=TRUE, stringsAsFactors = FALSE)
#load("/Users/shreyashekhar/Downloads/gtex_alldata.csv")


#THCA_tcga_raw formatting
colnames(THCA_tcga_raw) = gsub('[\\.]','-',colnames(THCA_tcga_raw))
#removing "scaled_estimate" and "transcript_id" columns
THCA_tcga_raw = subset(THCA_tcga_raw,select = THCA_tcga_raw[1,] %ni% c("scaled_estimate", "transcript_id"))
#finding and removing patients in npath_del, who have M1 status
coldels <- (grep(paste0( npath_del,collapse = "|"), colnames(THCA_tcga_raw)))
THCA_tcga_raw<- THCA_tcga_raw[,-coldels]

tcga_raw = THCA_tcga_raw
tcga_raw<- tcga_raw[-1,]
rownames(tcga_raw)<- tcga_raw$`Hybridization-REF`
tcga_raw<- tcga_raw[,-1]
tcga_raw<- data.frame(apply(tcga_raw, 2, as.numeric))
rawg = THCA_tcga_raw$`Hybridization-REF`[-1]
rownames(tcga_raw)<- rawg
tcga_raw_genes = sub("\\|.*","", rownames(tcga_raw))
tcga_raw<- cbind(gene = tcga_raw_genes, tcga_raw)

THCA_tcga<- immunesep(tcga_raw)



# rawg = THCA_tcga_raw$`Hybridization-REF`[-1]
# rownames(tcga_raw)<- rawg
# tcga_raw <-1000*tcga_raw
# tcga_raw<- log2(tcga_raw+1)
# tcga_raw<- cbind(gene = rownames(tcga_raw), tcga_raw)

T_vec = "CD3G, CD3D, TRA, CD6, CD5, NPDC1, CD28, CAMK4, GFI1, GATA3, SH2D1A, TRB, TNFRSF25, NK4, TACTILE, BCL11B, CD3E, INPP4B, MAL, NPDC1, ITM2A, ITK, LCK, NFATC3, RORA, MGC19764, TCF7, ZAP70, LEF1, SPOCK2, PRKCQ, SATB1, RASGRP1, LRIG1, DPP4, CD3Z, PDE4D, FYN, WWP1, LAT, DUSP16, KIAA0748, CDR2, STAT4, FLT3LG, IL6ST"
#use gsub to remove all white spaces in the string (sep can only be 1 char, which is ",")
T_vec = gsub(" ", "", T_vec, fixed = TRUE)
#print T_vec and paste output into TextCon.. reassign T_vec to the char vec (only sep is "," in text)
T_vec = scan(what=character(0), sep=",", file=textConnection(T_vec))
grows <- (grep(paste0( "^", T_vec,"$",  collapse = "|"), tcga_raw$gene))
tcga_rawcomp <-tcga_raw[grows,]
#ig_df =subset(ag_df,rownames(ag_df) %in% T_vec)
tcga_rawcomp$ic = "Tcell"




analyzeFeature2 = function(expr, feat, pos, neg, biostr, ctype) {
  s = apply(expr,1,sd)
  #p = apply(expr[s>1,],1,function(x) t.test(x[feat %in% pos],x[feat %in% neg])$statistic)
  #p = apply(expr[s>1,],1,function(x) t.test(x[str_detect(feat, pos)],x[str_detect(feat, neg)])$statistic)
  t.results = apply(expr[s>1,],1,function(x) t.test(x[grepl(paste(pos, collapse  = "|"), feat)],
                                                    x[grepl(paste(neg, collapse = "|"), feat)]))
  #p = apply(expr[s>1,],1,function(x) t.test(x[grepl(paste(pos, collapse  = "|"), feat)],
  #                                         x[grepl(paste(neg, collapse = "|"), feat)])$statistic)
  
  #p.value = apply(expr[s>1,],1,function(x) t.test(x[grepl(paste(pos, collapse  = "|"), feat)],
  #                                                x[grepl(paste(neg, collapse = "|"), feat)])$p.value)
  
  t.statistic = unlist(lapply(t.results,FUN=function(x) x$statistic)) 
  p.value = unlist(lapply(t.results,FUN=function(x) x$p.value))
  
  #cb(sort(p.value))
  
  p.adj = p.adjust(p.value, "fdr")
  
  # QQ plot
  pdf(paste0("qqplot_", biostr, "_", ctype, ".pdf"))
  o = -log10(sort(p.value,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(x= e, y= o, xlab = "Expected p-value", 
       ylab = "Observed p-value")
  dev.off()
  
  #ggplot() + geom_point(aes(x=p.adj, y= p.value)) + geom_abline(intercept = 0, slope = 1) +
  #  labs(x = "Observed p-value (p.value)", y = "Expected p-value (p.adjusted)") 
  #ggplot() +geom_point(aes(x= -log10(p.adj), y= -log10(p.value))) + geom_abline(intercept = 0, slope = 1)
  #ggsave(paste0("QQplot_", biostr, "_", ctype, ".pdf"))
  
  #cb(sort(p.adj))
  #cb(sort(t.statistic))
  #cb(sort(t.statistic, decreasing=TRUE))
  
  pl <- list()
  i = 1
  
  t.stats = data.frame(t.statistic)
  p.val = data.frame(p.value)
  p.adj_df = data.frame(p.adj)
  colnames(p.val) = 'p'
  p.val$t = t.stats$t.statistic
  #p.val$t = abs(t.stats$t.statistic)
  p.val$padj = p.adj_df$p.adj
  p.val = p.val[order(p.val$padj),]
  
  
  #temp_df <- data.frame(head(sort(p,decreasing = TRUE)))
  #qval_df <- data.frame(head(sort(p.adj)))
  #temp_df <- data.frame(head(sort(,decreasing = TRUE)))
  #temp_p_adj_df <- data.frame(head(sort(p.adj)))
  #temp_p_val_df <- data.frame(head(sort(p.value)))
  #temp_df <- data.frame(head(sort(p)))
  #colnames(temp_df) = 'x'
  #temp_df$y = rownames(temp_df)
  #temp_df$y=sapply(strsplit(rownames(temp_df),".",fixed = TRUE),"[[",1)
  #temp_df$z <- factor(temp_df$y, levels = temp_df$y[order(-temp_df$x)])
  temp_df = head(p.val)
  temp_df$r = rownames(temp_df)
  temp_df$z <- factor(temp_df$r, levels = temp_df$r[order(temp_df$padj)])
  #pl[[i]] = ggplot(temp_df,aes(x=z,y=x, fill = y)) +geom_bar(stat='identity') +
  pl[[i]] = ggplot(temp_df,aes(x=z,y=t, fill = r)) +geom_bar(stat='identity') +
    labs(x = "Genes", y = "T-statistic", fill = "Genes") +
    #labs(x = "Genes", y = "Q-value", fill = "Genes") +
    ggtitle("T-test Statistics") +
    #ggtitle("Q-values") +
    theme(axis.text.x=element_blank(), axis.title=element_text(size=9),
          #axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(size=10),
          legend.position="none")
  #ggsave(paste("ggplot_ttest_", biostr, ".pdf"))
  i = i+1
  
  #print the top genes save as pdf
  #top_genes = c(temp_df$y)
  #top_genes = c(temp_df$r)
  top_genes = c(temp_df$r)
  
  #pdf(paste("roc_plots_", biostr, ".pdf"))
  
  #plt <- ggplot()
  #genes <-list()
  auc_list <- as.numeric(list())
  auc_glist <- as.numeric(list())
  tmp_df <- data.frame(feat)
  colnames(tmp_df) = 'level'
  tmp_df$level0[grepl(paste(pos, collapse  = "|"), tmp_df$level)] = 'Positive'
  tmp_df$level0[grepl(paste(neg, collapse  = "|"), tmp_df$level)] = 'Negative'
  names(tmp_df) [ncol(tmp_df)] <- 'lvl0'
  
  #cb(tmp_df$lvl0)
  
  samples <- data.frame(tmp_df$lvl0)
  sample_counts <- table(samples, exclude=NULL)
  print(sample_counts)
  samp <- data.frame(sample_counts)
  group.colors <- c("Positive" = "aquamarine2", "Negative" = "deepskyblue1")
  #group.colors <- c('Positive' = "#300BFF", 'Negative' = "#CC6600", 
  #                                  'Indeterminate'="#9633FF", 
  #                                  'Equivocal' = "#E2FF33", 
  #                                  #' ' = "#000000",
  #                                  #' ' = "#FF0000",
  #                                  #' ' = "red",
  #                                  na.value = "#F0e664") 
  #E3DB71")
  pl[[i]] = ggplot(samp, aes(samples, Freq, fill = samples))+ geom_bar(stat="identity")+
    ggtitle(paste(biostr, "- samples")) +
    labs(y = "Num. of samples") +
    theme(axis.text.x=element_blank(), 
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          #axis.title=element_text(size=9),
          #legend.position="none",
          legend.position="bottom",
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_text(size=8), 
          legend.text = element_text(size=7),
          plot.title = element_text(size=10)) + 
    scale_fill_manual(values=group.colors, na.value = "yellow") #+
  #scale_color_discrete(na.value="yellow")
  #ggsave(paste("barplot_", biostr,".pdf"))
  i=i+1
  
  j=1
  #tmp_df$level1[tmp_df$level %in% neg] = 0
  #tmp_df$level1[tmp_df$level %in% pos] = 1
  tmp_df$level1[grepl(paste(pos, collapse  = "|"), tmp_df$level)] = 1
  tmp_df$level1[grepl(paste(neg, collapse  = "|"), tmp_df$level)] = 0
  
  names(tmp_df) [ncol(tmp_df)] <- 'lvl'
  
  #cb(tmp_df$lvl)
  
  for (g in top_genes) {
    #x = as.numeric(expr[g,])
    #a = roc(feat, x, levels = c("Positive", "Negative"))
    #plot(roc(feat, x, levels = c("Positive", "Negative")), main = list(paste(g, "~", echofn(her2), "\n(AUC: ", a$auc, ")"), cex = 1, col = "blue", font = 1.5))
    #print(g)
    #tmp_df <- data.frame(feat)
    #colnames(tmp_df) = 'level'
    tmp_df$gene = as.numeric(expr[g,])
    names(tmp_df) [ncol(tmp_df)] <- g
    tmp_df$tgene = as.numeric(expr[g,])
    #tmp_df$level1[tmp_df$level %in% neg] = 0
    #tmp_df$level1[tmp_df$level %in% pos] = 1
    #names(tmp_df) [ncol(tmp_df)] <- paste0(g, "-lvl")
    #a = roc(feat, x, levels = c("Positive", "Negative"))
    #genes[[i]] <- tmp_df
    
    #if (i==1) {
    #plt = ggplot() + geom_roc(data = tmp_df, aes(x=level1, y=gene, color = gene))
    #}
    #else {
    #plt = plt + geom_roc(data = tmp_df, aes(x=level1, y=gene, color = gene))
    #}
    
    #pl[[i]] <- ggplot(tmp_df, aes(d=level1, m=gene)) + geom_roc() 
    #+ labs(title = paste(g, "~", echofn(her2)), subtitle = paste("AUC: ", a$auc)) 
    #pl[[i]]<-  pl[[i]]    + style_roc(theme = theme_grey) +
    #theme(axis.text = element_text(colour = "blue")) +
    #ggtitle(paste(g, "~", echofn(her2), "(AUC: ", round(calc_auc(pl[[i]])$AUC, 4), ")"))  
    # + annotate("text", x = .75, y = .25, 
    #         label = "AUC") +
    #scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))
    
    plt <- ggplot(tmp_df, aes(d=lvl, m= tgene)) + geom_roc() 
    auc_list[j] <- round(calc_auc(plt)$AUC, 4)
    auc_glist[g] <- round(calc_auc(plt)$AUC, 4)
    if (auc_list[j] < 0.5) {
      auc_list[j] = 1- auc_list[j]
      #  tmp_df$gene = 
    }
    print(auc_list[j])
    j = j+1
    #plt<-  plt    + style_roc(theme = theme_grey) +
    #theme(axis.text = element_text(colour = "blue")) +
    #ggtitle(paste(g, "~", echofn(her2), "(AUC: ", round(calc_auc(pl[[i]])$AUC, 4), ")"))  
    
    #plt <- ggplot() + geom_roc(data = tmp_df, aes(x=level1, y=gene, color = gene)) + geom_roc() 
    #+ labs(title = paste(g, "~", echofn(her2))) 
    #subtitle = paste("AUC: ", a$auc)) + 
    #+ theme (plot.title = element_text(color="blue", size=10, face="bold.italic"))
    #plot.subtitle = element_text(color="blue", size=8, face="bold.italic"))
    #ggtitle(paste(top_genes[g], "~", echofn(her2), "\n(AUC: ", a$auc, ")")) 
    #i = i+1
  }
  #longtest <- melt_roc(tmp_df, "lvl", c(temp_df$y))
  longtest <- melt_roc(tmp_df, "lvl", c(temp_df$r))
  for (g in top_genes) {
    if (auc_glist[g] < 0.5) {
      longtest <- within(longtest, D[D == 1 & name %in% g ] <- 2)
      longtest <- within(longtest, D[D == 0 & name %in% g ] <- 1)
      longtest <- within(longtest, D[D == 2 & name %in% g ] <- 0)
    }
    print(auc_glist[g])
  }  
  pl[[i]] = ggplot(longtest, aes(d = D, m = M, color = name)) + 
    #ggplot(longtest, aes(x = D, y = M, label = c, colour = name)) +
    geom_roc(labels = FALSE) + style_roc(theme = theme_gray) +
    theme(axis.text = element_text(colour = "blue"), 
          axis.title=element_text(size=9),
          legend.position="none",
          plot.title = element_text(size=10)) +
    ggtitle(paste(biostr, "- ROC Curves")) + 
    #annotate("text", x = .75, y = .25, 
    #label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
    # scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))
    scale_x_continuous(breaks = seq(0, 1, by = .2)) 
  #ggsave(paste("rocCurves_", biostr, ".pdf"))
  i=i+1
  
  
  
  auc_df = data.frame(auc_list)
  #print (auc_list)
  #colnames(auc_df) = c(temp_df$y)
  #colnames(auc_df) = c(temp_df$r)
  #auc_tbl <- melt(auc_df)
  #auc_tbl$y = c(temp_df$y)
  #auc_tbl$y = c(temp_df$r)
  colnames(auc_df) = 'x'
  auc_df$y = c(temp_df$r)
  auc_df = auc_df[order(-auc_df$x),]
  auc_df$z <- factor(auc_df$y, levels = auc_df$y[order(-auc_df$x)])
  print(auc_df)
  #pl[[i]] = ggplot(auc_tbl, aes(variable, value, fill = y))+ geom_bar(stat="identity")+
  #pl[[i]] = ggplot(auc_df, aes(x=y, y=x, fill = y))+ geom_bar(stat="identity")+
  # pl[[i]] = ggplot(auc_df, aes(y, x, fill = y))+ geom_bar(stat="identity")+
  pl[[i]] = ggplot(auc_df, aes(z, x, fill = z))+ geom_bar(stat="identity")+
    ggtitle(paste(biostr, " AUC values")) + labs(x = "Genes", y = "AUC", fill = "Genes") +
    theme(axis.text.x=element_blank(), #axis.title=element_text(size=9), 
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          #axis.text=element_text(size=10),
          legend.position="bottom", 
          #legend.position="none",
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size=8), 
          legend.text = element_text(size=7),
          plot.title = element_text(size=10))
  #ggsave(paste("auc_barplot_", biostr, ".pdf"))
  
  i=i+1
  
  #temp3_df <- data.frame(head(sort(p,decreasing = TRUE), 3))
  #temp3_df <- data.frame(head(sort(p.adj), 3))
  #temp3_df <- data.frame(head(sort(p), 3))
  #colnames(temp3_df) = 'x'
  #temp3_df$y = rownames(temp3_df)
  
  temp3_df <- data.frame(head(auc_df, 3))
  print(temp3_df)
  #temp_df$z <- factor(temp_df$y, levels = temp_df$y[order(-temp_df$x)])
  #btest <- melt(tmp_df,id.vars = c("level","lvl"), measure.vars = c(temp3_df$y))
  btest <- melt(tmp_df,id.vars = c("lvl0","lvl"), measure.vars = c(temp3_df$y))
  pl[[i]] <- ggplot(btest, aes(x=factor(lvl0),y=value, fill = factor(lvl0))) +
    geom_boxplot(outlier.size=0,show.legend = F,lwd=0.1) + labs(title="Box plots for Top 3 Genes") +facet_wrap(~variable) +
    theme(axis.text.x=element_blank(),
          #axis.title=element_text(size=9),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.title = element_text(size=6), 
          legend.text = element_text(size=4),
          plot.title = element_text(size=10),
          legend.position="bottom") +
    scale_fill_manual(values=group.colors, na.value="yellow")
  i=i+1
  
  #plt<-  plt    + style_roc(theme = theme_grey) +
  #  theme(axis.text = element_text(colour = "blue")) +
  #  ggtitle(paste(g, "~", echofn(her2), "(AUC: ", round(calc_auc(plt$AUC, 4), ")")))  
  #grid.arrange(grobs=pl, ncol=3)
  #ml <- marrangeGrob(pl, nrow=2, ncol=2, top= paste(biostr, " - Analysis"))
  #leg_pl <- pl[[2]]
  #legend <- get_legend(leg_pl)
  
  lay <- rbind(c(2,1),
               c(3,4),
               c(5,5))
  ml = grid.arrange(grobs = pl, layout_matrix = lay, top= paste(biostr, " in ", ctype)) 
  #legend = legend)
  
  #gs <- arrangeGrob(grobs = pl, layout_matrix = lay)
  #ml <- marrangeGrob(grobs=pl, nrow=2, ncol=2, top= paste(biostr, "in ", ctype))
  ggsave(paste0("plots_", biostr, "_in_", ctype, ".pdf"), plot=ml, width = 6, height = 8, units = "in")
  #genes_df <- data.frame(genes)
  #plt = ggplot() + geom_roc(data = genes_df, aes(x=level1, y=gene, group = genes, color = genes))
  #print(plt)
  #do.call(grid.arrange, pl)
  #dev.off()
  
} 



mlasso = function(expr, feat, feat_str){
  data  <-expr
  
  data <- data.frame(t(data))
  
  #making data set
  feat_bin <- as.numeric( feat == "1" )
  
  data$feat_str <- feat_bin
  feat_target <- data$feat_str
  
  #training + testing sets
  feat_trsize<- floor(0.75 * nrow(data))
  feat_tr <- sample(seq_len(nrow(data)), size = feat_trsize)
  feat_train <- data[feat_tr, ]
  feat_test <- data[-feat_tr, ]
  
  #training model
  feat_responsecol <- which(colnames(feat_train) == "feat_str")
  feat_trainxdm <- data.matrix(feat_train[, -feat_responsecol])
  feat_lasso_model <- cv.glmnet(x = feat_trainxdm, y = feat_train$feat_str, alpha = 1)
  
  #test set predictions
  feat_testxdm <- data.matrix(feat_test[, -feat_responsecol])
  
  feat_predictions_lasso <- predict(feat_lasso_model, newx = feat_testxdm, type = "response", 
                                    s = "lambda.min")[, 1]
  
  #Model performance (AUC)
  auc(feat_test$feat_str, feat_predictions_lasso)
}

#pathologic M1 + N1b
analyzeFeature2(THCA_Clin.ge,path_df$path_clas, c("1"), c("0"), echofn(THCA-Metastasis), echofn(THCA-Pathological))
#pathologic M1 only
analyzeFeature2(THCA_Clin.ge,path_df$path_clas_m, c("1"), c("0"), echofn(THCA-Metastasis-M), echofn(THCA-Pathological))
#pathologic N1b only
analyzeFeature2(THCA_Clin.ge,path_df$path_clas_n, c("1"), c("0"), echofn(THCA-Metastasis-N), echofn(THCA-Pathological))
#Pathologic N1b without M1 samples included
analyzeFeature2(THCA_N.ge,npath_df$nclas, c("1"), c("0"), echofn(THCA-Met-N1b), echofn(THCA-Pathological))


library(gridExtra)

analyzeFeature3 = function(expr, feat, pos, neg, biostr, ctype, gene_vec) {
  s = apply(expr,1,sd)
  #p = apply(expr[s>1,],1,function(x) t.test(x[feat %in% pos],x[feat %in% neg])$statistic)
  #p = apply(expr[s>1,],1,function(x) t.test(x[str_detect(feat, pos)],x[str_detect(feat, neg)])$statistic)
  t.results = apply(expr,1,function(x) t.test(x[grepl(paste(pos, collapse  = "|"), feat)],
                                                    x[grepl(paste(neg, collapse = "|"), feat)]))
  #p = apply(expr[s>1,],1,function(x) t.test(x[grepl(paste(pos, collapse  = "|"), feat)],
  #                                         x[grepl(paste(neg, collapse = "|"), feat)])$statistic)
  
  #p.value = apply(expr[s>1,],1,function(x) t.test(x[grepl(paste(pos, collapse  = "|"), feat)],
  #                                                x[grepl(paste(neg, collapse = "|"), feat)])$p.value)
  
  t.statistic = unlist(lapply(t.results,FUN=function(x) x$statistic)) 
  p.value = unlist(lapply(t.results,FUN=function(x) x$p.value))
  
  #cb(sort(p.value))
  
  p.adj = p.adjust(p.value, "fdr")
  
  # QQ plot
  pdf(paste0("qqplot_", biostr, "_", ctype, ".pdf"))
  o = -log10(sort(p.value,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(x= e, y= o, xlab = "Expected p-value", 
       ylab = "Observed p-value")
  dev.off()
  
  #ggplot() + geom_point(aes(x=p.adj, y= p.value)) + geom_abline(intercept = 0, slope = 1) +
  #  labs(x = "Observed p-value (p.value)", y = "Expected p-value (p.adjusted)") 
  #ggplot() +geom_point(aes(x= -log10(p.adj), y= -log10(p.value))) + geom_abline(intercept = 0, slope = 1)
  #ggsave(paste0("QQplot_", biostr, "_", ctype, ".pdf"))
  
  #cb(sort(p.adj))
  #cb(sort(t.statistic))
  #cb(sort(t.statistic, decreasing=TRUE))
  
  pl <- list()
  i = 1
  
  t.stats = data.frame(t.statistic)
  p.val = data.frame(p.value)
  p.adj_df = data.frame(p.adj)
  colnames(p.val) = 'p'
  p.val$t = t.stats$t.statistic
  #p.val$t = abs(t.stats$t.statistic)
  p.val$padj = p.adj_df$p.adj
  #p.val = p.val[order(p.val$padj),]
  
  #CREATE DATA FRAME FOR ALL GENES-include t-stat & p-val
  Genes_Info = data.frame(t.stats$t.statistic, p.val$p)
  rownames(Genes_Info) = rownames(p.val)
  colnames(Genes_Info) = c("T-statistic", "P-value")
  #temp_df <- data.frame(head(sort(p,decreasing = TRUE)))
  #qval_df <- data.frame(head(sort(p.adj)))
  #temp_df <- data.frame(head(sort(,decreasing = TRUE)))
  #temp_p_adj_df <- data.frame(head(sort(p.adj)))
  #temp_p_val_df <- data.frame(head(sort(p.value)))
  #temp_df <- data.frame(head(sort(p)))
  #colnames(temp_df) = 'x'
  #temp_df$y = rownames(temp_df)
  #temp_df$y=sapply(strsplit(rownames(temp_df),".",fixed = TRUE),"[[",1)
  #temp_df$z <- factor(temp_df$y, levels = temp_df$y[order(-temp_df$x)])
  #temp_df = head(p.val)
  #temp_df$r = rownames(temp_df)
  #temp_df$z <- factor(temp_df$r, levels = temp_df$r[order(temp_df$padj)])
  #pl[[i]] = ggplot(temp_df,aes(x=z,y=x, fill = y)) +geom_bar(stat='identity') +
  #pl[[i]] = ggplot(temp_df,aes(x=z,y=t, fill = r)) +geom_bar(stat='identity') +
    #labs(x = "Genes", y = "T-statistic", fill = "Genes") +
    #labs(x = "Genes", y = "Q-value", fill = "Genes") +
    #ggtitle("T-test Statistics") +
    #ggtitle("Q-values") +
    #theme(axis.text.x=element_blank(), axis.title=element_text(size=9),
          #axis.title.x=element_blank(),
          #axis.ticks.x=element_blank(),
          #plot.title = element_text(size=10),
          #legend.position="none")
  #ggsave(paste("ggplot_ttest_", biostr, ".pdf"))
  #i = i+1
  
  #print the top genes save as pdf
  #top_genes = c(temp_df$y)
  #top_genes = c(temp_df$r)
  #top_genes = c(temp_df$r)
  all_genes = rownames(p.val)
  #pdf(paste("roc_plots_", biostr, ".pdf"))
  
  #plt <- ggplot()
  #genes <-list()
  auc_list <- as.numeric(list())
  auc_glist <- as.numeric(list())
  tmp_df <- data.frame(feat)
  colnames(tmp_df) = 'level'
  tmp_df$level0[grepl(paste(pos, collapse  = "|"), tmp_df$level)] = 'Positive'
  tmp_df$level0[grepl(paste(neg, collapse  = "|"), tmp_df$level)] = 'Negative'
  names(tmp_df) [ncol(tmp_df)] <- 'lvl0'
  
  #cb(tmp_df$lvl0)
  
  #samples <- data.frame(tmp_df$lvl0)
  #sample_counts <- table(samples, exclude=NULL)
  #print(sample_counts)
  #samp <- data.frame(sample_counts)
  #group.colors <- c("Positive" = "aquamarine2", "Negative" = "deepskyblue1")
  #group.colors <- c('Positive' = "#300BFF", 'Negative' = "#CC6600", 
  #                                  'Indeterminate'="#9633FF", 
  #                                  'Equivocal' = "#E2FF33", 
  #                                  #' ' = "#000000",
  #                                  #' ' = "#FF0000",
  #                                  #' ' = "red",
  #                                  na.value = "#F0e664") 
  #E3DB71")
  #pl[[i]] = ggplot(samp, aes(samples, Freq, fill = samples))+ geom_bar(stat="identity")+
    #ggtitle(paste(biostr, "- samples")) +
    #labs(y = "Num. of samples") +
    #theme(axis.text.x=element_blank(), 
         #axis.title.x=element_blank(),
          #axis.ticks.x=element_blank(),
          #axis.title=element_text(size=9),
          #legend.position="none",
          #legend.position="bottom",
          #legend.key.size = unit(0.3, "cm"),
          #legend.title = element_text(size=8), 
          #legend.text = element_text(size=7),
          #plot.title = element_text(size=10)) + 
    #scale_fill_manual(values=group.colors, na.value = "yellow") #+
  #scale_color_discrete(na.value="yellow")
  #ggsave(paste("barplot_", biostr,".pdf"))
  #i=i+1
  
  j=1
  #tmp_df$level1[tmp_df$level %in% neg] = 0
  #tmp_df$level1[tmp_df$level %in% pos] = 1
  tmp_df$level1[grepl(paste(pos, collapse  = "|"), tmp_df$level)] = 1
  tmp_df$level1[grepl(paste(neg, collapse  = "|"), tmp_df$level)] = 0
  
  names(tmp_df) [ncol(tmp_df)] <- 'lvl'
  
  #cb(tmp_df$lvl)
  
  for (g in all_genes) {
    #x = as.numeric(expr[g,])
    #a = roc(feat, x, levels = c("Positive", "Negative"))
    #plot(roc(feat, x, levels = c("Positive", "Negative")), main = list(paste(g, "~", echofn(her2), "\n(AUC: ", a$auc, ")"), cex = 1, col = "blue", font = 1.5))
    #print(g)
    #tmp_df <- data.frame(feat)
    #colnames(tmp_df) = 'level'
    tmp_df$gene = as.numeric(expr[g,])
    names(tmp_df) [ncol(tmp_df)] <- g
    tmp_df$tgene = as.numeric(expr[g,])
    #tmp_df$level1[tmp_df$level %in% neg] = 0
    #tmp_df$level1[tmp_df$level %in% pos] = 1
    #names(tmp_df) [ncol(tmp_df)] <- paste0(g, "-lvl")
    #a = roc(feat, x, levels = c("Positive", "Negative"))
    #genes[[i]] <- tmp_df
    
    #if (i==1) {
    #plt = ggplot() + geom_roc(data = tmp_df, aes(x=level1, y=gene, color = gene))
    #}
    #else {
    #plt = plt + geom_roc(data = tmp_df, aes(x=level1, y=gene, color = gene))
    #}
    
    #pl[[i]] <- ggplot(tmp_df, aes(d=level1, m=gene)) + geom_roc() 
    #+ labs(title = paste(g, "~", echofn(her2)), subtitle = paste("AUC: ", a$auc)) 
    #pl[[i]]<-  pl[[i]]    + style_roc(theme = theme_grey) +
    #theme(axis.text = element_text(colour = "blue")) +
    #ggtitle(paste(g, "~", echofn(her2), "(AUC: ", round(calc_auc(pl[[i]])$AUC, 4), ")"))  
    # + annotate("text", x = .75, y = .25, 
    #         label = "AUC") +
    #scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))
    
    plt <- ggplot(tmp_df, aes(d=lvl, m= tgene)) + geom_roc() 
    auc_list[j] <- round(calc_auc(plt)$AUC, 4)
    auc_glist[g] <- round(calc_auc(plt)$AUC, 4)
    if (auc_list[j] < 0.5) {
      auc_list[j] = 1- auc_list[j]
      #  tmp_df$gene = 
    }
    if (auc_glist[g] < 0.5) {
      auc_glist[g] = 1- auc_glist[g]
      #  tmp_df$gene = 
    }
    print(auc_list[j])
    j = j+1
    
    #plt<-  plt    + style_roc(theme = theme_grey) +
    #theme(axis.text = element_text(colour = "blue")) +
    #ggtitle(paste(g, "~", echofn(her2), "(AUC: ", round(calc_auc(pl[[i]])$AUC, 4), ")"))  
    
    #plt <- ggplot() + geom_roc(data = tmp_df, aes(x=level1, y=gene, color = gene)) + geom_roc() 
    #+ labs(title = paste(g, "~", echofn(her2))) 
    #subtitle = paste("AUC: ", a$auc)) + 
    #+ theme (plot.title = element_text(color="blue", size=10, face="bold.italic"))
    #plot.subtitle = element_text(color="blue", size=8, face="bold.italic"))
    #ggtitle(paste(top_genes[g], "~", echofn(her2), "\n(AUC: ", a$auc, ")")) 
    #i = i+1
  }
  #ADD AUCS TO Genes_Info
  auc_all = data.frame(auc_glist)
  Genes_Info$AUC = auc_all$auc_glist
  
  genes_diff = setdiff(gene_vec, rownames(expr))
  
  Genes_Info[genes_diff,] <- NA
  
  Genes_Info=Genes_Info[order(Genes_Info$`P-value`),]
  
  pdf(paste0("gene table_", biostr, "_in_", ctype, ".pdf"), height=100, width=8.5)
  grid.table(Genes_Info)
  dev.off()
   
  #longtest <- melt_roc(tmp_df, "lvl", c(temp_df$y))
  #longtest <- melt_roc(tmp_df, "lvl", c(temp_df$r))
  #longtest <- melt_roc(tmp_df, "lvl", c(rownames(p.val)))
  #for (g in top_genes) {
   # if (auc_glist[g] < 0.5) {
     # longtest <- within(longtest, D[D == 1 & name %in% g ] <- 2)
      #longtest <- within(longtest, D[D == 0 & name %in% g ] <- 1)
      #longtest <- within(longtest, D[D == 2 & name %in% g ] <- 0)
   # }
    #print(auc_glist[g])
 # }  
  #pl[[i]] = ggplot(longtest, aes(d = D, m = M, color = name)) + 
    #ggplot(longtest, aes(x = D, y = M, label = c, colour = name)) +
   # geom_roc(labels = FALSE) + style_roc(theme = theme_gray) +
    #theme(axis.text = element_text(colour = "blue"), 
      #    axis.title=element_text(size=9),
        #  legend.position="none",
        #  plot.title = element_text(size=10)) +
   #ggtitle(paste(biostr, "- ROC Curves")) + 
    #annotate("text", x = .75, y = .25, 
    #label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2))) +
    # scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))
    #scale_x_continuous(breaks = seq(0, 1, by = .2)) 
  #ggsave(paste("rocCurves_", biostr, ".pdf"))
  #i=i+1
  
  
  
  #auc_df = data.frame(auc_list)
  #print (auc_list)
  #colnames(auc_df) = c(temp_df$y)
  #colnames(auc_df) = c(temp_df$r)
  #auc_tbl <- melt(auc_df)
  #auc_tbl$y = c(temp_df$y)
  #auc_tbl$y = c(temp_df$r)
  #colnames(auc_df) = 'x'
  #auc_df$y = c(temp_df$r)
  #auc_df = auc_df[order(-auc_df$x),]
  #auc_df$z <- factor(auc_df$y, levels = auc_df$y[order(-auc_df$x)])
  #print(auc_df)
  #pl[[i]] = ggplot(auc_tbl, aes(variable, value, fill = y))+ geom_bar(stat="identity")+
  #pl[[i]] = ggplot(auc_df, aes(x=y, y=x, fill = y))+ geom_bar(stat="identity")+
  # pl[[i]] = ggplot(auc_df, aes(y, x, fill = y))+ geom_bar(stat="identity")+
  #pl[[i]] = ggplot(auc_df, aes(z, x, fill = z))+ geom_bar(stat="identity")+
    #ggtitle(paste(biostr, " AUC values")) + labs(x = "Genes", y = "AUC", fill = "Genes") +
    #theme(axis.text.x=element_blank(), #axis.title=element_text(size=9), 
          #axis.title.x=element_blank(),
         # axis.ticks.x=element_blank(),
          #axis.text=element_text(size=10),
          #legend.position="bottom", 
          #legend.position="none",
          #legend.key.size = unit(0.4, "cm"),
         # legend.title = element_text(size=8), 
         # legend.text = element_text(size=7),
         # plot.title = element_text(size=10))
  #ggsave(paste("auc_barplot_", biostr, ".pdf"))
  
  #i=i+1
  
  #temp3_df <- data.frame(head(sort(p,decreasing = TRUE), 3))
  #temp3_df <- data.frame(head(sort(p.adj), 3))
  #temp3_df <- data.frame(head(sort(p), 3))
  #colnames(temp3_df) = 'x'
  #temp3_df$y = rownames(temp3_df)
  
  #temp3_df <- data.frame(head(auc_df, 3))
  #print(temp3_df)
  #temp_df$z <- factor(temp_df$y, levels = temp_df$y[order(-temp_df$x)])
  #btest <- melt(tmp_df,id.vars = c("level","lvl"), measure.vars = c(temp3_df$y))
  #btest <- melt(tmp_df,id.vars = c("lvl0","lvl"), measure.vars = c(temp3_df$y))
  #pl[[i]] <- ggplot(btest, aes(x=factor(lvl0),y=value, fill = factor(lvl0))) +
   # geom_boxplot(outlier.size=0,show.legend = F,lwd=0.1) + labs(title="Box plots for Top 3 Genes") +facet_wrap(~variable) +
    #theme(axis.text.x=element_blank(),
          #axis.title=element_text(size=9),
         # axis.title.x=element_blank(),
         # axis.ticks.x=element_blank(),
         # legend.title = element_text(size=6), 
         # legend.text = element_text(size=4),
         # plot.title = element_text(size=10),
         # legend.position="bottom") +
   # scale_fill_manual(values=group.colors, na.value="yellow")
  #i=i+1
  
  #plt<-  plt    + style_roc(theme = theme_grey) +
  #  theme(axis.text = element_text(colour = "blue")) +
  #  ggtitle(paste(g, "~", echofn(her2), "(AUC: ", round(calc_auc(plt$AUC, 4), ")")))  
  #grid.arrange(grobs=pl, ncol=3)
  #ml <- marrangeGrob(pl, nrow=2, ncol=2, top= paste(biostr, " - Analysis"))
  #leg_pl <- pl[[2]]
  #legend <- get_legend(leg_pl)
  
  #lay <- rbind(c(2,1),
            #   c(3,4),
             #  c(5,5))
 # ml = grid.arrange(grobs = pl, layout_matrix = lay, top= paste(biostr, " in ", ctype)) 
  #legend = legend)
  
  #gs <- arrangeGrob(grobs = pl, layout_matrix = lay)
  #ml <- marrangeGrob(grobs=pl, nrow=2, ncol=2, top= paste(biostr, "in ", ctype))
 # ggsave(paste0("plots_", biostr, "_in_", ctype, ".pdf"), plot=ml, width = 6, height = 8, units = "in")
  #genes_df <- data.frame(genes)
  #plt = ggplot() + geom_roc(data = genes_df, aes(x=level1, y=gene, group = genes, color = genes))
  #print(plt)
  #do.call(grid.arrange, pl)
  #dev.off()
  
} 

#Filtering out immune cells function
#wont work for stf data because there gene= gene_id and because the gene name is not the start of string. 
#need to modify if want to apply to stf data
immunesep = function(ag_df){
  #set T_vec equal to vector of genes
  T_vec = "CD3G, CD3D, TRA, CD6, CD5, NPDC1, CD28, CAMK4, GFI1, GATA3, SH2D1A, TRB, TNFRSF25, NK4, TACTILE, BCL11B, CD3E, INPP4B, MAL, NPDC1, ITM2A, ITK, LCK, NFATC3, RORA, MGC19764, TCF7, ZAP70, LEF1, SPOCK2, PRKCQ, SATB1, RASGRP1, LRIG1, DPP4, CD3Z, PDE4D, FYN, WWP1, LAT, DUSP16, KIAA0748, CDR2, STAT4, FLT3LG, IL6ST"
  #use gsub to remove all white spaces in the string (sep can only be 1 char, which is ",")
  T_vec = gsub(" ", "", T_vec, fixed = TRUE)
  #print T_vec and paste output into TextCon.. reassign T_vec to the char vec (only sep is "," in text)
  T_vec = scan(what=character(0), sep=",", file=textConnection(T_vec))
  ig_rows <- (grep(paste0( "^", T_vec, "$",collapse = "|"), ag_df$gene))
  ig_df <-ag_df[ig_rows,]
  #ig_df =subset(ag_df,rownames(ag_df) %in% T_vec)
  ig_df$ic = "Tcell"
  
  TCD8_vec = "FCGBP, C1ORF21, PHEMX, KLRG1, ZNF145, ADRB2, DUSP2, IL2RB, CCL5, GZMC, TBX21, CCL4L, GZMH, PRF1, GNLY, CST7, GPR56, KLRC1, S100B, D12S2489E, CD8B1, CD8A"
  TCD8_vec = gsub(" ", "", TCD8_vec, fixed = TRUE)
  TCD8_vec = scan(what=character(0), sep=",", file=textConnection(TCD8_vec))
  ig_rows <- (grep(paste0( "^", TCD8_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  #ig_temp = subset(ag_df,rownames(ag_df) %in% TCD8_vec)
  ig_temp$ic = "T+CD8cells"
  ig_df <- rbind(ig_df, ig_temp)
  
  Treg_vec = "AKAP2, ANXA2, C8FW, CALM2, CDH13, ENTPD1, EPSTI1, F5, GBP2, GBP5, HS3ST3B1, HSPCA, ICA1, LAIR2, LGALS3, LOC51191, MKP-7, NINJ2, P53DINP1, PMAIP1, PTPLA, GPR2, HLA-DMA, HLA-DPA1, HLA-DPB1, HLA-DQB1, HLA-DRA, HLA-DRB4, HLA-DRB5, HPGD, PTT1, SAT, SHMT2, TMOD3, TNFRSF6, VIL2, LRRC32, CHAC1, EGR1, PYCR1, IL13, C3AR1, ZBED2, IL10, GAB3, NR4A1, IL1R2, ULBP1, IL17F, EIF4EBP1, RGS16, BCAT1, VEGFA, GIMAP1, GPR44, GPT2, PSPH, NLN, MS4A3, UCK2, TNIP3, TNFRSF9, IL7R, NR4A3, CSGALNA, CT1, SLC43A1, ZBTB32, SCFD2, HK2, C17ORF96, C1ORF186, LYZ, SFXN4, CYP1B1, SLC7A5, C17ORF58, SLC4A5, CD52, KLRB1, S100A4, POLA1, SLC7A1, ATPBD4, UBXN8, CASP8, FAIM3, IL12RB2, MPP7, ZNF831, BNIP3L, LRP8, UTRN, SOCS3, LY9, EVI2B, SFXN3, TNFRSF4, ING4, BLM, TNFSF11, SLFN5, CCR8, CD244, MYO1F, HOMER1, ATF3, AGPAT4, CTH, JAKMIP1, TFRC, RPL26L1, MICAL1, GPRASP1, TMX4, FYB, ZNF792, WDR4, HUWE1, CD83, TRIB2, LARP1B, SLC38A5, AUH, CTPS, CYTH1, APOL3, CHD6, IL1R1, LPCAT4, CITED2, DUSP4, SSBP2, LAYN, ZFAND2B, PIK3R1, GCNT4, PPP3CA, TIMP1, SIRPG, CCR3, PTP4A3, SERINC5, MAGOHB, GIMAP4, KLHL24, CCNG2, WDR70, MYBL1, CNPY4, TIMM8A, CCDC28A, MNAT1, MS4A2, GZMB, E2F5, HIBCH, C9ORF91, NHP2, CDC25B, VDR, CSF1, GPATCH4, NFKBID, TMCO7, EHBP1, FAM173B, FARS2, FLVCR2, NPM3, PNPO, P2RY8, TMEM63A, MXI1, MRTO4, SESN2, C1ORF38, STK38, SLC6A6, ALG14, LXN, TNFRSF18, ACSS1, UBA7, MTSS1, PSMA6, FAM8A1, RAB31, TNRC6B, SESN1, CEP68, LAPTM5, PLAGL2, LRRC33, EIF2B3, URB2, DDX21, ITGA4, RGS2, RORC, IFNAR2, CSDA, SORL1, LRRC37B, IL10RA, PRKACB, CBWD1, SMAD3, C1ORF96, FCER1A, HSP90AB1, EMG1"
  Treg_vec = gsub(" ", "", Treg_vec, fixed = TRUE)
  Treg_vec = scan(what=character(0), sep=",", file=textConnection(Treg_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% Treg_vec)
  ig_rows <- (grep(paste0( "^", Treg_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "Tregcells"
  ig_df <- rbind(ig_df, ig_temp)
  
  B_vec = "PSEN2, FZD5, PLEKHF2, SIDT2, MTPN, TUBB6, FLJ21127, SMC6L1, KIAA1026, NUP88, GNA12, RNTRE, GSTZ1, RFX5, XYLT1, GL012, EPB41L2, LOC51760, SCRN1, PTD008, CYB561D2, SCAP2, MMP11, ERDJ5, TTC7A, PMAIP1, GARNL4, TMEPAI, SH3BP5, EDD, UCP4, EVI5L, PIK3C2B, SMAD3, LHFPL2, LOC57228, GM2A, HRK, ZNF207, LOC92497, F5, DMXL1, CLIC4, PRCP, DKFZP434C0328, STRBP, PHF16, JUP, TEAD2, PPP3CA, CNR2, ATP5B, ABCA1, BMF, ZCCHC7, KIAA0274, SLC22A3, SPI1, AMFR, ANXA4, EBF, ITPR1, HIST1H2BK, COPS3, COL14A1, SAV1, APOBEC3B, CCNG2, STAG3, LOC348938, POLD4, TBC1D1, MK2S4, APG4A, DAPP1, TNFRSF18, HERPUD1, POU2F2, ADK, CKIP-1, ACTA2, PARP14, MTSS1, DTX1, DDR1, RRAS2, RNF141, SYPL, SSPN, CORO1C, NFKBIE, CHD7, SP140, WDR34, BTNL9, ATP6V0A1, UROS, HLA-DOA, C22ORF13, IFI27, JMJD2B, WEE1, ODC1, KIAA0746, WDR11, SYNGR2, FOXP1, BLR1, MYO1E, MYBL2, GYLTL1B, SETBP1, KLF1, INPPL1, LOC54103, ZNF154, MHC2TA, CD24, ARHGEF3, BIRC3, TRIM56, HLA-DPB1, UVRAG, CD38, PEA15, FLJ10697, TLR10, ARHGAP10, HLA-DRA, SHMT2, LGALS9, FLJ10853, CEBPB, PRICKLE1, LCP21, CR1, KLHL14, TLR7, C20ORF72, HLA-DQB2, CORO1A, PRKCE HECW2, FLJ25604, TRIM26, FBXO41, GCNT1, LRMP, UBE2J1, IGKV3D, PCCA, GLDC, CYSLTR1, BTLA, NET5, HLA-DPA1, BLNK, CD79B, TM4SF8, TFEB, CD11ORF24,ADRBK2, CDKN2A, RIPK2, STX7, NCF4, SRGAP2, SLC2A5, CTSZ, AIM2, TCF4, TPD52, MAP3K8, MOBKL2B, HLA-DOB, IFIT3, MGC24039, BRD4, RALGPS2, SEMA4B, BSG,IGKC, HLA-DRB6, FUBP1, UNC84B, IFNGR2, HSPA5, HLA-DRB5, MARCKS, KYNU, PACAP, RHOBTB2, FLJ12363, TMED8, FGD2, FBXO10, IL4R, CD1C, MGC50844,MRPL49, CTSH, LYN, WASPIP, C3ORF6, EGR1, IGKV15, SPAP1, CHERP, IGHG3, ADAM28, HLA-DMA, HLA-DMB, PALM2-AKAP2, CD86, LAF4, LOC283663, FREB,DKFZP586A0522, UREB1, IGLL1, HLA-DQA2, KIAA1219, PLCG2, MARCH-1, BCNP1, PNOC, CD20, PSCD1, BCL11A, VPS28, SWAP70, SYK, BCL7A, NAP1L, TEM6, BTK,RAB30, FCGR2B, STAT6, BRDG1, LOC201895, ITGB1, CD74, CORO2B, VPREB3, HSPC182, MGC15619, RHOH, IGHA, CYBASC3, KIAA0125, RAM2, CD83, FCRH1, FLII,SNX10, IGHG1, EPHX1, FLJ10979, GTPAP, IGLL3, E2F5, TNFRSF17, BLK, FLJ00332, SNX2, CD200, HLA-DQB1, TRIO, CYBB, PIK3AP1, IRF4, CD22, BACE2, IGJ, PAX5,RGS13, TCL1A, MEF2C, C13ORF18, POU2AF1, HHEX, IRF8, BANK1, OSBPL10, SLC2A1, NCF1, IGL@, IGLC2, SPAP1, IGHG3, HLA-DMA, HLA-DMB, CD79A, HSPA6, IGHD,KIAA0476, SLC7A7, NKG7, CD19, SAMD9, LY86, SPIB, NAPSB, RNASE6, IGHM, LY64, CD72, IRTA2, CD1D, HLA-DQA1, SAS, CTSH, LYN, WASPIP, C3ORF6, EGR1, FREB,DKFZP586A0522"
  B_vec = gsub(" ", "", B_vec, fixed = TRUE)
  B_vec = scan(what=character(0), sep=",", file=textConnection(B_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% B_vec)
  ig_rows <- (grep(paste0( "^", B_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "Bcells"
  ig_df <- rbind(ig_df, ig_temp)
  
  NK_vec = "IL18RAP, PRF1, GZMH, TKTL1"
  NK_vec = gsub(" ", "", NK_vec, fixed = TRUE)
  NK_vec = scan(what=character(0), sep=",", file=textConnection(NK_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% NK_vec)
  ig_rows <- (grep(paste0( "^", NK_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "NKcells"
  ig_df <- rbind(ig_df, ig_temp)
  
  ClasMon_vec = "S100A12, ALOX5AP, PAD14, NRG1, MCEMP1, THBS1, CRISPLD2, F13A1, MOSC1, CYP27A1, CD163, QPCT, ADAM19, ASGR2, RNASE4, ALDH1A1, EDG3, PROK2, CLEC4D, S100A9, OSM, CSPG2, GALS2, ANG, VNN2, EBI2, UTS2, RPPH1, MGST1, IL8, CCR2, SELL, DYSF, MT1F, CD14, RANSE6, CSF3R, STAB1, CYP1B1, ADHFE1, PLA2G7, MMP25, ZNF395, SLC2A3, EGR1, CMTM2, CD9"
  ClasMon_vec = gsub(" ", "", ClasMon_vec, fixed = TRUE)
  ClasMon_vec = scan(what=character(0), sep=",", file=textConnection(ClasMon_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% ClasMon_vec)
  ig_rows <- (grep(paste0( "^", ClasMon_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "ClasMon"
  ig_df <- rbind(ig_df, ig_temp)
  
  IntMon_vec = "GFRA2, NKG7, PLVAP, PLAC8, MARCKSL1, E2F2, G8P4, CLEC10A, SCD, COTL1, SLC29A1, DDIT4, TGM2, LILRA3, ATF5, GPA33, C1QC, EVA1, MERTK, MGLL, DDEFL1, MARCO, NR1H3, FBP1, ACP2, GBP1, GPBAR1, SASH1, OLFM1, TIMP1, HLA-DOA, CAMK1, POUFUT1, EPB4IL3, H19, ZNF703, SNX5, CLEC10, AK-ALPHA-1, DPB1, DHRS9, MTHFD2, RGL1, PRDM1, FADS1, SLC2A8, CSK, ISOC2, CD300C, FGD2"
  IntMon_vec = gsub(" ", "", IntMon_vec, fixed = TRUE)
  IntMon_vec = scan(what=character(0), sep=",", file=textConnection(IntMon_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% IntMon_vec)
  ig_rows <- (grep(paste0( "^", IntMon_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "IntMon"
  ig_df <- rbind(ig_df, ig_temp)
  
  NonclasMon_vec = "C1QA, HES4, CDKN1C, C1QB, RHOC, CLEC4F, ADA, TAGLN, TCF7L2, CD798, CTSL, SFTPD, ABI3, SETBP1, FGFRL1, CD798, SPRED1, SNFT, EVL, MTSS1, SH2D1B, CYFIP2, NELF, STS-1, SIGLEC10, CKS1B, PTP4A3, LYPD2, SIDT2, INSIG1, LTB, PAPSS2, ABCC3, CASP5, VOM1, HSPB1, GNGT2, HMOX1, RRAS, RNF122, CDH23, FMNL2, RGS12, SGGB3A1, FER1L3, CALML4, IFITM3, CKB, CEACAM3, IFITM1"
  NonclasMon_vec = gsub(" ", "", NonclasMon_vec, fixed = TRUE)
  NonclasMon_vec = scan(what=character(0), sep=",", file=textConnection(NonclasMon_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% NonclasMon_vec)
  ig_rows <- (grep(paste0( "^", NonclasMon_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "NonclasMon"
  ig_df <- rbind(ig_df, ig_temp)
  
  MatDC_vec = "ACCL5, CXCL10, CCR7, IL15, IF127, IF144L, IFIH1, IFIT1, MX1, ISG15, ISG20, IRF7, GBP4, DUSP5, NFKB1A, ATF3, TNFSF10, IL6, IL8, IL7R, CCL4, TNFA1P6, IFIT3, OASL, GBP1, HES4, CYP27B1, RIPK2, TNFRSF9, SOD2, CD38, CD44, CD80, CD83, CD86, INDO MT2A, TRAF1, GADD45B, MT1M, MTIP2, BIRC3, USP18, TUBB2A, CCL8, EBI3, IFITM1, MT1B, MT1E, MT1G, MT1H, GADD45A, CD200, LAMP3, RGS1, SAT1"
  MatDC_vec = gsub(" ", "", MatDC_vec, fixed = TRUE)
  MatDC_vec = scan(what=character(0), sep=",", file=textConnection(MatDC_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% MatDC_vec)
  ig_rows <- (grep(paste0( "^", MatDC_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "MatDC"
  ig_df <- rbind(ig_df, ig_temp)
  
  MonDC_vec = "CD1A, CD1B, CD1C, CD23A, MRC1, CD36, HLA-DQA, HLA-DP light chain, HLA-DG protein 41, RIL, RAP1GAP, CCND2, DUSP5, PPIC, STAC, PRKACB, TRIB2, SHB"
  MonDC_vec = gsub(" ", "", MonDC_vec, fixed = TRUE)
  MonDC_vec = scan(what=character(0), sep=",", file=textConnection(MonDC_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% MonDC_vec)
  ig_rows <- (grep(paste0( "^", MonDC_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "MonDC"
  ig_df <- rbind(ig_df, ig_temp)
  
  MatNeut_vec = "ALPL, IL8RB, FCGR3B, SEMA3C, HM74, SOD2, FCGR3A, IL-8, STHM, IL8RA, FCGR2A, CSF3R, NCF2, AOAH"
  MatNeut_vec = gsub(" ", "", MatNeut_vec, fixed = TRUE)
  MatNeut_vec = scan(what=character(0), sep=",", file=textConnection(MatNeut_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% MatNeut_vec)
  ig_rows <- (grep(paste0( "^", MatNeut_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "MatNeut"
  ig_df <- rbind(ig_df, ig_temp)
  
  ImNeut_vec = "AZU1, ELA2, BPI, LCN2, MPO, CTSG, MMP8, DEFA4, DEFA3, CAMP, X-CGD"
  ImNeut_vec = gsub(" ", "", ImNeut_vec, fixed = TRUE)
  ImNeut_vec = scan(what=character(0), sep=",", file=textConnection(ImNeut_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% ImNeut_vec)
  ig_rows <- (grep(paste0( "^", ImNeut_vec, "$",collapse = "|"), ag_df$gene))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "ImNeut"
  ig_df <- rbind(ig_df, ig_temp)
  
  return(ig_df)
}

#altered immunesep for stf data since it needs to be retrieved from gene_id(diff format)

immunesepstf = function(ag_df){
  #set T_vec equal to vector of genes
  T_vec = "CD3G, CD3D, TRA, CD6, CD5, NPDC1, CD28, CAMK4, GFI1, GATA3, SH2D1A, TRB, TNFRSF25, NK4, TACTILE, BCL11B, CD3E, INPP4B, MAL, NPDC1, ITM2A, ITK, LCK, NFATC3, RORA, MGC19764, TCF7, ZAP70, LEF1, SPOCK2, PRKCQ, SATB1, RASGRP1, LRIG1, DPP4, CD3Z, PDE4D, FYN, WWP1, LAT, DUSP16, KIAA0748, CDR2, STAT4, FLT3LG, IL6ST"
  #use gsub to remove all white spaces in the string (sep can only be 1 char, which is ",")
  T_vec = gsub(" ", "", T_vec, fixed = TRUE)
  #print T_vec and paste output into TextCon.. reassign T_vec to the char vec (only sep is "," in text)
  T_vec = scan(what=character(0), sep=",", file=textConnection(T_vec))
  ig_rows <- (grep(paste0( ".+_", T_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_df <-ag_df[ig_rows,]
  #ig_df =subset(ag_df,rownames(ag_df) %in% T_vec)
  ig_df$ic = "Tcell"
  
  TCD8_vec = "FCGBP, C1ORF21, PHEMX, KLRG1, ZNF145, ADRB2, DUSP2, IL2RB, CCL5, GZMC, TBX21, CCL4L, GZMH, PRF1, GNLY, CST7, GPR56, KLRC1, S100B, D12S2489E, CD8B1, CD8A"
  TCD8_vec = gsub(" ", "", TCD8_vec, fixed = TRUE)
  TCD8_vec = scan(what=character(0), sep=",", file=textConnection(TCD8_vec))
  ig_rows <- (grep(paste0( ".+_", TCD8_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  #ig_temp = subset(ag_df,rownames(ag_df) %in% TCD8_vec)
  ig_temp$ic = "T+CD8cells"
  ig_df <- rbind(ig_df, ig_temp)
  
  Treg_vec = "AKAP2, ANXA2, C8FW, CALM2, CDH13, ENTPD1, EPSTI1, F5, GBP2, GBP5, HS3ST3B1, HSPCA, ICA1, LAIR2, LGALS3, LOC51191, MKP-7, NINJ2, P53DINP1, PMAIP1, PTPLA, GPR2, HLA-DMA, HLA-DPA1, HLA-DPB1, HLA-DQB1, HLA-DRA, HLA-DRB4, HLA-DRB5, HPGD, PTT1, SAT, SHMT2, TMOD3, TNFRSF6, VIL2, LRRC32, CHAC1, EGR1, PYCR1, IL13, C3AR1, ZBED2, IL10, GAB3, NR4A1, IL1R2, ULBP1, IL17F, EIF4EBP1, RGS16, BCAT1, VEGFA, GIMAP1, GPR44, GPT2, PSPH, NLN, MS4A3, UCK2, TNIP3, TNFRSF9, IL7R, NR4A3, CSGALNA, CT1, SLC43A1, ZBTB32, SCFD2, HK2, C17ORF96, C1ORF186, LYZ, SFXN4, CYP1B1, SLC7A5, C17ORF58, SLC4A5, CD52, KLRB1, S100A4, POLA1, SLC7A1, ATPBD4, UBXN8, CASP8, FAIM3, IL12RB2, MPP7, ZNF831, BNIP3L, LRP8, UTRN, SOCS3, LY9, EVI2B, SFXN3, TNFRSF4, ING4, BLM, TNFSF11, SLFN5, CCR8, CD244, MYO1F, HOMER1, ATF3, AGPAT4, CTH, JAKMIP1, TFRC, RPL26L1, MICAL1, GPRASP1, TMX4, FYB, ZNF792, WDR4, HUWE1, CD83, TRIB2, LARP1B, SLC38A5, AUH, CTPS, CYTH1, APOL3, CHD6, IL1R1, LPCAT4, CITED2, DUSP4, SSBP2, LAYN, ZFAND2B, PIK3R1, GCNT4, PPP3CA, TIMP1, SIRPG, CCR3, PTP4A3, SERINC5, MAGOHB, GIMAP4, KLHL24, CCNG2, WDR70, MYBL1, CNPY4, TIMM8A, CCDC28A, MNAT1, MS4A2, GZMB, E2F5, HIBCH, C9ORF91, NHP2, CDC25B, VDR, CSF1, GPATCH4, NFKBID, TMCO7, EHBP1, FAM173B, FARS2, FLVCR2, NPM3, PNPO, P2RY8, TMEM63A, MXI1, MRTO4, SESN2, C1ORF38, STK38, SLC6A6, ALG14, LXN, TNFRSF18, ACSS1, UBA7, MTSS1, PSMA6, FAM8A1, RAB31, TNRC6B, SESN1, CEP68, LAPTM5, PLAGL2, LRRC33, EIF2B3, URB2, DDX21, ITGA4, RGS2, RORC, IFNAR2, CSDA, SORL1, LRRC37B, IL10RA, PRKACB, CBWD1, SMAD3, C1ORF96, FCER1A, HSP90AB1, EMG1"
  Treg_vec = gsub(" ", "", Treg_vec, fixed = TRUE)
  Treg_vec = scan(what=character(0), sep=",", file=textConnection(Treg_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% Treg_vec)
  ig_rows <- (grep(paste0( ".+_", Treg_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "Tregcells"
  ig_df <- rbind(ig_df, ig_temp)
  
  B_vec = "PSEN2, FZD5, PLEKHF2, SIDT2, MTPN, TUBB6, FLJ21127, SMC6L1, KIAA1026, NUP88, GNA12, RNTRE, GSTZ1, RFX5, XYLT1, GL012, EPB41L2, LOC51760, SCRN1, PTD008, CYB561D2, SCAP2, MMP11, ERDJ5, TTC7A, PMAIP1, GARNL4, TMEPAI, SH3BP5, EDD, UCP4, EVI5L, PIK3C2B, SMAD3, LHFPL2, LOC57228, GM2A, HRK, ZNF207, LOC92497, F5, DMXL1, CLIC4, PRCP, DKFZP434C0328, STRBP, PHF16, JUP, TEAD2, PPP3CA, CNR2, ATP5B, ABCA1, BMF, ZCCHC7, KIAA0274, SLC22A3, SPI1, AMFR, ANXA4, EBF, ITPR1, HIST1H2BK, COPS3, COL14A1, SAV1, APOBEC3B, CCNG2, STAG3, LOC348938, POLD4, TBC1D1, MK2S4, APG4A, DAPP1, TNFRSF18, HERPUD1, POU2F2, ADK, CKIP-1, ACTA2, PARP14, MTSS1, DTX1, DDR1, RRAS2, RNF141, SYPL, SSPN, CORO1C, NFKBIE, CHD7, SP140, WDR34, BTNL9, ATP6V0A1, UROS, HLA-DOA, C22ORF13, IFI27, JMJD2B, WEE1, ODC1, KIAA0746, WDR11, SYNGR2, FOXP1, BLR1, MYO1E, MYBL2, GYLTL1B, SETBP1, KLF1, INPPL1, LOC54103, ZNF154, MHC2TA, CD24, ARHGEF3, BIRC3, TRIM56, HLA-DPB1, UVRAG, CD38, PEA15, FLJ10697, TLR10, ARHGAP10, HLA-DRA, SHMT2, LGALS9, FLJ10853, CEBPB, PRICKLE1, LCP21, CR1, KLHL14, TLR7, C20ORF72, HLA-DQB2, CORO1A, PRKCE HECW2, FLJ25604, TRIM26, FBXO41, GCNT1, LRMP, UBE2J1, IGKV3D, PCCA, GLDC, CYSLTR1, BTLA, NET5, HLA-DPA1, BLNK, CD79B, TM4SF8, TFEB, CD11ORF24,ADRBK2, CDKN2A, RIPK2, STX7, NCF4, SRGAP2, SLC2A5, CTSZ, AIM2, TCF4, TPD52, MAP3K8, MOBKL2B, HLA-DOB, IFIT3, MGC24039, BRD4, RALGPS2, SEMA4B, BSG,IGKC, HLA-DRB6, FUBP1, UNC84B, IFNGR2, HSPA5, HLA-DRB5, MARCKS, KYNU, PACAP, RHOBTB2, FLJ12363, TMED8, FGD2, FBXO10, IL4R, CD1C, MGC50844,MRPL49, CTSH, LYN, WASPIP, C3ORF6, EGR1, IGKV15, SPAP1, CHERP, IGHG3, ADAM28, HLA-DMA, HLA-DMB, PALM2-AKAP2, CD86, LAF4, LOC283663, FREB,DKFZP586A0522, UREB1, IGLL1, HLA-DQA2, KIAA1219, PLCG2, MARCH-1, BCNP1, PNOC, CD20, PSCD1, BCL11A, VPS28, SWAP70, SYK, BCL7A, NAP1L, TEM6, BTK,RAB30, FCGR2B, STAT6, BRDG1, LOC201895, ITGB1, CD74, CORO2B, VPREB3, HSPC182, MGC15619, RHOH, IGHA, CYBASC3, KIAA0125, RAM2, CD83, FCRH1, FLII,SNX10, IGHG1, EPHX1, FLJ10979, GTPAP, IGLL3, E2F5, TNFRSF17, BLK, FLJ00332, SNX2, CD200, HLA-DQB1, TRIO, CYBB, PIK3AP1, IRF4, CD22, BACE2, IGJ, PAX5,RGS13, TCL1A, MEF2C, C13ORF18, POU2AF1, HHEX, IRF8, BANK1, OSBPL10, SLC2A1, NCF1, IGL@, IGLC2, SPAP1, IGHG3, HLA-DMA, HLA-DMB, CD79A, HSPA6, IGHD,KIAA0476, SLC7A7, NKG7, CD19, SAMD9, LY86, SPIB, NAPSB, RNASE6, IGHM, LY64, CD72, IRTA2, CD1D, HLA-DQA1, SAS, CTSH, LYN, WASPIP, C3ORF6, EGR1, FREB,DKFZP586A0522"
  B_vec = gsub(" ", "", B_vec, fixed = TRUE)
  B_vec = scan(what=character(0), sep=",", file=textConnection(B_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% B_vec)
  ig_rows <- (grep(paste0( ".+_", B_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "Bcells"
  ig_df <- rbind(ig_df, ig_temp)
  
  NK_vec = "IL18RAP, PRF1, GZMH, TKTL1"
  NK_vec = gsub(" ", "", NK_vec, fixed = TRUE)
  NK_vec = scan(what=character(0), sep=",", file=textConnection(NK_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% NK_vec)
  ig_rows <- (grep(paste0( ".+_", NK_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "NKcells"
  ig_df <- rbind(ig_df, ig_temp)
  
  ClasMon_vec = "S100A12, ALOX5AP, PAD14, NRG1, MCEMP1, THBS1, CRISPLD2, F13A1, MOSC1, CYP27A1, CD163, QPCT, ADAM19, ASGR2, RNASE4, ALDH1A1, EDG3, PROK2, CLEC4D, S100A9, OSM, CSPG2, GALS2, ANG, VNN2, EBI2, UTS2, RPPH1, MGST1, IL8, CCR2, SELL, DYSF, MT1F, CD14, RANSE6, CSF3R, STAB1, CYP1B1, ADHFE1, PLA2G7, MMP25, ZNF395, SLC2A3, EGR1, CMTM2, CD9"
  ClasMon_vec = gsub(" ", "", ClasMon_vec, fixed = TRUE)
  ClasMon_vec = scan(what=character(0), sep=",", file=textConnection(ClasMon_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% ClasMon_vec)
  ig_rows <- (grep(paste0( ".+_", ClasMon_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "ClasMon"
  ig_df <- rbind(ig_df, ig_temp)
  
  IntMon_vec = "GFRA2, NKG7, PLVAP, PLAC8, MARCKSL1, E2F2, G8P4, CLEC10A, SCD, COTL1, SLC29A1, DDIT4, TGM2, LILRA3, ATF5, GPA33, C1QC, EVA1, MERTK, MGLL, DDEFL1, MARCO, NR1H3, FBP1, ACP2, GBP1, GPBAR1, SASH1, OLFM1, TIMP1, HLA-DOA, CAMK1, POUFUT1, EPB4IL3, H19, ZNF703, SNX5, CLEC10, AK-ALPHA-1, DPB1, DHRS9, MTHFD2, RGL1, PRDM1, FADS1, SLC2A8, CSK, ISOC2, CD300C, FGD2"
  IntMon_vec = gsub(" ", "", IntMon_vec, fixed = TRUE)
  IntMon_vec = scan(what=character(0), sep=",", file=textConnection(IntMon_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% IntMon_vec)
  ig_rows <- (grep(paste0( ".+_", IntMon_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "IntMon"
  ig_df <- rbind(ig_df, ig_temp)
  
  NonclasMon_vec = "C1QA, HES4, CDKN1C, C1QB, RHOC, CLEC4F, ADA, TAGLN, TCF7L2, CD798, CTSL, SFTPD, ABI3, SETBP1, FGFRL1, CD798, SPRED1, SNFT, EVL, MTSS1, SH2D1B, CYFIP2, NELF, STS-1, SIGLEC10, CKS1B, PTP4A3, LYPD2, SIDT2, INSIG1, LTB, PAPSS2, ABCC3, CASP5, VOM1, HSPB1, GNGT2, HMOX1, RRAS, RNF122, CDH23, FMNL2, RGS12, SGGB3A1, FER1L3, CALML4, IFITM3, CKB, CEACAM3, IFITM1"
  NonclasMon_vec = gsub(" ", "", NonclasMon_vec, fixed = TRUE)
  NonclasMon_vec = scan(what=character(0), sep=",", file=textConnection(NonclasMon_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% NonclasMon_vec)
  ig_rows <- (grep(paste0( ".+_", NonclasMon_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "NonclasMon"
  ig_df <- rbind(ig_df, ig_temp)
  
  MatDC_vec = "ACCL5, CXCL10, CCR7, IL15, IF127, IF144L, IFIH1, IFIT1, MX1, ISG15, ISG20, IRF7, GBP4, DUSP5, NFKB1A, ATF3, TNFSF10, IL6, IL8, IL7R, CCL4, TNFA1P6, IFIT3, OASL, GBP1, HES4, CYP27B1, RIPK2, TNFRSF9, SOD2, CD38, CD44, CD80, CD83, CD86, INDO MT2A, TRAF1, GADD45B, MT1M, MTIP2, BIRC3, USP18, TUBB2A, CCL8, EBI3, IFITM1, MT1B, MT1E, MT1G, MT1H, GADD45A, CD200, LAMP3, RGS1, SAT1"
  MatDC_vec = gsub(" ", "", MatDC_vec, fixed = TRUE)
  MatDC_vec = scan(what=character(0), sep=",", file=textConnection(MatDC_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% MatDC_vec)
  ig_rows <- (grep(paste0( ".+_", MatDC_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "MatDC"
  ig_df <- rbind(ig_df, ig_temp)
  
  MonDC_vec = "CD1A, CD1B, CD1C, CD23A, MRC1, CD36, HLA-DQA, HLA-DP light chain, HLA-DG protein 41, RIL, RAP1GAP, CCND2, DUSP5, PPIC, STAC, PRKACB, TRIB2, SHB"
  MonDC_vec = gsub(" ", "", MonDC_vec, fixed = TRUE)
  MonDC_vec = scan(what=character(0), sep=",", file=textConnection(MonDC_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% MonDC_vec)
  ig_rows <- (grep(paste0( ".+_", MonDC_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "MonDC"
  ig_df <- rbind(ig_df, ig_temp)
  
  MatNeut_vec = "ALPL, IL8RB, FCGR3B, SEMA3C, HM74, SOD2, FCGR3A, IL-8, STHM, IL8RA, FCGR2A, CSF3R, NCF2, AOAH"
  MatNeut_vec = gsub(" ", "", MatNeut_vec, fixed = TRUE)
  MatNeut_vec = scan(what=character(0), sep=",", file=textConnection(MatNeut_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% MatNeut_vec)
  ig_rows <- (grep(paste0( ".+_", MatNeut_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "MatNeut"
  ig_df <- rbind(ig_df, ig_temp)
  
  ImNeut_vec = "AZU1, ELA2, BPI, LCN2, MPO, CTSG, MMP8, DEFA4, DEFA3, CAMP, X-CGD"
  ImNeut_vec = gsub(" ", "", ImNeut_vec, fixed = TRUE)
  ImNeut_vec = scan(what=character(0), sep=",", file=textConnection(ImNeut_vec))
  #ig_temp = subset(ag_df,rownames(ag_df) %in% ImNeut_vec)
  ig_rows <- (grep(paste0( ".+_", ImNeut_vec, "$",collapse = "|"), ag_df$gene_id))
  ig_temp <-ag_df[ig_rows,]
  ig_temp$ic = "ImNeut"
  ig_df <- rbind(ig_df, ig_temp)
  
  return(ig_df)
}
#TABLES FOR ENRICHED GENES using THCA_N.ge- deleted M1 samples
#Forming THCA_tcga_ic in between (under analyzefeature3 calls)

THCA_Clin_t.ge = data.frame(t(THCA_Clin.ge))
THCA_N_t.ge = data.frame(t(THCA_N.ge))

#T_rich
#set T_vec equal to vector of genes
#T_vec = "CD3G, CD3D, TRA, CD6, CD5, NPDC1, CD28, CAMK4, GFI1, GATA3, SH2D1A, TRB, TNFRSF25, NK4, TACTILE, BCL11B, CD3E, INPP4B, MAL, NPDC1, ITM2A, ITK, LCK, NFATC3, RORA, MGC19764, TCF7, ZAP70, LEF1, SPOCK2, PRKCQ, SATB1, RASGRP1, LRIG1, DPP4, CD3Z, PDE4D, FYN, WWP1, LAT, DUSP16, KIAA0748, CDR2, STAT4, FLT3LG, IL6ST"
#use gsub to remove all white spaces in the string (sep can only be 1 char, which is ",")
#T_vec = gsub(" ", "", T_vec, fixed = TRUE)
#print T_vec and paste output into TextCon.. reassign T_vec to the char vec (only sep is "," in text)
#T_vec = scan(what=character(0), sep=",", file=textConnection(T_vec)) 
#Make data frame of .ge of only genes in T_vec
#T_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% T_vec)

#analyzeFeature2(t(T_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-Tcells-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(T_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-Tcells-N1b), echofn(THCA-Pathological), T_vec)

#THCA_tcga_ic <- data.frame(t(T_rich.ge))
#THCA_tcga_ic$ic = "Tcell"

#TCD8_rich.ge
#TCD8_vec = "FCGBP, C1ORF21, PHEMX, KLRG1, ZNF145, ADRB2, DUSP2, IL2RB, CCL5, GZMC, TBX21, CCL4L, GZMH, PRF1, GNLY, CST7, GPR56, KLRC1, S100B, D12S2489E, CD8B1, CD8A"
#TCD8_vec = gsub(" ", "", TCD8_vec, fixed = TRUE)
#TCD8_vec = scan(what=character(0), sep=",", file=textConnection(TCD8_vec))
#TCD8_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% TCD8_vec)

#analyzeFeature2(t(TCD8_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-T&CD8cells-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(TCD8_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-T&CD8cells-N1b), echofn(THCA-Pathological), TCD8_vec)

#ic_temp <- data.frame(t(TCD8_rich.ge))
#ic_temp$ic = "T+CD8cells"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#Treg_rich.ge
#Treg_vec = "AKAP2, ANXA2, C8FW, CALM2, CDH13, ENTPD1, EPSTI1, F5, GBP2, GBP5, HS3ST3B1, HSPCA, ICA1, LAIR2, LGALS3, LOC51191, MKP-7, NINJ2, P53DINP1, PMAIP1, PTPLA, GPR2, HLA-DMA, HLA-DPA1, HLA-DPB1, HLA-DQB1, HLA-DRA, HLA-DRB4, HLA-DRB5, HPGD, PTT1, SAT, SHMT2, TMOD3, TNFRSF6, VIL2, LRRC32, CHAC1, EGR1, PYCR1, IL13, C3AR1, ZBED2, IL10, GAB3, NR4A1, IL1R2, ULBP1, IL17F, EIF4EBP1, RGS16, BCAT1, VEGFA, GIMAP1, GPR44, GPT2, PSPH, NLN, MS4A3, UCK2, TNIP3, TNFRSF9, IL7R, NR4A3, CSGALNA, CT1, SLC43A1, ZBTB32, SCFD2, HK2, C17ORF96, C1ORF186, LYZ, SFXN4, CYP1B1, SLC7A5, C17ORF58, SLC4A5, CD52, KLRB1, S100A4, POLA1, SLC7A1, ATPBD4, UBXN8, CASP8, FAIM3, IL12RB2, MPP7, ZNF831, BNIP3L, LRP8, UTRN, SOCS3, LY9, EVI2B, SFXN3, TNFRSF4, ING4, BLM, TNFSF11, SLFN5, CCR8, CD244, MYO1F, HOMER1, ATF3, AGPAT4, CTH, JAKMIP1, TFRC, RPL26L1, MICAL1, GPRASP1, TMX4, FYB, ZNF792, WDR4, HUWE1, CD83, TRIB2, LARP1B, SLC38A5, AUH, CTPS, CYTH1, APOL3, CHD6, IL1R1, LPCAT4, CITED2, DUSP4, SSBP2, LAYN, ZFAND2B, PIK3R1, GCNT4, PPP3CA, TIMP1, SIRPG, CCR3, PTP4A3, SERINC5, MAGOHB, GIMAP4, KLHL24, CCNG2, WDR70, MYBL1, CNPY4, TIMM8A, CCDC28A, MNAT1, MS4A2, GZMB, E2F5, HIBCH, C9ORF91, NHP2, CDC25B, VDR, CSF1, GPATCH4, NFKBID, TMCO7, EHBP1, FAM173B, FARS2, FLVCR2, NPM3, PNPO, P2RY8, TMEM63A, MXI1, MRTO4, SESN2, C1ORF38, STK38, SLC6A6, ALG14, LXN, TNFRSF18, ACSS1, UBA7, MTSS1, PSMA6, FAM8A1, RAB31, TNRC6B, SESN1, CEP68, LAPTM5, PLAGL2, LRRC33, EIF2B3, URB2, DDX21, ITGA4, RGS2, RORC, IFNAR2, CSDA, SORL1, LRRC37B, IL10RA, PRKACB, CBWD1, SMAD3, C1ORF96, FCER1A, HSP90AB1, EMG1"
#Treg_vec = gsub(" ", "", Treg_vec, fixed = TRUE)
#Treg_vec = scan(what=character(0), sep=",", file=textConnection(Treg_vec))
#Treg_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% Treg_vec)

#analyzeFeature2(t(Treg_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-Tregulatory-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(Treg_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-Tregulatory-N1b), echofn(THCA-Pathological), Treg_vec)

#ic_temp <- data.frame(t(Treg_rich.ge))
#ic_temp$ic = "Tregcells"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#B_rich.ge
#B_vec = "PSEN2, FZD5, PLEKHF2, SIDT2, MTPN, TUBB6, FLJ21127, SMC6L1, KIAA1026, NUP88, GNA12, RNTRE, GSTZ1, RFX5, XYLT1, GL012, EPB41L2, LOC51760, SCRN1, PTD008, CYB561D2, SCAP2, MMP11, ERDJ5, TTC7A, PMAIP1, GARNL4, TMEPAI, SH3BP5, EDD, UCP4, EVI5L, PIK3C2B, SMAD3, LHFPL2, LOC57228, GM2A, HRK, ZNF207, LOC92497, F5, DMXL1, CLIC4, PRCP, DKFZP434C0328, STRBP, PHF16, JUP, TEAD2, PPP3CA, CNR2, ATP5B, ABCA1, BMF, ZCCHC7, KIAA0274, SLC22A3, SPI1, AMFR, ANXA4, EBF, ITPR1, HIST1H2BK, COPS3, COL14A1, SAV1, APOBEC3B, CCNG2, STAG3, LOC348938, POLD4, TBC1D1, MK2S4, APG4A, DAPP1, TNFRSF18, HERPUD1, POU2F2, ADK, CKIP-1, ACTA2, PARP14, MTSS1, DTX1, DDR1, RRAS2, RNF141, SYPL, SSPN, CORO1C, NFKBIE, CHD7, SP140, WDR34, BTNL9, ATP6V0A1, UROS, HLA-DOA, C22ORF13, IFI27, JMJD2B, WEE1, ODC1, KIAA0746, WDR11, SYNGR2, FOXP1, BLR1, MYO1E, MYBL2, GYLTL1B, SETBP1, KLF1, INPPL1, LOC54103, ZNF154, MHC2TA, CD24, ARHGEF3, BIRC3, TRIM56, HLA-DPB1, UVRAG, CD38, PEA15, FLJ10697, TLR10, ARHGAP10, HLA-DRA, SHMT2, LGALS9, FLJ10853, CEBPB, PRICKLE1, LCP21, CR1, KLHL14, TLR7, C20ORF72, HLA-DQB2, CORO1A, PRKCE HECW2, FLJ25604, TRIM26, FBXO41, GCNT1, LRMP, UBE2J1, IGKV3D, PCCA, GLDC, CYSLTR1, BTLA, NET5, HLA-DPA1, BLNK, CD79B, TM4SF8, TFEB, CD11ORF24,ADRBK2, CDKN2A, RIPK2, STX7, NCF4, SRGAP2, SLC2A5, CTSZ, AIM2, TCF4, TPD52, MAP3K8, MOBKL2B, HLA-DOB, IFIT3, MGC24039, BRD4, RALGPS2, SEMA4B, BSG,IGKC, HLA-DRB6, FUBP1, UNC84B, IFNGR2, HSPA5, HLA-DRB5, MARCKS, KYNU, PACAP, RHOBTB2, FLJ12363, TMED8, FGD2, FBXO10, IL4R, CD1C, MGC50844,MRPL49, CTSH, LYN, WASPIP, C3ORF6, EGR1, IGKV15, SPAP1, CHERP, IGHG3, ADAM28, HLA-DMA, HLA-DMB, PALM2-AKAP2, CD86, LAF4, LOC283663, FREB,DKFZP586A0522, UREB1, IGLL1, HLA-DQA2, KIAA1219, PLCG2, MARCH-1, BCNP1, PNOC, CD20, PSCD1, BCL11A, VPS28, SWAP70, SYK, BCL7A, NAP1L, TEM6, BTK,RAB30, FCGR2B, STAT6, BRDG1, LOC201895, ITGB1, CD74, CORO2B, VPREB3, HSPC182, MGC15619, RHOH, IGHA, CYBASC3, KIAA0125, RAM2, CD83, FCRH1, FLII,SNX10, IGHG1, EPHX1, FLJ10979, GTPAP, IGLL3, E2F5, TNFRSF17, BLK, FLJ00332, SNX2, CD200, HLA-DQB1, TRIO, CYBB, PIK3AP1, IRF4, CD22, BACE2, IGJ, PAX5,RGS13, TCL1A, MEF2C, C13ORF18, POU2AF1, HHEX, IRF8, BANK1, OSBPL10, SLC2A1, NCF1, IGL@, IGLC2, SPAP1, IGHG3, HLA-DMA, HLA-DMB, CD79A, HSPA6, IGHD,KIAA0476, SLC7A7, NKG7, CD19, SAMD9, LY86, SPIB, NAPSB, RNASE6, IGHM, LY64, CD72, IRTA2, CD1D, HLA-DQA1, SAS, CTSH, LYN, WASPIP, C3ORF6, EGR1, FREB,DKFZP586A0522"
#B_vec = gsub(" ", "", B_vec, fixed = TRUE)
#B_vec = scan(what=character(0), sep=",", file=textConnection(B_vec))
#B_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% B_vec)

#analyzeFeature2(t(B_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-Bcells-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(B_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-Bcells-N1b), echofn(THCA-Pathological), B_vec)

#ic_temp <- data.frame(t(B_rich.ge))
#ic_temp$ic = "Bcells"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#NK_rich.ge
#NK_vec = "IL18RAP, PRF1, GZMH, TKTL1"
#NK_vec = gsub(" ", "", NK_vec, fixed = TRUE)
#NK_vec = scan(what=character(0), sep=",", file=textConnection(NK_vec))
#NK_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% NK_vec)

#analyzeFeature2(t(NK_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-NKcells-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(NK_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-NKcells-N1b), echofn(THCA-Pathological), NK_vec)

#ic_temp <- data.frame(t(NK_rich.ge))
#ic_temp$ic = "NKcells"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#ClasMon_rich.ge
#ClasMon_vec = "S100A12, ALOX5AP, PAD14, NRG1, MCEMP1, THBS1, CRISPLD2, F13A1, MOSC1, CYP27A1, CD163, QPCT, ADAM19, ASGR2, RNASE4, ALDH1A1, EDG3, PROK2, CLEC4D, S100A9, OSM, CSPG2, GALS2, ANG, VNN2, EBI2, UTS2, RPPH1, MGST1, IL8, CCR2, SELL, DYSF, MT1F, CD14, RANSE6, CSF3R, STAB1, CYP1B1, ADHFE1, PLA2G7, MMP25, ZNF395, SLC2A3, EGR1, CMTM2, CD9"
#ClasMon_vec = gsub(" ", "", ClasMon_vec, fixed = TRUE)
#ClasMon_vec = scan(what=character(0), sep=",", file=textConnection(ClasMon_vec))
#ClasMon_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% ClasMon_vec)

#analyzeFeature2(t(ClasMon_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-ClassicalMonocytes-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(ClasMon_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-ClassicalMonocytes-N1b), echofn(THCA-Pathological), ClasMon_vec)

#ic_temp <- data.frame(t(ClasMon_rich.ge))
#ic_temp$ic = "ClasMon"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#IntMon_rich.ge
#IntMon_vec = "GFRA2, NKG7, PLVAP, PLAC8, MARCKSL1, E2F2, G8P4, CLEC10A, SCD, COTL1, SLC29A1, DDIT4, TGM2, LILRA3, ATF5, GPA33, C1QC, EVA1, MERTK, MGLL, DDEFL1, MARCO, NR1H3, FBP1, ACP2, GBP1, GPBAR1, SASH1, OLFM1, TIMP1, HLA-DOA, CAMK1, POUFUT1, EPB4IL3, H19, ZNF703, SNX5, CLEC10, AK-ALPHA-1, DPB1, DHRS9, MTHFD2, RGL1, PRDM1, FADS1, SLC2A8, CSK, ISOC2, CD300C, FGD2"
#IntMon_vec = gsub(" ", "", IntMon_vec, fixed = TRUE)
#IntMon_vec = scan(what=character(0), sep=",", file=textConnection(IntMon_vec))
#IntMon_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% IntMon_vec)

#analyzeFeature2(t(IntMon_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-IntermediateMonocytes-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(IntMon_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-IntermediateMonocytes-N1b), echofn(THCA-Pathological), IntMon_vec)

#ic_temp <- data.frame(t(IntMon_rich.ge))
#ic_temp$ic = "IntMon"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#NonclasMon_rich.ge
#NonclasMon_vec = "C1QA, HES4, CDKN1C, C1QB, RHOC, CLEC4F, ADA, TAGLN, TCF7L2, CD798, CTSL, SFTPD, ABI3, SETBP1, FGFRL1, CD798, SPRED1, SNFT, EVL, MTSS1, SH2D1B, CYFIP2, NELF, STS-1, SIGLEC10, CKS1B, PTP4A3, LYPD2, SIDT2, INSIG1, LTB, PAPSS2, ABCC3, CASP5, VOM1, HSPB1, GNGT2, HMOX1, RRAS, RNF122, CDH23, FMNL2, RGS12, SGGB3A1, FER1L3, CALML4, IFITM3, CKB, CEACAM3, IFITM1"
#NonclasMon_vec = gsub(" ", "", NonclasMon_vec, fixed = TRUE)
#NonclasMon_vec = scan(what=character(0), sep=",", file=textConnection(NonclasMon_vec))
#NonclasMon_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% NonclasMon_vec)

#analyzeFeature2(t(NonclasMon_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-NonclassicalMonocytes-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(NonclasMon_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-NonclassicalMonocytes-N1b), echofn(THCA-Pathological), NonclasMon_vec)

#ic_temp <- data.frame(t(NonclasMon_rich.ge))
#ic_temp$ic = "NonclasMon"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#MatDC_rich.ge
#MatDC_vec = "ACCL5, CXCL10, CCR7, IL15, IF127, IF144L, IFIH1, IFIT1, MX1, ISG15, ISG20, IRF7, GBP4, DUSP5, NFKB1A, ATF3, TNFSF10, IL6, IL8, IL7R, CCL4, TNFA1P6, IFIT3, OASL, GBP1, HES4, CYP27B1, RIPK2, TNFRSF9, SOD2, CD38, CD44, CD80, CD83, CD86, INDO MT2A, TRAF1, GADD45B, MT1M, MTIP2, BIRC3, USP18, TUBB2A, CCL8, EBI3, IFITM1, MT1B, MT1E, MT1G, MT1H, GADD45A, CD200, LAMP3, RGS1, SAT1"
#MatDC_vec = gsub(" ", "", MatDC_vec, fixed = TRUE)
#MatDC_vec = scan(what=character(0), sep=",", file=textConnection(MatDC_vec))
#MatDC_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% MatDC_vec)

#analyzeFeature2(t(MatDC_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-MatureDC-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(MatDC_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-MatureDC-N1b), echofn(THCA-Pathological), MatDC_vec)

#ic_temp <- data.frame(t(MatDC_rich.ge))
#ic_temp$ic = "MatDC"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#MonDC_rich.ge
#MonDC_vec = "CD1A, CD1B, CD1C, CD23A, MRC1, CD36, HLA-DQA, HLA-DP light chain, HLA-DG protein 41, RIL, RAP1GAP, CCND2, DUSP5, PPIC, STAC, PRKACB, TRIB2, SHB"
#MonDC_vec = gsub(" ", "", MonDC_vec, fixed = TRUE)
#MonDC_vec = scan(what=character(0), sep=",", file=textConnection(MonDC_vec))
#MonDC_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% MonDC_vec)

#analyzeFeature2(t(MonDC_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-MonocyteDC-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(MonDC_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-MonocyteDC-N1b), echofn(THCA-Pathological), MonDC_vec)

#ic_temp <- data.frame(t(MonDC_rich.ge))
#ic_temp$ic = "MonDC"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#MatNeut_rich.ge
#MatNeut_vec = "ALPL, IL8RB, FCGR3B, SEMA3C, HM74, SOD2, FCGR3A, IL-8, STHM, IL8RA, FCGR2A, CSF3R, NCF2, AOAH"
#MatNeut_vec = gsub(" ", "", MatNeut_vec, fixed = TRUE)
#MatNeut_vec = scan(what=character(0), sep=",", file=textConnection(MatNeut_vec))
#MatNeut_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% MatNeut_vec)

#analyzeFeature2(t(MatNeut_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-MatureNeutrophils-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(MatNeut_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-MatureNeutrophils-N1b), echofn(THCA-Pathological), MatNeut_vec)

#ic_temp <- data.frame(t(MatNeut_rich.ge))
#ic_temp$ic = "MatNeut"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#ImNeut_rich.ge
#ImNeut_vec = "AZU1, ELA2, BPI, LCN2, MPO, CTSG, MMP8, DEFA4, DEFA3, CAMP, X-CGD"
#ImNeut_vec = gsub(" ", "", ImNeut_vec, fixed = TRUE)
#ImNeut_vec = scan(what=character(0), sep=",", file=textConnection(ImNeut_vec))
#ImNeut_rich.ge = subset(THCA_N_t.ge,select = names(THCA_N_t.ge) %in% ImNeut_vec)

#analyzeFeature2(t(ImNeut_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-ImmatureNeutrophils-N1b), echofn(THCA-Pathological))
#analyzeFeature3(t(ImNeut_rich.ge),npath_df$nclas, c("1"), c("0"), echofn(THCA-ImmatureNeutrophils-N1b), echofn(THCA-Pathological), ImNeut_vec)

#ic_temp <- data.frame(t(ImNeut_rich.ge))
#ic_temp$ic = "ImNeut"
#THCA_tcga_ic <- rbind(THCA_tcga_ic, ic_temp)

#Use immunesep to get THCA_tcga_ic
THCA_N_ic.ge<- cbind(gene = rownames(THCA_N.ge), THCA_N.ge)
THCA_tcga_ic <- immunesep(THCA_N_ic.ge)

#Concatenated into THCA_tcga_ic export into excel
#THCA_tcga_ic <- cbind(gene = rownames(THCA_tcga_ic), THCA_tcga_ic)

write.xlsx(THCA_tcga_ic, file = "/Users/shreyashekhar/GRIPS/THCA_tcga_ic.xlsx", sheetName="THCA_tcga_ic")


# #ic, gene
# tcga_ictemp = THCA_tcga_ic[,!colnames(THCA_tcga_ic) %in% c("ic", "gene")]
# THCA_tcga <- data.frame(rowMeans(tcga_ictemp))
# tcga_fkpm_stats <-lapply(rownames(tcga_ictemp), function(x) t.test(tcga_ictemp[x,]))
# tcga_fkpm_conf <-unlist(lapply(tcga_fkpm_stats,FUN=function(x) x$conf.int))
# THCA_tcga$cilow <-tcga_fkpm_conf[ c(T,F) ] 
# THCA_tcga$cihigh <-tcga_fkpm_conf[ c(F,T) ] 
# THCA_tcga$ic <- THCA_tcga_ic$ic
# THCA_tcga <- cbind(gene = THCA_tcga_ic$gene, THCA_tcga)
# write.xlsx(THCA_tcga, file = "/Users/shreyashekhar/GRIPS/THCA_tcga_fpkm.xlsx", sheetName="THCA_tcga_fpkm")


#ADDING ISOFORM TPMS FOR TCGA
library("Metrics")
library("dplyr")    
THCA_iso_tcga.ge = data.frame(t(THCA_iso_tcga.ge))
#HAVE TO RUN EVERYTHING IN ORDER BECAUSE THERE IS A TEMP VAR & BECAUSE TRANSPOSED

#C1orf43
#find kgids for all C1orf43 isoforms, and put into character vector (^ means beginning of string & $ means end- 
#to match exact string(gene) and not strings(genes) where C1orf43 is a substring)
C1orf43_genes = as.character(THCA_transtogene$X.kgID[grepl("^C1orf43$", THCA_transtogene$geneSymbol)])
#make dataframe with only rows of the isoform 
iso_temp = subset(THCA_iso_tcga.ge, select = C1orf43_genes)
#CHANGES THE CLASS OF COLUMNS FROM FACTOR TO NUMERIC!!!! OTHERWISE CANNOT USE ROWSUMS
iso_temp[,C1orf43_genes]<-lapply(C1orf43_genes, function(x) as.numeric(as.character(iso_temp[,x])))
#finds the class of each column (check to make sure they are numeric)
lapply(iso_temp, class)
#puts row sums (isoform sums) for C1orf43 into data frame
#to create new, use this
#THCA_TPM_tcga <- data.frame(rowSums(iso_temp))
#to only redo this column, use this
THCA_TPM_tcga$C1orf43 = rowSums(iso_temp)
#name column with C1orf43 isoform sums "C1orf43"
colnames(THCA_TPM_tcga) = "C1orf43"

#CHMP2A
CHMP2A_genes = as.character(THCA_transtogene$X.kgID[grepl("^CHMP2A$", THCA_transtogene$geneSymbol)])
#GPI
GPI_genes = as.character(THCA_transtogene$X.kgID[grepl("^GPI$", THCA_transtogene$geneSymbol)])
iso_temp = subset(THCA_iso_tcga.ge, select = GPI_genes)
iso_temp[,GPI_genes]<-lapply(GPI_genes, function(x) as.numeric(as.character(iso_temp[,x])))
THCA_TPM_tcga$GPI = rowSums(iso_temp)

#PSMB2
PSMB2_genes = as.character(THCA_transtogene$X.kgID[grepl("^PSMB2$", THCA_transtogene$geneSymbol)])

#PSMB4
PSMB4_genes = as.character(THCA_transtogene$X.kgID[grepl("^PSMB4$", THCA_transtogene$geneSymbol)])
iso_temp = subset(THCA_iso_tcga.ge, select = PSMB4_genes)
iso_temp[,PSMB4_genes]<-lapply(PSMB4_genes, function(x) as.numeric(as.character(iso_temp[,x])))
THCA_TPM_tcga$PSMB4 = rowSums(iso_temp)


#RAB7A
RAB7A_genes = as.character(THCA_transtogene$X.kgID[grepl("^RAB7A$", THCA_transtogene$geneSymbol)])
iso_temp = subset(THCA_iso_tcga.ge, select = RAB7A_genes)
iso_temp[,RAB7A_genes]<-lapply(RAB7A_genes, function(x) as.numeric(as.character(iso_temp[,x])))
THCA_TPM_tcga$RAB7A = rowSums(iso_temp)

#REEP5
REEP5_genes = as.character(THCA_transtogene$X.kgID[grepl("^REEP5$", THCA_transtogene$geneSymbol)])
iso_temp = subset(THCA_iso_tcga.ge, select =REEP5_genes)
iso_temp[,REEP5_genes]<-lapply(REEP5_genes, function(x) as.numeric(as.character(iso_temp[,x])))
THCA_TPM_tcga$REEP5 = rowSums(iso_temp)

#SNRPD3
SNRPD3_genes = as.character(THCA_transtogene$X.kgID[grepl("^SNRPD3$", THCA_transtogene$geneSymbol)])
iso_temp = subset(THCA_iso_tcga.ge, select =SNRPD3_genes)
iso_temp[,SNRPD3_genes]<-lapply(SNRPD3_genes, function(x) as.numeric(as.character(iso_temp[,x])))
THCA_TPM_tcga$SNRPD3 = rowSums(iso_temp)

#VCP
VCP_genes = as.character(THCA_transtogene$X.kgID[grepl("^VCP$", THCA_transtogene$geneSymbol)])
iso_temp = subset(THCA_iso_tcga.ge, select =VCP_genes)
iso_temp[,VCP_genes]<-lapply(VCP_genes, function(x) as.numeric(as.character(iso_temp[,x])))
THCA_TPM_tcga$VCP = rowSums(iso_temp)

#VPS29
VPS29_genes = as.character(THCA_transtogene$X.kgID[grepl("^VPS29$", THCA_transtogene$geneSymbol)])

#GAPDH-not found in tcga
GAPDH_genes = as.character(THCA_transtogene$X.kgID[grepl("^GAPDH$", THCA_transtogene$geneSymbol)])

#GUSB
GUSB_genes = as.character(THCA_transtogene$X.kgID[grepl("^GUSB$", THCA_transtogene$geneSymbol)])

#HMBS- not found in tcga
HMBS_genes = as.character(THCA_transtogene$X.kgID[grepl("^HMBS$", THCA_transtogene$geneSymbol)])

#HPRT1
HPRT1_genes = as.character(THCA_transtogene$X.kgID[grepl("HPRT1", THCA_transtogene$geneSymbol)])

#rRNA
rRNA_genes = as.character(THCA_transtogene$X.kgID[grepl("^rRNA$", THCA_transtogene$geneSymbol)])

#transpose after adding all housekeeping genes (found ones) and then do row stats
THCA_TPM_tcga = data.frame(t(THCA_TPM_tcga))

#add row of column means at the TOP of the dataframe
THCA_TPM_tcga <- rbind(colMeans(THCA_TPM_tcga), THCA_TPM_tcga)
#name the column added as "mean"
rownames(THCA_TPM_tcga)[1] <- "mean"
#filter out all columns but mean into df
tcga_tpm_df = THCA_TPM_tcga[!rownames(THCA_TPM_tcga) %in% "mean",]
#apply t-test to each column in the dataframe tcga_tpm_df
tcga_tpm_tests <-lapply(colnames(tcga_tpm_df), function(x) t.test(tcga_tpm_df[,x]))
#separate only the confidence interval values from the t-tests stored in tcga_tpm_tests
tpm_confvals <- unlist(lapply(tcga_tpm_tests,FUN=function(x) x$conf.int - x$estimate))
#take every other value in the list 
tpm_conf = tpm_confvals[ c(F,T) ] 
#add conf +- values into dataframe as CI+-
THCA_TPM_tcga <- rbind(tpm_conf, THCA_TPM_tcga)
rownames(THCA_TPM_tcga)[1] <- "CI +-"
#log2 + 1 transformation for the gene means
tpm_meantrans <-log2(subset(THCA_TPM_tcga, rownames(THCA_TPM_tcga) == "mean") + 1)
#append meantrans to df
THCA_TPM_tcga <- rbind(tpm_meantrans, THCA_TPM_tcga)
rownames(THCA_TPM_tcga)[1] <- "trans mean"
#log2 + 1 trans for ci vals

#gets all conf int values (both low & upper)
tpm_conflimits <-unlist(lapply(tcga_tpm_tests,FUN=function(x) x$conf.int))
#gets only lower conf int vals
tpm_conflow <-tpm_conflimits[ c(T,F) ]
#transforms the lower lims log2+1
tpm_conflowtrans <- log2(tpm_conflow + 1)
#subtracts trans lower lims from trans mean to get +- trans CI val
tpm_conftrans<- subset(THCA_TPM_tcga, rownames(THCA_TPM_tcga) == "trans mean")-tpm_conflowtrans
#append conftrans to df & name trans ci
THCA_TPM_tcga <- rbind(tpm_conftrans, THCA_TPM_tcga)
rownames(THCA_TPM_tcga)[1] <- "trans ci"


#rmse:
rmse(subset(THCA_TPM_tcga, rownames(THCA_TPM_tcga) == "trans mean1"), HK_Genes_Data$TPM)

#IMMUNE CATEGORIZED DATA FOR THCA NORMAL + TUMOR TISSUE
#NEEDS to be run in order
#T cells
# expr_rows <- (grep(paste0(".+_", T_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_norm_ic <-THCA_norm[expr_rows,]
# THCA_norm_ic$ic = "Tcell"
# 
# #CD8 + Tcells
# expr_rows <- (grep(paste0(".+_", TCD8_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "T+CD8cells"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #Tregcells
# expr_rows <- (grep(paste0(".+_", Treg_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "Tregcells"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #Bcells
# expr_rows <- (grep(paste0(".+_", B_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "Bcells"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #NKcells
# expr_rows <- (grep(paste0(".+_", NK_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "NKcells"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #ClasMon
# expr_rows <- (grep(paste0(".+_", ClasMon_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "ClasMon"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #IntMon
# expr_rows <- (grep(paste0(".+_", IntMon_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "IntMon"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #NonclasMon
# expr_rows <- (grep(paste0(".+_", NonclasMon_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "NonclasMon"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #MatDC
# expr_rows <- (grep(paste0(".+_", MatDC_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "MatDC"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #MonDC
# expr_rows <- (grep(paste0(".+_", MonDC_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "MonDC"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #MatNeut
# expr_rows <- (grep(paste0(".+_", MatNeut_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "MatNeut"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# #ImNeut
# expr_rows <- (grep(paste0(".+_", ImNeut_vec, "$",collapse = "|"), THCA_norm$gene_id))
# THCA_ictemp<-THCA_norm[expr_rows,]
# THCA_ictemp$ic = "ImNeut"
# THCA_norm_ic <- rbind(THCA_norm_ic, THCA_ictemp)
# 
# THCA_norm_ic <- THCA_norm_ic[,-2]
# #pdf(paste0("THCA_Normal_IC.pdf"), height=190, width=20)
# #grid.table(THCA_norm_ic)
# #dev.off()
# THCA_norm_ic$FPKM<-log2(THCA_norm_ic$FPKM + 1)
# THCA_norm_ic$FPKM_ci_lower_bound<-log2(THCA_norm_ic$FPKM_ci_lower_bound + 1)
# THCA_norm_ic$FPKM_ci_upper_bound<-log2(THCA_norm_ic$FPKM_ci_upper_bound + 1)
# #REMOVE gene PAR_Y_P2RY8
# THCA_norm_ic<- subset(THCA_norm_ic, !rownames(THCA_norm_ic) %in% "15228")
# THCA_norm_genes = unlist(lapply(THCA_norm_ic$gene_id,function(x) regmatches(x, regexpr("[A-Z,\\p{Pd},A-Z, 0-9]+$",x, perl=TRUE))))
# THCA_norm_ic <- cbind(gene = THCA_norm_genes, THCA_norm_ic)
# write.xlsx(THCA_norm_ic, file = "/Users/shreyashekhar/GRIPS/THCA_norm_ic.xlsx", sheetName="THCA_norm_ic")


THCA_norm <- immunesepstf(THCA_normex)
#REMOVE gene PAR_Y_P2RY8
THCA_norm<- subset(THCA_norm, !rownames(THCA_norm) %in% "15228")
#put gene names as first column
THCA_norm_genes = unlist(lapply(THCA_norm$gene_id,function(x) regmatches(x, regexpr("[A-Z,\\p{Pd},A-Z, 0-9]+$",x, perl=TRUE))))
THCA_norm <- cbind(gene = THCA_norm_genes, THCA_norm)



#TUMOR TISSUE

# #T cells
# expr_rows <- (grep(paste0(".+_", T_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_tum_ic <-THCA_tum[expr_rows,]
# THCA_tum_ic$ic = "Tcell"
# 
# #CD8 + Tcells
# expr_rows <- (grep(paste0(".+_", TCD8_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "T+CD8cells"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #Tregcells
# expr_rows <- (grep(paste0(".+_", Treg_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "Tregcells"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #Bcells
# expr_rows <- (grep(paste0(".+_", B_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "Bcells"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #NKcells
# expr_rows <- (grep(paste0(".+_", NK_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "NKcells"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #ClasMon
# expr_rows <- (grep(paste0(".+_", ClasMon_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "ClasMon"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #IntMon
# expr_rows <- (grep(paste0(".+_", IntMon_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "IntMon"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #NonclasMon
# expr_rows <- (grep(paste0(".+_", NonclasMon_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "NonclasMon"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #MatDC
# expr_rows <- (grep(paste0(".+_", MatDC_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "MatDC"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #MonDC
# expr_rows <- (grep(paste0(".+_", MonDC_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "MonDC"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #MatNeut
# expr_rows <- (grep(paste0(".+_", MatNeut_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "MatNeut"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# #ImNeut
# expr_rows <- (grep(paste0(".+_", ImNeut_vec, "$",collapse = "|"), THCA_tum$gene_id))
# THCA_ictemp<-THCA_tum[expr_rows,]
# THCA_ictemp$ic = "ImNeut"
# THCA_tum_ic <- rbind(THCA_tum_ic, THCA_ictemp)
# 
# THCA_tum_ic <- THCA_tum_ic[,-2]
# #pdf(paste0("THCA_Tumor_IC.pdf"), height=190, width=20)
# #grid.table(THCA_tum_ic)
# #dev.off()
# THCA_tum_ic$FPKM<-log2(THCA_tum_ic$FPKM + 1)
# THCA_tum_ic$FPKM_ci_lower_bound<-log2(THCA_tum_ic$FPKM_ci_lower_bound + 1)
# THCA_tum_ic$FPKM_ci_upper_bound<-log2(THCA_tum_ic$FPKM_ci_upper_bound + 1)
# THCA_tum_ic<- subset(THCA_tum_ic, !rownames(THCA_tum_ic) %in% "15228")
# THCA_tum_genes = unlist(lapply(THCA_tum_ic$gene_id,function(x) regmatches(x, regexpr("[A-Z,\\p{Pd},A-Z, 0-9]+$",x, perl=TRUE))))
# THCA_tum_ic <- cbind(gene = THCA_tum_genes, THCA_tum_ic)
# write.xlsx(THCA_tum_ic, file = "/Users/shreyashekhar/GRIPS/THCA_tum_ic.xlsx", sheetName="THCA_tum_ic")


THCA_tum <- immunesepstf(THCA_tumex)
#REMOVE gene PAR_Y_P2RY8
THCA_tum<- subset(THCA_tum, !rownames(THCA_tum) %in% "15228")
#put gene names as first column
THCA_tum_genes = unlist(lapply(THCA_tum$gene_id,function(x) regmatches(x, regexpr("[A-Z,\\p{Pd},A-Z, 0-9]+$",x, perl=TRUE))))
THCA_tum <- cbind(gene = THCA_tum_genes, THCA_tum)



#GTEX Manip
rownames(GTEX_Ann) <- GTEX_Ann$SAMPID
gtex_ids = rownames(subset(GTEX_Ann, SMTS =="Thyroid"))
GTEX_thyroid<- GTEX_Ann[rownames(GTEX_Ann) %in% gtex_ids,]

write.xlsx(GTEX_thyroid, file = "/Users/shreyashekhar/GRIPS/GTEX_thyroid.xlsx", sheetName="GTEX_thyroid")

#Patients with Normal + Autoimmune Disease
rownames(GTEX_thyroid_Grouped) <- GTEX_thyroid_Grouped$SAMPID
gtex_normids = rownames(subset(GTEX_thyroid_Grouped, SMUBRID =="N"))
gtex_autoids =rownames(subset(GTEX_thyroid_Grouped, SMUBRID =="A"))

GTEX_thyroid_norm<- GTEX_thyroid_Grouped[rownames(GTEX_thyroid_Grouped) %in% gtex_normids,]
GTEX_thyroid_auto<- GTEX_thyroid_Grouped[rownames(GTEX_thyroid_Grouped) %in% gtex_autoids,]

# #colSums GTEX_genes_data
# rownames(GTEX_genes_data)<- GTEX_genes_data$gene
# GTEX_genes_data<- GTEX_genes_data[,-1]
# 
# GTEX_sample_sums = data.frame(colSums(GTEX_genes_data))
# GTEX_sample_sums <- cbind(sample = rownames(GTEX_sample_sums), GTEX_sample_sums)
# write.xlsx(GTEX_sample_sums, file = "/Users/shreyashekhar/GRIPS/GTEX_sample_sums.xlsx", sheetName="GTEX_sample_sums")

#GTEX immune genes + effective length

#only run code below if ran colSums before, to reset the data to normal
#GTEX_genes_data <- cbind(gene = rownames(GTEX_genes_data), GTEX_genes_data)
#imgenes!!!
GTEX_imgenes <- immunesep(GTEX_genes_data)
#deriving gene names from the gene_ids in GTEX_imgenes
#gene.names = unlist(lapply(imgenes_stlength$gene_id,function(x) regmatches(x, regexpr("[A-Z,0-9]+$",x))))
# gene.names = unlist(lapply(imgenes_stlength$gene_id,function(x) regmatches(x, regexpr("[A-Z,\\p{Pd},A-Z, 0-9]+$",x, perl=TRUE))))
# imgenes_stlength$gene<- gene.names
# 
# #run after, make sure GTEX_imgenes is fresh out of immunesep (format needs to be exact)
# #GTEX_imgenes$gene_names <- GTEX_imgenes$gene
# GTEX_imgenes<- GTEX_imgenes %>% dplyr::left_join(dplyr::distinct(imgenes_stlength, gene, .keep_all = T))
# #GTEX_imgenes now has cols gene (names) and gene_id(ensg..+gene name) and ic(immune ctype) and eff length(of gene)
# #only one row(gene) has NA---> MCEMP1
# #one gene_id is also NA
# #replace gene MCEMP1 NA with value 3026
# GTEX_imgenes$eff_length[is.na(GTEX_imgenes$eff_length)] <- 3026
# write.xlsx(GTEX_imgenes, file = "/Users/shreyashekhar/GRIPS/GTEX_imgenes.xlsx", sheetName="GTEX_imgenes")
# 
# #GTEX data expected counts to FPKM 
# GTEX_imexpr = GTEX_imgenes[,!colnames(GTEX_imgenes) %in% c("ic", "gene", "gene_id", "eff_length")]
# #dividing each column of data by its respective colsum/1,000,000 (divide by per million for each sample)
# GTEX_fpkm<-sweep(GTEX_imexpr,MARGIN=2,FUN="/",STATS=colSums(GTEX_imexpr)/1000000)
# ##divides values in each row by respective effective length of gene/1,000
# GTEX_fpkm <-sweep(GTEX_fpkm,MARGIN=1,FUN="/",STATS=(GTEX_imgenes$eff_length)/1000)

#not rlly needed
GTEX_imexpr = GTEX_imgenes[,!colnames(GTEX_imgenes) %in% c("ic", "gene")]

# #############take from GTEX_imgenes

# #Getting gtex norm + auto
colnames(GTEX_imgenes) = gsub('[\\.]','-',colnames(GTEX_imgenes))
# #gtex normal patients- found 377/473 -- no gene names or ic
gtex_normex<- GTEX_imgenes[,colnames(GTEX_imgenes) %in% gtex_normids]
# #gtex autoimmune disease patients-found 28/40-- no gene names or ic
gtex_autoex<- GTEX_imgenes[,colnames(GTEX_imgenes) %in% gtex_autoids]

#GTEX norm&auto with gene names & ic
GTEX_norm <- gtex_normex
GTEX_norm <- cbind(gene = GTEX_imgenes$gene, GTEX_norm)
GTEX_norm$ic <- GTEX_imgenes$ic

GTEX_auto<- gtex_autoex
GTEX_auto <- cbind(gene = GTEX_imgenes$gene, GTEX_auto)
GTEX_auto$ic <- GTEX_imgenes$ic



##CREATE function which does this on its own- input: rows and expr, output- avg, cilow, cihigh
#attach ic and gene after !!

createstats = function(fpkm_df){
  fpkm_stats_df <- data.frame(rowMeans(fpkm_df))
  fpkm_stats <-lapply(rownames(fpkm_df), function(x) t.test(fpkm_df[x,]))
  fpkm_conf <-unlist(lapply(fpkm_stats,FUN=function(x) x$conf.int))
  
  fpkm_stats_df$cilow <-fpkm_conf[ c(T,F) ] 
  fpkm_stats_df$cihigh <-fpkm_conf[ c(F,T) ] 
  
  return (fpkm_stats_df)
  
}

#Creating gtex norm + auto datasets
GTEX_norm_stats<-createstats(gtex_normex)
GTEX_auto_stats<-createstats(gtex_autoex)

GTEX_norm_stats$ic <- GTEX_imgenes$ic
GTEX_norm_stats<- cbind(gene = GTEX_imgenes$gene, GTEX_norm_stats)
GTEX_norm_stats$rowMeans.fpkm_df.<-log2(GTEX_norm_stats$rowMeans.fpkm_df. + 1)
GTEX_norm_stats$cilow<-log2(GTEX_norm_stats$cilow + 1)
GTEX_norm_stats$cihigh<-log2(GTEX_norm_stats$cihigh + 1)

GTEX_auto_stats$ic <- GTEX_imgenes$ic
GTEX_auto_stats<- cbind(gene = GTEX_imgenes$gene, GTEX_auto_stats)
#setting na values to 0
GTEX_auto_stats$cihigh[is.na(GTEX_auto_stats$cihigh)]<- 0
GTEX_auto_stats$cilow[is.na(GTEX_auto_stats$cilow)]<- 0
#set cilow to 0 for IGLL1 to 0 since it is neg (wont transform)
GTEX_auto_stats[333, "cilow"]<- 0
GTEX_auto_stats$rowMeans.fpkm_df.<-log2(GTEX_auto_stats$rowMeans.fpkm_df. + 1)
GTEX_auto_stats$cilow<-log2(GTEX_auto_stats$cilow + 1)
GTEX_auto_stats$cihigh<-log2(GTEX_auto_stats$cihigh + 1)

#input = dataframe w/ all genes --fpkm, cilow, ,cihigh, ic
#output= dataframe w/ rownames = class + avg rank 
rankic = function(ge_stdf){
  ge_stdf$rank<- rank(-ge_stdf$FPKM, na.last = "keep", ties.method = "random")
  class_rank <- data.frame(mean(ge_stdf$rank[grepl("Tcell", ge_stdf$ic)]))
  rownames(class_rank)<- "Tcell"
  class_rank["T+CD8cells",] <- (mean(ge_stdf$rank[grepl("CD8", ge_stdf$ic)]))
  class_rank["Tregcells",] <- (mean(ge_stdf$rank[grepl("Tregcells", ge_stdf$ic)]))
  class_rank["Bcells",] <- (mean(ge_stdf$rank[grepl("Bcells", ge_stdf$ic)]))
  class_rank["NKcells",] <- (mean(ge_stdf$rank[grepl("NKcells", ge_stdf$ic)]))
  class_rank["ClasMon",] <- (mean(ge_stdf$rank[grepl("ClasMon", ge_stdf$ic)]))
  class_rank["IntMon",] <- (mean(ge_stdf$rank[grepl("IntMon", ge_stdf$ic)]))
  class_rank["NonclasMon",] <- (mean(ge_stdf$rank[grepl("NonclasMon", ge_stdf$ic)]))
  class_rank["MatDC",] <- (mean(ge_stdf$rank[grepl("MatDC", ge_stdf$ic)]))
  class_rank["MonDC",] <- (mean(ge_stdf$rank[grepl("MonDC", ge_stdf$ic)]))
  class_rank["MatNeut",] <- (mean(ge_stdf$rank[grepl("MatNeut", ge_stdf$ic)]))
  class_rank["ImNeut",] <- (mean(ge_stdf$rank[grepl("ImNeut", ge_stdf$ic)]))
  
  return (class_rank)
  
}


#Ranking between all 5 datasets
#GTEX_auto_stats, GTEX_norm_stats, THCA_norm_ic, THCA_tum_ic, THCA_tcga
#
tcga_genes <- data.frame(paste(THCA_tcga$gene,THCA_tcga$ic,sep="_"))
colnames(tcga_genes)<- "gene"
stf_genes <- data.frame(paste(THCA_norm_ic$gene,THCA_norm_ic$ic,sep="_"))
colnames(stf_genes)<- "gene"
gtex_genes <- data.frame(paste(GTEX_auto_stats$gene,GTEX_auto_stats$ic,sep="_"))
colnames(gtex_genes)<- "gene"

allgenes<- data.frame(tcga_genes, stringsAsFactors = FALSE)
g_temp <- data.frame(setdiff(stf_genes$gene, allgenes$gene), stringsAsFactors = FALSE)
colnames(g_temp)<- "gene"
allgenes <- rbind( allgenes,g_temp)
g_temp <- data.frame(setdiff(gtex_genes$gene, allgenes$gene), stringsAsFactors = FALSE)
colnames(g_temp)<- "gene"
allgenes<- rbind(allgenes,g_temp)

allgenes <-cSplit(allgenes, "gene", "_")

write.xlsx(allgenes, file = "/Users/shreyashekhar/GRIPS/allgenes.xlsx", sheetName="allgenes")


#creates rows with default rank 5 for each of 5 dataframes---NOT NEEDED
#allgenes$GTEX_norm_rank <-5
#allgenes$GTEX_auto_rank <-5
#allgenes$THCA_norm_rank <-5
#allgenes$THCA_tum_rank <-5
#allgenes$THCA_tcga_rank <-5

allex_temp<-lapply(allgenes$gene_1,FUN=function(x) GTEX_norm_stats$rowMeans.fpkm_df.[head(grep(x, GTEX_norm_stats$gene), 1)] )
allex_temp[allex_temp=="numeric(0)"]<- NA
allex<- data.frame(unlist(allex_temp))
colnames(allex)<- "gtex_norm"

allex_temp<-lapply(allgenes$gene_1,FUN=function(x) GTEX_auto_stats$rowMeans.fpkm_df.[head(grep(x, GTEX_auto_stats$gene), 1)] )
allex_temp[allex_temp=="numeric(0)"]<- NA
allex$gtex_auto <- as.numeric((unlist(allex_temp)))

allex_temp <- lapply(allgenes$gene_1,FUN=function(x) THCA_norm_ic$FPKM[head(grep(x, THCA_norm_ic$gene), 1)] )
allex_temp[allex_temp=="numeric(0)"]<- NA
allex$thca_norm <- as.numeric((unlist(allex_temp)))

allex_temp <- lapply(allgenes$gene_1,FUN=function(x) THCA_tum_ic$FPKM[head(grep(x, THCA_tum_ic$gene), 1)] )
allex_temp[allex_temp=="numeric(0)"]<- NA
allex$thca_tum <- as.numeric((unlist(allex_temp)))

allex_temp <- lapply(allgenes$gene_1,FUN=function(x) THCA_tcga$rowMeans.tcga_ictemp.[head(grep(x, THCA_tcga$gene), 1)] )
allex_temp[allex_temp=="numeric(0)"]<- NA
allex$thca_tcga <- as.numeric((unlist(allex_temp)))

#df with all genes and expression saved
allex_xlsx <-cbind(gene = allgenes$gene_1, allex)
write.xlsx(allex_xlsx, file = "/Users/shreyashekhar/GRIPS/allex.xlsx", sheetName="allex")


ex_rank<- lapply(rownames(allex), FUN=function(x) rank(-allex[x,], na.last = "keep"))
allex_ranks <- data.frame(t(sapply(ex_rank,c)))
allex_ranks[is.na(allex_ranks)] <- 5

#adding columns for gene & ic in the ranking df
allex_ranks<- cbind(gene = allgenes$gene_1, allex_ranks)
allex_ranks$ic <- allgenes$gene_2

#input-- df with columns for each of the dataframes used and ranks in each row
#also input name of df to get means for as a "" string
#output-- df with name of ics as rownames and mean ranks next to them 
meanrank = function(rankdf, tdf){
  
  class_rank <- data.frame(mean(rank_df[,tdf][grepl("Tcell", rank_df$ic)]))
rownames(class_rank)<- "Tcell"
class_rank["T+CD8cells",] <- (mean(rank_df[,tdf][grepl("CD8", rank_df$ic)]))
class_rank["Tregcells",] <- (mean(rank_df[,tdf][grepl("Tregcells", rank_df$ic)]))
class_rank["Bcells",] <- (mean(rank_df[,tdf][grepl("Bcells", rank_df$ic)]))
class_rank["NKcells",] <- (mean(rank_df[,tdf][grepl("NKcells", rank_df$ic)]))
class_rank["ClasMon",] <- (mean(rank_df[,tdf][grepl("ClasMon", rank_df$ic)]))
class_rank["IntMon",] <- (mean(rank_df[,tdf][grepl("IntMon", rank_df$ic)]))
class_rank["NonclasMon",] <- (mean(rank_df[,tdf][grepl("NonclasMon", rank_df$ic)]))
class_rank["MatDC",] <- (mean(rank_df[,tdf][grepl("MatDC", rank_df$ic)]))
class_rank["MonDC",] <- (mean(rank_df[,tdf][grepl("MonDC", rank_df$ic)]))
class_rank["MatNeut",] <- (mean(rank_df[,tdf][grepl("MatNeut", rank_df$ic)]))
class_rank["ImNeut",] <- (mean(rank_df[,tdf][grepl("ImNeut", rank_df$ic)]))

return (class_rank)

}

gtex_norm_ranks <- meanrank(allex_ranks, "gtex_norm")
write.xlsx(gtex_norm_ranks, file = "/Users/shreyashekhar/GRIPS/gtex_norm_ranks.xlsx", sheetName="gtex_norm_ranks", row.names = TRUE)

gtex_auto_ranks<-meanrank(allex_ranks, "gtex_auto")
write.xlsx(gtex_auto_ranks, file = "/Users/shreyashekhar/GRIPS/gtex_auto_ranks.xlsx", sheetName="gtex_auto_ranks", row.names = TRUE)

thca_norm_ranks<-meanrank(allex_ranks, "thca_norm")
write.xlsx(thca_norm_ranks, file = "/Users/shreyashekhar/GRIPS/thca_norm_ranks.xlsx", sheetName="thca_norm_ranks", row.names = TRUE)

thca_tum_ranks<-meanrank(allex_ranks, "thca_tum")
write.xlsx(thca_tum_ranks, file = "/Users/shreyashekhar/GRIPS/thca_tum_ranks.xlsx", sheetName="thca_tum_ranks", row.names = TRUE)

thca_tcga_ranks<-meanrank(allex_ranks, "thca_tcga")
write.xlsx(thca_tcga_ranks, file = "/Users/shreyashekhar/GRIPS/thca_tcga_ranks.xlsx", sheetName="thca_tcga_ranks", row.names = TRUE)

  
#EBSEQQQQQQQQQQQQQQQQ

#find intersecting genes for all datasets::
genes_mlist<- intersect(intersect(GTEX_norm$gene, THCA_norm$gene), THCA_tcga$gene)
  
#extract only intersecting genes from all and set equal  
GTEX_norm <-subset(GTEX_norm,GTEX_norm$gene %in% genes_mlist)
GTEX_auto <-subset(GTEX_auto,GTEX_auto$gene %in% genes_mlist)

THCA_norm<-subset(THCA_norm,THCA_norm$gene %in% genes_mlist)
THCA_tum<-subset(THCA_tum,THCA_tum$gene %in% genes_mlist)

THCA_tcga <-subset(THCA_tcga,THCA_tcga$gene %in% genes_mlist)

#combined dataset of all samples
GeneMat = cbind( GTEX_norm[,!colnames(GTEX_norm) %in% c("ic", "gene")],  GTEX_auto[,!colnames(GTEX_auto) %in% c("ic", "gene")],  THCA_norm[,colnames(THCA_norm) %in% c("expected_count")],  THCA_tum[,colnames(THCA_tum) %in% c("expected_count")] , THCA_tcga[,!colnames(THCA_tcga) %in% c("ic", "gene")] )
#need to find out how to input rownames??

 
  
#gtex norm vs thca tcga

gmat<- cbind( GTEX_norm[,!colnames(GTEX_norm) %in% c("ic", "gene")],THCA_tcga[,!colnames(THCA_tcga) %in% c("ic", "gene")])
gmat_sizes <- QuantileNorm(gmat, 0.75)

gmat <- as.matrix(gmat)
rownames(gmat) <- NULL
as.data.frame(gmat)
rownames(gmat) <- GTEX_norm$gene

Conditions = as.factor(append(rep("GTEX_norm", 377), rep("THCA_norm", 557)))

EB_norm_tcga=EBTest(Data=as.matrix(gmat),Conditions = as.factor(append(rep("GTEX_norm", 377), rep("THCA_tcga", 557))),sizeFactors= gmat_sizes, maxround=10)

#try with only t-cells

gmat<- cbind( GTEX_norm[,!colnames(GTEX_norm) %in% c("ic", "gene")],THCA_tcga[,!colnames(THCA_tcga) %in% c("ic", "gene")])

gtex_norm_t = GTEX_norm[rownames(subset(GTEX_norm, ic=="Tcell")),]
thca_tcga_t = THCA_tcga[rownames(subset(THCA_tcga, ic=="Tcell")),]
  
gmat_t<- cbind( gtex_norm_t[,!colnames(gtex_norm_t) %in% c("ic", "gene")], thca_tcga_t[,!colnames(thca_tcga_t) %in% c("ic", "gene")])


gmat_t <- as.matrix(gmat_t)
rownames(gmat_t) <- NULL
as.data.frame(gmat_t)
rownames(gmat_t) <- as.character(GTEX_norm$gene[grepl("Tcell", GTEX_norm$ic)])

EB_15_t =EBTest(Data=as.matrix(gmat_t),Conditions = as.factor(append(rep("GTEX_norm", 377), rep("THCA_tcga", 557))),sizeFactors= gmat_sizes, maxround=10)

EB_15_tres=GetDEResults(EB_15_t, FDR=0.05)

#GTEX DATA start point(removed empty row)
#QuantileNorm- 0.75 library sizes
#ALL library size factors in order gtex norm-377, gtex auto-28, thca norm-1, thca tum-1, and thca tcga-557
#removing row with "" as rowname from GTEX_grnes_data (its a null row + rowname)
GTEX_genes_data<- GTEX_genes_data[-24493,]
colnames(GTEX_genes_data) = gsub('[\\.]','-',colnames(GTEX_genes_data))
GTEX_normsizes <- QuantileNorm(GTEX_genes_data[,colnames(GTEX_genes_data) %in% gtex_normids], 0.75)
GTEX_autosizes  <- QuantileNorm(GTEX_genes_data[,colnames(GTEX_genes_data) %in% gtex_autoids], 0.75)
THCA_sizes <- QuantileNorm(data.frame(cbind(THCA_normex$expected_count, THCA_tumex$expected_count)), 0.75)

GTEX_allsizes <- QuantileNorm(GTEX_genes_data[,colnames(GTEX_genes_data) %ni% c("gene")], 0.75)

#separate data frame for gtex norm and gtex auto with ALL 24,492 genes & add gene names
GTEX_normag<-GTEX_genes_data[,colnames(GTEX_genes_data) %in% gtex_normids]
GTEX_autoag<- GTEX_genes_data[,colnames(GTEX_genes_data) %in% gtex_autoids]
  
GTEX_normag<- cbind(gene = GTEX_genes_data$gene, GTEX_normag)
GTEX_autoag<- cbind(gene = GTEX_genes_data$gene, GTEX_autoag)

#GTEX_all dataframe with extra "1" added to end of every column for difference when comparing with other gtex sets (norm)
GTEX_allag<-GTEX_genes_data[,colnames(GTEX_genes_data) %ni% c("gene")]
colnames(GTEX_allag) <- paste( colnames(GTEX_allag), "1", sep = "")
GTEX_allag<- cbind(gene = GTEX_genes_data$gene, GTEX_allag)

#remove the NA line of gene ATRIP in tcga_raw:
tcga_raw= tcga_raw[-which(rownames(tcga_raw) %in% c("ATRIP|84126")), ]

tcga_sizes <- QuantileNorm(tcga_raw[,!colnames(tcga_raw) %in% c("gene")], 0.75)

#add gene names to THCA_normex & THCA-tumex
thca_ag= unlist(lapply(THCA_normex$gene_id,function(x) regmatches(x, regexpr("[A-Z,\\p{Pd},A-Z,a-z, 0-9]+$",x, perl=TRUE))))

THCA_normex<- cbind(gene = thca_ag, THCA_normex)
THCA_tumex<- cbind(gene = thca_ag, THCA_tumex)


#can have all genes & ic in dfs
#input: two dataframes containing comparable data,cell type,  Conditions =, large sizefactors matrix w/ needed ones extracted
#output:: pdf with results-- alpha, beta, P over iterations, QQplot, DE genes with type

EB_dea = function(a_df, b_df,ctype, Conditions, libsizefactors ){
  
  #taking only cells of specific ic type
  a_df = a_df[rownames(subset(a_df, ic==ctype)),]
  b_df = b_df[rownames(subset(b_df, ic==ctype)),]
  
  #join dataframes by gene & ic with a_df ids first & b_df ids second
  c_df<- a_df %>% dplyr::left_join(dplyr::distinct(b_df, gene, .keep_all = T))
  #store gene names in a list
  gene_list<- c_df$gene
  #remove all columns but the ids (remove gene & ic)
  c_df <- data.frame(c_df[grepl(paste(c("GTEX", "TCGA", "expected"),collapse="|"), colnames(c_df))])
  #make rownames equal to genes list
  rownames(c_df)<- gene_list
  
  
  EB_out =EBTest(Data=as.matrix(c_df),Conditions = Conditions ,sizeFactors= libsizefactors, maxround=10)
  
  EBDERes=GetDEResults(EB_out, FDR=0.05)
  
  # pl <- list()
  # pl[[1]]= data.frame(EB_out$Alpha)
  # pl[[2]]= data.frame(EB_out$Beta)
  # pl[[3]]= data.frame(EB_out$P)
  # pl[[4]]= data.frame(EBDERes$DEfound)
  # pl[[5]]= QQP(EB_out)
  
  
  #pdf(paste0("DE-Analysis",strsplit(levels(Conditions), ","),strsplit(levels(Conditions), ",")[2]," for ", ctype, ".pdf"), height=20, width=8.5)
  # grid.table(data.frame(EB_out$Alpha))
  # grid.table(data.frame(EB_out$Beta))
  # grid.table(data.frame(EB_out$P))
  # grid.table(data.frame(EBDERes$DEfound))
  # QQP(EB_out)
  
 # dev.off()
  #lay <- rbind(c(1,2),
               #c(2),
               #c(3),
               #c(4),
              # c(5))
  #ml = marrangeGrob(grob = pl, nrow=5, ncol=1)
  #ggsave(filename = "DE-Analysis",strsplit(levels(Conditions), ","),strsplit(levels(Conditions), ",")[2]," for ", ctype, ".pdf" ,device = pdf, ml, height = 40, width =15, units = "in", dpi = 300)
   #lout = grid.arrange(grobs = pl, layout_matrix = lay, top= paste("DE-Analysis",strsplit(levels(Conditions), ","),strsplit(levels(Conditions), ",")[2]," for ", ctype, ".pdf")) 
  #ggsave(filename = "DE-Analysis",strsplit(levels(Conditions), ","),strsplit(levels(Conditions), ",")[2]," for ", ctype, ".pdf", device = pdf, plot=ml, width =10, height = 40, units = "in", dpi = 300)

  Table1 <- tableGrob(data.frame(EB_out$Alpha))
  Table2 <- tableGrob(data.frame(EB_out$Beta))
  Table3 <- tableGrob(data.frame(EB_out$P))
  Table4 <- tableGrob(data.frame(EBDERes$DEfound))
  
    
  pdf(paste0("DE-Analysis",strsplit(levels(Conditions), ","),strsplit(levels(Conditions), ",")[2]," for ", ctype, ".pdf"), height = 10, width = 5)
  grid.arrange(Table1, Table2, ncol = 1, nrow = 2)
  #grid::grid.newpage()
  grid.arrange(Table3, ncol = 1, nrow = 1)
  grid.arrange(Table4,  ncol = 1, nrow = 1)
  #grid::grid.newpage()
  QQP(EB_out)
  dev.off()
  
}


EB_dea(GTEX_norm, THCA_tcga, "Tcell", as.factor(append(rep("GTEX_norm", 377), rep("THCA_tcga", 557))), append(GTEX_normsizes, tcga_sizes))

EB_dea(GTEX_norm, THCA_tcga, "T+CD8cells", as.factor(append(rep("GTEX_norm", 377), rep("THCA_tcga", 557))), append(GTEX_normsizes, tcga_sizes))



#find all gene intersection to do transcript wide 
all_trome<- intersect(intersect(GTEX_autoall$gene, THCA_normex$gene), tcga_raw$gene)

#maybe make expression sets for all 5 dfs with all_trome genes
#make new function for all_trome genes and run it for gtexnorm & thcatcga
#look into DE genes

GTEX_normall <-subset(GTEX_normag,GTEX_normag$gene %in% all_trome)
GTEX_autoall<-subset(GTEX_autoag,GTEX_autoag$gene %in% all_trome)

THCA_normall<- subset(THCA_normex,THCA_normex$gene %in% all_trome)
THCA_tumall<-subset(THCA_tumex,THCA_tumex$gene %in% all_trome)

#remove duplicated gene name rows in THCA_normall & THCA_tumall

THCA_normall<- THCA_normall %>% distinct(gene, .keep_all = TRUE)
THCA_tumall<- THCA_tumall %>% distinct(gene, .keep_all = TRUE)
colnames(THCA_tumall)[3]<- "expected_count.1"

THCA_tcgall<-subset(tcga_raw,tcga_raw$gene %in% all_trome)

#EBseq function for all_trome genes

EB_deall = function(a_df, b_df, Conditions, libsizefactors ){
  
  #join dataframes by gene & ic with a_df ids first & b_df ids second
  c_df<- a_df %>% dplyr::left_join(dplyr::distinct(b_df, gene, .keep_all = T))
  #store gene names in a list
  gene_list<- c_df$gene
  #remove all columns but the ids (remove gene & ic)
  c_df <- data.frame(c_df[grepl(paste(c("GTEX", "TCGA", "expected"),collapse="|"), colnames(c_df))])
  #make rownames equal to genes list
  rownames(c_df)<- gene_list
  
  
  EB_out =EBTest(Data=as.matrix(c_df),Conditions = Conditions ,sizeFactors= libsizefactors, maxround=10)
  
  EBDERes=GetDEResults(EB_out, FDR=0.05)
  GeneFC=PostFC(EB_out)
  
  Table1 <- tableGrob(data.frame(EB_out$Alpha))
  Table2 <- tableGrob(data.frame(EB_out$Beta))
  Table3 <- tableGrob(data.frame(EB_out$P))
  Table4 <- tableGrob(data.frame(EBDERes$DEfound))
  
  maxrow = 45
  npages = ceiling(nrow(Table4)/maxrow);
  
  
  #probably move this under the for loop so DE_df can go in as table 4 instead
  pdf(paste0("DE-Analysis",strsplit(levels(Conditions), ","),strsplit(levels(Conditions), ",")[2]," alltrans", ".pdf"), height = 15, width = 5)
  grid.arrange(Table1, Table2, ncol = 1, nrow = 2)
  
  grid.arrange(Table3, ncol = 1, nrow = 1)
  for (i in 1:npages){
    #idx = seq(1+((i-1)*maxrow), i*maxrow)
    from_idx = 1+((i-1)*maxrow)
    to_idx = i*maxrow
    if (to_idx > nrow(Table4)) {
      to_idx = nrow(Table4);
    }
    #grid.newpage()
    #grid.table(Table4[idx, ])
    grid.arrange(Table4[from_idx:to_idx,],  ncol = 1, nrow = 1)
  }
  #grid::grid.newpage()
  QQP(EB_out)
  par(mar=c(1,2,1,2))
  PlotPostVsRawFC(EB_out,GeneFC)
  
  dev.off()
  
  DE_genes <- EBDERes$DEfound
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
  #return(DE_df)
  #resList <- list("DEgenes" = DE_df, "normdf"= data.frame(EB_out$DataNorm))
  #return(resList)
  return(DE_df)
  
}

EB_gtexnorm_tcga <- EB_deall(GTEX_normall, THCA_tcgall,as.factor(append(rep("GTEX_norm", 377), rep("THCA_tcga", 557))), append(GTEX_normsizes, tcga_sizes) )
EB_tnorm_tum <- EB_deall(THCA_normall, THCA_tumall,as.factor(append(rep("THCA_norm", 1), rep("THCA_tum", 1))), THCA_sizes )
EB_tcga_tum <- EB_deall(THCA_tcgall, THCA_tumall,as.factor(append(rep("THCA_tcga", 557), rep("THCA_tum", 1))), append(tcga_sizes, THCA_sizes[2]) )
EB_gtexnorm_tnorm <- EB_deall(GTEX_normall, THCA_normall,as.factor(append(rep("GTEX_norm", 377), rep("THCA_norm", 1))), append(GTEX_normsizes, THCA_sizes[1]) )

#GTEX comparisons done with ALL genes, since obtained from same experiment (GTEx) & have same genes
EB_gtexnorm_auto <- EB_deall(GTEX_normag, GTEX_autoag,as.factor(append(rep("GTEX_norm", 377), rep("GTEX_auto", 28))), append(GTEX_normsizes, GTEX_autosizes) )

#GTEX comparison of ALL samples with normal + all genes
EB_gtexall_gnorm <- EB_deall(GTEX_allag, GTEX_normag,as.factor(append(rep("GTEX_all", 446), rep("GTEX_norm", 377))), append(GTEX_allsizes, GTEX_normsizes) )

#edgeR TMM norm for all 5 dfs
#ALL library size factors in order gtex norm-377, gtex auto-28, thca norm-1, thca tum-1, and thca tcga-557
adf <- cbind(GTEX_normall, GTEX_autoall,THCA_normall,THCA_tumall, THCA_tcgall, by='gene')
rownames(adf)<- adf$gene
adf <- data.frame(adf[grepl(paste(c("GTEX", "TCGA", "expected"),collapse="|"), colnames(adf))])

adf_tmm <- calcNormFactors(adf , method = "TMM", refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,Acutoff=-1e10)

#TMM for only GTEX data involving all 24,492 genes
gtexdf<-cbind(GTEX_normag, GTEX_autoag, by = 'gene')
rownames(gtexdf)<- gtexdf$gene
gtexdf <- data.frame(gtexdf[grepl(paste(c("GTEX", "TCGA", "expected"),collapse="|"), colnames(gtexdf))])
gtex_tmm <- calcNormFactors(gtexdf , method = "TMM", refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,Acutoff=-1e10)


#TMM for all GTEX comparison with normal gtexnc(normcomp)df
gtexncdf<-cbind(GTEX_allag, GTEX_normag, by = 'gene')
rownames(gtexncdf)<- gtexncdf$gene
gtexncdf <- data.frame(gtexncdf[grepl(paste(c("GTEX", "TCGA", "expected"),collapse="|"), colnames(gtexncdf))])
gtexnc_tmm <- calcNormFactors(gtexncdf , method = "TMM", refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,Acutoff=-1e10)




#TRY DESEQ FOR ALL
adf_ceil <- ceiling(adf)

adf_pt <- data.frame(colnames(adf))
adf_pt$condition <- c(rep("GTEX_norm", 377),rep("GTEX_auto", 28), rep("THCA_norm", 1), rep("THCA_tum", 1), rep("THCA_tcga", 557))
adf_deseq<- DESeqDataSetFromMatrix(adf_ceil , adf_pt, design= ~condition)
sizeFactors(adf_deseq)<- adf_tmm


register(MulticoreParam(6))

adf_deseq <- DESeq(adf_deseq, parallel = TRUE, BPPARAM=MulticoreParam(6))
adf_deseqres <- results(adf_deseq)
gn_tcga <- results(adf_deseq, contrast=c("condition","GTEX_norm","THCA_tcga"))
tn_tt <- results(adf_deseq, contrast=c("condition","THCA_norm","THCA_tum"), parallel = TRUE)
tcga_tt <- results(adf_deseq, contrast=c("condition","THCA_tcga","THCA_tum"), parallel = TRUE)
gn_tn <- results(adf_deseq, contrast=c("condition","GTEX_norm","THCA_norm"), parallel = TRUE)

#gn_ga <- results(adf_deseq, contrast=c("condition","GTEX_norm","GTEX_auto"), parallel = TRUE)

#DESeq for GTEX data only wth all 24492 genes
gtexdf_ceil <- ceiling(gtexdf)

gtexdf_pt <- data.frame(colnames(gtexdf))
gtexdf_pt$condition <- c(rep("GTEX_norm", 377),rep("GTEX_auto", 28))
gtexdf_deseq<- DESeqDataSetFromMatrix(gtexdf_ceil , gtexdf_pt, design= ~condition)
sizeFactors(gtexdf_deseq)<- gtex_tmm

gtexdf_deseq <- DESeq(gtexdf_deseq, parallel = TRUE, BPPARAM=MulticoreParam(6))
gn_ga <- results(gtexdf_deseq, contrast=c("condition","GTEX_norm","GTEX_auto"), parallel = TRUE)

gtexdf_vst <- varianceStabilizingTransformation(gtexdf_deseq, blind=FALSE)

#DESeq for all GTEX with GTEX_norm
gtexncdf_ceil <- ceiling(gtexncdf)

gtexncdf_pt <- data.frame(colnames(gtexncdf))
gtexncdf_pt$condition <- c(rep("GTEX_all", 446),rep("GTEX_norm", 377))
gtexncdf_deseq<- DESeqDataSetFromMatrix(gtexncdf_ceil , gtexncdf_pt, design= ~condition)
sizeFactors(gtexncdf_deseq)<- gtexnc_tmm

gtexncdf_deseq <- DESeq(gtexncdf_deseq, parallel = TRUE, BPPARAM=MulticoreParam(6))
gall_gn <- results(gtexncdf_deseq, contrast=c("condition","GTEX_all","GTEX_norm"), parallel = TRUE)

gtexncdf_vst <- varianceStabilizingTransformation(gtexncdf_deseq, blind=FALSE)


#adf_deseqredo<- DESeqDataSetFromMatrix(adf_ceil , adf_pt, design= ~condition)
#adf_deseqvst <-varianceStabilizingTransformation(adf_deseqredo, blind=FALSE)
  
  
#load("/Users/shreyashekhar/Downloads/tn_tt.csv")
#load("/Users/shreyashekhar/Downloads/tcga_tt.csv")
#load("/Users/shreyashekhar/Downloads/gn_tn.csv")

#ALL res tables in order of p.adj values
gn_tcga_ord<- data.frame(gn_tcga[order(gn_tcga$padj),])
tn_tt_ord<-data.frame(tn_tt[order(tn_tt$padj),])
tcga_tt_ord<-data.frame(tcga_tt[order(tcga_tt$padj),])
gn_tn_ord<-data.frame(gn_tn[order(gn_tn$padj),])

gn_ga_ord <- data.frame(gn_ga[order(gn_ga$padj),])

gall_gn_ord<- data.frame(gall_gn[order(gall_gn$padj),])

#remove 26 genes with NA because of outliers/etc. 
#gn_tcga_ord<- na.omit(gn_tcga_ord)
#hist of log2foldchange
hist(gn_tcga_ord$log2FoldChange, breaks = 60, col = "red")

hist(gn_tcga_ord$pvalue, breaks = 50, col = "red")

hist(gn_tcga_ord$padj, breaks = 50, col = "red")
plotMA(adf_deseq)
plotMA(gn_tcga)



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


#input- DESeq results turned into a df object 
#output: all de genes + their immune type (if applicable) also, plots histogram + plotMA?? in a pdf

extractDEG = function(DEres, comp1, comp2){
  
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

gn_tcga_de <- extractDEG(gn_tcga_ord, "GTEX_norm", "THCA_tcga")
tn_tt_de <- extractDEG(tn_tt_ord, "THCA_norm", "THCA_tum")
tcga_tt_de <- extractDEG(tcga_tt_ord, "THCA_tcga", "THCA_tum")
gn_tn_de <- extractDEG(gn_tn_ord, "GTEX_norm", "THCA_norm")

gn_ga_de <- extractDEG(gn_ga_ord, "GTEX_norm", "GTEX_auto")

gall_gn_de<- extractDEG(gall_gn_ord, "GTEX_all", "GTEX_norm")

#input: output from extractDEG or EB_deall(the DE_df)(make sure for EB_deall dfs, the genes are the ROWNAMES of the df)
#output: only immune cells in the DE list
onlyimmune = function(DEG_df){
  
  imlabel <- isimmune(rownames(DEG_df))
  imlabel <- subset(imlabel, trimws(unlist.iclist.) !="")
  imlabel<- imlabel[order(imlabel$unlist.iclist.),]
  return(imlabel)
  
}

gn_tcga_imde<- onlyimmune(gn_tcga_de)
tn_tt_imde <- onlyimmune(tn_tt_de)
tcga_tt_imde<- onlyimmune(tcga_tt_de)
gn_tn_imde<- onlyimmune(gn_tn_de)

gn_ga_imde <- onlyimmune(gn_ga_de)

gall_gn_imde<- onlyimmune(gall_gn_de)

#relies on the fact that EB_list$DEgenes has at least 1 immune gene
# DEcomptable = function(EB_list, DES_df, numsamps1,samp1, samp2){
#   
#   norm_data<- EB_list$normdf
#   
#   EB_ig= EB_list$DEgenes
#   EB_ig<- subset(EB_ig, trimws(unlist.iclist.) !="")
#   
#   all_ig <- rbind(DES_df, EB_ig)
#   all_ig <- all_ig[order(all_ig$unlist.iclist.),]
#   
#   de_ig <- as.character(all_ig$DE_genes)
#   
#   #wont allow to take mean if only one column selected--
#   samp1means <- rowMeans(norm_data[rownames(norm_data) %in% de_ig, 1:numsamps1])
#   samp2means <- rowMeans(norm_data[rownames(norm_data) %in% de_ig, as.numeric(numsamps1+1):ncol(norm_data)])
#   
#   all_ig[,samp1] = samp1means
#   all_ig[,samp2] = samp2means
#   
#   #hkgenes = c("CHMP2A", "PSMB2", "RAB7A", "REEP5", "SNRPD3")
#   #hkmeans_df <- data.frame(hkgenes)
#   #hkmeans_df[,samp1]<- rowMeans(norm_data[rownames(norm_data) %in% hkgenes, 1:numsamps1])
#   #hkmeans_df[,samp2]<- rowMeans(norm_data[rownames(norm_data) %in% hkgenes, as.numeric(numsamps1+1):ncol(norm_data)])
#   
#   
# }





#_________________ Heatmaps______________________________________________________

#DESeq data objects for all four analyses

DEseqData = function(df1, df2, sdata, sfacts){
  
  combdf<- df1 %>% dplyr::left_join(dplyr::distinct(df2, gene, .keep_all = T))
  gene_list<- combdf$gene
  combdf <- data.frame(combdf[grepl(paste(c("GTEX", "TCGA", "expected"),collapse="|"), colnames(combdf))])
  rownames(combdf)<- gene_list
  
  ptdf <- data.frame(colnames(combdf))
  ptdf$condition <- sdata
  
  combdf <- ceiling(combdf)
  
  df_deseq<- DESeqDataSetFromMatrix(combdf , ptdf, design= ~condition)
  
  sizeFactors(df_deseq)<- sfacts
  
  return(df_deseq)
  
}


#getting all differentially expressed immune genes from both EBSeq and DESeq
#edited to eliminate redundancy
getDEIG = function(EB_df, imde_df){
  
  EB_ig<- subset(EB_df, trimws(unlist.iclist.) !="")
  
  
  all_ig <- rbind(imde_df, EB_ig)
  all_ig <- all_ig[order(all_ig$unlist.iclist.),]
  
  all_ig<- all_ig %>% distinct(DE_genes, .keep_all = TRUE)
  
  return(all_ig)
  
}

#getting all differentially expressed immune genes from both EBSeq and DESeq
DEIG_gntcga<- getDEIG(EB_gtexnorm_tcga, gn_tcga_imde)
DEIG_tntt<-getDEIG(EB_tnorm_tum, tn_tt_imde)
DEIG_tcgatt<-getDEIG(EB_tcga_tum, tcga_tt_imde)
DEIG_gntn<-getDEIG(EB_gtexnorm_tnorm, gn_tn_imde)

DEIG_gnga<- getDEIG(EB_gtexnorm_auto, gn_ga_imde)

DEIG_gallgn<- getDEIG(EB_gtexall_gnorm, gall_gn_imde)

#DDS_gn_tcga<- DEseqData(GTEX_normall, THCA_tcgall, c(rep("GTEX_norm", 377), rep("THCA_tcga", 557)), append(adf_tmm[1:377], adf_tmm[408:964]))
#RLD_gn_tcga <- rlog(DDS_gn_tcga)

#adf_rld <- rlog(adf_deseq, blind = FALSE)
#adf_vst <- getVarianceStabilizedData(adf_deseq, blind = FALSE)
#did the vst transform on the cluster and loaded the data here- is a DESeqTransform Object
#load("/Users/shreyashekhar/Downloads/vst_deseq.csv")

#GENE REFERENCE dataframes
adf_generef <-  data.frame(rownames(adf_deseqres))
colnames(adf_generef)<- "gene"

gtexdf_generef <- data.frame(rownames(gtexdf_deseq))
colnames(gtexdf_generef)<- "gene"



#NEED TO SPECIFY THE GENEREF DATAFRAME(GENES + ROWNUMS) IN INPUT

heatde = function(cond1df, cond2df,cond1, cond2, DEIG, vst_df, generef){
  ids <- append(colnames(cond1df), colnames(cond2df))
  ids = gsub('-','\\.',ids)
  
  rownums = unlist(lapply(DEIG$DE_genes,function(x) grep( paste0("^", x ,"$"), generef$gene)))
  
  
  mat<- assay(vst_df)[ as.numeric(rownums),colnames(assay(vst_df)) %in% ids ]
  
  mat <- mat - rowMeans(mat)
  #rownames(mat) <- DEIG$DE_genes
  rownames(mat)<- paste(DEIG$DE_genes,DEIG$unlist.iclist.)
  
  meta_df <- data.frame(colData(vst_df)[,c("condition")])
  colnames(meta_df)<- "condition"
  meta_df<- subset(meta_df,  meta_df$condition %in% c(cond1, cond2))
  rownames(meta_df)<- colnames(mat)
  
  #meta_df$condition <- gsub('GTEX_norm','Normal Tissue',meta_df$condition)
  #meta_df$condition <- gsub('GTEX_auto','Hashimoto\'s Tissue',meta_df$condition)
  
  
  pdf(paste0("pheatmap", cond1, " vs ", cond2, ".pdf" ), height = 30, width = 7, pointsize = 1)
  #pdf(paste0("pheatmap", cond1, " vs ", cond2, ".pdf" ), height = 25, width = 7)
  
  #pheatmap(mat, annotation_col=meta_df, cellheight = 10)
  #pheatmap(mat, annotation_col=meta_df, cellheight = 10)
  pheatmap(mat, annotation_col=meta_df, cellheight = 4, cluster_cols = FALSE)
  
  #to get table of only our data sample
  #if ncol of cond1df> 30 and ncol of cond2df <30(this one is our data)
  mat_df<- data.frame()
  
  #make sure the # to be less than is accurate- shouldn't include an unneeded dataset(GTEX_auto is <30, but shouldn't be made a table for)
  if(ncol(cond1df) >30  & ncol(cond2df) <5 ){
    mat_df<- data.frame(mat[, grepl("expected", colnames(mat))])
    colnames(mat_df)<- cond2
    
    Table1 <- tableGrob(data.frame(mat_df))
    grid.arrange(Table1, ncol = 1, nrow = 1)
    
  }
  
  dev.off()
  
  return(mat_df)
  
}

comptbl_gntcga<- heatde(GTEX_normall, THCA_tcgall, "GTEX_norm", "THCA_tcga", DEIG_gntcga, vst_deseq, adf_generef)
comptbl_tntt<- heatde(THCA_normall, THCA_tumall, "THCA_norm", "THCA_tum", DEIG_tntt, vst_deseq, adf_generef)
comptbl_tcgatt<- heatde(THCA_tcgall, THCA_tumall, "THCA_tcga", "THCA_tum", DEIG_tcgatt, vst_deseq, adf_generef)
comptbl_gntn<- heatde(GTEX_normall, THCA_normall, "GTEX_norm", "THCA_norm", DEIG_gntn, vst_deseq, adf_generef)

comptbl_gnga<- heatde(GTEX_normag, GTEX_autoag, "GTEX_norm", "GTEX_auto", DEIG_gnga, gtexdf_vst, gtexdf_generef)

#run on ALL de genes gnga with cell height 2- had to rename only gtex pdfs before running so they wouldn't be overwritten
comptbl_allgnga<- heatde(GTEX_normag, GTEX_autoag, "GTEX_norm", "GTEX_auto", isimmune(allde_gnga), gtexdf_vst, gtexdf_generef)

#run on top 100 genes with highest variance from de genes
gtexdf_devst <- assay(gtexdf_vst)[allde_gnga,]
gtex_topVar <- rownames(gtexdf_devst[head(order(-rowVars(gtexdf_devst)),100),])
#gtex_topVar <- rownames(gtexdf_devst[head(order(-rowVars(gtexdf_devst)),50),])

comptbl_topvargnga<- heatde(GTEX_normag, GTEX_autoag, "GTEX_norm", "GTEX_auto", isimmune(gtex_topVar), gtexdf_vst, gtexdf_generef)

#for only b cells
comptbl_bgnga<- heatde(GTEX_normag, GTEX_autoag, "GTEX_norm", "GTEX_auto", subset(DEIG_gnga, grepl("Bcells", DEIG_gnga$unlist.iclist.)), gtexdf_vst, gtexdf_generef)


#for only t cells
comptbl_tgnga<- heatde(GTEX_normag, GTEX_autoag, "GTEX_norm", "GTEX_auto", subset(DEIG_gnga, grepl("Tcell", DEIG_gnga$unlist.iclist.)), gtexdf_vst, gtexdf_generef)


#for only tregcells
comptbl_treggnga<- heatde(GTEX_normag, GTEX_autoag, "GTEX_norm", "GTEX_auto", subset(DEIG_gnga, grepl("Treg", DEIG_gnga$unlist.iclist.)), gtexdf_vst, gtexdf_generef)

#GTEX_all vs GTEX_norm
comptbl_gallgn<- heatde(GTEX_allag, GTEX_normag, "GTEX_all", "GTEX_norm", DEIG_gallgn, gtexncdf_vst, gtexdf_generef)
comptbl_allgallgn<- heatde(GTEX_allag, GTEX_normag, "GTEX_all", "GTEX_norm", isimmune(allde_gallgn), gtexncdf_vst, gtexdf_generef)


HKheat = function(cond1df, cond2df, cond1, cond2, vst_df, generef){
  
  ids <- append(colnames(cond1df), colnames(cond2df))
  ids = gsub('-','\\.',ids)
  
  hklist <- c("CHMP2A", "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", "SNRPD3", "VCP", "VPS29", "GUSB", "HMBS", "HPRT1")
  
  rownums = unlist(lapply(hklist,function(x) grep( paste0("^", x ,"$"), generef$gene)))
  
  mat<- assay(vst_df)[ as.numeric(rownums),colnames(assay(vst_df)) %in% ids ]
  
  mat <- mat - rowMeans(mat)
  rownames(mat) <- hklist
  
  meta_df <- data.frame(colData(vst_df)[,c("condition")])
  colnames(meta_df)<- "condition"
  meta_df<- subset(meta_df,  meta_df$condition %in% c(cond1, cond2))
  rownames(meta_df)<- colnames(mat)
  
  pdf(paste0("hkheatmap", cond1, " vs ", cond2, ".pdf" ), height = 10, width = 7)
  
  pheatmap(mat, annotation_col=meta_df,  cellheight = 10)
  
  mat_df<- data.frame()
  
  #make sure the # to be less than is accurate- shouldn't include an unneeded dataset(GTEX_auto is <30, but shouldn't be made a table for)
  if(ncol(cond1df) >30  & ncol(cond2df) <20 ){
    mat_df<- data.frame(mat[, grepl("expected", colnames(mat))])
    colnames(mat_df)<- cond2
    
    Table1 <- tableGrob(data.frame(mat_df))
    grid.arrange(Table1, ncol = 1, nrow = 1)
    
  }
  
  dev.off()
  
  return(mat_df)
  
}

hktbl_gntcga<- HKheat(GTEX_normall, THCA_tcgall, "GTEX_norm", "THCA_tcga", vst_deseq, adf_generef)
hktbl_tntt<- HKheat(THCA_normall, THCA_tumall, "THCA_norm", "THCA_tum", vst_deseq, adf_generef)
hktbl_tcgatt<- HKheat(THCA_tcgall, THCA_tumall, "THCA_tcga", "THCA_tum", vst_deseq, adf_generef)
hktbl_gntn<- HKheat(GTEX_normall, THCA_normall, "GTEX_norm", "THCA_norm", vst_deseq, adf_generef)

hktbl_gnga<- HKheat(GTEX_normall, GTEX_autoall, "GTEX_norm", "GTEX_auto", gtexdf_vst, gtexdf_generef)

#_________________________GO______________________________

#USE UNIQUE() OF THESE BECAUSE NAMES ARE DUPLICATED WITH EBSEQ + DESEQ
allde_gntcga<-unique(append(as.character(EB_gtexnorm_tcga$DE_genes), rownames(gn_tcga_de)))
allde_tntt<-unique(append(as.character(EB_tnorm_tum$DE_genes), rownames(tn_tt_de)))
allde_tcgatt<-unique(append(as.character(EB_tcga_tum$DE_genes), rownames(tcga_tt_de)))
allde_gntn<-unique(append(as.character(EB_gtexnorm_tnorm$DE_genes), rownames(gn_tn_de)))

allde_gnga <- unique(append(as.character(EB_gtexnorm_auto$DE_genes), rownames(gn_ga_de)))

allde_gallgn <- unique(append(as.character(EB_gtexall_gnorm$DE_genes), rownames(gall_gn_de)))

printDE = function(allde, exp){

   options(max.print=999999)
  allde<- unique(allde)
  # 
  # pdf(paste0("AllDE", exp, "Genes" ), height = 30, width = 10)
  # Table1 <- tableGrob(allde)
  # grid.arrange(Table1, ncol = 1, nrow = 1)
  # #print(allde)
  # 
  # dev.off()
  # 
  # options(max.print=999)
  
  lapply(allde, write, exp, ".txt", append=TRUE, ncolumns=1000)
  
}


#MAKE SURE TO DELETE FILE BEFORE RERUNNING- INSTEAD OF REWRITING FILE THIS DUPLICATES NAMES IN FILE
printDE(allde_gntcga, "gntcga")
printDE(allde_tntt, "tntt")
printDE(allde_tcgatt, "tcgatt")
printDE(allde_gntn, "gntn")

printDE(allde_gnga, "gnga")

printDE(allde_gallgn, "gallgn")


#GTEX norm and auto analysis extra plots___________________________________


#boxplots for cd19, cd27, cd86, ptpn22-- use vst data
gtexval_df<- data.frame(t(data.frame(assay(gtexdf_vst)[c("CD19", "CD27", "CD86", "PTPN22"),])))
gtexval_df$condition<- as.character(gtexdf_pt$condition)
gtexval_df<- melt(gtexval_df)



ggplot(gtexval_df, aes(x=condition,y=value, fill = condition)) +
  geom_boxplot(outlier.size=0,show.legend = F,lwd=0.1) + labs(title="Box plots for 4 Validated Genes") +facet_wrap(~variable) +
  theme(axis.text.x=element_blank(),
        #axis.title=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=6), 
        legend.text = element_text(size=4),
        plot.title = element_text(size=10),
        legend.position="bottom") +
  scale_fill_manual(values=wes_palette(n=2, name="Darjeeling1"), na.value="yellow")


pdf(paste0("GTEX_Validated_Genes.pdf" ), height = 5, width = 7)

ggplot(gtexval_df, aes(x=condition,y=value, fill = condition)) +
  geom_boxplot(outlier.size=0,show.legend = T,lwd=0.1) + labs(title="Box plots for 4 Validated Genes", y= "Mean exp", color = "Legend Title\n") +facet_wrap(~variable) +
  theme(axis.text.x=element_blank(),
        #axis.title=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=6), 
        legend.text = element_text(size=4),
        plot.title = element_text(size=10),
        legend.position="bottom") +
  scale_color_manual(labels = c("GTEX_norm", "GTEX_auto"), values=wes_palette(n=2, name="Darjeeling1"), na.value="yellow")

dev.off()

#maybe also do boxplots for the top de genes like cd79A, pou2af1, blk, cd19
gtextop_df<- data.frame(t(data.frame(assay(gtexdf_vst)[c("CD79A", "POU2AF1", "BLK", "CD19"),])))
gtextop_df$condition<- as.character(gtexdf_pt$condition)
gtextop_df<- melt(gtextop_df)


pdf(paste0("GTEX_Top_Genes.pdf" ), height = 5, width = 7)

ggplot(gtextop_df, aes(x=condition,y=value, fill = condition)) +
  geom_boxplot(outlier.size=0,show.legend = T,lwd=0.1) + labs(title="Box plots for 4 Top Genes", y= "Mean exp", color = "Legend Title\n") +facet_wrap(~variable) +
  theme(axis.text.x=element_blank(),
        #axis.title=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=6), 
        legend.text = element_text(size=4),
        plot.title = element_text(size=10),
        legend.position="bottom") +
  scale_color_manual(labels = c("GTEX_norm", "GTEX_auto"), values=wes_palette(n=2, name="Darjeeling1"), na.value="yellow")

dev.off()



#Corrplots?

#REACTOMEPA pathway enrichment
#Using clusterProfiler to convert to EntrezIds
#1180/1378 mapped
gnga_ids <- bitr(allde_gnga, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gnga_reactome <- enrichPathway(gene=gnga_ids$ENTREZID,pvalueCutoff=0.05, readable=T)

pdf(paste0("GTEX_ReactomePA.pdf" ), height = 5, width = 12)


barplot(gnga_reactome, showCategory=8)
#grid::grid.newpage()

dotplot(gnga_reactome, showCategory=15)

dev.off()
#grid::grid.newpage()

#clusterProfiler for GO - go enrichment

#GO over-representation test
#same as go enrichment

gnga_ego <- enrichGO(gene= gnga_ids$ENTREZID, OrgDb= org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable      = TRUE)

pdf(paste0("GTEX_GOplots.pdf" ), height = 5, width = 12)

barplot(gnga_ego, showCategory=8)

dotplot(gnga_ego)

emapplot(gnga_ego)

cnetplot(gnga_ego, categorySize="pvalue")

goplot(gnga_ego)

dev.off()




#convert DE gene symbols to GO IDs for plugging into revigo

ensembl<- useMart("ensembl",dataset="hsapiens_gene_ensembl")

gnga_goids<- getBM(attributes=c('hgnc_symbol', 'go_id'), filters = 'hgnc_symbol', values= allde_gnga, mart= ensembl)

printDE(as.character(gnga_goids$go_id), "gnga_goids")


#T-test to find top genes in the set of differentially expressed immune genes
gtexdf_icdevst<- data.frame(assay(gtexdf_vst)[as.character(DEIG_gnga$DE_genes),])
colnames(gtexdf_icdevst) = gsub('[\\.]','-',colnames(gtexdf_icdevst))


gtexdf_tresults = apply(gtexdf_icdevst,1,function(x) t.test(x[grepl(paste(gtex_normids, collapse  = "|"), colnames(gtexdf_icdevst))], x[grepl(paste(gtex_autoids, collapse = "|"), colnames(gtexdf_icdevst))]))
gtexdf_tstats = unlist(lapply(gtexdf_tresults,FUN=function(x) x$statistic))

gtexdf_tlow<- sort(gtexdf_tstats, decreasing = F)
gtexdf_thigh<- sort(gtexdf_tstats, decreasing = T)

#only for b cells
gtexdf_bdiff<- data.frame(assay(gtexdf_vst)[as.character(subset(DEIG_gnga$DE_genes, grepl("Bcells", DEIG_gnga$unlist.iclist.))),])
colnames(gtexdf_bdiff) = gsub('[\\.]','-',colnames(gtexdf_bdiff))

gtexdf_btresults = apply(gtexdf_bdiff,1,function(x) t.test(x[grepl(paste(gtex_normids, collapse  = "|"), colnames(gtexdf_bdiff))], x[grepl(paste(gtex_autoids, collapse = "|"), colnames(gtexdf_bdiff))]))

gtexdf_btstats = unlist(lapply(gtexdf_btresults,FUN=function(x) x$statistic))

gtexdf_bsort<- sort(gtexdf_btstats, decreasing = F)

#only for t cells
gtexdf_tdiff<- data.frame(assay(gtexdf_vst)[as.character(subset(DEIG_gnga$DE_genes, grepl("Tcell", DEIG_gnga$unlist.iclist.))),])
colnames(gtexdf_tdiff) = gsub('[\\.]','-',colnames(gtexdf_tdiff))

gtexdf_ttresults = apply(gtexdf_tdiff,1,function(x) t.test(x[grepl(paste(gtex_normids, collapse  = "|"), colnames(gtexdf_tdiff))], x[grepl(paste(gtex_autoids, collapse = "|"), colnames(gtexdf_tdiff))]))

gtexdf_ttstats = unlist(lapply(gtexdf_ttresults,FUN=function(x) x$statistic))

gtexdf_tsort<- sort(gtexdf_ttstats, decreasing = F)

#______________________________________________________________________
#GTEX BD COMP
colnames(gtex_alldata) = gsub('[\\.]','-',colnames(gtex_alldata))
rownames(gtex_alldata)<- gtex_alldata$gene
gtex_alldata <- gtex_alldata[,-1]
gtex_alldata <- data.frame(t(gtex_alldata))
gtex_alldata$SAMPID <- rownames(gtex_alldata)

#finding overlapping ids + joining by those ids
#when joined, type of tissue was added in SMTS col
gtex_allids = intersect(as.character(rownames(gtex_alldata)), as.character(GTEX_Ann$SAMPID))

GTEX_compAnn<-GTEX_Ann[GTEX_Ann$SAMPID %in% gtex_allids,]
GTEX_compAnn<- GTEX_compAnn[,colnames(GTEX_compAnn) %in% c("SAMPID", "SMTS")]

gtex_alldata<- gtex_alldata %>% dplyr::left_join(dplyr::distinct(GTEX_compAnn, SAMPID, .keep_all = T))
rownames(gtex_alldata)<- gtex_alldata$SAMPID


#ALL GTEX DATA ALL TISSUES!!!_____________________
#load("/Users/shreyashekhar/Downloads/gtex_atd.csv")
#remove NA gene row at bottom
gtex_alltissues<- gtex_alltissues[-24493,]

gtexat<- gtex_alltissues
#TMM for all GTEX tissues
rownames(gtexat)<- gtexat$gene
gtexat <- data.frame(gtexat[grepl(paste(c("GTEX"),collapse="|"), colnames(gtexat))])
gtexat_tmm <- calcNormFactors(gtexat , method = "TMM", refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,Acutoff=-1e10)

#DESEQ matrix from gtexat
gtexat_ceil <- ceiling(gtexat)

gtexat_pt <- data.frame(colnames(gtexat))
gtexat_pt$tissue <- GTEX_compAnn$SMTS
gtexat_deseq<- DESeqDataSetFromMatrix(gtexat_ceil , gtexat_pt, design= ~tissue)
sizeFactors(gtexat_deseq)<- gtexat_tmm
#saved for cluster vst
save(gtexat_deseq, file = "/Users/shreyashekhar/Downloads/gtexat_deseq.csv")

gtexat_nzd<- counts(gtexat_deseq, normalized=TRUE)

#--str
gtexat_nzd<- data.frame(gtexat_nzd)
colnames(gtexat_nzd) = gsub('[\\.]','-',colnames(gtexat_nzd))

#**Heatmaps for gtex comparisons

#Code for normalized(not vst) data

#create normalized matrix with "norm" thyroid data appended to end (added 1's to ids)
heat_nzd<- gtexat_nzd
compdf<- heat_nzd[,colnames(heat_nzd) %in% gtex_normids]
colnames(compdf) <- paste( colnames(compdf), "1", sep = "")
heat_nzd<-cbind(heat_nzd, compdf)

#find rownums of DEG
DEG<- isimmune(allde_gallgn)
rownums = unlist(lapply(DEG$DE_genes,function(x) grep( paste0("^", x ,"$"), gtexdf_generef$gene)))


hmat<- heat_nzd[as.numeric(rownums),]
hmat <- hmat - rowMeans(hmat)
rownames(hmat)<- paste(DEG$DE_genes,DEG$unlist.iclist.)

#added 377 "Thyroid_norm" tisnames for the tyroid norm data appended to end of normalized df(head_nzd)
tisname<- append(gtexat_pt$tissue, rep("Thyroid_norm", 377))
meta_heat<- data.frame(tisname)
colnames(meta_heat)<- "tissue"
rownames(meta_heat)<- colnames(hmat)



pdf(paste0("pheatmap_gtexthycomp.pdf" ), height = 30, width = 30, pointsize = 1)

pheatmap(hmat, annotation_col=meta_heat, cellheight = 4, cluster_cols = FALSE)

dev.off()

#MEANS OF EACH TISSUE TYPE INSTEAD...

#making df of exp for only DEGs with all ids
hm_nzd<- heat_nzd[as.numeric(rownums),]
#naming ids by the tissue type instead of gtex id
colnames(hm_nzd)<- meta_heat$tissue
#taking mean of all columns with same tissue name based on gene, putting in df
tis_means<- data.frame(sapply(split.default(hm_nzd, names(hm_nzd)), rowMeans))
#heatmap with differences from mean of genes
tm_df<- tis_means
tm_df <- tm_df - rowMeans(tm_df)
rownames(tm_df)<- paste(DEG$DE_genes,DEG$unlist.iclist.)


meta_tm<- data.frame(colnames(tm_df))
colnames(meta_tm)<- "tissue"
rownames(meta_tm)<- meta_tm$tissue

pdf(paste0("pheatmap_gtexthycomp_meandiffs.pdf" ), height = 30, width = 7, pointsize = 1)

pheatmap(tm_df, annotation_col=meta_tm, cellheight = 4, cluster_cols = FALSE)

dev.off()

#heatmap with means of genes
rownames(tis_means)<- paste(DEG$DE_genes,DEG$unlist.iclist.)

pdf(paste0("pheatmap_gtexthycomp_means.pdf" ), height = 30, width = 7, pointsize = 1)

pheatmap(tis_means, annotation_col=meta_tm, cellheight = 4, cluster_cols = FALSE)

dev.off()


#remove salivary gland and see
#tis_test<- tis_means[,-22]
#meta_test<- data.frame(meta_tm[-22,])
#colnames(meta_test)<- "tissue"
#rownames(meta_test)<- meta_test$tissue

#This is better(for vis), scales across rows(genes) with z-score normalization 
#(subtracts rowmean and divides by rowSD), so each exp value becomes a z score
pdf(paste0("pheatmap_gtexthycomp_meantest.pdf" ), height = 50, width = 7, pointsize = 1)

#pheatmap(tis_test, annotation_col=meta_test, cellheight = 4, cluster_cols = FALSE, scale = "row")
pheatmap(tis_means, annotation_col=meta_tm, cellheight = 7, cluster_cols = FALSE, scale = "row")


dev.off()


#check by manual z score normalization
tm_test<- tm_df/rowSds(as.matrix(tis_means))

pdf(paste0("pheatmap_gtexthycomp_meanz.pdf" ), height = 30, width = 7, pointsize = 1)

pheatmap(tm_test, annotation_col=meta_tm, cellheight = 4, cluster_cols = FALSE)


dev.off()


heatcomp = function(compids, DEG, vst_df, generef){
  
}

#GTEX ANNOTATION files with only 

GTEX_Sann<- GTEX_Ann[,colnames(GTEX_Ann) %in% c("SAMPID","SMPTHNTS", "SMTS")]


#Colon- 539
GTEX_Colon_Ann<- subset(GTEX_Sann, SMTS =="Colon")
write.xlsx(GTEX_Colon_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Colon_Ann.xlsx", sheetName="GTEX_Colon_Ann")

#Liver- 188 received
GTEX_Liver_Ann<- subset(GTEX_Sann, SMTS =="Liver")
write.xlsx(GTEX_Liver_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Liver_Ann.xlsx", sheetName="GTEX_Liver_Ann")


#Prostate- 159 received
GTEX_Prostate_Ann<- subset(GTEX_Sann, SMTS =="Prostate")
write.xlsx(GTEX_Prostate_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Prostate_Ann.xlsx", sheetName="GTEX_Prostate_Ann")


#Skin- 1362
GTEX_Skin_Ann<- subset(GTEX_Sann, SMTS =="Skin")
write.xlsx(GTEX_Skin_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Skin_Ann.xlsx", sheetName="GTEX_Skin_Ann")



#Stomach- 272 received
GTEX_Stomach_Ann<- subset(GTEX_Sann, SMTS =="Stomach")
write.xlsx(GTEX_Stomach_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Stomach_Ann.xlsx", sheetName="GTEX_Stomach_Ann")


#Heart- 727
GTEX_Heart_Ann<- subset(GTEX_Sann, SMTS =="Heart")


#Kidney- 50
GTEX_Kidney_Ann<- subset(GTEX_Sann, SMTS =="Kidney")
write.xlsx(GTEX_Kidney_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Kidney_Ann.xlsx", sheetName="GTEX_Kidney_Ann")


#Nerve- 498
GTEX_Nerve_Ann<- subset(GTEX_Sann, SMTS =="Nerve")

#Pancreas- 268
GTEX_Pancreas_Ann<- subset(GTEX_Sann, SMTS =="Pancreas")
write.xlsx(GTEX_Pancreas_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Pancreas_Ann.xlsx", sheetName="GTEX_Pancreas_Ann")


#Lung- 607
GTEX_Lung_Ann<- subset(GTEX_Sann, SMTS =="Lung")

#Blood- 2561
GTEX_Blood_Ann<- subset(GTEX_Sann, SMTS =="Blood")

#Bladder- 11
GTEX_Bladder_Ann<- subset(GTEX_Sann, SMTS =="Bladder")

#Ovary- 138
GTEX_Ovary_Ann<- subset(GTEX_Sann, SMTS =="Ovary")
write.xlsx(GTEX_Ovary_Ann, file = "/Users/shreyashekhar/GRIPS/GTEX_Ovary_Ann.xlsx", sheetName="GTEX_Ovary_Ann")


