## zzy 毕业设计 

# platform       x86_64-w64-mingw32               
# arch           x86_64                           
# os             mingw32                          
# crt            ucrt                             
# system         x86_64, mingw32                  
# status                                          
# major          4                                
# minor          3.3                              
# year           2024                             
# month          02                               
# day            29                               
# svn rev        86002                            
# language       R                                
# version.string R version 4.3.3 (2024-02-29 ucrt)
# nickname       Angel Food Cake 

# 设置根目录
setwd('/Users/ermine/Documents/毕设_华西做了一大半的/zzy/graduation_thesis')

# plot1 比较蛋白质组和磷酸化蛋白质组表达量差异
library(gcookbook) 
library(ggplot2)
library(reshape2)

# plot2 差异分析
# BiocManager::install("edgeR")
# BiocManager::install("EnhancedVolcano")
library(edgeR)
library(EnhancedVolcano)

# 名称转换
# BiocManager::install('clusterProfiler')
# BiocManager::install('org.Hs.eg.db')
library(clusterProfiler)
library(org.Hs.eg.db)

# 画韦恩图，看都上调，先上调再下调的基因
library(VennDiagram)

# 富集分析
library(WebGestaltR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(grid)
library(futile.logger)
library(VennDiagram)

# 分子网络绘制
# BiocManager::install(c("STRINGdb","igraph"),ask = F,update = F)
library(tidyverse)  # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(STRINGdb) 
library(igraph)

####################################  蛋白质组    ###############################################################################################
pro_file = "./protein_data_reviewed/human_pro/combined/txt/proteinGroups.txt"

full_data <- read.delim(pro_file,stringsAsFactors = FALSE)
colnames(full_data)

cols <- c('Reporter.intensity.1','Reporter.intensity.2','Reporter.intensity.3','Reporter.intensity.4','Reporter.intensity.5','Reporter.intensity.6')
pro_data <- full_data[cols]
rownames(pro_data) <- full_data$Protein.IDs
## 将集中的proteinID分散到新的行
pro_data_singleProteinID <- matrix( NA,ncol=ncol(pro_data) )
colnames(pro_data_singleProteinID) <- colnames(pro_data)
all_data_list <- list()
for(i in 1:nrow(pro_data)){
  names <- rownames(pro_data)[i]
  name_list <- strsplit(names, ';')[[1]]
  for (name in name_list) {
    all_data_list[[name]] <- rbind(all_data_list[[name]], pro_data[i, ])
  }
}
pro_data_singleProteinID <- do.call(rbind, all_data_list)
rownames(pro_data_singleProteinID) <- names(all_data_list)
## 只保留非CON__或REV__开头的
pro_data_singleProteinID_no_conrev <- pro_data_singleProteinID[!grepl("^CON__|^REV__", rownames(pro_data_singleProteinID)), ]
pro_data_singleProteinID_no_conrev <- pro_data_singleProteinID_no_conrev[1:3]
colnames(pro_data_singleProteinID_no_conrev) <- c('GFP','FHH_2_25d','FHH_5d')
# saveRDS(pro_data_singleProteinID_no_conrev,'./data/pro_data_singleProteinID_no_conrev.rds')

## 转换蛋白质名称，只保留能转换的
all_protein_list <- rownames(pro_data_singleProteinID_no_conrev)
gene_UNIPROT <- bitr(
  all_protein_list,
  fromType = "UNIPROT",
  toType = c("UNIPROT", "SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)
# 判断UNIPROT列中的行是否重复
duplicate_rows <- duplicated(gene_UNIPROT$UNIPROT)
# 使用逻辑索引选择非重复的行
gene_UNIPROT <- gene_UNIPROT[!duplicate_rows, ]

unique_protein_id <- unique(gene_UNIPROT[['UNIPROT']])
protein_transfered <- pro_data_singleProteinID_no_conrev[unique_protein_id,]

####################################  磷酸化蛋白质组    ##########################################################################################
phospro_file <- "./protein_data_reviewed/human_phospro/combined/txt/proteinGroups.txt"

full_data_phos <- read.delim(phospro_file,stringsAsFactors = FALSE)

cols <- c('Reporter.intensity.1','Reporter.intensity.2','Reporter.intensity.3','Reporter.intensity.4','Reporter.intensity.5','Reporter.intensity.6')
phospro_data <- full_data_phos[cols]
rownames(phospro_data) <- full_data_phos$Protein.IDs

phospro_data_singleProteinID <- matrix( NA,ncol=ncol(phospro_data) )
colnames(phospro_data_singleProteinID) <- colnames(phospro_data)
all_data_list <- list()
for(i in 1:nrow(phospro_data)){
  names <- rownames(phospro_data)[i]
  name_list <- strsplit(names, ';')[[1]]
  for (name in name_list) {
    all_data_list[[name]] <- rbind(all_data_list[[name]], phospro_data[i, ])
  }
}
phospro_data_singleProteinID <- do.call(rbind, all_data_list)
rownames(phospro_data_singleProteinID) <- names(all_data_list)
phospro_data_singleProteinID_no_conrev <- phospro_data_singleProteinID[!grepl("^CON__|^REV__", rownames(phospro_data_singleProteinID)), ]

phospro_data_singleProteinID_no_conrev <- phospro_data_singleProteinID_no_conrev[4:6]
colnames(phospro_data_singleProteinID_no_conrev) <- c('FHH_2_25d','FHH_5d','GFP')
# saveRDS(phospro_data_singleProteinID_no_conrev,'./data/phospro_data_singleProteinID_no_conrev.rds')
## 转换蛋白质名称，只保留能转换的
all_protein_list_phos <- rownames(phospro_data_singleProteinID_no_conrev)
gene_UNIPROT_phos <- bitr(
  all_protein_list_phos,
  fromType = "UNIPROT",
  toType = c("UNIPROT", "SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)
# 判断UNIPROT列中的行是否重复
duplicate_rows <- duplicated(gene_UNIPROT_phos$UNIPROT)
# 使用逻辑索引选择非重复的行
gene_UNIPROT_phos <- gene_UNIPROT_phos[!duplicate_rows, ]

unique_protein_id_phos <- unique(gene_UNIPROT_phos[['UNIPROT']])
phosprotein_transfered <- phospro_data_singleProteinID_no_conrev[unique_protein_id_phos,]

####################################  plot1:比较蛋白质组和磷酸化蛋白质组表达量差异    ################################################################
protein_data_for_plot <- pro_data_singleProteinID_no_conrev
phospro_data_for_plot <- phospro_data_singleProteinID_no_conrev

all_plot_list <- list()
for(condition in colnames(protein_data_for_plot)){
  intersect_rownames <- intersect(rownames(protein_data_for_plot),rownames(phospro_data_for_plot))
  phos_tmp <- phospro_data_for_plot[intersect_rownames,][[condition]]
  pro_tmp <- protein_data_for_plot[intersect_rownames,][[condition]]
  
  dataframe_for_plot <- data.frame(phosphoproteome=phos_tmp,proteome=pro_tmp)
  
  plot <- ggplot(dataframe_for_plot, aes(x=proteome, y=phosphoproteome)) +
    geom_point(size=0.2) + 
    ggtitle(condition) +
    xlim(0, 2500000) +
    ylim(0, 2500000) +
    stat_smooth(method=lm, level=0.99) #回归线
  
  all_plot_list[[condition]] <- plot
}
final_plot <- all_plot_list[['GFP']]+all_plot_list[['FHH_2_25d']]+all_plot_list[['FHH_5d']]
# ggsave('./要交的文件/results/plots/1蛋白质磷酸化蛋白质表达分布.pdf',final_plot,width = 13, height = 3,limitsize = FALSE)


####################################  转录组    ####################################################################################################
rna_file = "./rna/all.rawcount.tsv"
runtable_path = "./rna/SraRunTable.txt"

runTable <- read.table(runtable_path,sep=',')
# View(runTable)

rna_data <- read.table(rna_file, header=T, sep = '\t', fill=0)
rna_data <- data.frame(rna_data)
# 将第一列设置为行名
rownames(rna_data) <- rna_data[, 1]
# 删除第一列
rna_data <- rna_data[, -1]
new_order <- c("SRR14078360","SRR14078361","SRR14078362","SRR14078363","SRR14078364","SRR14078365","SRR14078366","SRR14078367")
rna_data <- rna_data[,new_order]
# 只要前3列，后面的是过表达PIM/GFP的
rna_data <- rna_data[,1:3]
colnames(rna_data) <- c('GFP','FHH_2_25d','FHH_5d')

# 转换id至ezid
all_rna_list <- rownames(rna_data)
gene_UNIPROT_rna <- bitr(
  all_rna_list,
  fromType = "SYMBOL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)
# 判断UNIPROT列中的行是否重复
duplicate_rows <- duplicated(gene_UNIPROT_rna$SYMBOL)
# 使用逻辑索引选择非重复的行
gene_UNIPROT_rna <- gene_UNIPROT_rna[!duplicate_rows, ]

rna_transfered <- rna_data[gene_UNIPROT_rna$SYMBOL,]


####################################  保存三个组学名称转换文件  ############################
name_transfer <- list()
name_transfer[['rna']] <- gene_UNIPROT_rna
name_transfer[['pro']] <- gene_UNIPROT[-656,]
name_transfer[['phospro']] <- gene_UNIPROT_phos[-486,]
# saveRDS(name_transfer,'./data/name_transfer.rds')
####################################  差异分析 #####################################################
# edger 函数
edger <- function(countmatrix){
  countMatrix <- countmatrix
  group <- factor(colnames(countMatrix))
  y <- DGEList(counts = countMatrix, group = group)
  exprSet  <- y
  
  bcv = 0.1#设置bcv为0.1
  et <- exactTest(exprSet, dispersion=bcv^2)
  DEG_edgeR=as.data.frame(topTags(et, n = nrow(exprSet$counts)))
  
  #筛选上下调，设定阈值
  fc_cutoff <- 1
  pvalue <- 0.05
  DEG_edgeR$regulated <- "normal"
  loc_up <- intersect(which(DEG_edgeR$logFC>log2(fc_cutoff)),
                      which(DEG_edgeR$PValue<pvalue))
  loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                        which(DEG_edgeR$PValue<pvalue))
  DEG_edgeR$regulated[loc_up] <- "up"
  DEG_edgeR$regulated[loc_down] <- "down" 
  table(DEG_edgeR$regulated)
  
  return(DEG_edgeR)
}

## 差异分析
condition_names <- c('GFP','FHH_2_25d','FHH_5d')
compare_condition_names <- c('F225_G225','F5_G225','F5_F225')
compare_result <- c('up','down','normal')


rna_data <- rna_transfered
pro_data <- protein_transfered
phospro_data <- phosprotein_transfered

data_for_edger <- list()

data_for_edger[['rna']][['F225_G225']] <- rna_data[c('FHH_2_25d','GFP')]
data_for_edger[['rna']][['F5_G225']] <- rna_data[c('FHH_5d','GFP')]
data_for_edger[['rna']][['F5_F225']] <- rna_data[c('FHH_5d','FHH_2_25d')]

data_for_edger[['pro']][['F225_G225']] <- pro_data[c('FHH_2_25d','GFP')]
data_for_edger[['pro']][['F5_G225']] <- pro_data[c('FHH_5d','GFP')]
data_for_edger[['pro']][['F5_F225']] <- pro_data[c('FHH_5d','FHH_2_25d')]

data_for_edger[['phospro']][['F225_G225']] <- phospro_data[c('FHH_2_25d','GFP')]
data_for_edger[['phospro']][['F5_G225']] <- phospro_data[c('FHH_5d','GFP')]
data_for_edger[['phospro']][['F5_F225']] <- phospro_data[c('FHH_5d','FHH_2_25d')]


edger_result_all <- list()
for(type in c('rna','pro','phospro')){
  for(compare in c('F225_G225','F5_G225','F5_F225')){
    edger_result_all[[type]][[compare]] <- edger(data_for_edger[[type]][[compare]])
  }
}

# saveRDS(edger_result_all,'./data/edger_result_all.rds')

####################################  表4:各组学上调下调基因数量 ###############################
edger_result_all <- readRDS('./data/edger_result_all.rds')
up_down_gene_number <- list()
for(type in c('rna','pro','phospro')){
  for(compare in c('F225_G225','F5_G225','F5_F225')){
    tmp <- edger_result_all[[type]][[compare]]
    up_down_gene_number[[type]][[compare]][['up']] <- nrow(tmp[tmp$regulated == 'up',])
    up_down_gene_number[[type]][[compare]][['down']] <- nrow(tmp[tmp$regulated == 'down',])
  }
}

# saveRDS(up_down_gene_number,'./data/up_down_gene_number.rds')

####################################  plot2:差异分析火山图 #####################################################
# 火山图函数
edger_result_all <- readRDS('./data/edger_result_all.rds')
EnhancedVolcano_plot <- function(name,data){
  plot <- EnhancedVolcano(
    data,
    x = "logFC",
    y = "FDR",
    lab = rownames(data),    # 基因名
    pCutoff = 0.05,      # p值截断值：水平线
    FCcutoff = 1,         # FC截断值：垂直线
    cutoffLineWidth = 0.7,      # 虚线粗细
    cutoffLineType = "twodash", # 线的类型
    xlim = c(-2.5, 2.5),  # x轴起始
    ylim = c(0, -log10(10e-9)), # y轴起始
    pointSize = 3,        # 点的大小
    labSize = 2.5,        # 标签大小
    xlab = bquote(~Log[2] ~ "fold change"),     # 此行为默认x轴名
    ylab = bquote(~-Log[10] ~ italic(p-value)), # 修改y轴名为斜体
    axisLabSize = 12,     # 坐标轴字体大小
    title = name,     # 修改标题
    titleLabSize = 16,    # 标题大小
    subtitle = bquote(italic("Volcano plot")),  # 修改子标题：斜体
    subtitleLabSize = 14, # 子标题大小
    legendLabSize = 11,   # legend大小
    col = c("grey20", "darkgreen", "royalblue", "red3"), # legend颜色
    colAlpha = 0.4,       # 透明度
    gridlines.major = FALSE,     # 背景网格
    gridlines.minor = FALSE
  )
  return(plot)
}

all_plot_list <- list()
for(type in c('rna','pro','phospro')){
  type_plot_list <- list()
  for(compare in c('F225_G225','F5_G225','F5_F225')){
    data <- edger_result_all[[type]][[compare]]
    name <- paste(type,compare,sep = '_')
    type_plot_list[[compare]] <- EnhancedVolcano_plot(name = name,data = data)
  }
  all_plot_list[[type]] <- type_plot_list[['F225_G225']]+type_plot_list[['F5_G225']]+type_plot_list[['F5_F225']]
}
final_volano_plot <- all_plot_list[['rna']]/all_plot_list[['pro']]/all_plot_list[['phospro']]

ggsave('./要交的文件/results/plots/2差异分析火山图.pdf',final_volano_plot,width = 20, height = 20,limitsize = FALSE)

####################################  plot3_1-3_4:都上调的基因和先上调再下调的基因 #####################################################
## 收集上调下调基因集
edger_result_all <- readRDS('./data/edger_result_all.rds')

up_lists <- list()
down_lists <- list()
for(type in names(edger_result_all)){
  data_this_type <- edger_result_all[[type]]
  up_lists[[type]][['F225_G225']] <- rownames(data_this_type$F225_G225[data_this_type$F225_G225[,'regulated']=='up',])
  up_lists[[type]][['F5_G225']] <- rownames(data_this_type$F5_G225[data_this_type$F5_G225[,'regulated']=='up',])
  up_lists[[type]][['F5_F225']] <- rownames(data_this_type$F5_F225[data_this_type$F5_F225[,'regulated']=='up',])
  
  
  down_lists[[type]][['F225_G225']] <- rownames(data_this_type$F225_G225[data_this_type$F225_G225[,'regulated']=='down',])
  down_lists[[type]][['F5_G225']] <- rownames(data_this_type$F5_G225[data_this_type$F5_G225[,'regulated']=='down',])
  down_lists[[type]][['F5_F225']] <- rownames(data_this_type$F5_F225[data_this_type$F5_F225[,'regulated']=='down',])
  
}
# saveRDS(up_lists,file = './data/up_regulates_all.rds')
# saveRDS(down_lists,file = './data/down_regulates_all.rds')


# 看总的上调下调基因个数有无差别，结果是没有差别
# rna_up_number 10405 rna_down_number 9843
# pro_up_number 1899 pro_down_number 1758
# phospro_up_number 631 phospro_down_number 668
rna_up_number <- length(unique(append(append(up_lists$rna$F225_G225,up_lists$rna$F5_G225),up_lists$rna$F5_F225)))
rna_down_number <- length(unique(append(append(down_lists$rna$F225_G225,down_lists$rna$F5_G225),down_lists$rna$F5_F225)))
pro_up_number <- length(unique(append(append(up_lists$pro$F225_G225,up_lists$pro$F5_G225),up_lists$pro$F5_F225)))
pro_down_number <- length(unique(append(append(down_lists$pro$F225_G225,down_lists$pro$F5_G225),down_lists$pro$F5_F225)))
phospro_up_number <- length(unique(append(append(up_lists$phospro$F225_G225,up_lists$phospro$F5_G225),up_lists$phospro$F5_F225)))
phospro_down_number <- length(unique(append(append(down_lists$phospro$F225_G225,down_lists$phospro$F5_G225),down_lists$phospro$F5_F225)))


## 韦恩图
venn_plot_2data <- function(data1,data2,sets_names,plot_name){
  sets <- list( set1 = data1, set2 = data2)
  names(sets) <- sets_names
  tmp_plot <- venn.diagram(
    x = sets, filename = NULL,
    col="white", 
    fill=c(colors()[616], colors()[38]),
    alpha=c(0.6, 0.6), 
    lwd=c(1, 1), 
    cex=2, 
    cat.dist=c(0.05, 0.05), 
    cat.pos=c(0, 0), 
    cat.cex=2,
    main = plot_name
  )
  return(tmp_plot)
}
venn_plot_3data <- function(data1,data2,data3,sets_names,plot_name){
  sets <- list( set1 = data1, set2 = data2, set3 = data3 )
  names(sets) <- sets_names
  tmp_plot <- venn.diagram(
    x = sets, filename = NULL,
    col="white", 
    fill=c(colors()[616], colors()[38], colors()[468]),
    alpha=c(0.6, 0.6, 0.6), 
    lwd=c(1, 1, 1), 
    cex=2, 
    cat.dist=c(0.05, 0.05, -0.45), 
    cat.pos=c(0, 0, 0), 
    cat.cex=2,
    main = plot_name
  )
  return(tmp_plot)
}

## 画关于compare条件的交集
set_names <- c('F225_G225','F5_G225','F5_F225')
up_plots <- list()
down_plots <- list()
for(type in c('rna','pro','phospro')){
  this_plot_name <- paste(type,' intersection',sep='')
  up_plots[[type]] <- venn_plot_3data(data1 = up_lists[[type]][['F225_G225']],
                                      data2 =  up_lists[[type]][['F5_G225']],
                                      data3 =  up_lists[[type]][['F5_F225']],
                                      sets_names = set_names,
                                      plot_name = this_plot_name)
  down_plots[[type]] <- venn_plot_3data(data1 = down_lists[[type]][['F225_G225']],
                                        data2 =  down_lists[[type]][['F5_G225']],
                                        data3 =  down_lists[[type]][['F5_F225']],
                                        sets_names = set_names,
                                        plot_name = this_plot_name)
}
up_plot_all <- cowplot::plot_grid(up_plots$rna)+cowplot::plot_grid(up_plots$pro)+cowplot::plot_grid(up_plots$phospro)
down_plot_all <- cowplot::plot_grid(down_plots$rna)+cowplot::plot_grid(down_plots$pro)+cowplot::plot_grid(down_plots$phospro)
# ggsave('./要交的文件/results/plots/3_3关于compare条件的上调交集.pdf',up_plot_all,width = 15, height = 6,limitsize = FALSE)
# ggsave('./要交的文件/results/plots/3_4关于compare条件的下调交集.pdf',down_plot_all,width = 15, height = 6,limitsize = FALSE)

## 画关于组学的交集
# 将pro和phospro转换为SYMBOL
pro_name_transfer <- gene_UNIPROT
phospro_name_transfer <- gene_UNIPROT_phos
up_plots_2 <- list()
down_plots_2 <- list()
set_names <- c('rna','pro','phospro')
for(compare in c('F225_G225','F5_G225','F5_F225')){
  this_plot_name <- paste(compare,' intersection',sep='')
  up_plots_2[[compare]] <- venn_plot_3data(data1 = up_lists[['rna']][[compare]],
                                      data2 =  pro_name_transfer[pro_name_transfer$UNIPROT %in% up_lists[['pro']][[compare]],]$SYMBOL,
                                      data3 =  phospro_name_transfer[phospro_name_transfer$UNIPROT %in% up_lists[['phospro']][[compare]],]$SYMBOL,
                                      sets_names = set_names,
                                      plot_name = this_plot_name)
  down_plots_2[[compare]] <- venn_plot_3data(data1 = down_lists[['rna']][[compare]],
                                           data2 =  pro_name_transfer[pro_name_transfer$UNIPROT %in% down_lists[['pro']][[compare]],]$SYMBOL,
                                           data3 =  phospro_name_transfer[phospro_name_transfer$UNIPROT %in% down_lists[['phospro']][[compare]],]$SYMBOL,
                                           sets_names = set_names,
                                           plot_name = this_plot_name)
}
up_plot_all_2 <- cowplot::plot_grid(up_plots_2$F225_G225)+cowplot::plot_grid(up_plots_2$F5_G225)+cowplot::plot_grid(up_plots_2$F5_F225)
down_plot_all_2 <- cowplot::plot_grid(down_plots_2$F225_G225)+cowplot::plot_grid(down_plots_2$F5_G225)+cowplot::plot_grid(down_plots_2$F5_F225)
# ggsave('./要交的文件/results/plots/3_1关于组学的上调交集.pdf',up_plot_all_2,width = 15, height = 6,limitsize = FALSE)
# ggsave('./要交的文件/results/plots/3_2关于组学的下调交集.pdf',down_plot_all_2,width = 15, height = 6,limitsize = FALSE)


####################################  plot3_5-3_6: 关于组学的先上调再下调，先下调再上调交集  ##########################################
## 画关于组学的先上调再下调，先下调再上调交集
set_names <- c('F225_G225','F5_F225')
up_down_plots <- list()
down_up_plots <- list()
for(type in c('rna','pro','phospro')){
  this_plot_name <- paste(type,' intersection',sep='')
  up_down_plots[[type]] <- venn_plot_2data(data1 = up_lists[[type]][['F225_G225']],
                                      data2 =  down_lists[[type]][['F5_F225']],
                                      sets_names = set_names,
                                      plot_name = this_plot_name)
  down_up_plots[[type]] <- venn_plot_2data(data1 = down_lists[[type]][['F225_G225']],
                                        data2 =  up_lists[[type]][['F5_F225']],
                                        sets_names = set_names,
                                        plot_name = this_plot_name)
}
up_down_plot_all <- cowplot::plot_grid(up_down_plots$rna)+cowplot::plot_grid(up_down_plots$pro)+cowplot::plot_grid(up_down_plots$phospro)
down_up_plot_all <- cowplot::plot_grid(down_up_plots$rna)+cowplot::plot_grid(down_up_plots$pro)+cowplot::plot_grid(down_up_plots$phospro)
# ggsave('./要交的文件/results/plots/3_5关于组学先上调再下调交集.pdf',up_down_plot_all,width = 15, height = 6,limitsize = FALSE)
# ggsave('./要交的文件/results/plots/3_6关于组学先下调再上调交集.pdf',down_up_plot_all,width = 15, height = 6,limitsize = FALSE)

####################################  准备给富集分析的geneList ###################################
all_gene_list_for_enrichment_analysis <- list()
# 关于组学的交集
for(type in c('rna','pro','phospro')){
  inter_up <- get.venn.partitions(list(data1 = up_lists[[type]][['F225_G225']],
                                  data2 =  up_lists[[type]][['F5_G225']],
                                  data3 =  up_lists[[type]][['F5_F225']]))
  all_gene_list_for_enrichment_analysis[[type]][['up']] <- as.character(unlist(inter_up[c(1,2,3,5),]$..values..))
  all_gene_list_for_enrichment_analysis[[type]][['all_up']] <- as.character(unlist(inter_up[c(1),]$..values..))
  
  
  inter_down <- get.venn.partitions(list(data1 = down_lists[[type]][['F225_G225']],
                                         data2 =  down_lists[[type]][['F5_G225']],
                                         data3 =  down_lists[[type]][['F5_F225']]))
  all_gene_list_for_enrichment_analysis[[type]][['down']] <- as.character(unlist(inter_down[c(1,2,3,5),]$..values..))
  all_gene_list_for_enrichment_analysis[[type]][['all_down']] <- as.character(unlist(inter_down[c(1),]$..values..))
  
  }
# 关于compare条件的交集
pro_name_transfer <- gene_UNIPROT
phospro_name_transfer <- gene_UNIPROT_phos
for(compare in c('F225_G225','F5_G225','F5_F225')){
  inter_up <- get.venn.partitions(list(data1 = up_lists[['rna']][[compare]],
                                       data2 =  pro_name_transfer[pro_name_transfer$UNIPROT %in% up_lists[['pro']][[compare]],]$SYMBOL,
                                       data3 =  phospro_name_transfer[phospro_name_transfer$UNIPROT %in% up_lists[['phospro']][[compare]],]$SYMBOL))
  all_gene_list_for_enrichment_analysis[[compare]][['up']] <- as.character(unlist(inter_up[c(1,2,3,4,5,6,7),]$..values..))
  # all_gene_list_for_enrichment_analysis[[compare]][['all_up']] <- as.character(unlist(inter_up[c(1),]$..values..))
  
  
  inter_down <- get.venn.partitions(list(data1 = down_lists[['rna']][[compare]],
                                         data2 =  pro_name_transfer[pro_name_transfer$UNIPROT %in% down_lists[['pro']][[compare]],]$SYMBOL,
                                         data3 =  phospro_name_transfer[phospro_name_transfer$UNIPROT %in% down_lists[['phospro']][[compare]],]$SYMBOL))
  all_gene_list_for_enrichment_analysis[[compare]][['down']] <- as.character(unlist(inter_down[c(1,2,3,4,5,6,7),]$..values..))
  # all_gene_list_for_enrichment_analysis[[compare]][['all_down']] <- as.character(unlist(inter_down[c(1),]$..values..))
  
}
# 上调下调
for(type in c('rna','pro','phospro')){
  inter_up_down <- get.venn.partitions(list(data1 = up_lists[[type]][['F225_G225']],
                                            data2 =  down_lists[[type]][['F5_F225']]))
  all_gene_list_for_enrichment_analysis[['up_down']][[type]] <- as.character(unlist(inter_up_down[c(1),]$..values..))

  
  inter_down_up <- get.venn.partitions(list(data1 = down_lists[[type]][['F225_G225']],
                                            data2 =  up_lists[[type]][['F5_F225']]))
  all_gene_list_for_enrichment_analysis[['down_up']][[type]] <- as.character(unlist(inter_down_up[c(1),]$..values..))
}
# saveRDS(all_gene_list_for_enrichment_analysis,'./data/all_gene_list_for_enrichment_analysis.rds')

####################################  富集分析    #####################################################
## 将上调、下调、先上调再下调的都富集到通路上
edger_result_all <- readRDS('./data/edger_result_all.rds')
all_gene_list_for_enrichment_analysis <- readRDS('./data/all_gene_list_for_enrichment_analysis.rds')

up_lists <- readRDS('./data/up_regulates_all.rds')
down_lists <- readRDS('./data/down_regulates_all.rds')

name_transfer <- readRDS('./data/name_transfer.rds')
rna_name_transfer <- name_transfer$rna
pro_name_transfer <- name_transfer$pro
phospro_name_transfer <- name_transfer$phospro

# ORA
ORA_enrichment_analysis <- function(gene_list,type,rna_name_transfer,pro_name_transfer,phospro_name_transfer){
  enrich_list <- list()
  
  gene_list <- gene_list
  type <- type
  rna_name_transfer<-rna_name_transfer
  pro_name_transfer<-pro_name_transfer
  phospro_name_transfer<-phospro_name_transfer
  
  if(type %in% c('rna','trans')){
    genes_symbol <- gene_list
    genes_id <- rna_name_transfer[rna_name_transfer$SYMBOL %in% gene_list,]$ENTREZID
  } else if(type == 'pro'){
    genes_symbol <- pro_name_transfer[pro_name_transfer$UNIPROT %in% gene_list,]$SYMBOL
    genes_id <- pro_name_transfer[pro_name_transfer$UNIPROT %in% gene_list,]$ENTREZID
  } else if(type == 'phospro'){
    genes_symbol <- phospro_name_transfer[phospro_name_transfer$UNIPROT %in% gene_list,]$SYMBOL
    genes_id <- phospro_name_transfer[phospro_name_transfer$UNIPROT %in% gene_list,]$ENTREZID
  }
  
  GO_result <- enrichGO(
    genes_symbol,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "MF",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    # universe,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  
  KEGG_result <- enrichKEGG(
    genes_id,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
  )
  
  enrich_list[['GO']] <- as.data.frame(GO_result) 
  enrich_list[['KEGG']] <- as.data.frame(KEGG_result) 
  
  return(enrich_list)
}

## GSEA(错的)
GSEA_enrichment_analysis <- function(dataframe,type,rna_name_transfer,pro_name_transfer,phospro_name_transfer){
  enrich_list <- list()
  
  dataframe <- dataframe
  type <- type
  rna_name_transfer<-rna_name_transfer
  pro_name_transfer<-pro_name_transfer
  phospro_name_transfer<-phospro_name_transfer
  
  if(type == 'rna'){
    name_transfer <- rna_name_transfer
  } else if(type == 'pro'){
    name_transfer<- pro_name_transfer
  } else if(type == 'phospro'){
    name_transfer <- phospro_name_transfer
  }
  
  geneList <- dataframe$logFC
  symbol <- name_transfer$SYMBOL
  names(geneList) <- symbol
  geneList=sort(geneList,decreasing = T) #从高到低排序
  
  gseaGO_result <- gseGO(
    geneList,
    ont = "ALL",#可以替换为BP CC MF，分别富集
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
  
  geneList <- dataframe$logFC
  ezid <- name_transfer$ENTREZID
  names(geneList) <- ezid
  geneList=sort(geneList,decreasing = T) #从高到低排序
  
  gseaKEGG_result<- gseKEGG(
    geneList,
    organism = "hsa",
    keyType = "kegg",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 1e-10,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
    use_internal_data = FALSE,
    seed = FALSE,
    by = "fgsea",
    scoreType = 'pos'   ### 不知道为啥选这个
  )
  
  enrich_list[['gseaGO']] <- as.data.frame(gseaGO_result)
  enrich_list[['gseaKEGG']] <- as.data.frame(gseaKEGG_result)
  
  return(enrich_list)
}

# ORA分析
ORA_result <- list()
for(name in c('rna','pro','phospro')){
  for(regulated_type in names(all_gene_list_for_enrichment_analysis[[name]])){
    gene_list <- all_gene_list_for_enrichment_analysis[[name]][[regulated_type]]
    type <- name
    ORA_result[[name]][[regulated_type]] <- ORA_enrichment_analysis(gene_list,type,rna_name_transfer,
                                                                    pro_name_transfer,phospro_name_transfer)
  }
}
for(name in c('F225_G225','F5_G225','F5_F225')){
  for(regulated_type in names(all_gene_list_for_enrichment_analysis[[name]])){
    gene_list <- all_gene_list_for_enrichment_analysis[[name]][[regulated_type]]
    type <- 'trans'
    ORA_result[[name]][[regulated_type]] <- ORA_enrichment_analysis(gene_list,type,rna_name_transfer,
                                                                    pro_name_transfer,phospro_name_transfer)
  }
}
for(name in c('up_down','down_up')){
  for(tome_type in names(all_gene_list_for_enrichment_analysis[[name]])){
    gene_list <- all_gene_list_for_enrichment_analysis[[name]][[tome_type]]
    type <- tome_type
    ORA_result[[name]][[tome_type]] <- ORA_enrichment_analysis(gene_list,type,rna_name_transfer,
                                                                    pro_name_transfer,phospro_name_transfer)
  }
}
# saveRDS(ORA_result,'./data/ORA_result.rds')

# GSEA分析
GSEA_result <- list()
for(type in c('rna','pro','phospro')){
  for(compare in c('F225_G225','F5_G225','F5_F225')){
    dataframe <- edger_result_all[[type]][[compare]]
    type <- type
    GSEA_result[[type]][[compare]] <- GSEA_enrichment_analysis(dataframe,type,rna_name_transfer,pro_name_transfer,phospro_name_transfer)
  }
}
# saveRDS(GSEA_result,'./data/GSEA_result.rds')

# 获取ORA分析结果中的通路列表
ORA_result_list <- list()
for(name in c('rna','pro','phospro')){
  for(regulated_type in c('up','all_up','down','all_down')){
    for(enrichment_way in c('GO','KEGG')){
      this_descriptions <- ORA_result[[name]][[regulated_type]][[enrichment_way]][['Description']]
      # this_dataframe_colname <- paste(name,regulated_type,enrichment_way,sep = '_')
      # this_dataframe <- data.frame(data = this_descriptions)
      # colnames(this_dataframe) <- this_dataframe_colname
      ORA_result_list[[name]][[regulated_type]][[enrichment_way]] <- this_descriptions
    }
  }
}
for (name in c('F225_G225','F5_G225','F5_F225')) {
  for(regulated_type in c('up','down')){
    for(enrichment_way in c('GO','KEGG')){
      this_descriptions <- ORA_result[[name]][[regulated_type]][[enrichment_way]][['Description']]
      # this_dataframe_colname <- paste(name,regulated_type,enrichment_way,sep = '_')
      # this_dataframe <- data.frame(data = this_descriptions)
      # colnames(this_dataframe) <- this_dataframe_colname
      ORA_result_list[[name]][[regulated_type]][[enrichment_way]] <- this_descriptions
    }
  }
}
for(name in c('up_down','down_up')){
  for(tome_name in c('rna','pro','phospro')){
    for(enrichment_way in c('GO','KEGG')){
      this_descriptions <- ORA_result[[name]][[tome_name]][[enrichment_way]][['Description']]
      # this_dataframe_colname <- paste(name,tome_name,enrichment_way,sep = '_')
      # this_dataframe <- data.frame(data = this_descriptions)
      # colnames(this_dataframe) <- this_dataframe_colname
      ORA_result_list[[name]][[tome_name]][[enrichment_way]] <- this_descriptions
    }
  }
}
# saveRDS(ORA_result_list,'./data/ORA_result_list.rds')


# 画富集分析条形图
ORA_result <- readRDS('./data/ORA_result.rds')
# 不同比对条件
up_plots_ORA_condition <- list()
down_plots_ORA_condition <- list()
for(name in c('F225_G225','F5_G225','F5_F225')){
  up_data <- rbind(ORA_result[[name]]$up$KEGG[c('Description','p.adjust')],
                ORA_result[[name]]$up$GO[c('Description','p.adjust')])
  up_data <- up_data[order(up_data$p.adjust),][1:20,]
  for(i in 1:nrow(up_data)){up_data[i,'Description'] <- strsplit(up_data[i,'Description'],',')[[1]][1]}
  up_plot <- up_data %>%
              # 通过reorder(site,-incid)来改变条形的顺序，按照incid值的大小排列条形
              ggplot(aes(x=reorder(Description,p.adjust), y=p.adjust)) +
              # 输出条形图，并把条形填充为蓝色
              geom_bar(stat="identity",fill="steelblue")+
              coord_flip()+
              labs(x = "Enriched Pathways")+
              ggtitle(name)
  
  down_data <- rbind(ORA_result[[name]]$down$KEGG[c('Description','p.adjust')],
                   ORA_result[[name]]$down$GO[c('Description','p.adjust')])
  down_data <- down_data[order(down_data$p.adjust),][1:20,]
  for(i in 1:nrow(down_data)){down_data[i,'Description'] <- strsplit(down_data[i,'Description'],',')[[1]][1]}
  
  down_plot <- down_data %>%
                # 通过reorder(site,-incid)来改变条形的顺序，按照incid值的大小排列条形
                ggplot(aes(x=reorder(Description,p.adjust), y=p.adjust)) +
                # 输出条形图，并把条形填充为蓝色
                geom_bar(stat="identity",fill="steelblue")+
                coord_flip()+
                labs(x = "Enriched Pathways")+
                ggtitle(name)
                
  
  up_plots_ORA_condition[[name]] <- up_plot
  down_plots_ORA_condition[[name]] <- down_plot
}
up_plots_ORA_condition_all <- up_plots_ORA_condition[['F225_G225']]+up_plots_ORA_condition[['F5_G225']]+up_plots_ORA_condition[['F5_F225']]
down_plots_ORA_condition_all <- down_plots_ORA_condition[['F225_G225']]+down_plots_ORA_condition[['F5_G225']]+down_plots_ORA_condition[['F5_F225']]

# ggsave('./要交的文件/results/plots/4_1_1上调的ORA.pdf',up_plots_ORA_condition_all,width = 16, height = 5,limitsize = FALSE)
# ggsave('./要交的文件/results/plots/4_2_1下调的ORA.pdf',down_plots_ORA_condition_all,width = 16, height = 5,limitsize = FALSE)


# 先上再下/先下再上
all_plot <- list()
for(name in c('up_down','down_up')){
  data <- rbind(ORA_result[[name]]$rna$KEGG[c('Description','p.adjust')],ORA_result[[name]]$rna$GO[c('Description','p.adjust')],
                        ORA_result[[name]]$pro$KEGG[c('Description','p.adjust')],ORA_result[[name]]$pro$GO[c('Description','p.adjust')],
                        ORA_result[[name]]$phospro$KEGG[c('Description','p.adjust')],ORA_result[[name]]$phospro$GO[c('Description','p.adjust')]
                        )
  data <- data[order(data$p.adjust),][1:20,]
  for(i in 1:nrow(data)){data[i,'Description'] <- strsplit(data[i,'Description'],',')[[1]][1]}
  plot <- data %>%
    # 通过reorder(site,-incid)来改变条形的顺序，按照incid值的大小排列条形
    ggplot(aes(x=reorder(Description,p.adjust), y=p.adjust)) +
    # 输出条形图，并把条形填充为蓝色
    geom_bar(stat="identity",fill="steelblue")+
    coord_flip()+
    labs(x = "Enriched Pathways")+
    ggtitle(name)
  
  all_plot[[name]] <- plot
}
all_plots <- all_plot[['up_down']]+all_plot[['down_up']]

# ggsave('./要交的文件/results/plots/4_3先上后下先下后上的ORA.pdf',all_plots,width = 16, height = 5,limitsize = FALSE)


####################################  用网页做GSEA #########################
edger_result_all <- readRDS('./data/edger_result_all.rds')
path <- './data/data_for_GSEA_web/'

# 全部的
for(type in c('rna','pro','phospro')){
  for(compare in c('F225_G225','F5_G225','F5_F225')){
    data <- edger_result_all[[type]][[compare]]
    dataframe <- data.frame(genes = rownames(data), logFC = data$logFC)
    data_path <- paste(path,type,'_',compare,'.rnk',sep = '')
    write.table(dataframe, file = data_path, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  }
}

# 按照compare条件来
name_transfer <- readRDS('./data/name_transfer.rds')
rna_name_transfer <- name_transfer$rna
pro_name_transfer <- name_transfer$pro
phospro_name_transfer <- name_transfer$phospro


pro_gsea_dataframe <- list()
phospro_gsea_dataframe <- list()
compare_gsea <- list()

for(compare in c('F225_G225','F5_G225','F5_F225')){
  rna_data <- edger_result_all[['rna']][[compare]]
  rna_dataframe <- data.frame(genes = rownames(rna_data), logFC = rna_data$logFC)
  
  pro_data <- edger_result_all[['pro']][[compare]]
  order <- rownames(pro_data)
  sorted_index <- order(match(pro_name_transfer$UNIPROT, order))
  sorted_pro_transfer <- pro_name_transfer[sorted_index, ]
  pro_dataframe <- data.frame(genes = sorted_pro_transfer$SYMBOL, logFC = pro_data$logFC)
  pro_gsea_dataframe[[compare]] <- pro_dataframe
  
  phospro_data <- edger_result_all[['phospro']][[compare]]
  order <- rownames(phospro_data)
  sorted_index <- order(match(phospro_name_transfer$UNIPROT, order))
  sorted_phospro_transfer <- phospro_name_transfer[sorted_index, ]
  phospro_dataframe <- data.frame(genes = sorted_phospro_transfer$SYMBOL, logFC = phospro_data$logFC)
  phospro_gsea_dataframe[[compare]] <- phospro_dataframe
  
  compare_gsea[[compare]] <- rbind(rna_dataframe,rbind(pro_dataframe,phospro_dataframe))
  
  data_path <- paste(path,compare,'.rnk',sep = '')
  # write.table(compare_gsea[[compare]], file = data_path, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}

####################################  分子网络图数据准备 ###############################
ORA_result <- readRDS('./data/ORA_result.rds')

data1 <- ORA_result$F225_G225$up$GO
net1 <- data1[data1$Description %in% c('protein kinase A binding','structural constituent of cytoskeleton'),][,c('Description','geneID')]
data2 <- ORA_result$F225_G225$down$GO
net2 <- data2[data2$Description %in% c('cholesterol transfer activity','chemokine activity'),][,c('Description','geneID')]

data3 <- ORA_result$F225_G225$down$KEGG
net3 <- data3[data3$Description %in% c('Antigen processing and presentation'),][,c('Description','geneID')]
net3 <- separate_rows(net3,geneID,sep = '/')
gene <-  bitr(
  separate_rows(net3,geneID,sep = '/')$geneID,
  fromType = "ENTREZID",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)['SYMBOL']
net3 <- data.frame(Description='Antigen processing and presentation',geneID=gene)
colnames(net3) <- c('Description','geneID')

data4 <- ORA_result$down_up$rna$KEGG
net4 <- data4[data4$Description %in% c('Viral myocarditis'),][,c('Description','geneID')]
net4 <- separate_rows(net4,geneID,sep = '/')
gene <-  bitr(
  separate_rows(net4,geneID,sep = '/')$geneID,
  fromType = "ENTREZID",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)['SYMBOL']
net4 <- data.frame(Description='Viral myocarditis',geneID=gene)
colnames(net4) <- c('Description','geneID')

data5 <- ORA_result$down_up$phospro$KEGG
net5 <- data5[data5$Description %in% c('Allograft rejection'),][,c('Description','geneID')]
net5 <- separate_rows(net5,geneID,sep = '/')
gene <-  bitr(
  separate_rows(net5,geneID,sep = '/')$geneID,
  fromType = "ENTREZID",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)['SYMBOL']
net5 <- data.frame(Description='Allograft rejection',geneID=gene)
colnames(net5) <- c('Description','geneID')

data <- rbind(net1,net2)
netdata <- separate_rows(data,geneID,sep = '/')
netdata <- rbind(netdata,net3,net4,net5)

colnames(netdata) <- c('source','target')
# write.csv(netdata,file = './data/csv_for_cytoscape.csv',quote = F)


