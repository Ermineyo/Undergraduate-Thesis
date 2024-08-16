# BiocManager::install('WebGestaltR')
# BiocManager::install('RDAVIDWebService ')
# BiocManager::install("enrichplot")
# BiocManager::install("VennDiagram")

library(WebGestaltR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(grid)
library(futile.logger)
library(VennDiagram)

edger_result_all <- readRDS('D:/zzy/graduation_thesis/data/edger_result_all.rds')

up_lists <- readRDS('D:/zzy/graduation_thesis/data/up_regulates_all.rds')
down_lists <- readRDS('D:/zzy/graduation_thesis/data/down_regulates_all.rds')












setwd('D:/zzy/graduation_thesis/enrichment_analysis')
# 富集分析
# GO
# obo:GO术语，gaf:从基因到GO术语的映射关系
go_obo_path <-  'go.obo.txt'
goa_human_gaf_path <- 'goa_human.gaf'
# KEGG
kegg_linksgenes_pathway_path <- 'linksgenes_pathway.list'
kegg_linksgenes_uniprot_path <- 'linksgenes_uniprot.list'
kegg_pathway_path <- 'pathway.list'


## 读取差异分析结果
setwd('D:/zzy/graduation_thesis/enrichment_analysis/new_data')
rna_folder <- "./rna"
pro_folder <- "./pro"
rna_files <- list.files(rna_folder, full.names = TRUE)
rna_data <- lapply(rna_files, readRDS)
pro_files <- list.files(pro_folder, full.names = TRUE)
pro_data <- lapply(pro_files, readRDS)
dif_ana_results <- list()
for (file in rna_files) {
  filename <- tools::file_path_sans_ext(basename(file))
  dif_ana_results[['rna']][[filename]] <- readRDS(file)
}

for (file in pro_files) {
  filename <- tools::file_path_sans_ext(basename(file))
  dif_ana_results[['pro']][[filename]] <- readRDS(file)
}
rm(file,filename,pro_files,rna_files,rna_folder,pro_folder,pro_data,rna_data)

# 构建GO分析数据库
# 构建数据库/背景文件
# 使用不同计算方法富集分析
## 绘图，使用差异表达基因
## 富集分析四种算法：ORA（david）/GSEA/通路拓扑(sp那个只有kegg)/xxx


# 准备WebGestaltR的数据
diff_gene_list <- function(diff_gene){
  path <- 'D:/zzy/graduation_thesis/enrichment_analysis/new_data/'
  result_list <- list()
  
  # lis <- strsplit(name,'_')
  # new_name <- paste(lis[[1]][1],lis[[1]][2],lis[[1]][3],sep='_')
  
  gene_diff <- diff_gene
  gene_diff <- gene_diff[order(gene_diff$FDR, gene_diff$logFC, decreasing = c(FALSE, TRUE)), ]
  gene_diff[which(gene_diff$logFC >= 2 & gene_diff$FDR < 0.05),'sig'] <- 'up'
  gene_diff[which(gene_diff$logFC <= -2 & gene_diff$FDR < 0.05),'sig'] <- 'down'
  gene_diff[which(abs(gene_diff$logFC) <= 1 | gene_diff$FDR >= 0.01),'sig'] <- 'none'
  gene_diff_select <- subset(gene_diff, sig %in% c('up', 'down'))
  # write.table(gene_diff_select, file = 'control_treat.glmQLFit.select.txt', sep = '\t', col.names = NA, quote = FALSE)
  
  #根据 up 和 down 分开输出
  gene_diff_up <- subset(gene_diff, sig == 'up')
  gene_diff_down <- subset(gene_diff, sig == 'down')
  # write.table(gene_diff_up, file = 'control_treat.glmQLFit.up.txt', sep = '\t', col.names = NA, quote = FALSE)
  # write.table(gene_diff_down, file = 'control_treat.glmQLFit.down.txt', sep = '\t', col.names = NA, quote = FALSE)
  
  ### KEGG/GO 富集分析
  # 输入：基因列表
  ## WebGestaltR
  up_gene_list <- rownames(gene_diff_up)
  down_gene_list <- rownames(gene_diff_down)
  gene_table_for_gsea <- data.frame(Gene =rownames(gene_diff),Score =gene_diff[['logFC']])
  
  result_list[['up_dataframe']] <- gene_diff_up
  result_list[['down_dataframe']] <- gene_diff_down
  result_list[['all_dataframe']] <- gene_diff
  
  return(result_list)
  # writeLines(rownames(diff_gene), paste(path,new_name,'_all_gene_list.txt',sep = '') )
  # writeLines(up_gene_list, paste(path,new_name,'_up_gene_list.txt',sep = '') )
  # writeLines(down_gene_list, paste(path,new_name,'_down_gene_list.txt',sep = '') )
  # write.table(gene_table_for_gsea, file = paste(path,new_name,'_for_gsea.rnk',sep = ''), sep = '\t', quote = FALSE,row.names = FALSE, col.names = FALSE)
}

list_for_enrichment_analysis <- list()

for(type in names(dif_ana_results)){
  tmp_list <- list()
  for(filename in names(dif_ana_results[[type]])){
    file <- dif_ana_results[[type]][[filename]]
    tmp_list[[filename]] <- diff_gene_list(file)
  }
  list_for_enrichment_analysis[[type]] <- tmp_list
}
names(list_for_enrichment_analysis$rna) <- c('F225_F5','F225_G225','F5_G225')
names(list_for_enrichment_analysis$pro) <- c('F225_F5','F225_G225','F5_G225')

# 富集分析
# ORA
enrich_analysis_ORA <- function(gene_data_frame, type){
  enrich_list <- list()
  
  if (type == 'rna') {
    genes <- rownames(gene_data_frame)
    keyType <- 'SYMBOL'
    gene_trans <- bitr(
      rownames(gene_data_frame),
      fromType = "SYMBOL",
      toType = c("SYMBOL", "ENTREZID"),
      OrgDb = org.Hs.eg.db
    )
    gene_trans <- gene_trans[!duplicated(gene_trans[, 'SYMBOL']), ]
    gene_id <- gene_trans[['ENTREZID']]
    
  } else if (type == 'pro') {
    genes <- rownames(gene_data_frame)
    keyType <- 'UNIPROT'
    gene_trans <-  bitr(
      rownames(gene_data_frame),
      fromType = "UNIPROT",
      toType = c("UNIPROT","SYMBOL", "ENTREZID"),
      OrgDb = org.Hs.eg.db
    )
    gene_trans <- gene_trans[!duplicated(gene_trans[, 'UNIPROT']), ]
    gene_id <- gene_trans[['ENTREZID']]
  }
  
  
  ## ORA
  GO_result <- enrichGO(
    genes,
    OrgDb = org.Hs.eg.db,
    keyType = keyType,
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
    gene_id,
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
# GSEA
enrich_analysis_GSEA <- function(gene_data_frame,type){
  enrich_list <- list()
  
  if(type=='rna'){
    gene_id_type <- 'SYMBOL'
  } else if(type=='pro'){
    gene_id_type <- 'UNIPROT'
    }
  
  gene_data_frame[[gene_id_type]] <- rownames(gene_data_frame)
  gene <- bitr(rownames(gene_data_frame),fromType=gene_id_type,toType="ENTREZID",OrgDb="org.Hs.eg.db")
  if(type=='rna'){
    gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
    gene_df <- data.frame(logFC=gene_data_frame$logFC, SYMBOL = gene_data_frame[[gene_id_type]]) 
  } else if(type=='pro'){
    gene <- dplyr::distinct(gene,UNIPROT,.keep_all=TRUE)
    gene_df <- data.frame(logFC=gene_data_frame$logFC, UNIPROT = gene_data_frame[[gene_id_type]]) 
  }
  gene_df <- merge(gene_df,gene,by=gene_id_type)
  geneList<-gene_df$logFC 
  names(geneList)=as.integer(gene_df$ENTREZID)  #使用转换好的ID
  geneList=sort(geneList,decreasing = T) #从高到低排序
  
  ## GSEA
  gseaGO_result <- gseGO(
    geneList,
    ont = "ALL",#可以替换为BP CC MF，分别富集
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
  
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

enrich_analysis_result <- list()
for(type in names(list_for_enrichment_analysis)){
  this_type_list <- list()
  for(condition in names(list_for_enrichment_analysis[[type]])){
    this_condition <- list_for_enrichment_analysis[[type]][[condition]]
    this_condition_list <- list()
    this_condition_list[['enrich_up']] <- enrich_analysis_ORA(gene_data_frame=this_condition$up_dataframe,type=type)
    this_condition_list[['enrich_down']] <- enrich_analysis_ORA(gene_data_frame=this_condition$down_dataframe,type=type)
    this_condition_list[['enrich_gsea']] <- enrich_analysis_GSEA(gene_data_frame=this_condition$all_dataframe,type=type)
    this_type_list[[condition]] <- this_condition_list
  }
  enrich_analysis_result[[type]] <- this_type_list
}
# saveRDS(enrich_analysis_result,file ='D:/zzy/graduation_thesis/enrichment_analysis/enrich_analysis_result.rds' )


### 看三种condition中上调下调的基因蛋白质富集通路的交集
enrich_analysis_result <- readRDS('D:/zzy/graduation_thesis/enrichment_analysis/enrich_analysis_result.rds')

integrate_pathway_list <- list()
for(condition in names(enrich_analysis_result$rna)){
  condition_up <- c()
  condition_down <- c()
  for(type in c('rna','pro')){
    up_GO <- enrich_analysis_result[[type]][[condition]]$enrich_up$GO$ID
    up_KEGG <- enrich_analysis_result[[type]][[condition]]$enrich_up$KEGG$ID
    down_GO <- enrich_analysis_result[[type]][[condition]]$enrich_down$GO$ID
    down_KEGG <- enrich_analysis_result[[type]][[condition]]$enrich_down$KEGG$ID
    tmp_up <- append(up_GO,up_KEGG)
    tmp_down <- append(down_GO, down_KEGG)
    condition_up <- append(condition_up,tmp_up)
    condition_down <- append(condition_down,tmp_down)
  }
  
  integrate_pathway_list[[condition]][['up']] <- unique(condition_up)
  integrate_pathway_list[[condition]][['down']] <- unique(condition_down)
}

F225_F5_up <- integrate_pathway_list$F225_F5$up
F225_G225_up <- integrate_pathway_list$F225_G225$up
F5_G225_up <- integrate_pathway_list$F5_G225$up


F225_F5_down <- integrate_pathway_list$F225_F5$down
F225_G225_down <- integrate_pathway_list$F225_G225$down
F5_G225_down <- integrate_pathway_list$F5_G225$down


venn_list_up <- list(F225_F5_up = F225_F5_up, F225_G225_up = F225_G225_up, F5_G225_up = F5_G225_up)
venn_list_down <- list(F225_F5_down = F225_F5_down, F225_G225_down = F225_G225_down, F5_G225_down = F5_G225_down)

setwd('D:/zzy/graduation_thesis/enrichment_analysis/')
up_path <- 'up.png'
down_path <- 'down.png'


venn.diagram(venn_list_up, filename = up_path, imagetype = 'png', 
             fill = c('red', 'blue', 'green'), alpha = 0.50, 
             cat.col = c('red', 'blue', 'green'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('red', 'blue', 'green'), cex = 1.5, fontfamily = 'serif')

venn.diagram(venn_list_down, filename = down_path, imagetype = 'png', 
             fill = c('red', 'blue', 'green'), alpha = 0.50, 
             cat.col = c('red', 'blue', 'green'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('red', 'blue', 'green'), cex = 1.5, fontfamily = 'serif')


# 结果保存到表格
inter_up <- get.venn.partitions(venn_list_up)
up_inter_pathways <- as.character(unlist(inter_up[c(1,2,3,5),]$..values..))
# writeLines(up_inter_pathways,file=)

inter_down <- get.venn.partitions(venn_list_down)
down_inter_pathways <- as.character(unlist(inter_down[c(1,2,3,5),]$..values..))

# write.table(inter[-c(5, 6)], 'venn_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)


########################################
# 
WebGestaltR(
  enrichMethod = "NTA",
  organism = "hsapiens",
  enrichDatabase = 'pathway_KEGG',
  interestGene = genes,
  interestGeneType = 'genesymbol',
  collapseMethod = "mean",
)



ridgeplot(Go_gseresult,10) #输出前十个结果







PPI <- getPPI(
  ID = 1,
  taxID = "auto",
  required_score = NULL,
  network_type = "functional",
  add_nodes = 0,
  show_query_node_labels = 0,
  output = "igraph"
)





GO_result_table <- setReadable(GO_result, OrgDb = org.Hs.eg.db)
# write.csv(summary(ego_ALL),"ALL-enrich.csv",row.names =FALSE)

##可视化--点图
dot <- dotplot(GO_result,title="EnrichmentGO_MF_dot")#点图，按富集的数从大到小的
##可视化--条形图
bar <- barplot(KEGG_result, showCategory=20,title="EnrichmentGO_MF")#条状图，按p从小到大排，绘制前20个Term






## convert gene ID to Symbol
edox <- setReadable(GO_result, 'org.Hs.eg.db', 'UNIPROT')
p1 <- cnetplot(edox, foldChange=test_genes)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=test_genes, max.overlaps=100000)
p3 <- cnetplot(edox, foldChange=test_genes, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))


edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')












