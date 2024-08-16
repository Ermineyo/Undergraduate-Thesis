test <- rbind(ORA_result$F225_G225$up$KEGG[c('Description','p.adjust')],
              ORA_result$F225_G225$up$GO[c('Description','p.adjust')])
test <- test[order(-test$p.adjust),]
test <- test[1:10,]

# 绘制横过来的 barplot
barplot(test$p.adjust, names.arg = test$Description, horiz = TRUE, las = 1, col = "skyblue", 
        main = "Top 10 Enriched Pathways", xlab = "Adjusted P-value", ylab = "Pathway")








## 

##  分析：看上调的通路，作图
## 下一步：看整合的情况，两个生长过程中都上升的通路，

####################################  分子网络绘制   #####################################################


string_db <- STRINGdb$new( version="11.5", #数据库版本。截止2022.5.24最新为11.5
                           species=9606,   #人9606，小鼠10090 
                           score_threshold=700)  #蛋白互作的得分 默认400, 低150，高700，极高900

log2FC_cutoff = log2(2)
pvalue_cutoff = 0.05
padj_cutoff = 0.05



test <- edger_result_all$rna$F225_G225

test$Gene <- rownames(test)

# 选择需要的列，并按照要求的格式重命名列名
test <- test[, c("Gene", "logFC")]

# 将数据框写入文件，以制表符分隔
write.table(test, file = "output.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



need_deg <- test[c('logFC','PValue','FDR')]

colnames(need_deg) <- c('log2FC','pvalue','padj')
need_deg$gene <- rownames(need_deg)
if(T){  
  gene_up=need_deg[with(need_deg,log2FC>log2FC_cutoff & pvalue<pvalue_cutoff & padj<padj_cutoff),]
  gene_down=need_deg[with(need_deg,log2FC < -log2FC_cutoff & pvalue<pvalue_cutoff & padj<padj_cutoff),]
  gene_diff=need_deg[with(need_deg,abs(log2FC)>log2FC_cutoff & pvalue<pvalue_cutoff & padj<padj_cutoff),]
}
dat <- gene_diff ##这里选取前100显著基因用于后续分析
# write.table(rownames(dat),'gene_diff100.txt',row.names = F,col.names = F,quote = F) #字符不要带引号 

dat_map <- string_db$map(my_data_frame=dat, 
                         my_data_frame_id_col_names="gene", #使用gene symbol或ENTREZID都可
                         removeUnmappedRows = TRUE )

hits <- dat_map$STRING_id 
data_links <- data_mapped$STRING_id[1:100] %>% string_db$get_interactions()

## PPI
png("string_PPI.png",units="in",width = 10,height = 10, res=400)
string_db$plot_network(hits)
dev.off()
## PPI_halo  #给PPI添加上下调信息
# filter by p-value and add a color column(i.e.green for down and red for up genes)
dat_map_color <- string_db$add_diff_exp_color(subset(dat_map, pvalue<0.01),
                                              logFcColStr="log2FC" )
payload_id <- string_db$post_payload(dat_map_color$STRING_id,
                                     colors=dat_map_color$color)
png("string_PPI_halo.png",units="in",width = 10,height = 10, res=400)
string_db$plot_network(hits, payload_id=payload_id )
dev.off()



up_lists <- readRDS('./data/up_regulates_all.rds')
down_lists <- readRDS('./data/down_regulates_all.rds')


PPI <- getPPI(
  ID = 1,
  taxID = "auto",
  required_score = NULL,
  network_type = "functional",
  add_nodes = 0,
  show_query_node_labels = 0,
  output = "igraph"
)

up_test <- up_lists$rna$F225_G225
test_data <- data.frame(node1 = up_lists$rna$F225_G225,
                        log2FC = edger_result_all$rna$F225_G225[edger_result_all$rna$F225_G225$regulated == 'up',]$logFC)

write.csv(test_data,'node_data.csv',row.names = F,quote = F) #字符不要带引号 





### 转换名称
writeLines(rownames(pro_data),'./data/all_protein_ids.txt')


kegg <- read.csv('./enrichment_analysis/linksgenes_uniprot.list',sep='\t',header = FALSE)
# 假设您有一个包含 UniProt 名称列表的数据框，名为 uniprot_list_df
uniprot_list <- rownames(pro_data)
# 创建一个空的数据框来存储结果
result_df <- data.frame(uniprot = uniprot_list, gene = NA)
# 创建保存未匹配
protein_id_not_matched <- list()
# 遍历 UniProt 名称列表，查找对应的基因并保存到结果数据框中
for (i in 1:length(uniprot_list)) {
  uniprot <- uniprot_list[i]
  # 在 kegg 数据中查找与当前 UniProt 名称匹配的行
  matching_row <- kegg[kegg$V2 == uniprot, ]
  if (nrow(matching_row) > 0) {
    # 如果找到匹配的行，则将对应的基因保存到结果数据框中
    result_df[i, "gene"] <- matching_row$V1
  }
  else{
    protein_id_not_matched <- append(protein_id_not_matched,uniprot)
  }
}

# 查看结果数据框
print(result_df)




## clusterprofiler 名称转换
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

all_protein_list <- readLines('./data/all_protein_ids.txt')

gene_UNIPROT <- bitr(
  all_protein_list,
  fromType = "UNIPROT",
  toType = c("UNIPROT", "SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)

diff <- setdiff(all_protein_list, gene_UNIPROT[['UNIPROT']])
# writeLines(diff,'./data/not_contain_ids.txt' )

idmapping <- read.table('./data/idmapping_2024_04_23.tsv', header = TRUE, sep = "\t")

em <- idmapping[['to']]
gene_UNIPROT <- bitr(
  em,
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)















setwd('.')
pro_file = "./protein_data_reviewed/human_pro/combined/txt/proteinGroups.txt"

full_data <- read.delim(pro_file,stringsAsFactors = FALSE)
colnames(full_data)

cols <- c('Reporter.intensity.1','Reporter.intensity.2','Reporter.intensity.3','Reporter.intensity.4','Reporter.intensity.5','Reporter.intensity.6')
pro_data <- full_data[cols]
rownames(pro_data) <- full_data$Gene.names
pro_data_singleProteinID <- matrix( NA,ncol=ncol(pro_data) )
colnames(pro_data_singleProteinID) <- colnames(pro_data)
all_data_list <- list()
for(i in 1:nrow(pro_data)){
  names <- full_data$Gene.names[i]
  name_list <- strsplit(names, ';')[[1]]
  for (name in name_list) {
    all_data_list[[name]] <- rbind(all_data_list[[name]], pro_data[i, ])
  }
}
pro_data_singleProteinID <- do.call(rbind, all_data_list)
rownames(pro_data_singleProteinID) <- names(all_data_list)


ids <- c("Q92928", "Q5VTU8", "Q8NHW5", "Q15695", "A6NIZ1", "Q5T1J5", "Q63ZY6", "Q58FG0", "Q6ZUV0", "A0A0B4J2D5",
         "Q6DN03", "Q58FF8", "O60361", "O75343", "Q8NFI4", "Q6DRA6", "Q9NQ39", "A6NEC2", "Q9BYX7", "P01893",
         "Q9NQA3", "A8MWD9", "P0C7P4", "Q309B1", "Q9BZK3", "A6NF01", "Q14568", "Q5VTE0", "Q6ZSR9", "A8MWX3",
         "Q56UQ5", "Q9BZ68", "A4QPH2", "Q6VEQ5", "A6NMY6", "Q75LS8", "Q6NVV1", "Q8N8J0", "Q5JNZ5", "P0DSN6",
         "A6NKH3", "C9JQL5", "Q58FG1", "Q6NSW5", "Q6A1A2", "Q58FF3", "A8MUU1", "C4AMC7", "P48741", "Q58FF6",
         "P59074", "Q6ZT62", "Q5T035", "P0CG22", "Q6P474", "Q9NSQ0", "O94854", "A6NGW2", "Q8IVE0", "Q9UKY3",
         "Q8IZP2", "O43423")






setwd('.')
pro_file = "./protein_data_reviewed/human_pro/combined/txt/proteinGroups.txt"

full_data <- read.delim(pro_file,stringsAsFactors = FALSE)
colnames(full_data)

cols <- c('Reporter.intensity.1','Reporter.intensity.2','Reporter.intensity.3','Reporter.intensity.4','Reporter.intensity.5','Reporter.intensity.6')
pro_data <- full_data[cols]
gene_names <- full_data$Gene.names
gene_name_list <- c()
for(x in gene_names){
  this_names <- strsplit(x, ';')[[1]]
  for (name in this_names) {
    gene_name_list <- append(gene_name_list,name)
  }
}

gene_UNIPROT <- bitr(
  gene_name_list,
  fromType = "SYMBOL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)

diff_genes <- setdiff(gene_name_list,gene_UNIPROT[['SYMBOL']])
length(diff_genes)
