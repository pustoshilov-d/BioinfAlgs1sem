if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")
install.packages("plyr")
install.packages("dplyr")
install.packages("ggplot2")
BiocManager::install("biomaRt")
BiocManager::install("GO.db")
BiocManager::install("org.Hs.eg.db")

library("org.Hs.eg.db")
library("GO.db")
library("ggplot2")
library("limma")
library("biomaRt")
library("plyr")
library("edgeR")
library("dplyr")

#read count table

origin_table <- read.table("data/E-GEOD-47718-raw-counts.tsv", sep="\t", header=TRUE)
design_table <- read.table("data/E-GEOD-47718-experiment-design.tsv", sep="\t", header=TRUE)

origin_table
design_table

counts_table <- origin_table
col_for_research = "Sample.Characteristic.clinical.information."



design_table$short_name <- rep(NaN, nrow(design_table))
length(colnames(counts_table))
colnames(design_table)




diseases_list <- rep(0, length(unique(design_table[[col_for_research]])))
diseases_list

names(diseases_list) <- unique(design_table[[col_for_research]])
names(diseases_list)



cols_names <- colnames(counts_table)
new_names <- c()

for (col in cols_names){
    if (length(design_table[design_table$Run == col,1]) != 0) {
      disease <- design_table[design_table$Run == col, col_for_research]
      diseases_list[disease] <- diseases_list[disease] + 1
      new_name <- paste(toupper(substr(disease,1,1)), toString(diseases_list[disease]),sep = ".")
      new_names <-  append(new_names, new_name)
      design_table[design_table$Run == col, "short_name"] <- new_name
    } else {
      new_names <-  append(new_names, col)
    }
}

new_names
diseases_list
names(counts_table) <- new_names

rownames(counts_table) <- counts_table$Gene.ID

counts_table <- subset(counts_table, select = c(-Gene.ID, -Gene.Name))
counts_table <- counts_table[ order(rownames(counts_table)), ]
counts_table



diseases_short <- unique(substring(colnames(counts_table),1,1))
diseases_short

#"S" "C" "O" "P" "U" "N"

what_to_compare <- c("N", "S")


#  



cur_cols <- colnames(counts_table)[(substring(colnames(counts_table),1,1) == what_to_compare[1] 
                                    | substring(colnames(counts_table),1,1) == what_to_compare[2])]


cur_table <- subset(counts_table, select = c(cur_cols))
cur_table


cur_lables <- substring(colnames(cur_table),1,1)
cur_lables

#
pca <- prcomp(t(cur_table))
pc <- round(pca$sdev^2/sum(pca$sdev^2)*100,1)

barplot(pc, xlab = 'PC number', ylab = '% of cariation PC accounts for')

#
pca.data <- data.frame(Type=cur_lables, PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3])
ggplot(pca.data, aes(x=PC1, y=PC2)) +
  geom_point(aes(col=Type)) +
  labs(x="PC1", y="PC2", col="Type") + 
  theme(legend.position = "bottom")

#
dist_data <- dist(t(cur_table), method = "euclidean")
clust1 <-hclust(dist_data)
plot(clust1)

#
plotMDS(cur_table)
#plotMDS(counts_table)


#
#suppressPackageStartupMessages(library(org.Hs.eg.db))

mart <- useMart("ENSEMBL_MART_ENSEMBL")
#listDatasets(mart = mart)

mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("ensembl_gene_id", "description", "chromosome_name",
                "start_position", "end_position", "gene_biotype",
                "strand", "hgnc_symbol", "entrezgene_id")

filters <- "ensembl_gene_id"
genes_annotation <- getBM(attributes = attributes, filters = filters, values = rownames(cur_table), mart = mart)
View(genes_annotation)

#table(rownames(cur_table) == genes_annotation$ensembl_gene_id)



paste("Intersect: ", toString(length(intersect(rownames(cur_table), genes_annotation$ensembl_gene_id))))
paste("In annotation: ", toString(length(genes_annotation$ensembl_gene_id)))
paste("In data: ", toString(length(rownames(cur_table))))

#row.names(cur_table)[!(row.names(cur_table) %in% genes_annotation$ensembl_gene_id)]
#getBM(attributes = attributes, filters = filters, values = "ENSG00000268568", mart = mart)
# 


cur_table_cut <- subset(cur_table, row.names(cur_table) %in% genes_annotation$ensembl_gene_id)
annotation_cut <- subset(genes_annotation, genes_annotation$ensembl_gene_id %in% row.names(cur_table_cut))
paste("In annotation after cut: ", toString(length(annotation_cut$ensembl_gene_id)))

annotation_cut <- annotation_cut %>% distinct(ensembl_gene_id, .keep_all = TRUE)
paste("In annotation after drop dupl: ", toString(length(annotation_cut$ensembl_gene_id)))
paste("In data after cut: ", toString(length(rownames(cur_table_cut))))



# анализ диф.экспрессии
# лист
DGE <- DGEList(counts = cur_table, group = cur_lables)

# нормализация
DGE <- calcNormFactors(DGE, method = "TMM")

# дисперсия
DGE <- estimateCommonDisp(DGE)
DGE <- estimateTagwiseDisp(DGE)

# тест на диф.экспрессию
DE_res <- exactTest(DGE)

topTags(DE_res)
DE_res_top_table <- topTags(DE_res)$table
View(DE_res_top_table)



DGE_cut <- DGEList(counts = cur_table_cut, genes = annotation_cut, group = cur_lables)
DGE_cut <- calcNormFactors(DGE_cut, method = "TMM")
DGE_cut <- estimateCommonDisp(DGE_cut)
DGE_cut <- estimateTagwiseDisp(DGE_cut)
DE_res_cut <- exactTest(DGE_cut)

topTags(DE_res_cut)
DE_res_top_table_cut <- topTags(DE_res_cut)$table
View(DE_res_top_table_cut)

rownames(DE_res_top_table_cut) == rownames(DE_res_top_table)

DE_res_cut


# подтягиваем генную онтологию 
go <- goana.DGEExact(de=DE_res_cut, geneid = DE_res_cut$genes$entrezgene_id, species="Hs")
topGO(go, number = 10)

# подтягиваем метаболические пути
kegg <- kegga.DGEExact(de=DE_res_cut, geneid = DE_res_cut$genes$entrezgene_id, species="Hs")
topKEGG(kegg, number = 10)
