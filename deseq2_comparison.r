suppressPackageStartupMessages(library("DESeq2"))               
suppressPackageStartupMessages(library("BiocParallel"))         
suppressPackageStartupMessages(library("gprofiler2"))
suppressPackageStartupMessages(library("rlist"))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(DOSE))
register(MulticoreParam(48))                                    

print("Loaded Libraries")

workingDirectory <- "" #CHANGE THIS VARIABLE TO THE DESIRED WORKING DIRECTORY

### Read in count matrix from featureCounts output without column 1 to 5
countData <- read.table(paste0(workingDirectory, "counts/all_samples_counts.txt"), header = TRUE, row.names = 1, stringsAsFactors = FALSE)[, -c(1:5)]

### Rename columns to "sample_1", "sample_2", ..., "sample_43"
colnames(countData) <- paste0("sample_", 1:43)

### Read in sample information
sampleTable <- read.table(paste0(workingDirectory, "sample_information.csv"), header = TRUE, sep=";")

### Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleTable, design = ~ Group)

print("Running dds file")

### Run DESeq2
dds <- DESeq(dds, quiet = TRUE, parallel=TRUE, BPPARAM=MulticoreParam(48))

### Write normalized counts to file
write.csv(as.data.frame(counts(dds, normalized=TRUE)), file=paste0(workingDirectory, "counts/all_samples_normalized_counts.csv"))

## Create folders for pairwise comparisons of treatments and write results to file
create folder for each comparison
for (i in 1:15) {
  for (j in (i+1):16) {
    dir.create(paste0(workingDirectory, "deseq2/", levels(dds$Group)[i], "_vs_", levels(dds$Group)[j]))
  }
}

print("start writing results")

# write results to file including MA plot for each comparison
for (i in 1:15) {
  for (j in (i+1):16) {
    res <- results(dds, parallel=TRUE, BPPARAM=MulticoreParam(48), contrast = c("Group", levels(dds$Group)[i], levels(dds$Group)[j]))
    resOrdered <- res[order(res$padj),]
    write.csv(as.data.frame(resOrdered), file=paste0(workingDirectory, "deseq2/", levels(dds$Group)[i], "_vs_", levels(dds$Group)[j], "/result.csv"))
    pdf(paste0(workingDirectory, "deseq2/", levels(dds$Group)[i], "_vs_", levels(dds$Group)[j], "/MA_plot.pdf"))
    plotMA(res, ylim=c(-2,2))
    dev.off()
    res_sig = subset(resOrdered, padj < 0.05)
    res_sig_up = subset(res_sig, log2FoldChange > 0)
    res_sig_down = subset(res_sig, log2FoldChange < 0)
    goc = cbind("GO:BP", "GO:MF", "GO:CC")
    for (category in goc)
    {
      multi_gp = gost(list("up-regulated" = row.names(res_sig_up), "down-regulated" = row.names(res_sig_down)), organism = "athaliana" , multi_query = FALSE, evcodes = TRUE, sources = category)
      if (length(multi_gp$result) > 0) 
      {
        gp_mod = multi_gp$result[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
        gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
        gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
        names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
        gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
        write.csv(as.data.frame(gp_mod), file=paste0(workingDirectory, "deseq2/", levels(dds$Group)[i], "_vs_", levels(dds$Group)[j], paste0("/",gsub(":","_",category),"_enrichment.csv")))
        gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)
        pdf(paste0(workingDirectory, "deseq2/", levels(dds$Group)[i], "_vs_", levels(dds$Group)[j], paste0("/",gsub(":","_",category),"_enrichment.pdf")))
        print(enrichplot::dotplot(gp_mod_cluster, x = "Cluster", showCategory = 15, by = "GeneRatio", font.size = 8, label_format = 50, title = paste0("Top 15 enriched ", category," terms down-/up-regulated genes")))
        dev.off()
      }
    }
    print(paste0(levels(dds$Group)[i], "_vs_", levels(dds$Group)[j]))
  }
}

# goc = cbind("GO:BP", "GO:MF", "GO:CC")
# for (category in goc)
# {
#   multi_gp = gost(list("up-regulated" = row.names(res_sig_up), "down-regulated" = row.names(res_sig_down)), organism = "athaliana" , multi_query = FALSE, evcodes = TRUE, sources = category)
#   gp_mod = multi_gp$result[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
#   gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
#   gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
#   names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", "query_size", "Count", "term_size", "effective_domain_size", "geneID", "GeneRatio", "BgRatio")
#   gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
#   write.csv(as.data.frame(gp_mod), file=paste0(workingDirectory, "deseq2/", levels(dds$Group)[i], "_vs_", levels(dds$Group)[j], paste0("/",gsub(":","_",category),"_enrichment.csv")))
#   gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)
#   pdf(paste0(workingDirectory, "deseq2/", levels(dds$Group)[i], "_vs_", levels(dds$Group)[j], paste0("/",gsub(":","_",category),"_enrichment.pdf")))
#   enrichplot::dotplot(gp_mod_cluster, x = "Cluster", showCategory = 15, by = "GeneRatio", font.size = 8, label_format = 50, title = paste0("Top 15 enriched ", gsub(":","_",category)," terms down-/up-regulated genes"))
#   dev.off()
# }


# # gp_up = gost(row.names(res_sig_up), organism = "athaliana")
# # gp_up$result <- list.remove(gp_up$result, "parents")
# # gp_down = gost(row.names(res_sig_down), organism = "athaliana") 
# # gp_down$result <- list.remove(gp_down$result, "parents")
# # write.csv(as.data.frame(gp_up$result), file=paste0(workingDirectory, "deseq2/", levels(dds$Group)[i], "_vs_", levels(dds$Group)[j], "/go_enrichment_upregulated_genes.csv"))
# # write.csv(as.data.frame(gp_down$result), file=paste0(workingDirectory, "deseq2/", levels(dds$Group)[i], "_vs_", levels(dds$Group)[j], "/go_enrichment_downregulated_genes.csv"))
# # print(i)
# # print(j)

#       write.csv(as.data.frame(gp_mod), file="/data/fass2/bic_projects/SE20231212_167/deseq2/CLLF_1_vs_CLLF_4/GO_BP_enrichment.csv")
# gp_mod <- read.csv("/data/fass2/bic_projects/SE20231212_167/deseq2/CLLF_1_vs_CLLF_4/GO_BP_enrichment.csv")
# gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)


# list_of_files <- list.files(paste0(workingDirectory,"deseq2"), pattern = "GO_BP_enrichment.csv", recursive = T)
# for (file in list_of_files)
# {
#   gp_mod <- read.csv(paste0(workingDirectory,"deseq2/",file))
#   gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)
#   pdf(paste0(workingDirectory, "deseq2/", gsub("/GO_BP_enrichment.csv","",file),"/GO_BP_enrichment.pdf"))
#   print(enrichplot::dotplot(gp_mod_cluster, x = "Cluster", showCategory = 15, by = "GeneRatio", font.size = 8, label_format = 50, title = "Top 15 enriched GO:BP terms down-/up-regulated genes"))
#   dev.off()
# }