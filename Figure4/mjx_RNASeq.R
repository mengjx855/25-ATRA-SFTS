#### Jinxin Meng, 20240614, 20251122 ####
setwd('F:/project/20240821_SFTSV_ATRA/git/Figure4/')
pacman::p_load(tidyverse, ggpubr, openxlsx)

#### data ####
group_color <- c('mock' = '#8da0cb', 'ATRA' = '#66c2a5', 
                 'SFTSV' = '#fc8d62', 'ATRA_SFTSV' = '#e78ac3')

group <- read.xlsx('../data/RNASeq_THP-1.xlsx', sheet = 'group')
profile <- read.xlsx('../data/RNASeq_THP-1.xlsx', sheet = 'profile') %>% 
  column_to_rownames('name')
profile <- profile[rowSums(profile > 2) > 3, ]

gene_info <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(profile), 
                                   keytype = 'ENSEMBL', columns = c('ENTREZID', 'SYMBOL')) 

#### Figure 4a ####
metadata <- data.frame(row.names = colnames(profile)) %>%
  mutate(group = group$group[match(rownames(.), group$sample)],
         group = factor(group))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = profile, colData = metadata, design = ~ group) # 构建dds对象
des <- DESeq2::DESeq(dds) # deseqs差异分析

comparisons <- list( # 构建分组比较对, 'group'是group表的第二列名称
  c('group', 'ATRA', 'mock'),
  c('group', 'SFTSV', 'mock'),
  c('group', 'ATRA_SFTSV', 'mock'),
  c('group', 'SFTSV', 'ATRA'),
  c('group', 'ATRA_SFTSV', 'ATRA'),
  c('group', 'ATRA_SFTSV', 'SFTSV')
)

difference <- map(
  comparisons, ~ 
    data.frame(DESeq2::results(des, contrast = .x)) %>%
    rownames_to_column('name') %>% 
    mutate(enriched = ifelse(log2FoldChange > 1 & padj < 0.05, .x[2], 
                             ifelse(log2FoldChange < -1 & padj < 0.05, .x[3], 'none')))
  ) %>% 
  set_names(map_vec(comparisons, ~ paste0(.x[2:3], collapse = '_vs_')))

saveRDS(difference, 'difference.rds')
write.xlsx(difference, 'difference.xlsx')

set.seed(2025)

genes <- c('CYP27A1','DBI','PPARG','FADS2','CPT1A','ACSL1','ACSL4','ACOX1','ANGPTL4','MMP1')

plot_data <- diff[['SFTSV_vs_mock']] %>% 
  left_join(gene_info, by = c('name' = 'ENSEMBL')) %>% 
  mutate(.padj = -log10(padj),
         .label = ifelse(SYMBOL %in% genes, SYMBOL, NA)) %>% 
  filter(abs(log2FoldChange) < 6, .padj < 190)

plot_label <- filter(plot_data, !is.na(.label) & enriched != 'none')

plot_data <- rbind(
  plot_label, 
  filter(plot_data, !name %in% plot_label$name) %>% 
    slice_sample(prop = .2))
  
ggscatter(plot_data, 'log2FoldChange', '.padj', color = 'enriched', 
          xlab = 'log2FoldChange', ylab = '-log10 FDR',
          palette = c(group_color, 'none' = 'grey'), size = 1) +
  ggrepel::geom_label_repel(aes(log2FoldChange, .padj, label = .label), plot_label, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  theme(aspect.ratio = 1)
ggsave('volcano.pdf', width = 6, height = 5)

#### Figure S4f ####
BiocParallel::register(BiocParallel::SerialParam())

gseGO <- map(
  difference, ~ 
    pull(.x, log2FoldChange, name = name) %>%
    sort(decreasing = T) %>% 
    clusterProfiler::gseGO(ont = 'ALL', OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                           keyType = 'ENSEMBL', pvalueCutoff = 1) %>% 
    data.frame())
write.xlsx(gseGO, 'gseGO.xlsx')

gseKEGG <- map(
  difference, ~
    left_join(.x, gene_info, by = c('name' = 'ENSEMBL')) %>%
    filter(!is.na(ENTREZID)) %>% 
    group_by(ENTREZID) %>% 
    slice_head() %>% 
    pull(log2FoldChange, name = ENTREZID) %>%
    sort(decreasing = T) %>% 
    clusterProfiler::gseKEGG(organism = 'hsa', pvalueCutoff = 1) %>% 
    data.frame())
write.xlsx(gseKEGG, 'gseKEGG.xlsx')

data <- read.xlsx('gsea_custom.xlsx') %>% 
  mutate(comparison = factor(comparison, c('ATRA_SFTSV_vs_SFTSV', 'SFTSV_vs_mock')))

level <- filter(data, comparison == 'ATRA_SFTSV_vs_SFTSV') %>% 
  arrange((NES)) %>% 
  pull(Description)

data <- mutate(data, Description = factor(Description, rev(level)))

ggdotchart(data, 'Description', 'NES', color = 'comparison', rotate = T,
          add = 'segments', add.params = list(color = 'comparison', size = 1), sorting = 'de',
          size = 4, xlab = '', ylab = 'NES-value', width = .7, palette = c('#66c2a5','#fc8d62')) +
  geom_hline(yintercept = 0, linetype = 'longdash')
ggsave('gsea.dotchart.pdf', width = 6.5, height = 5)

#### path展示 ####
library(pathview)
library(org.Hs.eg.db)

kegg_info <- download_KEGG('hsa')
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = c('ENSEMBL', 'ENTREZID'))
diff <- readRDS('diff.rds')[c('ATRA_vs_mock', 'SFTSV_vs_mock', 'ATRA_SFTSV_vs_SFTSV')]
terms <- c('MAPK signaling pathway', 'Retinol metabolism', 'NF-kappa B signaling pathway',
           'Toll-like receptor signaling pathway','B cell receptor signaling pathway',
           'JAK-STAT signaling pathway','T cell receptor signaling pathway',
           'NOD-like receptor signaling pathway', 'RIG-I-like receptor signaling pathway')
terms <- c('PPAR signaling pathway')

# fmt-1
map(terms, \(x)
    map(names(diff), \(y)
        {
          suffix <- paste0(gsub(' ', '_', x), '.', y)
          path_id <- kegg_info$KEGGPATHID2NAME$from[kegg_info$KEGGPATHID2NAME$to == x]
          data <- kegg_info$KEGGPATHID2EXTID %>% 
            filter(from == path_id) %>% 
            left_join(gene_info, by = c('to' = 'ENTREZID')) %>% 
            left_join(diff[[y]], by = c('ENSEMBL' = 'name')) %>% 
            filter(!is.na(pvalue)) %>% 
            filter(enriched != 'none') %>% 
            mutate(log2FoldChange = ifelse(log2FoldChange > 0, 1, -1))
          gene <- structure(data$log2FoldChange, names = data$to)
          try(pathview(gene.data = gene, pathway.id = path_id, species = 'hsa', 
                       out.suffix = paste0('fmt_1.', suffix), low = '#abd9e9', high = '#fdae61'))
    } ) )

# fmt-2
map(terms, \(x)
    map(names(diff), \(y)
        {
          suffix <- paste0(gsub(' ', '_', x), '.', y)
          path_id <- kegg_info$KEGGPATHID2NAME$from[kegg_info$KEGGPATHID2NAME$to == x]
          data <- kegg_info$KEGGPATHID2EXTID %>% 
            filter(from == path_id) %>% 
            left_join(gene_info, by = c('to' = 'ENTREZID')) %>% 
            left_join(diff[[y]], by = c('ENSEMBL' = 'name')) %>% 
            filter(!is.na(log2FoldChange))
          gene <- structure(data$log2FoldChange, names = data$to)
          try(pathview(gene.data = gene, pathway.id = path_id, species = 'hsa', 
                       out.suffix = paste0('fmt_2.', suffix), low = '#abd9e9', high = '#fdae61'))
    } ) )

#### MAPK gene heatmap ####
library(ComplexHeatmap)

fpkm <- read.delim('../data/profile_fpkm.txt')
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), 
                                   columns = c('ENSEMBL', 'ENTREZID', 'SYMBOL'))

# N01592 GF-RTK-RAS-ERK signaling pathway
# N00186 IL1-IL1R-p38 signaling pathway
# N00188 IL1-IL1R-JNK signaling pathway

genes <- read.delim('/database/KEGG/v20230401/network_hsa.txt', header = F, 
                    col.names = c('ENTREZID', 'NE')) %>% 
  mutate(ENTREZID = gsub('hsa:', '', ENTREZID),
         NE = gsub('ne:', '', NE)) %>% 
  filter(NE %in% c('N01592')) %>% 
  mutate(SYMBOL = gene_info$SYMBOL[match(ENTREZID, gene_info$ENTREZID)])

plot_data <- filter(fpkm, name %in% genes$SYMBOL) %>% 
  column_to_rownames('name')

col_data <- group %>% 
  mutate(group = factor(group, names(gp_col))) %>% 
  arrange(group)
col_annotation <- column_to_rownames(group, 'sample')
annotation_colors <- list(group = gp_col)
col_split <- factor(col_data$group)

pdf('MAPK_ab_N01592.heatmap.pdf', width = 6, height = 9)
pheatmap(plot_data, scale = 'row',
         cluster_cols = F,
         show_row_dend = F, show_column_dend = F,
         color = colorRampPalette(c('#307cc0', 'white', '#e43589'))(100),
         cellwidth = 14, cellheight = 8,
         fontsize_col = 12, fontsize_row = 8,
         border_color = 'white',
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         column_split = col_split,
         column_gap = unit(0, 'mm'),
         heatmap_legend_param = list(border = 'black', title = 'Scale FPKM'))
dev.off()


#### Custom heatmap ####
fpkm <- read.delim('../data/profile_fpkm.txt')
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), 
                                   columns = c('ENSEMBL', 'ENTREZID', 'SYMBOL'))

genes <- c('FOSB','FOSL2-AS1','FOS','FOSL1','FOSL2','JUNB','JUND','JUN',
           'MAPK8','MAPK9','MAPK10','MAPK11','MAPK12','MAPK13','MAP2K3','MAP2K7',
           'IKBKB','IKBKE','NLRP3','NFKB1','NFKB2','NFKBIZ','RELA','RELB','RIPK2',
           'IL6','IL10','TNF','NOD2','IL1A', 'IL1B','STAT3','TAB2', 'CASP1',
           'PYCARD','NLRP1','GSDMD')

plot_data <- filter(fpkm, name %in% genes) %>% 
  column_to_rownames('name')

col_data <- group %>% 
  mutate(group = factor(group, names(gp_col))) %>% 
  arrange(group)
col_annotation <- column_to_rownames(group, 'sample')
annotation_colors <- list(group = gp_col)
col_split <- factor(col_data$group)

pdf('MAPK_ab_Custom.heatmap.pdf', width = 6, height = 7)
pheatmap(plot_data, scale = 'row',
         cluster_cols = F,
         treeheight_row = 20,
         show_column_dend = F,
         color = colorRampPalette(c('#307cc0', 'white', '#e43589'))(100),
         cellwidth = 14, cellheight = 10,
         fontsize_col = 12, fontsize_row = 8,
         border_color = 'white',
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         column_split = col_split,
         column_gap = unit(0, 'mm'),
         heatmap_legend_param = list(border = 'black', title = 'Scale FPKM'))
dev.off()

