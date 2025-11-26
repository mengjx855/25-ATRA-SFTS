#### Jinxin Meng, 20241113, 20251122 ####
setwd('F:/project/20240821_SFTSV_ATRA/git/Figure1/')
pacman::p_load(tidyverse, ggpubr, openxlsx)
source('../scripts/calcu_difference.R')

#### info ####
group_color <- c('SF' = '#fc8d62', 'NonSF' = '#66c2a5')

profile <- read.xlsx('../data/MBX_data.xlsx', sheet = 'profile', check.names = F) %>% 
  column_to_rownames('name') %>% 
  na.omit() %>% 
  apply(2, \(x) x / sum(x) * 100) %>% 
  data.frame(check.names = F)

metabolite_info <- read.xlsx('../data/MBX_data.xlsx', sheet = 'info')

group <- read.xlsx('../data/MBX_data.xlsx', sheet = 'group') %>% 
  select(sample, group = group2) %>% 
  mutate(group = factor(group, names(group_color)))

#### Figure S1a ####
plsda <- ropls::opls(x = data.frame(t(profile), check.names = F), orthoI = 0, predI = 2)

pca_point <- data.frame(plsda@scoreMN) %>% 
  rownames_to_column('sample') %>% 
  left_join(group, by = 'sample')
pca_label <- paste0(c('PC1', 'PC2'), ' (', plsda@modelDF$R2X * 100, '%)')

ggscatter(pca_point, 'p1', 'p2', color = 'group', palette = group_color, size = 1,
          xlab = pca_label[1], ylab = pca_label[2]) +
  stat_ellipse(aes(color = group, fill = group), geom = 'polygon', level = 0.9, 
               alpha = .05, lty = 2, lwd = .3, show.legend = F) +
  theme(aspect.ratio = .9)
ggsave('PLS.pdf', width = 5, height = 5)

#### Figure S1b ####
oplsda <- ropls::opls(x = t(profile), pull(group, group),
                      orthoI = NA, predI = 1)
  
difference <- difference_analysis(profile, group, comparison = c('SF', 'NonSF'), 
                                  method = 'wilcox')
  
result <- data.frame(vip = oplsda@vipVn) %>% 
  rownames_to_column('cpd_id') %>% 
  left_join(difference, ., by = c('name' = 'cpd_id')) %>% 
  left_join(metabolite_info, ., by = c('cpd_id' = 'name'))

write.xlsx(result, 'difference.xlsx')

plot_data <- result %>% 
  filter(!is.na(log2FC)) %>% 
  select(cpd_id, log2FC, padj, vip) %>% 
  mutate(.padj = -log10(padj),
         sig = ifelse(padj < 0.05 | vip > 1, 'yes', 'no'),
         enriched = ifelse(sig == 'yes' & log2FC > 0, 'SF', 
                           ifelse(sig == 'yes' & log2FC < 0, 'NonSF', 'none'))) %>% 
  within(., {
    .padj[.padj > 25] = 25
    log2FC[log2FC > 2] = 2
    log2FC[log2FC < -2] = -2
  } )

ggscatter(plot_data, 'log2FC', '.padj', color = 'enriched', 
          xlab = 'log2FoldChange', ylab = '-log10 PDR',
          palette = c(group_color, 'none' = 'grey'), size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  theme(aspect.ratio = 1)
ggsave('volcano.SF_vs_NonSF.pdf', width = 3.6, height = 4)

#### Figure 1a ####
bg_data <- read.delim('../data/cpd2path_enrichment.tsv')
kegg_info <- read.delim('../data/path_description')

difference <- read.xlsx('difference.xlsx') %>% 
  filter(FC > 0.5 & (vip > 2 | padj < 0.001)) %>% 
  filter(KEGG_ID != '')

eKEGG <- clusterProfiler::enricher(
  gene = unique(difference$KEGG), TERM2GENE = bg_data, minGSSize = 1, 
  pvalueCutoff = 1, qvalueCutoff = 1) %>% 
  data.frame()
write.xlsx(eKEGG, 'eKEGG.xlsx')

plot_data <- eKEGG %>% 
  mutate(name = stringr::str_extract(ID, 'map\\d+'), .before = 1) %>% 
  left_join(select(kegg_info, name, category), by = 'name') %>% 
  filter(category == 'Metabolism' & 
           name %in% c('map00310','map00380','map00360','map00300','map01200','map00650','map00830','map00220')) %>% 
  select(name = Description, GeneRatio, pvalue) %>% 
  mutate(GeneRatio = as.numeric(stringr::str_split_i(GeneRatio, '/', 1)) /
           as.numeric(stringr::str_split_i(GeneRatio, '/', 2))) %>% 
  mutate(pvalue = round(pvalue, 3),
         GeneRatio = GeneRatio + rnorm(nrow(.), mean = 0, sd = .005)) %>% 
  arrange(desc(GeneRatio))

ggscatter(plot_data, 'name', 'GeneRatio', fill = 'pvalue', rotate = T, size = 5,
          shape = 21, legend = 'right', xlab = '', ylab = 'Gene Ratio') +
  scale_fill_viridis_c(begin = .6)
ggsave('eKEGG.pdf', width = 5.6, height = 4)
