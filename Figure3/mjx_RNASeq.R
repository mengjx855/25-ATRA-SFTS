#### Jinxin Meng, 20240614, 20251122 ####
setwd('F:/project/20240821_SFTSV_ATRA/git/Figure3/')
pacman::p_load(tidyverse, ggpubr, openxlsx)

#### Figure 3e ####
group_color <- c('mock' = '#8da0cb', 'ATRA' = '#66c2a5', 
                 'SFTSV' = '#fc8d62', 'ATRA_SFTSV' = '#e78ac3')

group <- read.xlsx('../data/RNASeq_THP-1.xlsx', sheet = 'group')
fpkm <- read.xlsx('../data/RNASeq_THP-1.xlsx', sheet = 'fpkm')

genes <- c('FOSB','FOSL2-AS1','FOS','FOSL1','FOSL2','JUNB','JUND','JUN',
           'MAPK8','MAPK9','MAPK10','MAPK11','MAPK12','MAPK13','MAP2K3','MAP2K7',
           'IKBKB','IKBKE','NLRP3','NFKB1','NFKB2','NFKBIZ','RELA','RELB','RIPK2',
           'IL6','IL10','TNF','NOD2','IL1A', 'IL1B','STAT3','TAB2', 'CASP1',
           'PYCARD','NLRP1','GSDMD')

plot_data <- filter(fpkm, name %in% genes) %>% 
  column_to_rownames('name')

col_data <- group %>% 
  mutate(group = factor(group, names(group_color))) %>% 
  arrange(group)
col_annotation <- column_to_rownames(group, 'sample')
annotation_colors <- list(group = group_color)
col_split <- factor(col_data$group)

pdf('MAPK_ab_Custom.heatmap.pdf', width = 6, height = 7)
pheatmap(plot_data, scale = 'row',
         cluster_cols = F,
         treeheight_row = 20,
         show_column_dend = F,
         color = colorRampPalette(c('#307cc0', 'white', '#e43589'))(100),
         cellwidth = 14, cellheight = 10,
         show_colnames = F, fontsize_row = 8,
         fontface_row = 'italic',
         border_color = 'white',
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         column_split = col_split,
         column_gap = unit(0, 'mm'),
         heatmap_legend_param = list(border = 'black', title = 'Scale FPKM'))
dev.off()
