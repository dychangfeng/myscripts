rna.tissue= read.delim('~/Downloads/rna_tissue.tsv', header = T)
library(dplyr)
library(tidyr)
rna.df = rna.tissue%>%select(-Unit)%>%
  group_by(Gene.name)%>%spread(key = Sample, value = Value )%>%ungroup()
write.table(rna.df, '~/Desktop/rna_tissue.txt', sep = '\t', row.names = F)
##-----------------------------------------------------------------------------
GSE45332_RNA_norm = read.delim('Downloads/GSE45332_GEO_RNA_normalization.txt', header = T)
GSE45332_transcripts = read.delim('Downloads/GSE45332_GEO_transcription_data.txt', header = T)
##--------------------------------------------------------------------
RNA.cellline = read.delim('~/Desktop/RNA_cell_line/rna_celline.tsv', header = T)
##-------------------------
colon.cancer = read.delim('~/Downloads/GSE90830_cl_gene_counts.table.txt', header = T)
