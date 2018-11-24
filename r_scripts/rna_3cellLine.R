library(readr)
library(tidyr)
library(dplyr)
library(stringi)
library(stringr)
library(reshape2)
rna_3cell = read.delim('~/Desktop/RNA_cell_line/u87_hela_hepG2.txt', sep  = '\t', header  = T)
colnames(rna_3cell) = c('ensemblID', 'geneSymbol', 'cellLine', 'TPM', 'unit')
rna_3a = rna_3cell%>%spread(cellLine, TPM)
rna_3 = rna_3a%>%select( `U-87 MG`, HeLa, `Hep G2`)
rna_3 = rna_3[complete.cases(rna_3),]
rna_3 = data.frame(t(rna_3))
rownames(rna_3) = c('U87', 'HeLa', 'HepG2')
p_code = read.delim('~/Desktop/RNA_cell_line/code_raw.txt', sep= '\t', header =  T)
p_code  = data.frame(t(p_code))
colnames(p_code) = c('SV40','TATA_Box','SPHII','CoreC','PVUII')

rna_3_code = cbind(rna_3, p_code)

pred = colnames(rna_3)[17:20]
fla <- paste(" cbind(SV40, TATA_Box, SPHII, CoreC, PVUII) ~", paste(pred, collapse=" + "))
as.formula(fla) ## change the form to formular
nf_lm = lm(as.formula(fla), data = rna_3_code)
nf_lm$coefficients
nf_lm$coefficients[!is.na((nf_lm$coefficients))]
lm_fit = lm.fit(as.matrix(rna_3), as.matrix(p_code))
lm_fit$coefficients[apply(lm_fit$coefficients, 2, FUN = function(x){complete.cases(x)}),]
