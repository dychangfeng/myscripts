library(DESeq2)
work_dir = '/Users/Yun/Documents/bacteria_G4/D_thermus/D_r2/'
counts.matrix <- read.delim('/Users/Yun/Documents/bacteria_G4/D_thermus/D_r2/RNAseq2/GSM2572202_Z_0608Z0225-TR1.txt', row.names=1)
## convert float to int
counts.matrix2 = apply(counts.matrix, 2, as.integer)
counts.matrix2 = counts.matrix2[complete.cases(counts.matrix2),] ## complete cases label rows with NAs
rownames(counts.matrix2) = rownames(counts.matrix[complete.cases(counts.matrix),]) ## reomve rownames with NAs
coldata = read.delim('/Users/Yun/Documents/bacteria_G4/D_thermus/D_r2/RNAseq2/coldata2.txt', row.names = 1)
dataset <- DESeqDataSetFromMatrix(countData = counts.matrix2, colData = coldata, design = ~ condition)
nrow(dataset)
dds <- dataset
results.dataset <- DESeq(dds)
res<- results(results.dataset)
sum(res$padj < 0.05, na.rm=TRUE) #5% false pos acceptable
sum(res$log2FoldChange[res$padj<0.05] > 0.5 ,na.rm = T) # 133 has more than 1 fold change
sum(res$log2FoldChange[res$padj<0.05] < -0.5, na.rm = T ) # 187 has more than -1 fold change
## --------plot data -------------
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
install.packages('calibrate')
library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=rownames(counts.matrix), cex=.8))
## up regulation
up = subset(res, padj<.05 & log2FoldChange>0.25)
down = subset(res, padj<.05 & log2FoldChange < -0.25)
write.table(down, file = '~/Documents/bacteria_G4/D_thermus/D_radiodurans/RNA_seq/down_regulate_quater_fold.txt', sep = '\t',
            fileEncoding = "")
write.table(up, file = '~/Documents/bacteria_G4/D_thermus/D_radiodurans/RNA_seq/up_regulate_quater_fold.txt', sep = '\t')
##---------------RNA seq3--------------------
RNAseq3 <- read.delim('/Users/Yun/Documents/bacteria_G4/D_thermus/D_r2/RNAseq3/GSE9636_OxyRmergeddata.txt', stringsAsFactors = FALSE)
colnames(RNAseq3) = c('genes', 'gene name', 'M', 'p.value')
RNAseq3 = RNAseq3[,1:4]
RNAseq3 = RNAseq3[complete.cases(RNAseq3),]
RNAseq3 = RNAseq3[RNAseq3$p.value<0.05,]
up_30nt = read.delim(paste(work_dir, 'G4s_genes_list_within_30nt.txt', sep = ''), header = FALSE)
colnames(up_30nt) = 'genes'
merge(up_30nt, RNAseq3)
