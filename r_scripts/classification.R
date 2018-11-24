library(devtools)
library(ggbiplot) ## to make a pretty ggbiplot
library(ggplot2)
library(tidyr)
library(mosaic) 
library(stringi)
library(readr)
library(stringr)
library(pvclust) ## calculate p values for clustering
library(reshape2)
## scree plot
PVEplot <- function(x){qplot(c(1:10), x[1:10]) + 
    geom_line() + 
    xlab("Principal Component") + 
    ylab("PVE") +
    ggtitle("Scree Plot") +
    ylim(0, 1)}
## biplot
my_theme = theme(legend.direction = 'horizontal', legend.position = 'top',
                 axis.text = element_text(size = 16),
                 axis.title = element_text(size = 18,face = "bold"),
                 plot.title = element_text(size = 20, face = 'bold'),
                 strip.text.x = element_text(size=16, face="bold"),
                 strip.background = element_rect(colour="burlywood", fill="burlywood"),
                 legend.text = element_text(size = 16, face = 'bold' ),
                 legend.title = element_text(size=16, face='bold'))
gbiplot <- function(pca.a, groups){ggbiplot(pca.a, obs.scale = 1, var.scale = 1, 
                                          ellipse = TRUE, groups = groups,
                                          circle = TRUE, inherit.aes = FALSE) + scale_color_discrete(name = '') 
}

## function to calculate pve
pve = function(pca.b) {pca.b$sdev^2/sum(pca.b$sdev^2)}
save_path = '/Users/Yun/Documents/bacteria_G4/D_thermus/stats_fig/'
##---------------------------Dt---------------------------------
dt = read.csv('/Users/Yun/Documents/bacteria_G4/D_thermus/D_t_H_cluster.csv', stringsAsFactors = FALSE)

dt = dt[, 2:25] # remove the index
dt = dt[!duplicated(dt$Name),]
dt_g4_raw = dt[,5:24]
norm = function(x){
  x/sum(x)
}
dt_g4_sd = data.frame(t(apply(dt_g4_raw, 1, norm)))
m_dt = t(apply(dt_g4_raw, 1, norm))
pca.dt = prcomp(dt_g4_sd,scale = TRUE)
pca.dt$sdev
scores.dt = data.frame(pca.dt$x) ## all PCs
pve.dt = pve(pca.dt)
PVEplot(pve.dt)
ggsave('Scree_plot_dt.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
##-----------boxplot-------------
dt.t = data.frame(t(dt_g4_sd))
colnames(dt.t) = dt$GCA_access
dt.t['Distance_to_TSS'] = as.factor(seq(-300, 270, 30))
dt.t = melt(dt.t, id.vars = 'Distance_to_TSS', variable.name = 'GCA_access')
dt.t = merge(dt.t, dt[c('GCA_access', 'genus', 'order')], by = 'GCA_access')
ggplot(dt.t, aes(Distance_to_TSS, value)) + geom_violin(aes(fill = order)) + facet_wrap( ~ order, ncol = 1) +
  my_theme + labs(x='Distance to TSS (bp)', y = 'Fraction of PQSs around TSS')
ggsave('violin plot of DT order around TSS2.svg', path = save_path, scale = 1, width = 5,
       height = 4, units='in')
ggplot(dt.t, aes(Distance_to_TSS, y = value)) + 
      geom_boxplot(aes(fill = order), alpha = 0.8) + my_theme + ylim(0, 0.3) +
      labs(x = 'Distance to TSS', y = 'PQS fraction around TSS')


ggsave('boxplot of DT genus around TSS.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
## biplot(pca.dt)
dt.bi = gbiplot(pca.dt, dt$genus) + my_theme
dt.bi
ggsave('PCA of DT.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
pc12 = as.matrix(dt_g4_sd)%*%pca.dt$rotation[,1:2]
row.names(pc12) = dt$genus
pc12 = data.frame(pc12)
## ------plot on the first two PCs----------
ggplot(pc12, aes(PC1, PC2)) + 
  modelr::geom_ref_line(h = 0) +
  modelr::geom_ref_line(v = 0) +
  geom_text(aes(label = dt$genus, color = dt$genus), size = 4) +
  xlab("First Principal Component") + 
  ylab("Second Principal Component") + 
  ggtitle("First Two Principal Components of PQS distribution in D.thermus")
ggsave('First_two_pca_of_DT.png', path =  save_path,scale = 1, width = 10,
       height = 8, units='in')
##------------------dt classification--------------------
library(cluster)    # clustering algorithms
library(factoextra) 
library(ggdendro)
dt_g4_sd = data.frame(t(apply(dt_g4_raw, 1, norm)))
rownames(dt_g4_sd) = dt$Name
dt.dist = dist(dt_g4_sd)
dt.hc = hclust(dt.dist, method = 'ward.D2') ## the same result as kmeans
plot(dt.hc, cex = 0.6, hang = -1, labels = dt$genus)
dt.cl = cutree(dt.hc, k = 2)
dt$cluster = dt.cl
## use ggdendrogram
ggdendrogram(dt.hc, rotate = TRUE, size = 2) 
dhc <- as.dendrogram(dt.hc)
# Rectangular lines
ddata <- dendro_data(dhc, type = "rectangle")
# set the genus along with the label to set color
DT.genus = gsub(' .+', '', ddata$labels$label) 
p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_text(data = ddata$labels, aes(x, y, label = label, color = DT.genus, fontface = 'bold'),
            hjust = 1, angle = 0, size = 4) + 
  coord_flip() + 
  ylim(-0.6, 0.75) + 
  labs(x = 'Species', y = 'Distance') + 
  my_theme
p
## fill color with text
## get the order from the genus, map back. cannot use dt$order because the the names was rearranged when clustering
DT.order = gsub('.+era$', 'Deinococcales',gsub('.+ccus$', 'Deinococcales',gsub('.+mus$', 'Thermales', DT.genus)))
p2 <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_label(data = ddata$labels, aes(x, y, label = label, fill = DT.order), 
             colour = 'white', hjust = 1, angle = 0, size = 3, check_overlap = TRUE) + 
  coord_flip() + 
  ylim(-0.3, 0.75) + 
  labs(x = 'Species', y = 'Distance') + 
  my_theme
p2
ggsave('dendrogram_DT_fill.svg', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
### caculate the p value
dt.pv = pvclust(dt_g4_sd, method.hclust = 'ward.D2', method.dist = 'correlation', nboot = 100)
plot(dt.pv, hang = -1, cex = 0.5)
pvrect(dt.pv, alpha = 0.95 )



##-------------------all bacteria--------------
bacteria = read_csv('/Users/Yun/Documents/bacteria_G4/D_thermus/bacteria_all_G4_vector_.csv' )
## remove rows that have 0 G4s around TSS
bacteria = bacteria[!apply(bacteria[,6:25], 1, function(x){sum(x)==0}),]
bacteria_sd = data.frame(t(apply(bacteria[,6:25], 1, norm)))
# center data in the columns
bacteria_sd = data.frame(apply(bacteria_sd, 2, scale))
pca.b = prcomp(bacteria_sd)
pca.b$sdev
pve.b = pve(pca.b)
PVEplot(pve.b)
gbiplot(pca.b, bacteria$Group) + my_theme


##--------------------Acidobacteria-----------
pca.acido = prcomp(bacteria_sd[bacteria$Group=='Acidobacteria',], scale. = T)
pve.acido = pve(pca.acido)
PVEplot(pve.acido)
gbiplot(pca.acido, bacteria[bacteria$Group=='Acidobacteria',]$order)
##--------------------proteobacteria-------
pca.proteo = prcomp(bacteria_sd[bacteria$Group=='Proteobacteria',], scale. = T)
pve.proteo = pve(pca.proteo)
PVEplot(pve.proteo)
gbiplot(pca.proteo, bacteria[bacteria$Group=='Proteobacteria',]$SubGroup) + my_theme
##-------------box plot in proteobacteria-------------
colnames(bacteria_sd) = factor(seq(-285,285,30))
rownames(bacteria_sd) = bacteria$GCA_access
bacteria_sd['GCA_access'] = bacteria$GCA_access
melt.bacteria = melt(bacteria_sd, id.vars = 'GCA_access',variable.name = 'Distance_to_TSS')
melt.bacteria = merge(melt.bacteria, bacteria[,c('GCA_access', 'order','family', 'Group', 'SubGroup')], by = 'GCA_access')
## remove NA
melt.bacteria = melt.bacteria[!is.na(melt.bacteria$order),]


plot.subgroup = function(x, col = 'SubGroup'){
  ggplot(melt.bacteria[melt.bacteria[col]==x,], aes(factor(Distance_to_TSS), value)) + 
    geom_boxplot(aes(fill = col)) + facet_grid(order ~ .) + my_theme + 
  labs(x='Distance to TSS', y = 'Fraction of PQS around TSS') + ggtitle(i)
  
  ggsave(paste(gsub("[[:punct:]]", " ", x),'.png'), path = save_path, scale = 1, width = 10,
         height = 8, units='in')
}
for (i in unique(bacteria[bacteria$Group == 'Terrabacteria group',]$SubGroup)){
  plot.subgroup(i)}
for (i in unique(bacteria[bacteria$Group == 'Terrabacteria group',]$order)){
  plot.subgroup(i, col='order')}
for (i in unique(bacteria[bacteria$Group == 'Proteobacteria',]$SubGroup)){
  plot.subgroup(i)
}

for (i in unique(bacteria[bacteria$SubGroup == 'Alphaproteobacteria',]$order)){
  plot.subgroup(i, col = 'order')
}

for (i in unique(bacteria[bacteria$SubGroup == 'Gammaproteobacteria',]$order)){
  plot.subgroup(i, col = 'order')
}

##-------------gamma--------
table(proteo$SubGroup)
pca.ga = prcomp(bacteria_sd[bacteria$SubGroup=='Gammaproteobacteria',], scale. = T)
pve.ga = pve(pca.ga)
PVEplot(pve.ga)
gbiplot(pca.ga, bacteria[bacteria$SubGroup=='Gammaproteobacteria',]$order)
# Xanthomonadales has more G4 in x9, 10, 11
##--------------delta---------
pca.de = prcomp(bacteria_sd[bacteria$SubGroup=='delta/epsilon subdivisions',],scale. = T)
pve.de = pve(pca.de)
PVEplot(pve.de)
## one family in delta stands out a lot: Myxococcales
gbiplot(pca.de, bacteria[bacteria$SubGroup=='delta/epsilon subdivisions',]$order)

## --------------Terrabacteria-----------
pca.te = prcomp(bacteria_sd[bacteria$Group=='Terrabacteria group',],scale. = T)
pve.te = pve(pca.te)
PVEplot(pve.te)
gbiplot(pca.te, bacteria[bacteria$Group=='Terrabacteria group', ]$SubGroup) + my_theme
##--------------actino--------
pca.ac = prcomp(bacteria_sd[bacteria$SubGroup=='Actinobacteria',], scale. = T)
pve.ac = pve(pca.ac)
PVEplot(pve.ac)
## streptomyce stand out
gbiplot(pca.ac, bacteria[bacteria$SubGroup=='Actinobacteria',]$family) + my_theme
## ----------------cyanobacteria------------------------
pca.cyano = prcomp(bacteria_sd[bacteria$SubGroup=='Cyanobacteria/Melainabacteria group',])
pve.cyano = pve(pca.cyano)
PVEplot(pve.cyano)
gbiplot(pca.cyano, groups = bacteria[bacteria$SubGroup=='Cyanobacteria/Melainabacteria group',]$order)

##----------------all cyano-----------------------------
cyano_all = read.csv('/Users/Yun/Documents/bacteria_G4/cyno_bacteria/cyano_G4_vector.txt', stringsAsFactors = FALSE)
cyano_all = cyano_all[!apply(cyano_all[,4:23], 1, function(x){sum(x)==0}),]
cyano_all = cyano_all[!(cyano_all$order  %in% c('', 'Not assigned')), ]
cyano_all_scale = data.frame(t(apply(cyano_all[,4:23], 1, norm)))
pca.cyano2 = prcomp(cyano_all_scale, scale. = T)
pve.cyano.all = pve(pca.cyano2)
PVEplot(pve.cyano.all)
gbiplot(pca.cyano2, cyano_all$order) + my_theme
ggsave('PCA of the cyano bacteria.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
## plot
cyano_t = data.frame(t(cyano_all_scale))
colnames(cyano_t) = cyano_all$GCA_access
cyano_t['Distance_to_TSS'] = as.factor(seq(-285, 285, 30))
df = melt(cyano_t, id.vars = 'Distance_to_TSS', variable.name = 'GCA_access')
df = merge(df, cyano_all[c('GCA_access', 'family', 'order')], by = 'GCA_access')
ggplot(df, aes(x=Distance_to_TSS, y = value)) + geom_boxplot(aes(fill=order)) 
ggplot(df, aes(x=Distance_to_TSS, y = value)) + geom_boxplot(aes(fill = order)) + 
  facet_grid(order~.) + ylim(0, 0.4) +
  my_theme + 
  labs(x= 'Distance to TSS', y = 'G4 fraction around TSS')
ggsave('boxplot of cyano bacteria around TSS.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
ggplot(df, aes(x=Distance_to_TSS, y = value)) + geom_violin(aes(fill = order)) + 
  facet_grid(order~.) + 
  my_theme + 
  labs(x= 'Distance to TSS', y = 'G4 fraction around TSS')
ggsave('violin plot of cyano bacteria around TSS.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
ggplot(df[df$order=='Nostocales',], aes(x=Distance_to_TSS, y = value)) + geom_boxplot(aes(fill = family)) + 
  facet_grid(family~.) + 
  my_theme + 
  labs(x= 'Distance to TSS', y = 'G4 fraction around TSS')
ggsave('boxplot of Nostocales order around TSS.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
ggplot(df[df$order=='Chroococcales',], aes(x=Distance_to_TSS, y = value)) + geom_boxplot(aes(fill = family)) + 
  facet_grid(family~.) + 
  my_theme + 
  labs(x= 'Distance to TSS', y = 'G4 fraction around TSS')
ggsave('boxplot of Chroococcales order around TSS.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
## Stigonematales, highly enriched
ggplot(df[df$order=='Stigonematales',], aes(x=Distance_to_TSS, y = value)) + geom_boxplot(aes(fill = family)) + 
  facet_grid(family~.) + 
  my_theme + 
  labs(x= 'Distance to TSS', y = 'G4 fraction around TSS')
ggsave('boxplot of Stigonematales order around TSS.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
##------strep--------------------------------------------------------------- 
strep = t(bacteria_sd[bacteria$family=='Streptomycetaceae',])
rownames(strep) = as.factor(seq(-285, 285, 30))
strep = t(strep)
strep = melt(strep)
ggplot(strep, aes(x=as.factor(Var2), y = value)) + geom_violin(colour = 'blue') + 
  my_theme + labs(x = 'Distance to TSS', y = 'Fractions of PQS around TSS') + 
  ggtitle('Distribution of PQS around TSS in Streptomycetaceae')
ggsave('Distribution of PQS around TSS in Streptomycetaceae.png',path = save_path, scale = 1, width = 10,
       height = 8, units='in')

