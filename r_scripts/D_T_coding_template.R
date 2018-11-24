library(ggplot2)
library(tidyr)
library(mosaic) 
library(stringi)
library(readr)
library(stringr)
library(reshape2)
save_path = '/Users/Yun/Documents/bacteria_G4/D_thermus/stats_fig/'
my_theme = theme(legend.direction = 'horizontal', legend.position = 'top',
                 axis.text = element_text(size = 14),
                 axis.title = element_text(size = 18,face = "bold"),
                 plot.title = element_text(size = 18, face = 'bold'),
                 strip.text.x = element_text(size=16, face="bold"),
                 strip.background = element_rect(colour="burlywood", fill="burlywood"),
                 legend.text = element_text(size = 16, face = 'bold' ),
                 legend.title = element_text(size=16, face='bold'))
my_theme_2col = theme(legend.direction = 'horizontal', legend.position = 'top',
                 axis.text = element_text(size = 8),
                 axis.title = element_text(size = 18,face = "bold"),
                 plot.title = element_text(size = 20, face = 'bold'),
                 strip.text.x = element_text(size=12, face="bold"),
                 strip.background = element_rect(colour="burlywood", fill="burlywood"),
                 legend.text = element_text(size = 16, face = 'bold' ),
                 legend.title = element_text(size=16, face='bold'))
## ------------data processing------------------
dt_code = read_table2('/Users/Yun/Documents/bacteria_G4/D_thermus/spare_G4_code_distribution.txt')
dt_template = read_table2('/Users/Yun/Documents/bacteria_G4/D_thermus/spare_G4_template_distribution.txt')
dt_spare = read_table2('/Users/Yun/Documents/bacteria_G4/D_thermus/spare_G4_distribution.txt')
dt_template_all = read_table2('~/Documents/bacteria_G4/D_thermus/all_G4_template_distribution.txt')
dt_info = read.table('/Users/Yun/Documents/bacteria_G4/D_thermus/G4_in_D_thermus_final_summary.txt', sep = '\t')
dt_info = dt_info[,c('V1','V14', 'V16', 'V17')]
colnames(dt_info) = c('GCA_access', 'Name', 'Genus', 'Order')
dt_info = dt_info[2:50,]
dt_code = dt_code[, 1:21] # remove the index and other 
dt_template = dt_template[,1:21]
dt_template_all= dt_template_all[,1:21]

dt_spare = dt_spare[,1:21]
dt_code = merge(dt_info, dt_code)
dt_template = merge(dt_info, dt_template)
dt_template_all = merge(dt_info, dt_template_all)
dt_spare=merge(dt_info, dt_spare)
dt_code_raw = dt_code[,5:24]
dt_template_raw = dt_template[,5:24]
dt_template_all_raw = dt_template_all[,5:24]
dt_spare_raw=dt_spare[,5:24]
norm = function(x){
  x/sum(x)
}
dt_code_sd = data.frame(t(apply(dt_code_raw, 1, norm)))
dt_template_sd = data.frame(t(apply(dt_template_raw, 1, norm)))
dt_spare_sd=data.frame(t(apply(dt_spare_raw, 1, norm)))
dt_template_all_sd=data.frame(t(apply(dt_template_all_raw, 1, norm)))

##-------------------box plot------------------
dt.code.t = data.frame(t(dt_code_raw))
## use dt_code_raw for the raw counts visulization otherwise use the dt_code_sd for standarized data
dt.template.t = data.frame(t(dt_template_raw))

colnames(dt.code.t) = dt_code$GCA_access
colnames(dt.template.t) = dt_template$GCA_access

dt.code.t['Distance_to_TSS'] = factor(seq(-300,270, 30))
dt.template.t['Distance_to_TSS'] = factor(seq(-300,270, 30))

dt.code.t = melt(dt.code.t, id.vars = 'Distance_to_TSS',variable.name = "GCA_access")
dt.template.t = melt(dt.template.t, id.vars = 'Distance_to_TSS', variable.name = 'GCA_access')

dt.code.t = merge(dt_info, dt.code.t)
dt.template.t = merge(dt_info, dt.template.t)


dt.spare.t = data.frame(t(dt_spare_sd))
colnames(dt.spare.t) = dt_spare$GCA_access
dt.spare.t['Distance_to_TSS'] = factor(seq(-300, 270, 30))
dt.spare.t = melt(dt.spare.t, id.vars = 'Distance_to_TSS', variable.name = 'GCA_access')
dt.spare.t = merge(dt_info, dt.spare.t)


ggplot(dt.spare.t, aes(Distance_to_TSS, value)) + geom_violin(aes(fill = Order)) +
  my_theme + ylim(0,0.4)+labs(x='Distance to TSS CODING', y = 'Fraction of PQSs around TSS ')
ggsave('boxplot of DT order around TSS_CODING.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
ggplot(dt.template.t, aes(Distance_to_TSS, value)) + geom_boxplot(aes(fill = Order)) +
  my_theme + ylim(0,0.4) + labs(x='Distance to TSS TEMPLATE', y = 'Fraction of PQS around TSS ')
ggsave('boxplot of DT order around TSS_TEMPLATE.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')

ggplot() + geom_boxplot(data = dt.template.t[dt.template.t$Order == 'Deinococcales', ], aes(factor(Distance_to_TSS), value), fill = 'brown1') + 
  geom_boxplot(data = dt.code.t[dt.code.t$Order =='Deinococcales', ], aes(factor(Distance_to_TSS), value), fill = 'blue') +
  labs(x='Distance to TSS', y = 'Fraction of PQS around TSS ') + my_theme + 
  ylim(0, 0.4) 
ggsave('boxplot of DT order around TSS_TEMPLATE_coding_together.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
##-------------------combine coding and template ------------
dt.code.t$label = replicate(dim(dt.code.t)[1], 'coding')
dt.template.t$label = replicate(dim(dt.template.t)[1], 'template')
dt_label = rbind(dt.code.t, dt.template.t)
ggplot(dt_label) + geom_boxplot(aes(factor(Distance_to_TSS), value, fill= label), size = 0.2, outlier.shape = 19, outlier.size = 0.5,
                                outlier.stroke = 0.2, outlier.alpha = NULL) + 
  facet_grid(Order~.) + my_theme + labs(x='Distance to TSS (bp)', y = 'Fraction of PQSs')
ggsave('boxplot of DT order around TSS_TEMPLATE_coding_sd_spare_G4_raw.svg', path = save_path, scale = 1, width = 7,
       height = 3, units='in')
##------------------t.test on the template strand on spare G4s and all G4s on the central bin-------
#Only look at the Deinococcales 
t.test(dt_template_all_sd[dt_template_all$Order=='Deinococcales',]$G4_counts_at_tss, dt_template_sd[dt_template_all$Order=='Deinococcales',]$G4_counts_at_tss)
##-------------------dt,all------------------

dt_summary = read_table2('~/Documents/bacteria_G4/D_thermus/G4_in_D_thermus_final_summary.txt')
ggplot() + 
##-------------------cyano bacteria, coding strand-------------------
cyano_all.code = read_table2('/Users/Yun/Documents/bacteria_G4/cyno_bacteria/G4_coding_tss_300.txt')
cyano_all.code = cyano_all.code[!apply(cyano_all.code[,2:21], 1, function(x){sum(x)==0}),]
cyano_all.code = cyano_all.code[!(cyano_all.code$order  %in% c('', 'Not', NA)), ]
cyano_all_scale.code = data.frame(t(apply(cyano_all.code[,2:21], 1, norm)))

## plot
cyano_t_c = data.frame(t(cyano_all_scale.code))
colnames(cyano_t_c) = cyano_all.code$GCA_access
cyano_t_c['Distance_to_TSS'] = as.factor(seq(-285, 285, 30))
df.c = melt(cyano_t_c, id.vars = 'Distance_to_TSS', variable.name = 'GCA_access')
df.c = merge(df.c, cyano_all.code[c('GCA_access', 'family', 'order')], by = 'GCA_access')

ggplot(df.c, aes(x=Distance_to_TSS, y = value)) + geom_boxplot(aes(fill = order)) + 
  facet_grid(order~.) + ylim(0, 0.4) +
  my_theme + 
  labs(x= 'Distance to TSS CODING', y = 'G4 fraction around TSS')
ggsave('boxplot of cyano bacteria around TSS CODING.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
##------------------cyano bacteria stemplate strand-----------------
cyano_all.t = read_table2('/Users/Yun/Documents/bacteria_G4/cyno_bacteria/G4_template_tss_300.txt')
cyano_all.t = cyano_all.t[!apply(cyano_all.t[,2:21], 1, function(x){sum(x)==0}),]
cyano_all.t = cyano_all.t[!(cyano_all.t$order  %in% c('', 'Not', NA)), ]
cyano_all_scale.t = data.frame(t(apply(cyano_all.t[,2:21], 1, norm)))


##-----------combine coding and template---------------
df.c['label'] = replicate(dim(df.c)[1], 'coding')
df.t['label'] = replicate(dim(df.t)[1], 'tempalte')
df_cyano = rbind(df.c, df.t)
## plot------

cyano_t_t = data.frame(t(cyano_all_scale.t))
colnames(cyano_t_t) = cyano_all.t$GCA_access
cyano_t_t['Distance_to_TSS'] = as.factor(seq(-285, 285, 30))
df.t = melt(cyano_t_t, id.vars = 'Distance_to_TSS', variable.name = 'GCA_access')
df.t = merge(df.t, cyano_all.t[c('GCA_access', 'family', 'order')], by = 'GCA_access')

ggplot(df.t, aes(x=Distance_to_TSS, y = value)) + geom_boxplot(aes(fill = order)) + 
  facet_grid(order~.) + ylim(0, 0.4) +
  my_theme + 
  labs(x= 'Distance to TSS TEMPLATE', y = 'G4 fraction around TSS')
ggsave('boxplot of cyano bacteria around TSS TEMPLATE.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
ggplot(df_cyano[df_cyano$order!="Chamaesiphonales",]) + geom_boxplot(aes(factor(Distance_to_TSS), value,fill = label), alpha=0.95)+
  my_theme_2col+ facet_wrap(~order, ncol = 2) + scale_y_continuous(limits = c(0,0.35)) + 
  labs(x='Distance to TSS', y = 'Fraction of PQS around TSS')
