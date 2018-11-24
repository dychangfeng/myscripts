library(devtools)
library(ggplot2)
library(tidyr)
library(mosaic) 
library(stringi)
library(readr)
library(stringr)
library(reshape2)
my_theme = theme(legend.direction = 'horizontal', legend.position = 'top',
                 axis.text = element_text(size = 6),
                 axis.title = element_text(size = 20,face = "bold"),
                 plot.title = element_text(size = 20, face = 'bold'),
                 strip.text.x = element_text(size=16, face="bold"),
                 strip.background = element_rect(colour="burlywood", fill="burlywood"),
                 legend.text = element_text(size = 16, face = 'bold' ),
                 legend.title = element_text(size=16, face='bold'))
norm = function(x){
  x/sum(x)
}
save_path = '/Users/Yun/Documents/bacteria_G4/stat_figures/'
bacteria = read_csv('/Users/Yun/Documents/bacteria_G4/D_thermus/bacteria_all_G4_vector_.csv' )

check_na = function(x){
  !any(is.na(x))
}
## remove rows that have 0 G4s around TSS
bacteria = bacteria[!apply(bacteria[,6:25], 1, function(x){sum(x)==0}),]
bacteria_sd = data.frame(t(apply(bacteria[,6:25], 1, norm)))
colnames(bacteria_sd) = factor(seq(-285,285,30))
rownames(bacteria_sd) = bacteria$GCA_access
bacteria_sd['GCA_access'] = bacteria$GCA_access
melt.bacteria = melt(bacteria_sd, id.vars = 'GCA_access',variable.name = 'Distance_to_TSS')
melt.bacteria = merge(melt.bacteria, bacteria[,c('GCA_access', 'order','family', 'Group', 'SubGroup')], by = 'GCA_access')
melt.bacteria = melt.bacteria[complete.cases(melt.bacteria),] ## remove rows with NAs
##----------coding-----------------------
coding = read_tsv('/Users/Yun/Documents/bacteria_G4/bacteria_representive/G4_counts_closest_tss_300_code.txt')
coding = coding[,1:21] ## only keep the vector
coding = coding[!apply(coding[,2:21],1, function(x){sum(x)==0}),]
coding_sd = data.frame(t(apply(coding[,2:21], 1, norm)))
colnames(coding_sd) = factor(seq(-285,285,30))
rownames(coding_sd) = coding$GCA_access
coding_sd['GCA_access'] = coding$GCA_access
melt.coding = melt(coding_sd, id.vars = 'GCA_access',variable.name = 'Distance_to_TSS')
melt.coding = merge(melt.coding, bacteria[,c('GCA_access', 'order','family', 'Group', 'SubGroup')], by = 'GCA_access')
melt.coding = melt.coding[complete.cases(melt.coding),] ## remove rows with NAs
##----------template------------------
template = read_tsv('/Users/Yun/Documents/bacteria_G4/bacteria_representive/G4_counts_closest_tss_300_template.txt')
template = template[,1:21] ## only keep the vector
template = template[!apply(template[,2:21],1, function(x){sum(x)==0}),]
template_sd = data.frame(t(apply(template[,2:21], 1, norm)))
colnames(template_sd) = factor(seq(-285,285,30))
rownames(template_sd) = template$GCA_access
template_sd['GCA_access'] = template$GCA_access
melt_template = melt(template_sd, id.vars = 'GCA_access',variable.name = 'Distance_to_TSS')
melt_template = merge(melt_template, bacteria[,c('GCA_access', 'order','family', 'Group', 'SubGroup')], by = 'GCA_access')
melt_template = melt_template[complete.cases(melt_template),] ## remove rows with NAs
##----------------------combine coding and template data frames---------------------
#melt.bacteria$label = replicate(dim(melt.bacteria)[1], 'all')
melt.coding$label = replicate(dim(melt.coding)[1], 'coding')
melt_template$label = replicate(dim(melt_template)[1], 'template')
bacteria_lable = rbind(melt.coding, melt_template)
##-----------------functions to get data and plot data---------------
#function to get the specific order of data from the all the bacteria
# return the list of three data frames, all, coding, template
get_specific_data = function(df,df_coding,df_template, specific = 'order' ,t= 'Pseudomonadales' ){
  pse = df[df[,specific] == t,][!is.na(df[,specific]),]
  pse_coding = df_coding[df_coding[,specific] == t,][!is.na(df_coding[,specific]),]
  pse_template = df_template[df_template[,specific] == t,][!is.na(df_template[,specific]),]
  return(list(pse, pse_coding, pse_template))
  ## return the list of targeted data frame 
}

bac_box = function(data){ggplot(data, aes(factor(Distance_to_TSS), value)) + geom_violin(aes(fill = order)) +
    my_theme + labs(x='Distance to TSS', y = 'Fraction of PQS around TSS')}
bac_box_family = function(data){ggplot(data, aes(factor(Distance_to_TSS), value)) + geom_boxplot(aes(fill = family)) +
    my_theme + labs(x='Distance to TSS', y = 'Fraction of PQS around TSS')}
code_temp_plot = function(x, my_col = 'order'){
  ggplot(bacteria_lable[bacteria_lable[,my_col]==x,]) + geom_boxplot(aes(factor(Distance_to_TSS), value, fill = label), alpha=0.95)+
    my_theme+ facet_wrap(~order, ncol = 3) + scale_y_continuous(limits = c(0,0.4)) + 
    labs(x='Distance to TSS', y = 'Fraction of PQS around TSS')
  ggsave(paste(gsub("[[:punct:]]", " ", x),'.png'), path = save_path, scale = 1, width = 10,
         height = 8, units='in')
}
##----------------------Xanthomonadales--------------------------------------
xan_list = get_specific_data(melt.bacteria, melt.coding, melt_template, t = 'Xanthomonadales')
xan = xan_list[[1]]
xan_coding = xan_list[[2]]
xan_template = xan_list[[3]]
bac_box(xan) 
bac_box(xan_coding)
bac_box(xan_template)
bac_box_family(xan_template)

ggsave('Xanthomonadales_coding_template.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in') 
save_path = '/Users/Yun/Documents/bacteria_G4/Xanthomonas_campestris/figures/'
ggsave('Xanthomonadales.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
code_temp_plot('Xanthomonadales')
subgroups_care = unique(bacteria[bacteria$Group %in% c("Terrabacteria group","Proteobacteria","FCB group"),]$SubGroup)

##------using for loop to plot all the orders in the major SugGroups 

for (i in subgroups_care){
  #ERROR HANDLING
  possibleError <- tryCatch({
    code_temp_plot(i, my_col = 'SubGroup')
  }
  ,
  error=function(e) {
    e
    print(paste("Oops! --> Error in Loop ",i,sep = ""))
  }
  ) #catch the possible error and then use the next to skip the error
  if(inherits(possibleError, "error")) next ## skip the ones that has errors in them
  print(paste("  End Loop ",i,sep = ""))
}


##---------------Pseudomonadales--------------

pse_list = get_specific_data(melt.bacteria, melt.coding, melt_template, t_order = 'Pseudomonadales')
pse = pse_list[[1]]
pse_coding = pse_list[[2]]
pse_template = pse_list[[3]]

bac_box(pse)
bac_box(pse_coding)
bac_box(pse_template)
bac_box_family(pse_coding)
bac_box_family(pse_template)
bac_box_family(pse_template[pse_template$family=="Pseudomonadaceae",])
code_temp_plot(pse_coding, pse_template)
##---------------------Streptomycetaceae----------------------
strep_list = get_specific_data(melt.bacteria, melt.coding, melt_template, specific = 'family', t = 'Streptomycetaceae')
strep = strep_list[[1]]
strep_coding = strep_list[[2]]
strep_template = strep_list[[3]]
bac_box_family(strep)
bac_box_family(strep_coding)
bac_box_family(strep_template)
code_temp_plot(strep_coding,strep_template)
