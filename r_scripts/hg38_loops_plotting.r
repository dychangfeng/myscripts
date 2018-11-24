install.packages('stringr')
library(stringr)
hg38_loop <- read.delim('~/ncbi_genomes/hg38_all_g4_loops.csv', sep = ',')
library('tidyr')
library('dplyr')
library('ggplot2')
hg38_loop <- hg38_loop%>%select(-X)
t <- hg38_loop%>%group_by(g_runs)%>%summarise(count = n())%>%arrange(desc(count))%>%
  filter(count>1000) ##popular g_runs, filter g_rus that has more than 1000 counts
#-----------------deal with G4 without spare tire and extra G------------------
hg38_loop_3333 <- hg38_loop%>%filter(g_runs == '3:3:3:3:')%>%select(-g_runs, -loop4, -loop5)
hg38_3333_group1 <- hg38_loop_3333%>%group_by(loop1)%>%summarise()
hg38_3333_tb <- data.frame(unclass(summary(hg38_loop_3333, 100)), 
                            check.names = FALSE, stringsAsFactors = FALSE)
colnames(hg38_3333_tb) <- c('loop1', 'loop2', 'loop3')

##------separate three data frames
loop1_3333 <- hg38_3333_tb%>%separate(loop1, into = c('loop1_bases', 'loop1_counts'), 
                                      sep= ':', extra ='merge')%>%select(loop1_bases,loop1_counts)
loop2_3333 <- hg38_3333_tb%>%separate(loop2, into = c('loop2_bases', 'loop2_counts'), 
                                      sep= ':', extra ='merge')%>%select(loop2_bases,loop2_counts)
loop3_3333 <- hg38_3333_tb%>%separate(loop3, into = c('loop3_bases', 'loop3_counts'), 
                                      sep= ':', extra ='merge')%>%select(loop3_bases,loop3_counts)
loop1_3333$loop1_counts <- as.integer(loop1_3333$loop1_counts)
loop2_3333$loop2_counts <- as.integer(loop2_3333$loop2_counts)
loop3_3333$loop3_counts <- as.integer(loop3_3333$loop3_counts)
##----------------remove extra space in the bases, merge to one data frame----------
loop1_3333$loop1_bases <-  sapply(loop1_3333$loop1_bases, function(x) str_replace_all(x, fixed(' '), ''))
loop2_3333$loop2_bases <-  sapply(loop2_3333$loop2_bases, function(x) str_replace_all(x, fixed(' '), ''))
loop3_3333$loop3_bases <-  sapply(loop3_3333$loop3_bases, function(x) str_replace_all(x, fixed(' '), ''))
loop1_2_3333 <- merge.data.frame(loop1_3333, loop2_3333, by.x = 'loop1_bases', by.y = 'loop2_bases',
                 all.x = TRUE, all.y = TRUE)
  
loop1_2_3_3333 <- merge.data.frame(loop1_2_3333, loop3_3333, by.x = 'loop1_bases', by.y = 'loop3_bases',
                                 all.x = TRUE, all.y = TRUE)
loop1_2_3_3333 <- loop1_2_3_3333[order(loop1_2_3_3333$loop1_counts, decreasing = TRUE),]%>%
                 mutate
loop1_2_3_3333_rmNA <- loop1_2_3_3333[!(apply(loop1_2_3_3333, 1, function(x) any(is.na(x)))),]%>%
                       mutate(all_loop = loop1_counts + loop2_counts + loop3_counts)
ggplot(loop1_2_3_3333_rmNA[2:10,]) + geom_point(aes(x=loop1_bases, y = all_loop))

