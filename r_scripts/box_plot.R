library(ggplot2)

my_theme <-   theme(axis.text = element_text(size = 16),
                    axis.title = element_text(size = 18,face = "bold"),
                    plot.title = element_text(size = 16, face = 'bold'),
                    strip.text.x = element_text(size=14, face="bold"),
                    strip.background = element_rect(colour="burlywood", fill="burlywood"),
                    legend.text = element_text(size = 16, face = 'bold' ),
                    legend.title = element_text(size=16, face='bold'))

ogg1_point5=c(332,347,319,347,343,302,357,290,318,300)
wt_point5=c(28,35,32,32,35,24,34,30,25,24)
ogg1=data.frame(matrix(ogg1_point5, nrow = 10, ncol = 1))
wt=data.frame(matrix(wt_point5, nrow=10, ncol = 1))
colnames(ogg1)='Counts'
ogg1$name1 = replicate(10, 'ogg1')
colnames(wt)='Counts'
wt$name1 = replicate(10,'wt')

mm=rbind(ogg1, wt)

ggplot()+geom_boxplot(data = mm, aes(name1, counts)) + geom_point(aes(x='ogg1', y=378), color='red',size=2) + 
  geom_point(aes(x='wt', y=31), color='red', size=3) + my_theme + labs(x='Sample', y='Counts')
