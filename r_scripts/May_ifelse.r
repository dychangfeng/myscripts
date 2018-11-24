sdata <- read.csv("~/Desktop/sdata2.csv")
sdata_test <- sdata[1:10, 1:12]
sdata_test%>%select(-X)%>%spread(ragender,r2socwk)
sdata44<- transform(sdata_test, r2socwk = ifelse(ragender==2,r2socwk,s2socwk), 
                    s2socwk = ifelse(ragender == 1, r2socwk,s2socwk ))
