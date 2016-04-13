rm(list=ls())
gc()
graphics.off()
clear
library(ggplot2)
library(reshape2)

setwd("/Users/liurui/Desktop/hmm_bioinformatics codes/for lung injury")
ch.dia<-read.table('bubble4_data.txt')
ch.dia2<-read.table('bubble4_data_time2.txt')
ch.dia<-rbind(ch.dia,ch.dia2)
colnames(ch.dia)=c('cluster','gene.num','pvalue','fold.change')
#ch1<-ch.dia[order(ch.dia$gene.num,decreasing=T),]
ch1<-ch.dia
ch1$pvalue=log(ch1$pvalue,5)

for(i in 1:2){#####color level definition: 1-non chronic disease; 2-out of top20; 3-in top20   
  b=rep(2,nrow(ch1))
  b[14:26]=1
  b=as.factor(b);
  if(i==1){ch1$cost.lv=b;}else ch1$CNT.lv=b
}

col=c("green","deeppink","red")
rgb=col2rgb(col)
rgb=rbind(rgb/255,c(0.5,0.5,1))
col=rgb(red=rgb[1,],green=rgb[2,],blue=rgb[3,],alpha=rgb[4,])


ce<-log(ch1$gene.num,12)*8
plot(fold.change~pvalue,data=ch1,pch=20,xlim=c(-6.5,-2.5),ylim=c(0.5,3.8),
     cex=ce,cex.lab=1.2,cex.axis=1.2,fg="white",
     col=col[ch1$cost.lv],xlab="P-value(log5)",ylab="Fold change",font.lab=2,font.axis=2
)
box(lwd=3)
text(ch1$pvalue,ch1$fold.change,ch1$cluster,cex=0.8)