rm(list=ls())
gc()
setwd('C:\\Users\\starml\\Desktop\\EM\\HMM\\lung injury')
control.data<-read.table('lung_injury_control.txt',header=T)
case.data<-read.table('lung_injury_case.txt',header=T)
control.scaled<-data.frame(t(apply(control.data[,-1],1,scale)))
rownames(control.scaled)<-control.data$Symbol
case.scaled<-data.frame(t(apply(case.data[,-1],1,scale)))
rownames(case.scaled)<-case.data$Symbol

fun.pvalue<-function(arg)  {
  col.num<-length(arg)
  d1<-arg[1:(col.num/2)]
  d2<-arg[(col.num/2+1):col.num]
  tmp<-t.test(d1,d2)
  return(tmp$p.value)
}

temp.data<-cbind(control.scaled[,25:30],case.scaled[,25:30])
pv<-apply(temp.data,1,fun.pvalue)
pved<-sort(pv,index.return=T)

thresh.turn=0.05
turn.num<-0
for(i in 1:500)  {
  ps<-t.test(case.scaled[pved$ix[i],7:18],case.scaled[pved$ix[i],31:42])
  if(ps$p.value<=thresh.turn)  {
    turn.num<-turn.num+1
  }
}
print(turn.num)
