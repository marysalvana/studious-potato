rm(list=ls())


getReference<-function(n,p,estimate,maxTimeLag)
{
  dataRef=NULL;dataX=NULL;
  
  #####--Generate data--#####
  
  for(times in 1:10)
  {
    Y=simulateEsimatedSepData(n,p,maxTimeLag,estimate$covar_L,estimate$Sigma11_L,estimate$mu);
    Y=matrix(Y,n);
    g=get_g(Y,n);
    dataRef=cbind(dataRef,g);
  }
  ct=dim(dataRef)[2];
  gg=fbplot(dataRef,plot=F);
  nonOut=setdiff(1:ct,gg$outpoint);
  dataRef=dataRef[,nonOut];
  
  for(times in 1:5)
  {
    Y=simulateEsimatedSepData(n,p,maxTimeLag,estimate$covar_L,estimate$Sigma11_L,estimate$mu);
    Y=matrix(Y,n);
    g=get_g(Y,n);
    dataX=cbind(dataX,g);
  }
  ct=dim(dataX)[2];
  gg=fbplot(dataX,plot=F);
  nonOut=setdiff(1:ct,gg$outpoint);
  dataX=dataX[,nonOut];
  
  return(list(dataX=dataX,dataRef=dataRef));
}

generateData4Crit<-function(n,p,estimate,critReplicate,maxTimeLag,crit.path)
{
  for(times in 1:critReplicate)
  {
    Y=simulateEsimatedSepData(n,p,maxTimeLag,estimate$covar_L,estimate$Sigma11_L,estimate$mu);
    save(Y,file=paste(crit.path,"/",times,".Rdata",sep=''));
  }
}

library(fda)
library(fields)

source('./function/function.R')
source('./function/tool.R')

blocksize=20;

# load data Y:
#             dim: n * p
#             n: number of spatial locations
#             p: number of temporal points

g=get_g(Y,n);

estimate=symmCovEstimate(Y,n,blocksize);

refSet=getReference(n,p,estimate,blocksize);
dataX=refSet$dataX;
dataRef=refSet$dataRef;

rank.value=rankTest(dataX,g,dataRef);

print(rank.value);

############ compute critical value ###############
crit.data.path='critData/';
crit.result.path='critResult/';
critReplicate=100;
generateData4Crit(n,p,estimate,critReplicate,blocksize,crit.data.path);

for(i in 1:critReplicate)
{
  load(paste(crit.data.path,"/",i,".Rdata",sep=''));
  Y=matrix(Y,n);
  g=get_g(Y,n);
  
  estimate=symmCovEstimate(Y,n,blocksize);
  
  refSet=getReference(n,p,estimate,blocksize);
  dataX=refSet$dataX;
  dataRef=refSet$dataRef;
  
  rank.value=rankTest(dataX,g,dataRef);
  
  write.table(rank.value,file=paste(crit.result.path,"/",i,".csv",sep=''),col.names = F, row.names = F);
}
