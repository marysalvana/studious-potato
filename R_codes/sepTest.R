rm(list=ls())

getReference<-function(n,p,estimate,maxTimeLag)
{
  dataRef=NULL;dataX=NULL;
  
  #####--Generate data--#####
  
  for(times in 1:10)
  {
    Y=simulateEsimatedSepData(n,p,maxTimeLag,estimate$covar_L,estimate$Sigma11_L,estimate$mu);
    Y=matrix(Y,n);
    f=get_f(Y,n);
    dataRef=cbind(dataRef,f);
  }
  ct=dim(dataRef)[2];
  ff=fbplot(dataRef,plot=F);
  nonOut=setdiff(1:ct,ff$outpoint);
  dataRef=dataRef[,nonOut];
  
  for(times in 1:5)
  {
    Y=simulateEsimatedSepData(n,p,maxTimeLag,estimate$covar_L,estimate$Sigma11_L,estimate$mu);
    Y=matrix(Y,n);
    f=get_f(Y,n);
    dataX=cbind(dataX,f);
  }
  ct=dim(dataX)[2];
  ff=fbplot(dataX,plot=F);
  nonOut=setdiff(1:ct,ff$outpoint);
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

f=get_f(Y,n);

estimate=sepCovEstimate(Y,n,blocksize);

refSet=getReference(n,p,estimate,blocksize);
dataX=refSet$dataX;
dataRef=refSet$dataRef;

rank.value=rankTest(dataX,f,dataRef);
print(rank.value)

############ compute critical value ###############
crit.data.path='critData/';
crit.result.path='critResult/';
critReplicate=100;
generateData4Crit(n,p,estimate,critReplicate,blocksize,crit.data.path);

for(i in 1:critReplicate)
{
  load(paste(crit.data.path,"/",i,".Rdata",sep=''));
  Y=matrix(Y,n);
  f=get_f(Y,n);
  
  estimate=sepCovEstimate(Y,n,blocksize);
  
  refSet=getReference(n,p,estimate,blocksize);
  dataX=refSet$dataX;
  dataRef=refSet$dataRef;
  
  rank.value=rankTest(dataX,f,dataRef);
  
  write.table(rank.value,file=paste(crit.result.path,"/",i,".csv",sep=''),col.names = F, row.names = F);
}
