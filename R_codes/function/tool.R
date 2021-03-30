library(fields)

simulateSepData<-function(size.dim,p,blocksize,beta)
{  
  n=size.dim*size.dim;
  data=integer(n*p);
  Coor=matrix(0,n,2);
  Coor[,1]=rep(1:size.dim,size.dim);
  Coor[,2]=kronecker(1:size.dim,rep(1,size.dim));
  
  distMatrix=rdist(Coor);
  
  ##-------- covariance matrix
  Sigma11=matrix(0,n*blocksize,n*blocksize);
  for(space1 in 1:n)
  {
    for(space2 in 1:n)
    {
      hlag=distMatrix[space1,space2];
      for(time1 in 1:blocksize)
      {
        for(time2 in 1:blocksize)
        {
          ulag=time2-time1;
          index1=(time1-1)*n+space1;
          index2=(time2-1)*n+space2;      
          Sigma11[index1,index2]=Gneiting(h=hlag/size.dim,u=ulag/5,a=1,c=1,beta=beta)
        }
      }
    }
  }
  
  Sigma12=matrix(0,n*blocksize,n*blocksize);
  for(space1 in 1:n)
  {
    for(space2 in 1:n)
    {
      hlag=distMatrix[space1,space2];
      for(time1 in 1:blocksize)
      {
        for(time2 in 1:blocksize+blocksize)
        {
          ulag=time2-time1;
          index1=(time1-1)*n+space1;
          index2=(time2-1-blocksize)*n+space2;      
          Sigma12[index1,index2]=Gneiting(h=hlag/size.dim,u=ulag/5,a=1,c=1,beta=beta)
        }
      }
    }
  }
  
  
  mu=Sigma12%*%solve(Sigma11);
  covar=Sigma11-mu%*%t(Sigma12);
  L=chol(covar);
  
  L11=chol(Sigma11);
  Z=rnorm(n*blocksize);
  data[1:(n*blocksize)]=t(L11)%*%Z;
  
  for(block in 0:(floor(p/blocksize)-2))
  {
    Z=rnorm(n*blocksize);
    tmp=t(L)%*%Z;
    data[(block+1)*blocksize*n+1:(n*blocksize)]=tmp+mu%*%data[block*blocksize*n+1:(n*blocksize)];
  }
  
  return(data);
}

simulateDataCH<-function(size.dim,p,blocksize,beta)
{  
  data=integer(n*floor(p/blocksize)*blocksize);
  
  if(beta==0) func=CnHsep else func=CnH;
  
  n=size.dim*size.dim;
  data=integer(n*p);
  Coor=matrix(0,n,2);
  Coor[,1]=rep(1:size.dim,size.dim);
  Coor[,2]=kronecker(1:size.dim,rep(1,size.dim));
  
  distMatrix=rdist(Coor);
  
  ##-------- covariance matrix
  Sigma11=matrix(0,n*blocksize,n*blocksize);
  for(space1 in 1:n)
  {
    for(space2 in 1:n)
    {
      hlag=distMatrix[space1,space2];
      for(time1 in 1:blocksize)
      {
        for(time2 in 1:blocksize)
        {
          ulag=time2-time1;
          index1=(time1-1)*n+space1;
          index2=(time2-1)*n+space2;      
          Sigma11[index1,index2]=func(h=hlag/size.dim,u=ulag/2)
        }
      }
    }
  }
  
  Sigma12=matrix(0,n*blocksize,n*blocksize);
  for(space1 in 1:n)
  {
    for(space2 in 1:n)
    {
      hlag=distMatrix[space1,space2];
      for(time1 in 1:blocksize)
      {
        for(time2 in 1:blocksize+blocksize)
        {
          ulag=time2-time1;
          index1=(time1-1)*n+space1;
          index2=(time2-1-blocksize)*n+space2;      
          Sigma12[index1,index2]=func(h=hlag/size.dim,u=ulag/2)
        }
      }
    }
  }
  
  mu=Sigma12%*%solve(Sigma11);
  covar=Sigma11-mu%*%t(Sigma12);
  L=chol(covar);
  
  L11=chol(Sigma11);
  Z=rnorm(n*blocksize);
  data[1:(n*blocksize)]=t(L11)%*%Z;
  
  for(block in 0:(floor(p/blocksize)-2))
  {
    Z=rnorm(n*blocksize);
    tmp=t(L)%*%Z;
    data[(block+1)*blocksize*n+1:(n*blocksize)]=tmp+mu%*%data[block*blocksize*n+1:(n*blocksize)];
  }
  
  return(data);
}


simulateSymmData<-function(size.dim,p,blocksize,lambda)
{  
  n=size.dim*size.dim;
  data=integer(n*p);
  Coor=matrix(0,n,2);
  Coor[,1]=kronecker(1:size.dim-1,rep(1,size.dim));
  Coor[,2]=rep(1:size.dim-1,size.dim);
  
  distMatrix=rdist(Coor);
  
  ##-------- covariance matrix
  Sigma11=matrix(0,n*blocksize,n*blocksize);
  for(space1 in 1:n)
  {
    for(space2 in 1:n)
    {
      hlag=distMatrix[space1,space2];
      for(time1 in 1:blocksize)
      {
        for(time2 in 1:blocksize)
        {
          ulag=time2-time1;
          h1=Coor[space2,1]-Coor[space1,1];
          index1=(time1-1)*n+space1;
          index2=(time2-1)*n+space2;      
          Sigma11[index1,index2]=GneitingSymm(h=hlag/size.dim,u=ulag/5,h1=h1/size.dim,lambda=lambda)
        }
      }
    }
  }
  
  Sigma12=matrix(0,n*blocksize,n*blocksize);
  for(space1 in 1:n)
  {
    for(space2 in 1:n)
    {
      hlag=distMatrix[space1,space2];
      for(time1 in 1:blocksize)
      {
        for(time2 in 1:blocksize+blocksize)
        {
          ulag=time2-time1;
          h1=Coor[space2,1]-Coor[space1,1];
          index1=(time1-1)*n+space1;
          index2=(time2-1-blocksize)*n+space2;      
          Sigma12[index1,index2]=GneitingSymm(h=hlag/size.dim,u=ulag/5,h1=h1/size.dim,lambda=lambda)
        }
      }
    }
  }
  
  mu=Sigma12%*%solve(Sigma11);
  covar=Sigma11-mu%*%t(Sigma12);
  L=chol(covar);
  
  L11=chol(Sigma11);
  Z=rnorm(n*blocksize);
  data[1:(n*blocksize)]=t(L11)%*%Z;
  
  for(block in 0:(floor(p/blocksize)-2))
  {
    Z=rnorm(n*blocksize);
    tmp=t(L)%*%Z;
    data[(block+1)*blocksize*n+1:(n*blocksize)]=tmp+mu%*%data[block*blocksize*n+1:(n*blocksize)];
  }
  
  return(data);
}

rankTest<-function(dataX,dataY,dataRef)
{
  n=dim(dataX)[2];
  m=dim(dataY)[2];
  r=dim(dataRef)[2];
  order=integer(n+m);
  for(i in 1:m)
  {
    sample=cbind(dataRef,dataY[,i]);
    result=fbplot(sample,plot=F,method="Both");
    order[i]=sum(result$depth[1:r]<=result$depth[r+1])
  }
  for(i in 1:n)
  {
    sample=cbind(dataRef,dataX[,i]);
    result=fbplot(sample,plot=F,method="Both");
    order[i+m]=sum(result$depth[1:r]<=result$depth[r+1])
  }
  rk=(rank(order)-1)/(n+m-1);
  W=mean(rk[1:m]);
}

get_f<-function(Y,n,lag.max=5)
{
  f=matrix(0,lag.max,n*(n-1));
  ct=0;
  for(space1 in 1:(n-1))
  {
    for(space2 in (space1+1):n)
    {
      ct=ct+1;
      tSeries1=Y[space1,];
      tSeries2=Y[space2,];
      
      res=ccf(tSeries1,tSeries2,type="covariance",lag.max=lag.max,plot=F);
      
      chu=res$acf[(lag.max+1):(2*lag.max+1)];
      part1=chu[-1]/chu[1];
      RevPart1=res$acf[lag.max:1]/chu[1];
      
      
      res=acf(tSeries1,lag.max=lag.max,type="covariance",plot=F);
      c0u=res$acf[1:(lag.max+1)];
      c0u=c0u[-1]/c0u[1];
      
      res=acf(tSeries2,lag.max=lag.max,type="covariance",plot=F);
      c0u_=res$acf[1:(lag.max+1)];
      c0u_=c0u_[-1]/c0u_[1];
      
      part2=(c0u+c0u_)/2;
      
      f[,ct]=part1-part2;
      ct=ct+1;
      f[,ct]=RevPart1-part2;
    }
  }
  f=f[,1:ct];
  ff=fbplot(f,plot=F);
  nonOut=setdiff(1:ct,ff$outpoint);
  f=f[,nonOut];
  return(f);
}

sepCovEstimate<-function(Y,n,blocksize)
{
  #####--Estimate Sigma11 Sigm12--#####
  sigma=var(as.vector(Y));
  spatio=integer(n);
  temporal=integer(blocksize*2);
  
  for(sp in 1:n)
  {
    tSeries=Y[sp,];
    res=acf(tSeries,lag.max=blocksize*2-1,type="covariance",plot=F,na.action = na.pass);
    temporal=temporal+as.numeric(res$acf);
  }
  temporal=temporal/n;
  spatio=var(t(Y),na.rm=T);
  
  temp=matrix(0,blocksize,blocksize);
  for(i in 2:blocksize-1)
    for(j in (i+1):blocksize)
      temp[i,j]=temporal[j-i+1];
  temp=temp+t(temp);
  diag(temp)=temporal[1];
  
  Sigma11=kronecker(temp,spatio)/sigma;
  
  temp=matrix(0,blocksize,blocksize);
  for(i in 1:blocksize)
    for(j in 1:blocksize)
      temp[i,j]=temporal[j-i+1+blocksize];
  
  Sigma12=kronecker(temp,spatio)/sigma;
  
  diag(Sigma11)=diag(Sigma11)+1e-4;
  mu=Sigma12%*%solve(Sigma11);
  covar=Sigma11-mu%*%t(Sigma12);
  
  covar_L=chol(covar);
  Sigma11_L=chol(Sigma11);
  
  return(list(mu=mu,covar_L=covar_L,Sigma11_L=Sigma11_L));
}

simulateEsimatedSepData<-function(n,p,maxTimeLag,covar_L,Sigma11_L,mu)
{  
  data=integer(n*floor(p/maxTimeLag)*maxTimeLag);
  
  Z=rnorm(n*maxTimeLag);
  data[1:(n*maxTimeLag)]=t(Sigma11_L)%*%Z;
  
  for(block in 0:(floor(p/maxTimeLag)-2))
  {
    Z=rnorm(n*maxTimeLag);
    tmp=t(covar_L)%*%Z;
    data[(block+1)*maxTimeLag*n+1:(n*maxTimeLag)]=tmp+mu%*%data[block*maxTimeLag*n+1:(n*maxTimeLag)];
  }
  
  return(data);
}









symmCovEstimate<-function(Y,n,blocksize)
{
  #####--Estimate Sigma11 Sigm12--#####
  Sigma11=matrix(0,n*blocksize,n*blocksize);
  for(space1 in 1:n)
  {
    for(space2 in space1:n)
    {
      tSeries1=Y[space1,];
      tSeries2=Y[space2,];
      res=ccf(tSeries1,tSeries2,lag.max=blocksize-1,type="covariance",plot=F);
      for(ulag in 1:blocksize-1)
        Sigma11[space1,ulag*n+space2]=(res[ulag]$acf[1]+res[-ulag]$acf[1])/2;
    }
  }
  for(n.id in 1:blocksize-1)
  {
    Sigma11[1:n,1:n+n.id*n]=Sigma11[1:n,1:n+n.id*n]+t(Sigma11[1:n,1:n+n.id*n]);
    diag(Sigma11[1:n,1:n+n.id*n])=diag(Sigma11[1:n,1:n+n.id*n])/2;
  }
  for(row in 2:blocksize-1)
    for(col in 1:blocksize-1)
    {
      if(col<row) 
        Sigma11[row*n+1:n,col*n+1:n]=Sigma11[col*n+1:n,row*n+1:n] else
          Sigma11[row*n+1:n,col*n+1:n]=Sigma11[(row-1)*n+1:n,(col-1)*n+1:n]
    }
  
  
  Sigma12=matrix(0,n*blocksize,n*blocksize);
  for(space1 in 1:n)
  {
    for(space2 in space1:n)
    {
      tSeries1=Y[space1,];
      tSeries2=Y[space2,];
      res=ccf(tSeries1,tSeries2,lag.max=blocksize+blocksize-1,type="covariance",plot=F);
      for(ulag in 1:blocksize-1+blocksize)
        Sigma12[space1,(ulag-blocksize)*n+space2]=(res[ulag]$acf[1]+res[-ulag]$acf[1])/2;
      for(ulag in 2:blocksize-1)
        Sigma12[(blocksize-ulag)*n+space1,space2]=(res[ulag]$acf[1]+res[-ulag]$acf[1])/2;
    }
  }
  
  for(n.id in 1:blocksize-1)
  {
    Sigma12[1:n,1:n+n.id*n]=Sigma12[1:n,1:n+n.id*n]+t(Sigma12[1:n,1:n+n.id*n]);
    diag(Sigma12[1:n,1:n+n.id*n])=diag(Sigma12[1:n,1:n+n.id*n])/2;
  }
  
  for(n.id in 2:blocksize-1)
  {
    Sigma12[n.id*n+1:n,1:n]=Sigma12[n.id*n+1:n,1:n]+t(Sigma12[n.id*n+1:n,1:n]);
    diag(Sigma12[n.id*n+1:n,1:n])=diag(Sigma12[n.id*n+1:n,1:n])/2;
  }
  
  for(row in 2:blocksize-1)
    for(col in 2:blocksize-1)
      Sigma12[row*n+1:n,col*n+1:n]=Sigma12[(row-1)*n+1:n,(col-1)*n+1:n]
  
  
  diag(Sigma11)=diag(Sigma11)+1e-4;
  mu=Sigma12%*%solve(Sigma11);
  covar=Sigma11-mu%*%t(Sigma12);
  
  covar_L=chol(covar);
  Sigma11_L=chol(Sigma11);
  
  return(list(mu=mu,covar_L=covar_L,Sigma11_L=Sigma11_L));
}


get_g<-function(Y,n,lag.max=5)
{
  g=matrix(0,lag.max,n*(n-1));
  ct=0;
  for(space1 in 0:(n-2))
  {
    for(space2 in (space1+1):(n-1))
    {
      ct=ct+1;
      tSeries1=Y[space1+1,];
      tSeries2=Y[space2+1,];
      
      res=ccf(tSeries2,tSeries1,lag.max=lag.max,type="covariance",plot=F);
      g[,ct]=res$acf[(lag.max+2):(2*lag.max+1)]-res$acf[lag.max:1];
    }
  }
  
  g=g[,1:ct];
  gg=fbplot(g,plot=F);
  nonOut=setdiff(1:ct,gg$outpoint);
  g=g[,nonOut];
  return(g);
}

