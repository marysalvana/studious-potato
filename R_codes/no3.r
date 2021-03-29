source("cls.r")
library(mapproj)


###################################################
#### Read data
###################################################
lon<-read.table("loc.txt")[,1]#Longitude
lat<-read.table("loc.txt")[,2]#Latitude
A<-read.table("train.id.txt")[,]#training site number
XX.new<-read.table("no3.txt")[,]#observation
n<-79#number of sites
kk<-20 #number of year
XX.new<-matrix(XX.new,kk,n)

####################################################
#### Map locations on [0,1]^2
####################################################
loc1.proj<-mapproject(lon,lat,proj="bonne",param=30)
loc1<-cbind(loc1.proj$x,loc1.proj$y)
max.range<-max((range(loc1[,1])[2]-range(loc1[,1])[1]),(range(loc1[,2])[2]-range(loc1[,2])[1]))*1.05
x.l<-(range(loc1[,1])[1])
y.l<-(range(loc1[,2])[1])
loc1[,1]<-(loc1[,1]-1.2*x.l)/max.range
loc1[,2]<-(loc1[,2]-1.005*y.l)/max.range

dist1<-as.matrix(dist(loc1))

####################################################
#### sample covariance
####################################################
sample.cov<-function(XX.obs){
n.obs<-dim(XX.obs)[2]
Corr.obs<-matrix(NA,n.obs,n.obs)
n.pair<-NULL
for(i in 1:n.obs){
	for(j in i:n.obs){
		XX.temp<-as.matrix(XX.obs[,c(i,j)])
		XX.temp<-matrix(XX.temp[is.na(XX.temp[,1])==0&is.na(XX.temp[,2])==0,],ncol=2)
		if(dim(XX.temp)[1]==1)XX.temp<-matrix(NA,1,2)
		Corr.temp<-t(XX.temp)%*%XX.temp
		Corr.temp<-diag(diag(Corr.temp)^-0.5,dim(Corr.temp)[1],dim(Corr.temp)[1])%*%Corr.temp%*%diag(diag(Corr.temp)^-0.5,dim(Corr.temp)[1],dim(Corr.temp)[1])
		Corr.obs[c(i,j),c(i,j)]<-Corr.temp
		n.pair<-c(n.pair,dim(XX.temp)[1])
	}
}
n.each<-apply(1-is.na(XX.obs),2,sum)
XX.sqr<-XX.obs^2
XX.sqr[is.na(XX.sqr)==1]<-0
Var.em<-(apply(XX.sqr,2,sum)/n.each)
Corr.em<-Corr.obs
Corr.em[is.na(Corr.em)==1]<-0
Sig.em<-diag(Var.em^0.5)%*%Corr.em%*%diag(Var.em^0.5)
Sig.em[Sig.em==0]<-NA
return(list(Sig.em=Sig.em,n.pair=n.pair) )
}


####################################################
#### function for finding alpha_l in stationary process W_l
####################################################
find.alpha<-function(Adist,ASig.em,p){
A.length<-dim(Adist)[1]
alpha.all<-NULL
for(i in 1:A.length){
	for(j in i:A.length){
	
		if (is.na(ASig.em[i,j])==0 & ASig.em[i,j]>0) al.temp<-(-Adist[i,j]/log(ASig.em[i,j]/sqrt(ASig.em[i,i]*ASig.em[j,j])))
    else  al.temp<-NA
		alpha.all<-c(alpha.all,al.temp)
		
	}
}
alpha.all<-alpha.all[is.na(alpha.all)==0]
alpha.set<-as.numeric(quantile(alpha.all,p))
return(alpha.set)
}



#####################################################
#### Gaussian kernal
#####################################################
sd.set<-c(1:4)/8
gg<-NULL
loc.mu<-NULL

for(i in 1:length(sd.set)){
vari<-sd.set[i]^2
loc.mui<-NULL
mu.m<-sd.set[i]
hex.nx<-floor((1+1.5*mu.m)/mu.m)
hex.ny<-floor((1+mu.m*sqrt(3)/2)/(mu.m*sqrt(3)/2))
for(k in -1:(hex.ny)){
	for(j in -1:(hex.nx-k%%2)){
	loc.jk<-c(mu.m*(j+(k%%2)/2),mu.m*k*sqrt(3)/2)
	loc.mui<-rbind(loc.mui,c(sd.set[i],loc.jk))
	temp.g<-diag((loc1[,1:2]-rep(loc.jk,each=n))%*%t(loc1[,1:2]-rep(loc.jk,each=n)))
 	temp.g<-exp(-0.5*temp.g/vari)
	gg<-rbind(gg,temp.g)
	}
	}
	loc.mu<-rbind(loc.mu,loc.mui)
}


####################################################
#### Gaussian kernal matrix for Lasso
####################################################
AM<-gg[,A]
k<-dim(AM)[1]
a<-dim(AM)[2]
si.fix<-NULL
si.temp1<-matrix(as.numeric(outer(c(1:a),c(1:a),"==")),a,a)
si.temp2<-NULL
for(s in 1:a){
	si.temp2<-c(si.temp2,si.temp1[s:a,s])

}
si.fix<-cbind(si.fix,si.temp2)

for(i in 1:k){
	si.temp1<-outer(AM[i,],AM[i,],"*")
	si.temp2<-NULL
for(s in 1:a){
	si.temp2<-c(si.temp2,si.temp1[s:a,s])
}
si.fix<-cbind(si.fix,si.temp2)
}
si.dim<-dim(si.fix)
si.fix<-matrix(as.numeric(si.fix),si.dim[1],si.dim[2])



##################################################
#### kriging
##################################################
mu.loc<-rep(0,n)
var.loc<-rep(0,n)
for(i in 1:n){
temp<-XX.new[,i]
na.id<-is.na(temp)
mu.loc[i]<-mean(temp[na.id==0])
var.loc[i]<-var(temp[na.id==0])
}

lm.var<-lm(log(var.loc[A])~log(mu.loc[A]))
n.na<-1-is.na(XX.new[,A])
n.obs.loc<-apply(n.na,2,sum)
ff<-exp(lm.var$coefficients[1])*(mu.loc[A]^(lm.var$coefficients[2]))/n.obs.loc

fun1<-function(para,input){
	cat(para,"\n")
	n<-length(input$y.mean)
	Sigma<- para[2]*exp(-input$dist.matrix/para[1])+diag(input$y.var)
	likeli<- log(det(Sigma))+ matrix(input$y.mean,1,n)%*%solve(Sigma)%*%matrix(input$y.mean,n,1)
	likeli
}

Adist<-dist1[A,A]
test<-optim( c(0.01,0.1), fun1, method="L-BFGS-B",lower=c(0,0),upper=c(Inf,Inf),input=list(y.mean=mu.loc[A],y.var=ff,dist.matrix=Adist))
k.s<-test$par[2]
k.r<-test$par[1]
cov.exp<-k.s*exp(-dist1/k.r)
predict.y<-cov.exp[,A]%*%solve(cov.exp[A,A]+diag(ff))%*%mu.loc[A]


###################################################
#### finding  alpha_l for W_l
###################################################
XX.obs<-XX.new[,A]-rep(predict.y[A],each=kk)
Sig.em.all<-sample.cov(XX.obs)
Sig.em<-Sig.em.all$Sig.em
n.pair<-Sig.em.all$n.pair
p<-c(1:20)/20-1/40
alpha.set<-find.alpha(dist1[A,A],Sig.em,p)


###################################################
#### matrix for lasso
###################################################
si<-NULL
for(i in 1:length(alpha.set)){
	exp.temp<-NULL
	for(s in 1:length(A)){
			temp<-exp(-abs(Adist[s:a,s])/alpha.set[i])
			exp.temp<-c(exp.temp,temp)
		}
	si<-cbind(si,exp.temp)
	}
si.dim<-dim(si)
si<-matrix(as.numeric(si),si.dim[1],si.dim[2])
si<-cbind(si.fix,si)



AY.lars<-NULL
loc.y<-NULL
for(i in 1:length(A)){
		AY.lars<-c(AY.lars,Sig.em[(i:a),i])
		loc.y<-rbind(loc.y,cbind(rep(i,(a-i+1)),c(i:a)))
	}
id.train<-c(1:length(AY.lars))[is.na(AY.lars)==0]



###############################################
#### CLS
###############################################
cat(date())
si.re<-si[id.train,]
BY.lars<-cbind(loc.y[id.train,],AY.lars[id.train])
test1<-larspositive(si.re,AY.lars[id.train],type="lasso")
temp1<-cv.lars.positive(si.re,BY.lars)
s.temp<-temp1$fraction[which.min(temp1$cv)]
pre1<-predict.larspositive(test1,s=s.temp,type="coefficient", mode="fraction")$coefficients
cat(date())

#The estimates might be different from those reported in the paper due to the CV method 
#randomly partition the points into K-folds. 
#If you want to get the same estimates as those in our paper, 
#please set s.temp<-0.5050505.


pre1[1]  #parameter estimate:error
M.n<-dim(gg)[1]
est.basis<-cbind(loc.mu,pre1[2:(M.n+1)])
est.basis[est.basis[,4]>0,]   #parameter estimate:basis function
est.stat<-cbind(alpha.set,pre1[-c(1:(M.n+1))])
est.stat[est.stat[,2]>0,] #parameter estimate: stationary component

####est. Sigma by CLS
M<-gg
M.n<-dim(gg)[1]
Sig.g<-(t(M)%*%diag(pre1[(2:(M.n+1))])%*%M)+diag(pre1[1],n)

Sig.exp<-matrix(0,n,n)
for(i in 1:length(alpha.set)){
	Sig.expt<-exp(-dist1/alpha.set[i])*pre1[-c(1:(M.n+1))][i]
	Sig.exp<-Sig.exp+Sig.expt
  }
Sig.CLS<-Sig.g+Sig.exp


###############################################
#### CWLS
###############################################

#compute BPHI
S<-Sig.CLS[A,A]
S.tensor<-tensor(S,S)

n.Knn<-dim(S)[1]
Knn<-Knn.func(n.Knn)

PHI.sec<-NULL       #indicator,only S_ij,j>i are used
ii.index<-NULL
  for (i in 1:n.Knn){
	 temp<-c(i:n.Knn)+n.Knn*(i-1)
	 PHI.sec<-c(PHI.sec,temp)
	 ii.index<-c(ii.index,temp[1])
  }
  

BPHI<-(Knn%*%S.tensor)
BPHI<-BPHI[PHI.sec,PHI.sec][id.train,id.train]+S.tensor[PHI.sec,PHI.sec][id.train,id.train]
p.n<-length(PHI.sec)
n.ijkl<-n.ijkl.fast(XX.new[,A]) # the number of the times that sites (s_i,s_j,s_k,s_l) all have observations,i,j,k,l=1,2,...,n
T.ijkl<-matrix(n.ijkl,p.n,p.n,byrow=T)/outer(n.pair,n.pair)
BPHI<-BPHI*T.ijkl[id.train,id.train]
R.m12<-diag(diag(BPHI)^(-0.5))



#estimation by positive lasso
si.re<-si[id.train,]
AY.re<-AY.lars[id.train]
BY.lars<-cbind(loc.y[id.train,],AY.re)
test1<-larspositive(si.re,AY.re,BPHI,is.CWLS=TRUE,type="lasso")
temp1<-cv.lars.positive(si.re,BY.lars,BPHI,is.CWLS=TRUE)
s.temp<-temp1$fraction[which.min(temp1$cv)]
pre1<-predict.larspositive(test1,s=s.temp,type="coefficient", mode="fraction")$coefficients

#The estimates might be different from those reported in the paper due to the CV method 
#randomly partition the points into K-folds. 
#If you want to get the same estimates as those in our paper, 
#please set s.temp<-0.7171717.

pre1[1]  #parameter estimate:error
M.n<-dim(gg)[1]
est.basis<-cbind(loc.mu,pre1[2:(M.n+1)])
est.basis[est.basis[,4]>0,]   #parameter estimate:basis function
est.stat<-cbind(alpha.set,pre1[-c(1:(M.n+1))])
est.stat[est.stat[,2]>0,] #parameter estimate: stationary component




####est. Sigma by CWLS
M<-gg
M.n<-dim(gg)[1]
Sig.g<-(t(M)%*%diag(pre1[(2:(M.n+1))])%*%M)+diag(pre1[1],n)

Sig.exp<-matrix(0,n,n)
for(i in 1:length(alpha.set)){
	Sig.expt<-exp(-dist1/alpha.set[i])*pre1[-c(1:(M.n+1))][i]
	Sig.exp<-Sig.exp+Sig.expt
  }
Sig.CWLS<-Sig.g+Sig.exp







