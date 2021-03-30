Cesare<-function(h,u,k2,k1=1,k3=1)
{
  part1=k1*exp(-h)*exp(-u);
  part2=k2*exp(-h);
  part3=k3*exp(-u)
  return(part1+part2+part3)/3;
}


CnH<-function(h,u)
{
  part1=abs(u)+1;
  part2=part1^2+h^2;
  return(part1/(part2^1.5));
}

CnHsep<-function(h,u)
{
  part1=(abs(u)+1)^2;
  part2=(1+h^2)^1.5;
  return(1/(part1*part2));
}


Gneiting<-function(h,u,a=0.5,c=1,beta)
{
  part1=a*abs(u)+1;
  denom=part1^(beta/2);
  part2=exp(-c*h/denom);
  return(part2/part1);
}

GneitingSymm<-function(h,u,h1,lambda)
{
  part1=1/(abs(u)+1);
  part2=exp(-h);

  term1=part1*part2;
  term1=term1*(1-lambda);
  
  part3=abs(h1-u)/2;
  part3=1-part3;
  part3=max(part3,0);
  term2=part3*lambda;
  
  return(term1+term2);
}
