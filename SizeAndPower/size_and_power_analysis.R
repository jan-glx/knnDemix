# Size and power behavior check--------------
library(data.table)
library(ggplot2)
library(knnDemix)

once <- function(n0=100,nm=100,n1=50,k=1,p=1,d=10,B=1000) {
  a<-replicate(B,{
    X0 <- matrix(rnorm(n0*p),ncol=p)
    Xm <- matrix(c(rnorm((nm-n1)*p),rnorm(n1*p)+d),ncol=p, byrow = TRUE)

    c(null=mixture.test(X0, Xm,alpha=(nm-n1)/nm, calc.CI=FALSE)$p.value,
      eins=mixture.test(X0, Xm,alpha=1, calc.CI=FALSE)$p.value)
  })
  a<-rowMeans(a<0.1)
  data.table(variable=c("size","power"),value=a, n0, nm, d=d, alpha=(nm-n1)/nm, k, p,id=runif(1))
}
#----------------
dt<-rbindlist(lapply(rep(c(1,4,16,64),3),function(k) once(k=k)))
ggplot(dt,aes(x=k,y=value,color=variable))+geom_point()+coord_cartesian(ylim=0:1)+geom_hline(yintercept=0.1,linetype="dashed")
ggsave("SizeAndPower/over_k.pdf")
#----------------
dt<-rbindlist(lapply(rep(c(0,0.5,1,1.5,2,3,4,8,16),3),function(d) once(d=d)))
ggplot(dt,aes(x=d,y=value,color=variable))+geom_point()+coord_cartesian(ylim=0:1)+geom_hline(yintercept=0.1,linetype="dashed")
ggsave("SizeAndPower/over_d.pdf")
#----------------
dt<-rbindlist(lapply(rep(2^(4:10),3),function(N) once(n0=N,n1=N/2,nm=N)))
ggplot(dt,aes(x=n0,y=value,color=variable))+geom_point()+coord_cartesian(ylim=0:1)+geom_hline(yintercept=0.1,linetype="dashed")
ggsave("SizeAndPower/over_N.pdf")
#----------------
dt<-rbindlist(lapply(rep(100*(1:9),3),function(N) once(n0=1000,n1=N,nm=1000)))
ggplot(dt,aes(x=alpha,y=value,color=variable))+geom_point()+coord_cartesian(ylim=0:1)+geom_hline(yintercept=0.1,linetype="dashed")
ggsave("SizeAndPower/over_alpha.pdf")
