macc <-
function(dat,model.type=c("single","multilevel","twolevel"),method=c("HL","TS","HL-TS"),
               delta=NULL,interval=c(-0.90,0.90),tol=10e-4,max.itr=500,conf.level=0.95,
               optimizer=c("optimx","bobyqa","Nelder_Mead"),mix.pkg=c("nlme","lme4"),
               random.indep=TRUE,random.var.equal=FALSE,u.int=FALSE,Sigma.update=TRUE,
               var.constraint=TRUE,random.var.update=TRUE,logLik.type=c("logLik","HL"),
               error.indep=TRUE,error.var.equal=FALSE,
               sens.plot=FALSE,sens.interval=seq(-1,1,by=0.01),legend.pos="topright",
               xlab=expression(delta),ylab=expression(hat(AB)),
               cex.lab=1,cex.axis=1,lgd.cex=1,lgd.pt.cex=1,plot.delta0=TRUE,...)
{
  if(model.type[1]=="single")
  {
    if(is.null(delta))
    {
      delta<-0
    }
    run.time<-system.time(re<-cma.uni.delta(dat,delta=delta,conf.level=conf.level))
    if(sens.plot)
    {
      re.sens<-cma.uni.sens(dat,delta=sens.interval,conf.level=conf.level)
      cma.uni.plot(re.sens,re,delta=NULL,legend.pos=legend.pos,xlab=xlab,ylab=ylab,
                   cex.lab=cex.lab,cex.axis=cex.axis,lgd.cex=lgd.cex,lgd.pt.cex=1,plot.delta0=plot.delta0,...)
    }
  }else
    if(model.type[1]=="multilevel")
    {
      if(is.null(delta))
      {
        if(method[1]=="TS")
        {
          run.time1<-system.time(re1<-optimize(cma.uni.mix.dhl,interval=interval,dat=dat,tol=tol,max.itr=0,optimizer=optimizer,
                                               mix.pkg=mix.pkg,random.indep=random.indep,random.var.equal=random.var.equal,
                                               u.int=u.int,Sigma.update=Sigma.update,logLik.type=logLik.type,maximum=TRUE))
          run.time2<-system.time(re<-cma.uni.mix(dat,delta=re1$maximum,conf.level=conf.level,optimizer=optimizer,mix.pkg=mix.pkg,
                                                 random.indep=random.indep,random.var.equal=random.var.equal,u.int=u.int))
          
          run.time<-run.time1+run.time2
        }else
        {
          run.time1<-system.time(re1<-optimize(cma.uni.mix.dhl,interval=interval,dat=dat,tol=tol,max.itr=max.itr,optimizer=optimizer,
                                               mix.pkg=mix.pkg,random.indep=random.indep,random.var.equal=random.var.equal,
                                               u.int=u.int,Sigma.update=Sigma.update,var.constraint=var.constraint,
                                               random.var.update=random.var.update,logLik.type=logLik.type,maximum=TRUE))
          if(method[1]=="HL-TS")
          {
            run.time2<-system.time(re<-cma.uni.mix(dat,delta=re1$maximum,conf.level=conf.level,optimizer=optimizer,mix.pkg=mix.pkg,
                                                   random.indep=random.indep,random.var.equal=random.var.equal,u.int=u.int))
          }
          if(method[1]=="HL")
          {
            run.time2<-system.time(re<-cma.uni.mix.hl(dat,delta=re1$maximum,tol=tol,max.itr=max.itr,alpha=1-conf.level,
                                                      random.indep=random.indep,optimizer=optimizer,mix.pkg=mix.pkg,
                                                      random.var.equal=random.var.equal,u.int=u.int,
                                                      Sigma.update=Sigma.update,var.constraint=var.constraint,
                                                      random.var.update=random.var.update))
          }
          
          run.time<-run.time1+run.time2
        }
      }else
      {
        if(method[1]=="TS")
        {
          run.time<-system.time(re<-cma.uni.mix(dat,delta=delta,conf.level=conf.level,optimizer=optimizer,mix.pkg=mix.pkg,
                                                random.indep=random.indep,random.var.equal=random.var.equal,u.int=u.int))
        }
        if(method[1]=="HL")
        {
          run.time<-system.time(re<-cma.uni.mix.hl(dat,delta=delta,tol=tol,max.itr=max.itr,alpha=1-conf.level,random.indep=random.indep,
                                                   optimizer=optimizer,mix.pkg=mix.pkg,random.var.equal=random.var.equal,u.int=u.int,
                                                   Sigma.update=Sigma.update,var.constraint=var.constraint,
                                                   random.var.update=random.var.update))
        }
      }
    }else
      if(model.type[1]=="twolevel")
      {
        if(is.null(delta))
        {
          if(method[1]=="TS")
          {
            run.time1<-system.time(re1<-optimize(cma.delta.lm.HL,interval=interval,dat=dat,max.itr=0,tol=tol,
                                                 error.indep=error.indep,error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                 var.constraint=var.constraint,maximum=TRUE))
            run.time2<-system.time(re<-cma.delta.lm(dat,delta=re1$maximum,max.itr=0,tol=tol,error.indep=error.indep,
                                                    error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                    var.constraint=var.constraint))
            
            run.time<-run.time1+run.time2
          }else
          {
            run.time1<-system.time(re1<-optimize(cma.delta.lm.HL,interval=interval,dat=dat,max.itr=max.itr,tol=tol,
                                                 error.indep=error.indep,error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                 var.constraint=var.constraint,maximum=TRUE))
            if(method[1]=="HL")
            {
              run.time2<-system.time(re<-cma.delta.lm(dat,delta=re1$maximum,max.itr=max.itr,tol=tol,error.indep=error.indep,
                                                      error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                      var.constraint=var.constraint))
            }
            if(method[1]=="HL-TS")
            {
              run.time2<-system.time(re<-cma.delta.lm(dat,delta=re1$maximum,max.itr=0,tol=tol,error.indep=error.indep,
                                                      error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                      var.constraint=var.constraint))
            }
            
            run.time<-run.time1+run.time2
          }
        }else
        {
          if(method[1]=="TS")
          {
            run.time<-system.time(re<-cma.delta.lm(dat,delta=delta,max.itr=0,tol=tol,error.indep=error.indep,
                                                   error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                   var.constraint=var.constraint))
          }
          if(method[1]=="HL")
          {
            run.time<-system.time(re<-cma.delta.lm(dat,delta=delta,max.itr=max.itr,tol=tol,error.indep=error.indep,
                                                   error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                   var.constraint=var.constraint))
          }
        }
      }
  
  re$time<-run.time
  
  return(re)
  if(as.numeric(re$Var.comp[1])<1e-5)
  {
    warning("The variance of A's random effect is less than 1e-5.")
  }
  return(re)
  if(as.numeric(re$Var.comp[2])<1e-5)
  {
    warning("The variance of C's random effect is less than 1e-5.")
  }
  return(re)
  if(as.numeric(re$Var.comp[3])<1e-5)
  {
    warning("The variance of B's random effect is less than 1e-5.")
  }
}
