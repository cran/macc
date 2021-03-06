cma.uni.mix.dhl <-
function(dat,delta,tol=10e-4,max.itr=50,alpha=0.05,optimizer=c("bobyqa","Nelder_Mead","optimx"),
                          mix.pkg=c("lme4","nlme"),random.indep=TRUE,random.var.equal=TRUE,u.int=FALSE,
                          Sigma.update=FALSE,var.constraint=FALSE,random.var.update=TRUE,logLik.type=c("logLik","HL"))
{
  if(max.itr==0)
  {
    if(logLik.type[1]=="HL")
    {
      return(cma.uni.mix(dat,delta,conf.level=1-alpha,optimizer=optimizer,mix.pkg=mix.pkg,random.indep=random.indep,
                         random.var.equal=random.var.equal,u.int=u.int)$HL)
    }else
    {
      return(cma.uni.mix(dat,delta,conf.level=1-alpha,optimizer=optimizer,mix.pkg=mix.pkg,random.indep=random.indep,
                         random.var.equal=random.var.equal,u.int=u.int)$logLik)
    }
  }else
  {
    return(cma.uni.mix.hl(dat,delta=delta,tol=tol,max.itr=max.itr,alpha=alpha,random.indep=random.indep,
                          optimizer=optimizer,mix.pkg=mix.pkg,random.var.equal=random.var.equal,u.int=u.int,
                          Sigma.update=Sigma.update,var.constraint=var.constraint,random.var.update=random.var.update)$HL)
  }
}
