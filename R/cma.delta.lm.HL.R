cma.delta.lm.HL <-
function(dat,delta,max.itr=500,tol=1e-4,error.indep=FALSE,error.var.equal=FALSE,Sigma.update=FALSE,
                          var.constraint=FALSE)
{
  re<-cma.delta.lm(dat,delta,max.itr,tol,error.indep,error.var.equal,Sigma.update,var.constraint)
#   if(re$convergence==TRUE)
#   {
#     return(re$HL)
#   }else
#   {
#     return(-10^8)
#   }
  return(re$HL)
}
