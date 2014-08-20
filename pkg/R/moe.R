moeridge <- function(y,x,lambda.ridge=0.1*log(nrow(x)),
              initial=list(alpha=rep(0.5,(ncol(x)+1)),
                            beta1=rep(0,ncol(x)),
                            beta2=rep(0,ncol(x)),
                            sigma=0.88),
                conv.eps=1e-8, maxiter=list(total=2500,em=15))
{
NCOMP=2


constant=c(nrow(x),ncol(x),2,3,
           (ncol(x)+1),NCOMP,maxiter$total,maxiter$em)

result=.C("Ridge_MOE_EM",PACKAGE="mixvarselect",as.double(y),as.double(as.vector(x)),as.integer(constant), 
as.double(initial$alpha[1]), as.double(initial$alpha[2:(ncol(x)+1)]), 
as.double(initial$beta1),as.double(initial$beta2), as.double(initial$sigma), 
  as.double(lambda.ridge), 
as.double(conv.eps), alpha=as.double(as.double(rep(0,(NCOMP-1)*ncol(x)))),
beta=as.double(rep(0,NCOMP*ncol(x))))

return(list(alpha=matrix(result$alpha,nrow=1),beta=matrix(result$beta,nrow=2)))
}

