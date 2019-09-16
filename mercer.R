# CAUTION: logit.est and var.est must be provided for admin1 in alphabetic order
# otherwise the graph file will not fit.
mercer_u1m = function(logit.est, var.est, graph.path){
  
  # our weighted estimates represent the outcome 
  # assumed to be normally distributed with the true mean and
  # fixed variance equal to var.est
  y = logit.est
  
  # our model includes an intercept, iid space and ICAR component
  # use hyperpriors as in Mercer et al. 
  formula = y ~ f(idx, model="iid", 
                  hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
                f(idx2, model="besag", graph="Kenyaadm1.graph", 
                  hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma")))
  
  ## generate dataset
  data = list(y=y, scale=1/var.est, idx=1:length(y), idx2=1:length(y))
  
  # be careful to fix the likelihood variance, since we have observation-specific
  # variances we need to use the argument scale and fix the overall precision to 1.
  result = inla(formula, family="gaussian", data=data, 
                control.family=list(hyper=list(prec=list(initial=0, fixed=TRUE))), 
                control.predictor=list(compute=TRUE, link=1, quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975)),
                scale = scale)
}

# modified version of the Mercer et al. model using the BYM2 model with joint PC prior
mercer_u1m2 = function(logit.est, var.est, graph.path, plotPriorPost=FALSE, previousResult=NULL, 
                       doValidation=FALSE, predCountyI=NULL, strictPrior=FALSE) {
  # remove observation at the given county if necessary for cross-validation
  if(!is.null(predCountyI))
    logit.est[predCountyI] = NA
    
  # our weighted estimates represent the outcome 
  # assumed to be normally distributed with the true mean and
  # fixed variance equal to var.est
  y = logit.est
  
  # our model includes an intercept, iid space and ICAR component
  # use BYM2 hyperpriors based on INLA rw1 and bym2 recommendations (this is 
  # probability scale rather than full real line, so shrink sigma slightly)
  if(!strictPrior)
    hyperList = list(param=c(1, 0.01), prior="pc.prec")
  else
    hyperList = list(param=c(.15, 0.01), prior="pc.prec")
  formula = y ~ f(idx, model="bym2", graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
                  hyper=list(prec=hyperList, 
                             phi=list(param=c(0.5, 2/3), prior="pc")))
  
  if(doValidation) {
    control.inla = list(strategy="laplace", int.strategy="grid", diff.logdens=4, npoints=21) 
  }
  else {
    control.inla = inla.set.control.inla.default()
  }
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousResult)) {
    # initialize the fitting process based on a previous optimum
    # modeControl$result = previousResult
    modeControl$theta = previousResult$mode$theta
    modeControl$x = previousResult$mode$x
    modeControl$restart = TRUE
  }
  
  ## generate dataset
  data = list(y=y, scale=1/var.est, idx=1:length(y))
  
  # be careful to fix the likelihood variance, since we have observation-specific
  # variances we need to use the argument scale and fix the overall precision to 1.
  result = inla(formula, family="gaussian", data=data, 
                control.family=list(hyper=list(prec=list(initial=0, fixed=TRUE))), 
                control.predictor=list(compute=TRUE, link=1),
                scale = scale, quantiles=c(0.1, 0.5, 0.9), 
                control.mode=modeControl, 
                control.compute=list(cpo=doValidation, dic=doValidation, waic=doValidation), 
                control.inla=control.inla)
  
  if(plotPriorPost) {
    maxX = max(result$marginals.hyperpar[[1]][,1]) * 1.1
    xs = seq(0, maxX, l=500)[-1]
    ys = inla.pc.dprec(xs, 1, 0.01)
    maxY = max(max(ys), max(result$marginals.hyperpar[[1]][,2]))
    plot(xs, ys, type="l", col="blue", xlab="Mercer (BYM2) Precision", ylab="Density", xlim=c(0, maxX), 
         ylim=c(0, maxY), main="BYM2 Precision Prior vs Posterior")
    lines(result$marginals.hyperpar[[1]], col="red")
    legend("topright", c("Prior", "Posterior"), col=c("blue", "red"), lty=1)
  }
  
  result
}