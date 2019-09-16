# library(haven)
library(survey)

# -- a few helpful functions -- #
expit<-function(x){
  exp(x)/(1+exp(x))
}

logit<-function(x){
  log(x/(1-x))
}

# function to extend dataset to binary form 
extendData <- function(clustDatRow, v001, divideWeight=TRUE, useNumChildrenDied=TRUE){
  
  # add extra columns for ageMonth, ageGrpD, v001, v002
  if(useNumChildrenDied) {
    clustDatRow$n = clustDatRow$numChildren
    clustDatRow$y = clustDatRow$died
  }
  n = clustDatRow$n
  # tmp = data.frame(clustDatRow[c(1, 6:16)])
  tmp = data.frame(clustDatRow[c(1, c(4, 6:ncol(clustDatRow)))])
  tmp$v001 = v001

  ageMonth = rep(0, n)
  ageGrpD = rep("[0,1)", n)
  v001 = rep(v001, n)
  # there is only one child and one mother per household.
  # All 25 households are sampled
  v002 = 1:n
  
  y = c(rep(0,n-clustDatRow$y), rep(1, clustDatRow$y))
  if(clustDatRow["urban"][1,1]){
    urbanRural = rep("urban", n)
  } else {
    urbanRural = rep("rural", n)
  }
  # admin1 = rep(clustDatRow$admin1, n)

  res = merge(data.frame(y, ageMonth, ageGrpD, v001, v002, urbanRural), tmp, by="v001")
  
  # the below line was commented out since each cluster only has one type of admin and urban level. 
  # The equivalent line has been added into the parent function
  # res$regionRural <- with(res, interaction(admin1, urbanRural), drop=TRUE)
  
  if(divideWeight)
    res$samplingWeight = res$samplingWeight / n
  return(res)
}

extendDataDat <- function(clustDatRow, v001, divideWeight=TRUE){
  
  # add extra columns for ageMonth, ageGrpD, v001, v002
  n = clustDatRow$n
  # the only things we need are admin1 and sampling weight, but we must get rid of 
  # urban, y, and the number of women since those will be recalculated
  # tmp = data.frame(clustDatRow[c(1, 6:16)])
  # tmp = data.frame(clustDatRow[c(1, c(4, 6:ncol(clustDatRow)))])
  tmp = data.frame(clustDatRow[c(1, c(4, 6:(ncol(clustDatRow) - 2)))])
  tmp$v001 = v001
  
  ageMonth = rep(0, n)
  # ageGrpD = rep("[0,1)", n)
  v001 = rep(v001, n)
  # there is only one child and one mother per household.
  # All 25 households are sampled
  v002 = 1:n
  
  y = c(rep(0,n-clustDatRow$y), rep(1, clustDatRow$y))
  if(clustDatRow["urban"][1,1]){
    urbanRural = rep("urban", n)
  } else {
    urbanRural = rep("rural", n)
  }
  # admin1 = rep(clustDatRow$admin1, n)
  
  # res = merge(data.frame(y, ageMonth, ageGrpD, v001, v002, urbanRural), tmp, by="v001")
  res = merge(data.frame(y, ageMonth, v001, v002, urbanRural), tmp, by="v001")
  
  # the below line was commented out since each cluster only has one type of admin and urban level. 
  # The equivalent line has been added into the parent function
  # res$regionRural <- with(res, interaction(admin1, urbanRural), drop=TRUE)
  
  if(divideWeight)
    res$samplingWeight = res$samplingWeight / n
  return(res)
}


# - a function that reads in a glm or svyglm - #
# - object and returns the estimate and SE - #
# - specifics in the supplementary materials - #
## This function takes care of the delta method
## to calculate the variance of u5m as a function
## of the age specific hazards, \beta_a .

get.est<-function(glm.ob){

  beta<-summary(glm.ob)$coef[,1]
 
  est <-expit(beta)
  var.est <- vcov(glm.ob)[1,1]
  
  # compute 80% CI intervals
  lower <- logit(est)+qnorm(c(0.9))*sqrt(var.est)
  upper <- logit(est)+qnorm(c(0.1))*sqrt(var.est)
  return(c(est,lower, upper,logit(est),var.est))
}

# -- a function to subset the design based on a region and time period -- #
# -- and then run the svyglm function and return the get.est() results -- #

## First line in function allows you to subset your data and ALSO the specified
## svydesign object into area (usually v024 variable in DHS) 
## and time (per5 is a variable we construct for the 5-year periods in the Stata step)
## Second line fits the survey-weighted glm

region.time.HT<-function(dataobj, svydesign, area){
  
  tmp<-subset(svydesign, (admin1==area))
  
  tt2 <- tryCatch(glmob<-svyglm(y.x~1,
                                design=tmp,family=quasibinomial, maxit=50), 
                  error=function(e) e, warning=function(w) w)
  
  if(is(tt2, "warning")){
    if(grepl("agegroups", tt2)){
      res <- get.est(glmob)
      res = c(res, 2)
    } else {
      res = c(rep(NA, 5), 3)
    }
    return(res)
  }
  if(is(tt2,"error")){
    res = c(rep(NA, 5), 1)
    return(res)
  } else {
    res <- get.est(glmob)
    res = c(res, 0)
    return(res)
  }
}

region.time.HTDat<-function(dataobj, svydesign, area, nationalEstimate){
  
  if(!nationalEstimate) {
    
    tmp<-subset(svydesign, (admin1==area))
    
    tt2 <- tryCatch(glmob<-svyglm(y~1,
                                  design=tmp,family=quasibinomial, maxit=50), 
                    error=function(e) e, warning=function(w) w)
  } else {
    thisUrban = area == 1
    tmp<-subset(svydesign, (urban==thisUrban))
    tt2 <- tryCatch(glmob<-svyglm(y~1,
                                  design=tmp,family=quasibinomial, maxit=50), 
                    error=function(e) e, warning=function(w) w)
  }
  
  if(is(tt2, "warning")){
    if(grepl("agegroups", tt2)){
      res <- get.est(glmob)
      res = c(res, 2)
    } else {
      res = c(rep(NA, 5), 3)
    }
    return(res)
  }
  if(is(tt2,"error")){
    res = c(rep(NA, 5), 1)
    return(res)
  } else {
    res <- get.est(glmob)
    res = c(res, 0)
    return(res)
  }
}



defineSurvey <- function(dat_obj, stratVar, useSamplingWeights=TRUE){
  
  options(survey.lonely.psu="adjust")
  
  # --- setting up a place to store results --- #
  regions <- sort(unique(dat_obj$admin1))
  regions_num  <- 1:length(regions)
  
  results<-data.frame(admin1=rep(regions,each=1))
  results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
  results$converge <- NA
  
  if(useSamplingWeights){
    dat_obj$wt <- dat_obj$samplingWeight
  } else {
    dat_obj$wt <- NULL
  }

  if(is.null(stratVar)){
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata =NULL,
                              weights=NULL, data=dat_obj)
  } else {
    ## not in all surveys does v022 contain the correct sampling strata
    ## Thus, the correct vector has to be provided externally
    dat_obj$strat <- stratVar
  
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata=~strat, nest=T, 
                              weights=~wt, data=dat_obj)
  }
  
  for(i in 1:nrow(results)){
    results[i, 2:7] <- region.time.HT(dataobj=dat_obj, svydesign=my.svydesign, 
                                      area=results$admin1[i])
  }
  return(results)
}

defineSurveyDat <- function(dat_obj, stratVar, useSamplingWeights=TRUE, nationalEstimate=FALSE, 
                             getContrast=nationalEstimate){
  
  options(survey.lonely.psu="adjust")
  
  # --- setting up a place to store results --- #
  regions <- sort(unique(dat_obj$admin1))
  regions_num  <- 1:length(regions)
  
  if(!nationalEstimate) {
    results<-data.frame(admin1=rep(regions,each=1))
    results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
    results$converge <- NA
  }
  else {
    results<-data.frame(urban=c(TRUE, FALSE))
    results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
    results$converge <- NA
  }
  
  if(useSamplingWeights){
    dat_obj$wt <- dat_obj$samplingWeight
  } else {
    dat_obj$wt <- NULL
  }
  
  if(is.null(stratVar)){
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata =NULL,
                              weights=NULL, data=dat_obj)
  } else {
    ## not in all surveys does v022 contain the correct sampling strata
    ## Thus, the correct vector has to be provided externally
    dat_obj$strat <- stratVar
    
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata=~strat, nest=T, 
                              weights=~wt, data=dat_obj)
  }
  
  for(i in 1:nrow(results)){
    if(!nationalEstimate) {
      results[i, 2:7] <- region.time.HTDat(dataobj=dat_obj, svydesign=my.svydesign, 
                                            area=results$admin1[i], nationalEstimate=nationalEstimate)
    }
    else {
      results[i, 2:7] <- region.time.HTDat(dataobj=dat_obj, svydesign=my.svydesign, 
                                            area=i, nationalEstimate=nationalEstimate)
    }
  }
  
  if(getContrast) {
    # out = svyby(~y, by = ~urban, design = svydesign, svymean)
    glmob<-svyglm(y~urban,
                  design=my.svydesign,family=quasibinomial, maxit=50)
    
    # get contrast mean and variance
    est = glmob$coefficients[2]
    urbanVar = vcov(glmob)[2,2]
    
    # get confidence interval
    lower = est + qnorm(0.025, sd=sqrt(urbanVar))
    upper = est + qnorm(0.975, sd=sqrt(urbanVar))
    contrastStats = list(est=est, sd=sqrt(urbanVar), lower95=lower, upper95=upper)
    return(list(results=results, contrastStats=contrastStats))
  } else {
    return(results)
  }
  
}

# Set dat_obj$admin1 to be something else for different kinds of aggregations
run_naive <- function(dat_obj){
  regions <- sort(unique(dat_obj$admin1))
  regions_num  <- 1:length(regions)
  
  results<-data.frame(admin1=rep(regions,each=1))
  results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
  results$converge <- NA
  
  for(i in 1:nrow(results)){
    my.glm <- glm(y.x~1, family=binomial, 
                  data=dat_obj, 
                  subset = admin1 == results$admin1[i] ) 
    # newdat = dat_obj[dat_obj$admin1==results$admin1[i], ]
    # my.glm2 <- glm(y.x~1, family=binomial, 
    #               data=newdat) 
    
    results[i, 2:7] <- c(get.est(my.glm),0)
  }
  return(results)
}

# running the analysis for the actual mortality dataset is slightly different
run_naiveDat <- function(dat_obj){
  regions <- sort(unique(dat_obj$admin1))
  regions_num  <- 1:length(regions)
  
  results<-data.frame(admin1=rep(regions,each=1))
  results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
  results$converge <- NA
  
  for(i in 1:nrow(results)){
    my.glm <- glm(y~1, family=binomial, 
                  data=dat_obj, 
                  subset = admin1 == results$admin1[i] ) 
    # newdat = dat_obj[dat_obj$admin1==results$admin1[i], ]
    # my.glm2 <- glm(y.x~1, family=binomial, 
    #               data=newdat) 
    
    results[i, 2:7] <- c(get.est(my.glm),0)
  }
  return(results)
}

