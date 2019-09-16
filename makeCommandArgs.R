# make the command arguments file for runAllSPDE.R
makeSpdeCommandArgs = function(tausqVec=c(0, 0.1^2), gammaVec=c(0, -1), margVarVec=c(0, 0.15^2), 
                               includeClustEffectVec=c(FALSE, TRUE), urbanEffectVec=c(FALSE, TRUE)) {
  spdeCommandArgs = list()
  i = 1
  for(i1 in 1:length(tausqVec)) {
    tausq=tausqVec[i1]
    
    for(i2 in 1:length(gammaVec)) {
      gamma = gammaVec[i2]
      
      for(i3 in 1:length(margVarVec)) {
        margVar = margVarVec[i3]
        
        # check to make sure this dataset even exists
        included = c(margVar != 0, gamma != 0, tausq != 0)
        if(!included[1]) {
          if(any(included[2:3]))
            next
        }
        if(!included[2]) {
          if(included[3])
            next
        }
        
        for(i4 in 1:length(includeClustEffectVec)) {
          includeClustEffect = includeClustEffectVec[i4]
          
          for(i5 in 1:length(urbanEffectVec)) {
            urbanEffect = urbanEffectVec[i5]
            
            # don't bother with the test argument, since most of the time here is spent in the aggregation 
            # over pixels and enumeration areas rather than on the fitting
            spdeCommandArgs[[i]] = list(tausq=tausq, gamma=gamma, margVar=margVar, 
                                        includeClustEffect=includeClustEffect, urbanEffect=urbanEffect)
            i = i + 1
          }
        }
      }
    }
  }
  
  save(spdeCommandArgs, file="spdeCommandArgs.RData")
}

# tausq=0.1^2, big=TRUE, test
# make the command arguments file for runDirectAll.R
makeDirectCommandArgs = function(tausqVec=c(0, 0.1^2), gammaVec=c(0, -1), margVarVec=c(0, 0.15^2), 
                                 bigVec=c(FALSE, TRUE), testVec=c(FALSE)) {
  directCommandArgs = list()
  i = 1
  for(i1 in 1:length(tausqVec)) {
    tausq=tausqVec[i1]
    
    for(i2 in 1:length(gammaVec)) {
      gamma = gammaVec[i2]
      
      for(i3 in 1:length(margVarVec)) {
        margVar = margVarVec[i3]
        
        # check to make sure this dataset even exists
        included = c(margVar != 0, gamma != 0, tausq != 0)
        if(!included[1]) {
          if(any(included[2:3]))
            next
        }
        if(!included[2]) {
          if(included[3])
            next
        }
        
        for(i4 in 1:length(bigVec)) {
          big = bigVec[i4]
          
          for(i5 in 1:length(testVec)) {
            test = testVec[i5]
            
            directCommandArgs[[i]] = list(tausq=tausq, gamma=gamma, margVar=margVar, 
                                          big=big, test=test)
            i = i + 1
          }
        }
      }
    }
  }
  
  save(directCommandArgs, file="directCommandArgs.RData")
}

# make the command arguments file for runMercerAll.R
makeMercerCommandArgs = function(tausqVec=c(0, 0.1^2), gammaVec=c(0, -1), margVarVec=c(0, 0.15^2), 
                                 testVec=c(FALSE)) {
  mercerCommandArgs = list()
  i = 1
  for(i1 in 1:length(tausqVec)) {
    tausq=tausqVec[i1]
    
    for(i2 in 1:length(gammaVec)) {
      gamma = gammaVec[i2]
      
      for(i3 in 1:length(margVarVec)) {
        margVar = margVarVec[i3]
        
        # check to make sure this dataset even exists
        included = c(margVar != 0, gamma != 0, tausq != 0)
        if(!included[1]) {
          if(any(included[2:3]))
            next
        }
        if(!included[2]) {
          if(included[3])
            next
        }
        
        for(i4 in 1:length(testVec)) {
          test = testVec[i4]
          
          mercerCommandArgs[[i]] = list(tausq=tausq, gamma=gamma, margVar=margVar, test=test)
          i = i + 1
        }
      }
    }
  }
  
  save(mercerCommandArgs, file="mercerCommandArgs.RData")
}

# make the command arguments file for runBYM2All.R
makeBYM2CommandArgs = function(tausqVec=c(0, 0.1^2), gammaVec=c(0, -1), margVarVec=c(0, 0.15^2), 
                               testVec=c(FALSE), includeUrbanRuralVec=c(FALSE, TRUE), 
                               includeClusterVec=c(FALSE, TRUE), aggregateByPopulationVec=c(FALSE, TRUE)) {
  bym2CommandArgs = list()
  i = 1
  for(i1 in 1:length(tausqVec)) {
    tausq=tausqVec[i1]
    
    for(i2 in 1:length(gammaVec)) {
      gamma = gammaVec[i2]
      
      for(i3 in 1:length(margVarVec)) {
        margVar = margVarVec[i3]
        
        # check to make sure this dataset even exists
        included = c(margVar != 0, gamma != 0, tausq != 0)
        if(!included[1]) {
          if(any(included[2:3]))
            next
        }
        if(!included[2]) {
          if(included[3])
            next
        }
        
        for(i4 in 1:length(testVec)) {
          test = testVec[i4]
          
          for(i5 in 1:length(includeUrbanRuralVec)) {
            includeUrbanRural = includeUrbanRuralVec[i5]
            
            for(i6 in 1:length(includeClusterVec)) {
              includeCluster = includeClusterVec[i6]
              
              for(i7 in 1:length(aggregateByPopulationVec)) {
                aggregateByPopulation = aggregateByPopulationVec[i7]
                
                bym2CommandArgs[[i]] = list(tausq=tausq, gamma=gamma, margVar=margVar, test=test, 
                                            includeUrbanRural=includeUrbanRural, includeCluster=includeCluster, 
                                            aggregateByPopulation=aggregateByPopulation)
                i = i + 1
              }
            }
          }
        }
      }
    }
  }
  
  save(bym2CommandArgs, file="bym2CommandArgs.RData")
}

makeCompareModelArgs = function(tausqVec=c(0, 0.1^2), gammaVec=c(0, -1), margVarVec=c(0, 0.15^2), 
                                resultTypeVec=c("county"), testVec=c(FALSE), 
                                samplingVec=c("SRS", "oversamp"), recomputeTruth=TRUE, modelsIList=list(1:2, 1:21), 
                                produceFigures=FALSE, bigVec=c(FALSE, TRUE), printIEvery=50, 
                                maxDataSets=NULL, nsim=10, saveResults=TRUE, loadResults=FALSE, 
                                xtable.args=list(digits=c(0, 2, 2, 2, 2, 1, 2), display=rep("f", 7), auto=TRUE), 
                                tableFormat=c("2", "1"), colScale=c(10^4, 10^5, 100^2, 10^3, 100, 100), 
                                colUnits=c(" ($\\times 10^{-4}$)", " ($\\times 10^{-5}$)", " ($\\times 10^{-4}$)", 
                                           " ($\\times 10^{-3}$)", " ($\\times 10^{-2}$)", " ($\\times 10^{-2}$)"), 
                                colDigits=c(2, 2, 2, 2, 1, 2)) {
  
  compareModelCommandArgs = list()
  i = 1
  
  # include modelsIList on the outside so that we can run a sequence of indices all at once 
  # corresponding to a given set of models
  for(i7 in 1:length(modelsIList)) {
    modelsI = modelsIList[[i7]]
    
    for(i1 in 1:length(tausqVec)) {
      tausq=tausqVec[i1]
      
      for(i2 in 1:length(gammaVec)) {
        gamma = gammaVec[i2]
        
        for(i3 in 1:length(margVarVec)) {
          margVar = margVarVec[i3]
          
          # check to make sure this dataset even exists
          included = c(margVar != 0, gamma != 0, tausq != 0)
          if(!included[1]) {
            if(any(included[2:3]))
              next
          }
          if(!included[2]) {
            if(included[3])
              next
          }
          
          for(i4 in 1:length(testVec)) {
            test = testVec[i4]
            
            for(i5 in 1:length(resultTypeVec)) {
              resultType = resultTypeVec[i5]
              
              for(i6 in 1:length(samplingVec)) {
                sampling = samplingVec[i6]
                
                for(i8 in 1:length(bigVec)) {
                  big = bigVec[i8]
                  
                  # only use the big data sets for the naive and direct models
                  if(!identical(modelsI, 1:2) && big)
                    next
                  if(identical(modelsI, 1:2) && !big)
                    next
                  
                  compareModelCommandArgs[[i]] = list(tausq=tausq, gamma=gamma, margVar=margVar, test=test, 
                                                      resultType=resultType, sampling=sampling, 
                                                      recomputeTruth=recomputeTruth, modelsI=modelsI, 
                                                      produceFigures=produceFigures, big=big, printIEvery=printIEvery, 
                                                      maxDataSets=maxDataSets, nsim=nsim, saveResults=saveResults, loadResults=loadResults, 
                                                      xtable.args=xtable.args, tableFormat=tableFormat, colScale=colScale, 
                                                      colUnits=colUnits, colDigits=colDigits)
                  i = i + 1
                }
              }
            }
          }
        }
      }
    }
  }
  
  save(compareModelCommandArgs, file="compareModelCommandArgs.RData")
}