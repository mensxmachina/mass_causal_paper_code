################## read working directory from arguments
args<-commandArgs(TRUE)
if(length(args) == 0){
  ###############
  ### ADJUST ####
  ###############
  config.file <- "configCytoEuler_inhibitor.R"
}else{
  print(args)
  numberOfArgs <- length(args)
  wd <- sub("-","",args[numberOfArgs])
  cat(paste("Setting working directory to", wd, "\n"))
  setwd(wd)
  config.file <- sub("-","",args[numberOfArgs-1])
  cat(paste("Reading inputs from ", config.file, "\n"))
}

################## source config file
source(config.file, echo = TRUE)
require(backShift)
require(sfsmisc)
require(parallel)
require(analyzeMC)

if(onEuler){
  numCores <- strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
}else if(localMC){
  numCores <- detectCores()
}
cat(paste("Number of cores used: ", numCores))

set.seed(myseed)

combinationsInhAct <- expand.grid(allActs, allInhibitors)
combinationsInhAct <- lapply(as.list(1:dim(combinationsInhAct)[1]), function(x) combinationsInhAct[x[1],])
combinationsInhAct <- combinationsInhAct[1:2]

results <-
    mclapply(combinationsInhAct, function(actInh){
        runCCNcytometryPar(dataPath = celltypeFile,
                          selInhibitor = actInh$Var2,
                          selDosage = dosage,
                          selDonor = donor,
                          seed = myseed,
                          CCNOptions = CCNOptions,
                          runID = paste(actInh$Var2, actInh$Var1,
                                        round(as.numeric(Sys.time()),0),
                                        sep = "_"),
                          saveResultsPath = saveResultsPath,
                          verbose = TRUE,
                          #removeActs = NULL,
                          keepActs = actInh$Var1,
                          removeVars = removeVars
                          # keepVars = NULL,
          )

    },
    mc.cores = numCores)
saveRDS(object = results, file = paste(saveResultsPath, celltype, "_res.rds", sep = ""))
createPlotsCytoInhibitor(results, verbose = TRUE)

