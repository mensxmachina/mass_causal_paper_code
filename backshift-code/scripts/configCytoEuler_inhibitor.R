onEuler <- FALSE
localMC <- TRUE

### setting
###############
### ADJUST ####
###############
celltype <- "cd8"
celltypeFile <-
  if(onEuler){
    "/cluster/scratch/heinzec/cytometry/BM/inhibitor-data/cd8+_ inhibitor_BM.rds"
  }else{
    "~/Code/flow-cytometry/_rds/BM/inhibitor-data/cd8+_ inhibitor_BM.rds"
  }


### select inhibitor
allInhibitors <-  c("Akti", "BTKi", "Crassin", "Dasatinib", "GDC-0941",
                    "Go69", "H89","IKKi","Imatinib","Jak1i","Jak2i","Jak3i",
                    "Lcki","Lestaurtinib","PP2","Rapamycin","Ruxolitinib",
                    "SB202","Sorafenib","SP6","Staurosporine","Streptonigrin",
                    "Sunitinib","Syki","Tofacitinib","U0126","VX680")

### select dosage
dosage <- c("D0", "D1", "D2", "D3", "D4", "D5", "D6", "D7")

### select donor
donor <- "M1"

### select vars / remove vars
allVarNames <-  c("Time","Cell_length","CD3","CD45","BC1","BC2","pNFkB",
                  "pp38","CD4","BC3","CD20","CD33","pStat5","CD123",
                  "pAkt","pStat1","pSHP2","pZap70","pStat3","BC4","CD14",
                  "pSlp76","BC5","pBtk","pPlcg2","pErk","BC6","pLat",
                  "IgM","pS6","HLA-DR","BC7","CD7","DNA-1","DNA-2")
removeVars <- c("Time", "EventNum", "DNA-2", "Cell_length", "IgM", "HLA-DR")
removeCDVars <- allVarNames[grep("CD", allVarNames)]
removeBCVars <- allVarNames[grep("BC", allVarNames)]
removeVars <- c(removeVars, removeBCVars, removeCDVars)

### backShift settings
myseed <- 2
EV <- 5
nsim <- 100
nrep <- 100
thres.pe <- 0.25
thres <- 0.75

CCNOptions <- list()
CCNOptions$method <- "backShift"
CCNOptions$options <- list(ev = EV, nsim = nsim, threshold = thres, covariance=TRUE)
CCNOptions$nrep <- nrep
CCNOptions$thres.pe <- thres.pe

### save results
###############
### ADJUST ####
###############
saveResultsPath <-
  if(onEuler){
    "/cluster/home/heinzec/causality/computationalcausality/packages/nonlinearICP/scripts/cytometry-results/inhibitor/cd8/"
  }else{
    "~/Code/flow-cytometry/test/"
  }


### select activators / remove activators
removeActs <- ".."
allActs <-  c("BCR","GCSF","GMCSF","IFNa","IFNg","IL12","IL2","IL3","LPS",
              "PMA-IONO","pVO4","Ref")
