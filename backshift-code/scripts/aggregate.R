require(analyzeMC)
require(backShift)
require(reshape2)
require(Hmisc)
###############
### ADJUST ####
###############
resultPath <- "/Users/heinzec/Code/flow-cytometry/test/"
thresConv <- 75
saveResultsPath <- "/Users/heinzec/Code/flow-cytometry/test/"
takeout <- "cd8_res.rds"
folderName <- "effects"
woInhibitors <- c("")

allInhibitors <-  c("Akti", "BTKi", "Crassin", "Dasatinib", "GDC-0941",
                    "Go69", "H89","IKKi","Imatinib","Jak1i","Jak2i","Jak3i",
                    "Lcki","Lestaurtinib","PP2","Rapamycin","Ruxolitinib",
                    "SB202","Sorafenib","SP6","Staurosporine","Streptonigrin",
                    "Sunitinib","Syki","Tofacitinib","U0126","VX680")
allActs <-  c("BCR","GCSF","GMCSF","IFNa","IFNg","IL12","IL2","IL3","LPS",
              "PMA-IONO","pVO4","Ref")

files <- aggregateResults(resultPath, thresConv, takeout, woInhibitors,
                          allInhibitors, allActs, saveResultsPath, folderName)
