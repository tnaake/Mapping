library (MetCirc)
setwd("/home/thomas/Documents/University/Master/MScArbeit/MSMS/")

##xsetMSMS <- xcmsSet(file = "CDFs/MS", method="centWave", ppm=20, snthresh=10, peakwidth=c(5,18))
##classes <- c("attMSWOS72", "obtMSWOS72", "clevMSWOS72", "quadMSWOS72", "x1027MSWOS72", "x57126MSWOS72")
## classes <- c("att20evWOS72", "att30evWOS72", "att40evWOS72", "attMSWOS72",
##             "obt20evWOS72", "obt30evWOS72", "obt40evWOS72", "obtMSWOS72",
##              "clev20evWOS72", "clev30evWOS72", "clev40evWOS72", "clevMSWOS72",
##              "quad20evWOS72", "quad30evWOS72", "quad40evWOS72", "quadMSWOS72",
##              "x102720evWOS72", "x102730evWOS72", "x102740evWOS72", "x1027MSWOS72",
##              "x5712620evWOS72", "x5712630evWOS72", "x5712640evWOS72", "x57126MSWOS72",
##              "att20evCon", "att30evCon", "att40evCon", "attMSCon",
##              "att20evMJ72", "att30evMJ72", "att40evMJ72", "attMSMJ72",
##              "obt20evCon", "obt30evCon", "obt40evCon", "obtMSCon",
##              "obt20evMJ72", "obt30evMJ72", "obt40evMJ72", "obtMSMJ72",
##              "cle20evCon", "cle30evCon", "cle40evCon", "cleMSCon",
##              "cle20evMJ72", "cle30evMJ72", "cle40evMJ72", "cleMSMJ72",
##              "quad20evCon", "quad30evCon", "quad40evCon", "quadMSCon",
##              "quad20evMJ72", "quad30evMJ72", "quad40evMJ72", "quadMSMJ72",
##              "x102720evCon", "x102730evCon", "x102740evCon", "x1027MSCon",
##              "x102720evMJ72", "x102730evMJ72", "x102740evMJ72", "x1027MSMJ72",
##              "x5712620evCon", "x5712630evCon", "x5712640evCon", "x57126MSCon",
##              "x5712620evMJ72", "x5712630evMJ72", "x5712640evMJ72", "x57126MSMJ72")
## add classes
##sampclass(xsetMSMS) <- classes
## 
##xset2MSMS <- group(xsetMSMS, method="density", minfrac=0.5, minsamp=1, bw=10, mzwid=0.01)
##xset3MSMS <- retcor(xset2MSMS, family= "s", plottype= "m", missing=1, extra=1, span=1)
##xset4MSMS <- group(xset3MSMS, method="density", mzwid=0.01, minfrac=0.5, 
##               minsamp=1, bw=5)
##xset5MSMS <- fillPeaks(xset4MSMS, method = "chrom")
##save(xsetMSMS, xset2MSMS, xset3MSMS, xset4MSMS, file = "MSMS_xcms.RData")
load("MSMS_xcms.RData")

## CAMERA
load("MSMS_CAMERA.RData")
##anMSMS <- xsAnnotate(xset5MSMS)
##anFMSMS <- groupFWHM(anMSMS, perfwhm = 0.6)
##anIMSMS <- findIsotopes(anFMSMS, mzabs=0.01)
##anICMSMS <- groupCorr(anIMSMS, cor_eic_th=0.75, graphMethod = "lpc")
##anFAMSMS <- findAdducts(anICMSMS, polarity="positive")
#peaklistMSMS <- getPeaklist(anFAMSMS)

##anMSMS2 <- xsAnnotate(xset4MSMS)
##anFMSMS2 <- groupFWHM(anMSMS2, perfwhm = 0.6)
##anIMSMS2 <- findIsotopes(anFMSMS2, mzabs = 0.01)
##anICMSMS2 <- groupCorr(anIMSMS2, cor_eic_th=0.75, graphMethod = "lpc")
##anFAMSMS2 <- findAdducts(anICMSMS2, polarity="positive")
peaklistMSMS2 <- getPeaklist(anFAMSMS2)
##save(anMSMS2, anFMSMS2, anIMSMS2, anICMSMS2, anFAMSMS2, peaklistMSMS2, file = "MSMS_CAMERA.RData")


## re-create xcms combined peaklist for WOS and MeJA (shell)
setwd("/home/thomas/Documents/University/Master/MScArbeit/Metabolic_profiling/WOS_MeJA_0h_72h/")
load("./metabolicProfiling.RData")
##xset <- xcmsSet(method="centWave", ppm=20, snthresh=10, peakwidth=c(5,18), BPPARAM = MulticoreParam(workers = 6))
##classes <- c(rep("att0", 5),  rep("att72WOS", 5),
##    rep("obt0", 5), rep("obt72WOS", 5),
##    rep("clev0", 5),  rep("clev72WOS", 5),
##    rep("quad0", 5), rep("quad72WOS", 5),
##    rep("x10270", 5),  rep("x1072WOS", 5),
##    rep("x570", 5), rep("x5772WOS", 5),
##    rep("att0", 5),  rep("att72MJ", 5),
##    rep("obt0", 5), rep("obt72MJ", 5),
##    rep("clev0", 5),  rep("clev72MJ", 5),
##    rep("quad0", 5), rep("quad72MJ", 5),
##    rep("x10270", 5),  rep("x1072MJ", 5),
##    rep("x570", 5), rep("x5772MJ", 5))
##sampclass(xset) <- classes
##xset2 <- group(xset, method="density", minfrac=0.5, minsamp=2, bw=10, mzwid=0.05) ## mzwid = 0.01
##xset3 <- retcor(xset2, family= "s", plottype= "m", missing=1, extra=1, span=1)
##xset4 <- group(xset3, method="density", bw=10, mzwid=0.05, minfrac=0.5, minsamp=2)
##xset5 <- fillPeaks(xset4, method = "chrom")
##cleanParallel()
##save(xset, xset2, xset3, xset4, xset5, file = "./metabolicProfiling.RData")

load("./CAMERA_complete.RData")
##an <- xsAnnotate(xset5)
##anF <- groupFWHM(an, perfwhm = 0.6)
##anI <- findIsotopes(anF, mzabs=0.01)
##anIC <- groupCorr(anI, cor_eic_th=0.75, graphMethod = "lpc")
##anFA <- findAdducts(anIC, polarity="positive")
peaklist <- getPeaklist(anFA)
# ## get colnames of samples (min and max ind)
colMin <- which(colnames(peaklist) == "X001_BA4_01_28001")
colMax <- which(colnames(peaklist) == "X180_RD2_01_28200")
cols <- colMin:colMax
peaklistConc <- peaklist
peaklistConc[, cols] <- apply(peaklist[, cols], 2, FUN = function(x) (x / quantile(x, 0.75)))
peaklist[, cols] <- apply(peaklist[, cols], 2, FUN = function(x) (x / quantile(x, 0.75) + 1))
# ## without peak filling (to get number of compounds)
# an2 <- xsAnnotate(xset4)
# anF2 <- groupFWHM(an2, perfwhm = 0.6)
# anI2 <- findIsotopes(anF2, mzabs = 0.01)
# anIC2 <- groupCorr(anI2, cor_eic_th = 0.75, graphMethod = "lpc")
# anFA2 <- findAdducts(anIC2, polarity = "positive")
peaklist2 <- getPeaklist(anFA2)
peaklist2[, cols] <- apply(peaklist2[, cols], 2, FUN = function(x) (x / quantile(x, 0.75, na.rm = TRUE) + 1))
## write pcgroups of peaklist to peaklist2
peaklist2[, "pcgroup"] <- peaklist[, "pcgroup"]
# save(an, anF, anI, anIC, anFA, peaklist, an2, anF2, anI2, anIC2, anFA2, peaklist2, file = "./CAMERA_complete.RData")

MSMS <- read.csv("../../MSMS/idmsms_3_3_0.8_minCor_rmv50.csv")
MSMS <- MSMS[,c(2,3,4,8)]


## remove these entries from the MSMS which do not have the precursor ions in the fragmentation
## start 
uniqPrecMZRTPC <- unique(as.character(MSMS[, "precursor"]))
PrecMZRTPC <- as.character(MSMS[, "precursor"])
 
MSMS_mod <- cbind(MSMS, "check" = numeric(dim(MSMS)[1]))
for (i in 1:length(uniqPrecMZRTPC)) {
    mzPC <- as.numeric(cutUniquePreMZ(uniqPrecMZRTPC[i], splitPattern = "_", splitInd = 1))
    inds <- which(PrecMZRTPC == uniqPrecMZRTPC[i])
    mzGROUPS <- MSMS[inds, "mz"]
    if( any(abs(mzGROUPS - mzPC) < 0.02) ) MSMS_mod[inds, "check"] <- TRUE
}
## remove lines which have check == 0
MSMS <- MSMS[MSMS_mod[, "check"] == 1, ]
## end remove


## rename column inten to intensity
colnames(MSMS)[which(colnames(MSMS) == "inten")] <- "intensity"
##colnames(MSMS)[which(colnames(MSMS) == "precursor")] <- "pcgroup_precursorMZ"
# 
# ## truncate MSMS (after 40 min there is only isocratic 15% A and 85% B)
# ## truncate MSMS after 41 min = 60*41 s = 2460 s
# MSMS <- MSMS[which(as.numeric(MSMS[, "rt"]) < 2460), ]

precursorMZ <- unlist(lapply(strsplit(as.character(MSMS[,"precursor"]), "_"), "[", 1))
precursorMZ <- as.numeric(precursorMZ)
precursorMZ_unique <- unique(precursorMZ)
precursorRT <- unlist(lapply(strsplit(as.character(MSMS[,"precursor"]), "_"), "[", 2))
precursorRT <- as.numeric(precursorRT)
precursorRT_unique <- unique(precursorRT)
uniqueMZRTPC <- as.character(unique(MSMS[, "precursor"]))
precursorRT_minute <- precursorRT / 60 - 1 ## minimum is now 0.027 minutes (in the gradient phase)
gradientMSMS <- numeric(length = length(precursorRT_minute))
gradientMSMS[which(precursorRT_minute <= 0)] <- 0.9
gradientMSMS[which(precursorRT_minute > 0)] <- 0.90 - 0.01923077 * precursorRT_minute[which(precursorRT_minute > 0)]
gradientMSMS[which(precursorRT_minute >= 39)] <- 0.15 ## between 1min and 40min there is gradient phase

## gradient = 0.90 - 0.0727 * minute --> minute = (0.9 - gradient) / 0.0727
## +1  since 1 min isocratic 90% A

##gradient_mapped_minute <- (0.9 - gradient ) / 0.07272727272727 + 1
##gradient_mapped <- gradient_mapped_minute * 60

precursorMZ[1]
## profiling data for W+OS
##load("/home/thomas/Documents/University/Master/MScArbeit/Metabolic_profiling/WOS_0h_72h/metabolicProfiling.RData")
##load("/home/thomas/Documents/University/Master/MScArbeit/Metabolic_profiling/WOS_0h_72h/CAMERA_complete.RData")
peaklist <- getPeaklist(anFA)
gradientProfiling <- numeric(length = length(peaklist[,"rt"]))
precursorRT_minute_profiling <- peaklist[,"rt"] / 60 - 1 ## 1 min isocratic
gradientProfiling[which(precursorRT_minute_profiling <= 0)] <- 0.9
gradientProfiling[which(precursorRT_minute_profiling > 0)] <- 0.9 - 0.0727272727272727 * precursorRT_minute_profiling[which(precursorRT_minute_profiling > 0)]
## from minute 12 (11) isocratic 0.1
gradientProfiling[which(precursorRT_minute_profiling > 11)] <- 0.1


## profiling data for MJ
#load("/home/thomas/Documents/University/Master/MScArbeit/Metabolic_profiling/MeJa_0h_72h/metabolicProfilingMJ.RData")
#load("/home/thomas/Documents/University/Master/MScArbeit/Metabolic_profiling/MeJa_0h_72h/CAMERA_completeMJ.RData")
# peaklistMJ <- getPeaklist(anFAMJ)
# gradientProfilingMJ <- numeric(length = length(peaklistMJ[,"rt"]))
# precursorRT_minute_profilingMJ <- peaklistMJ[,"rt"] / 60 - 1 ## 1 min isocratic
# gradientProfilingMJ[which(precursorRT_minute_profilingMJ <= 0)] <- 0.9
# gradientProfilingMJ[which(precursorRT_minute_profilingMJ > 0)] <- 0.9 - 0.0727272727272727 * precursorRT_minute_profilingMJ[which(precursorRT_minute_profilingMJ > 0)]
# ## from minute 12 (11) isocratic 0.1
# gradientProfilingMJ[which(precursorRT_minute_profilingMJ > 11)] <- 0.1

## prepare MSMS
## add column gradientWOS, gradientMJ and gradientMSMS
MSMS_mod <- cbind(MSMS, gradientWOS = numeric(dim(MSMS)[1]))
MSMS_mod <- cbind(MSMS_mod, gradientMJ = numeric(dim(MSMS)[1]))
MSMS_mod <- cbind(MSMS_mod, gradientMSMS = numeric(dim(MSMS)[1]))
MSMS_mod[, "gradientMSMS"] <- gradientMSMS
## add column mzWOS,rtWOS that is the mz and retention time of mapped features of peaklistWOS
MSMS_mod <- cbind(MSMS_mod, mzWOS = numeric(dim(MSMS)[1]))
MSMS_mod <- cbind(MSMS_mod, rtWOS = numeric(dim(MSMS)[1]))
## add column mzMJ,rtMJ that is the mz and retention time of mapped features of peaklistMJ
MSMS_mod <- cbind(MSMS_mod, mzMJ = numeric(dim(MSMS)[1]))
MSMS_mod <- cbind(MSMS_mod, rtMJ = numeric(dim(MSMS)[1]))
## add columns that have number of biological replicates that synthesise compound: 
## e.g. NattWOS72, NobtWOS72, ... NattMJ0, ... NattMJ72, ...
mm <- matrix(0, nrow = dim(MSMS)[1], ncol = 18)
colnames(mm) <-  c("att0", "obt0", "clev0", "quad0", "x10270", "x571260",
  "att72WOS", "obt72WOS", "clev72WOS", "quad72WOS", "x102772WOS", "x5712672WOS",
  "att72MJ", "obt72MJ", "clev72MJ", "quad72MJ", "x102772MJ", "x5712672MJ")
MSMS_mod <- cbind(MSMS_mod, mm)
## add column: mapped1WOS and mapped1MJ, i.e. was this feature mapped in the first round = mapped by deviance?
## add column: mapped2WOS and mapped2MJ, i.e. was this feature mapped in the second round = mapped by interval?
MSMS_mod <- cbind(MSMS_mod, mapped1WOS = numeric(dim(MSMS)[1]))
MSMS_mod <- cbind(MSMS_mod, mapped1MJ = numeric(dim(MSMS)[1]))
MSMS_mod <- cbind(MSMS_mod, mapped2WOS = numeric(dim(MSMS)[1]))
MSMS_mod <- cbind(MSMS_mod, mapped2MJ = numeric(dim(MSMS)[1]))
MSMS_mod[, "mapped1WOS"] <- factor(x = MSMS_mod[, "mapped1WOS"], levels = c(0,1))
MSMS_mod[, "mapped1MJ"] <- factor(x = MSMS_mod[, "mapped1MJ"], levels = c(0,1))
MSMS_mod[, "mapped2WOS"] <- factor(x = MSMS_mod[, "mapped2WOS"], levels = c(0,1))
MSMS_mod[, "mapped2MJ"] <- factor(x = MSMS_mod[, "mapped2MJ"], levels = c(0,1))


## W+OS
devWOS <- numeric(length(precursorMZ))
devMJ <- numeric(length(precursorMZ))

gradient <- 0.1 ## tolerated deviance in gradient, define greater deviance 
## since we use a "combined" peaklist of WOS and MeJA runs


## round 1: define deviance for gradient and deviance for mz and map based on
## these criteria
for (i in 1:length(precursorMZ)) {
    
    ## WOS+MJ
    ## shrink space of possible mapped features by gradient deviance
    ind <- which(abs(gradientMSMS[i] - gradientProfiling) <= gradient)
    ## get feature with minimum deviance to mz
    ind_minMZ <- ind[which.min(abs(precursorMZ[i] - peaklist[ind, "mz"] ))]
    minInPeaklist <- peaklist[ind_minMZ,]
    
    if (abs(minInPeaklist["mz"] - precursorMZ[i]) <= 0.01) { ## tolerated m/z deviance is 0.01
        ## set mapped1 to 1
        MSMS_mod[i, "mapped1WOS"] <- 1
        MSMS_mod[i, "mapped1MJ"] <- 1
        ## write gradient of profiling to column gradientWOS
        MSMS_mod[i, "gradientWOS"] <- MSMS_mod[i, "gradientMJ"] <- gradientProfiling[ind_minMZ]
        
        ## write mz of profiling to column mzWOS and mzMJ
        MSMS_mod[i, "mzWOS"] <- MSMS_mod[i, "mzMJ"] <- peaklist[ind_minMZ, "mz"]
        ## write retention time of profiling to column rtWOS and rtMJ
        MSMS_mod[i, "rtWOS"] <- MSMS_mod[i, "rtMJ"] <- peaklist[ind_minMZ, "rt"]
        
        ## write numbers of biological replicates where feature is present: C
        MSMS_mod[i, "att0"] <- peaklist[ind_minMZ, "att0"]
        MSMS_mod[i, "obt0"] <- peaklist[ind_minMZ, "obt0"]
        MSMS_mod[i, "clev0"] <- peaklist[ind_minMZ, "clev0"]
        MSMS_mod[i, "quad0"] <- peaklist[ind_minMZ, "quad0"]
        MSMS_mod[i, "x10270"] <- peaklist[ind_minMZ, "x10270"]
        MSMS_mod[i, "x571260"] <- peaklist[ind_minMZ, "x570"]
        
        ## write numbers of biological replicates where feature is present: WOS
        MSMS_mod[i, "att72WOS"] <- peaklist[ind_minMZ, "att72WOS"]
        MSMS_mod[i, "obt72WOS"] <- peaklist[ind_minMZ, "obt72WOS"]
        MSMS_mod[i, "clev72WOS"] <- peaklist[ind_minMZ, "clev72WOS"]
        MSMS_mod[i, "quad72WOS"] <- peaklist[ind_minMZ, "quad72WOS"]
        MSMS_mod[i, "x102772WOS"] <- peaklist[ind_minMZ, "x1072WOS"]
        MSMS_mod[i, "x5712672WOS"] <- peaklist[ind_minMZ, "x5772WOS"]
        
        ## write numbers of biological replicates where feature is present: MJ
        MSMS_mod[i, "att72MJ"] <- peaklist[ind_minMZ, "att72MJ"]
        MSMS_mod[i, "obt72MJ"] <- peaklist[ind_minMZ, "obt72MJ"]
        MSMS_mod[i, "clev72MJ"] <- peaklist[ind_minMZ, "clev72MJ"]
        MSMS_mod[i, "quad72MJ"] <- peaklist[ind_minMZ, "quad72MJ"]
        MSMS_mod[i, "x102772MJ"] <- peaklist[ind_minMZ, "x1072MJ"]
        MSMS_mod[i, "x5712672MJ"] <- peaklist[ind_minMZ, "x5772MJ"]
    }
}

# for (i in 1:length(precursorMZ)) {
#     
#     ## WOS 
#     ## shrink space of possible mapped features by gradient deviance
#     indWOS <- which(abs(gradientMSMS[i] - gradientProfiling) <= gradient)
#     ## get feature with minimum deviance to mz
#     indWOS_minMZ <- indWOS[which.min(abs(precursorMZ[i] - peaklist[indWOS, "mz"] ))]
#     minInPeaklistWOS <- peaklist[indWOS_minMZ,]
#     
#     if (abs(minInPeaklistWOS["mz"] - precursorMZ[i]) <= 0.01) { ## tolerated m/z deviance is 0.01
#         ## set mapped1 to 1
#         MSMS_mod[i, "mapped1WOS"] <- 1
#         ## write gradient of profiling to column gradientWOS
#         MSMS_mod[i, "gradientWOS"] <- gradientProfiling[indWOS_minMZ]
# 
#         ## write mz of profiling to column mzWOS
#         MSMS_mod[i, "mzWOS"] <- peaklist[indWOS_minMZ, "mz"]
#         ## write retention time of profiling to column rtWOS
#         MSMS_mod[i, "rtWOS"] <- peaklist[indWOS_minMZ, "rt"]
#         
#         ## write numbers of biological replicates where feature is present
#         MSMS_mod[i, "att72WOS"] <- peaklist[indWOS_minMZ, "att72"]
#         MSMS_mod[i, "obt72WOS"] <- peaklist[indWOS_minMZ, "obt72"]
#         MSMS_mod[i, "clev72WOS"] <- peaklist[indWOS_minMZ, "clev72"]
#         MSMS_mod[i, "quad72WOS"] <- peaklist[indWOS_minMZ, "quad72"]
#         MSMS_mod[i, "x102772WOS"] <- peaklist[indWOS_minMZ, "x1072"]
#         MSMS_mod[i, "x5712672WOS"] <- peaklist[indWOS_minMZ, "x5772"]
#     }
#         
#     ## MJ
#     ## shrink space of possible mapped features by gradient deviance
#     indMJ <- which(abs(gradientMSMS[i] - gradientProfilingMJ) <= gradient)
#     ## get feature with minimum deviance to mz
#     indMJ_minMZ <- indMJ[which.min(abs(precursorMZ[i] - peaklistMJ[indMJ, "mz"] ))]
#     minInPeaklistMJ <- peaklistMJ[indMJ_minMZ,]
#     
#     if (abs(peaklistMJ[indMJ_minMZ, "mz"] - precursorMZ[i]) <= 0.01) { ## tolerated m/z deviance is 0.01
#         ## set mapped1 to 1
#         MSMS_mod[i, "mapped1MJ"] <- 1
#         ## write gradient of profiling to column gradientMJ
#         MSMS_mod[i, "gradientMJ"] <- gradientProfilingMJ[indMJ_minMZ]
# 
#         ## write mz of profiling to column mzMJ
#         MSMS_mod[i, "mzMJ"] <- peaklistMJ[indMJ_minMZ, "mz"]
#         ## write retention time of profiling to column rtWOS
#         MSMS_mod[i, "rtMJ"] <- peaklistMJ[indMJ_minMZ, "rt"]
#         
#         ## write numbers of biological replicates where feature is present
#         MSMS_mod[i, "att0MJ"] <- peaklistMJ[indMJ_minMZ, "att0"]
#         MSMS_mod[i, "obt0MJ"] <- peaklistMJ[indMJ_minMZ, "obt0"]
#         MSMS_mod[i, "clev0MJ"] <- peaklistMJ[indMJ_minMZ, "clev0h"]
#         MSMS_mod[i, "quad0MJ"] <- peaklistMJ[indMJ_minMZ, "quad0h"]
#         MSMS_mod[i, "x10270MJ"] <- peaklistMJ[indMJ_minMZ, "x10270"]
#         MSMS_mod[i, "x571260MJ"] <- peaklistMJ[indMJ_minMZ, "x570"]
#         
#         MSMS_mod[i, "att72MJ"] <- peaklistMJ[indMJ_minMZ, "att72"]
#         MSMS_mod[i, "obt72MJ"] <- peaklistMJ[indMJ_minMZ, "obt72"]
#         MSMS_mod[i, "clev72MJ"] <- peaklistMJ[indMJ_minMZ, "clev72"]
#         MSMS_mod[i, "quad72MJ"] <- peaklistMJ[indMJ_minMZ, "quad72"]
#         MSMS_mod[i, "x102772MJ"] <- peaklistMJ[indMJ_minMZ, "x1072"]
#         MSMS_mod[i, "x5712672MJ"] <- peaklistMJ[indMJ_minMZ, "x5772"]
#     }
# }
## end of round 1

## round 2: use results from round 1 and define a retention time window 
## between already mapped features, check then in this window if other (not yet
## mapped) features can be mapped

for (i in 1:length(precursorMZ)) {
    
    ## for WOS
    ##mappedWOS <- which(MSMS_mod[, "mapped1WOS"] == 1 ) 
    indmapped <- which(MSMS_mod[, "mapped1WOS"] == 1)
    mappedGradients <- MSMS_mod[indmapped, "gradientMSMS"]
    mappedGradients_uni <- unique(mappedGradients)
    
    if (MSMS_mod[i,"mapped1WOS"] == 0) {
        
        ## get feature that has gradient of rank +-3 to calculated one
        devGradient <- MSMS_mod[i, "gradientMSMS"] - mappedGradients_uni
        ## upper and lower feature which is +-3
        devGradient_ind_u <- which( devGradient == sort(devGradient[devGradient > 0])[3] )
        devGradient_ind_l <- which( devGradient == sort(devGradient[devGradient < 0])[3] )
        
        ## retrieve respective mapped MSMS feature with lower and higher retention time
        ## that will use as a lower and upper bound for search space
        lower <- MSMS_mod[intersect(which(MSMS_mod[, "gradientMSMS"] ==  mappedGradients_uni[devGradient_ind_l]),indmapped), ]
        upper <- MSMS_mod[intersect(which(MSMS_mod[, "gradientMSMS"] ==  mappedGradients_uni[devGradient_ind_u]), indmapped), ]
        
        ## implement a rule for boundary values
        if (dim(upper)[1] == 0) {
            upper <- lower
            upper[,"rtWOS"] <- lower[,"rtWOS"] + 10
            upper[,"rtMJ"] <- lower[,"rtMJ"] + 10
        }
        if (dim(lower)[1] == 0) {
            lower <- upper
            lower[,"rtWOS"] <- upper[,"rtWOS"] - 10
            lower[,"rtMJ"] <- upper[,"rtMJ"] - 10
        }
        
        
        
        upperRT <- max(unique(upper[, "rtWOS"]))
        lowerRT <- min(unique(lower[, "rtWOS"]))
        
        ## implement a rule that there is a certain range of 20s to look into
        ## when 
        if (upperRT - lowerRT < 20) {
            upperRT <- upperRT + 10
            lowerRT <- lowerRT - 10
        }
        
        ind_tr <- intersect(which(peaklist[, "rt"] <= upperRT), which(peaklist[, "rt"] >= lowerRT))
        
        peaklist_tr <- peaklist[ind_tr, ]
        
        ind_mapped <- which.min(abs(peaklist_tr[, "mz"] - precursorMZ[i]))
        mapped <- peaklist_tr[ind_mapped, ]
        
        if (abs(mapped[, "mz"] - precursorMZ[i]) <= 0.008) {
            
            ind_minMZ <- ind_tr[ind_mapped]
            
            ## set mapped2 to 1
            MSMS_mod[i, "mapped2WOS"] <- MSMS_mod[i, "mapped2MJ"] <- 1
            ## write gradient of profiling to column gradientWOS
            MSMS_mod[i, "gradientWOS"] <- MSMS_mod[i, "gradientMJ"] <- gradientProfiling[ind_minMZ]
            
            ## write mz of profiling to column mzWOS and mzMJ
            MSMS_mod[i, "mzWOS"] <- MSMS_mod[i, "mzMJ"] <- peaklist[ind_minMZ, "mz"]
            ## write retention time of profiling to column rtWOS and rtMJ
            MSMS_mod[i, "rtWOS"] <- MSMS_mod[i, "rtMJ"] <- peaklist[ind_minMZ, "rt"]
            
            ## write numbers of biological replicates where feature is present: C
            MSMS_mod[i, "att0"] <- peaklist[ind_minMZ, "att0"]
            MSMS_mod[i, "obt0"] <- peaklist[ind_minMZ, "obt0"]
            MSMS_mod[i, "clev0"] <- peaklist[ind_minMZ, "clev0"]
            MSMS_mod[i, "quad0"] <- peaklist[ind_minMZ, "quad0"]
            MSMS_mod[i, "x10270"] <- peaklist[ind_minMZ, "x10270"]
            MSMS_mod[i, "x571260"] <- peaklist[ind_minMZ, "x570"]
            
            ## write numbers of biological replicates where feature is present: WOS
            MSMS_mod[i, "att72WOS"] <- peaklist[ind_minMZ, "att72WOS"]
            MSMS_mod[i, "obt72WOS"] <- peaklist[ind_minMZ, "obt72WOS"]
            MSMS_mod[i, "clev72WOS"] <- peaklist[ind_minMZ, "clev72WOS"]
            MSMS_mod[i, "quad72WOS"] <- peaklist[ind_minMZ, "quad72WOS"]
            MSMS_mod[i, "x102772WOS"] <- peaklist[ind_minMZ, "x1072WOS"]
            MSMS_mod[i, "x5712672WOS"] <- peaklist[ind_minMZ, "x5772WOS"]
            
            ## write numbers of biological replicates where feature is present: MJ
            MSMS_mod[i, "att72MJ"] <- peaklist[ind_minMZ, "att72MJ"]
            MSMS_mod[i, "obt72MJ"] <- peaklist[ind_minMZ, "obt72MJ"]
            MSMS_mod[i, "clev72MJ"] <- peaklist[ind_minMZ, "clev72MJ"]
            MSMS_mod[i, "quad72MJ"] <- peaklist[ind_minMZ, "quad72MJ"]
            MSMS_mod[i, "x102772MJ"] <- peaklist[ind_minMZ, "x1072MJ"]
            MSMS_mod[i, "x5712672MJ"] <- peaklist[ind_minMZ, "x5772MJ"]
            
        }
    }
}


# for (i in 1:length(precursorMZ)) {
# 
#     ## for WOS
#     ##mappedWOS <- which(MSMS_mod[, "mapped1WOS"] == 1 ) 
#     indmapped <- which(MSMS_mod[, "mapped1WOS"] == 1)
#     mappedGradients <- MSMS_mod[indmapped, "gradientMSMS"]
#     mappedGradients_uni <- unique(mappedGradients)
#     
#     if (MSMS_mod[i,"mapped1WOS"] == 0) {
#         
#         ## get feature that has gradient of rank +-3 to calculated one
#         devGradient <- MSMS_mod[i, "gradientMSMS"] - mappedGradients_uni
#         ## upper and lower feature which is +-3
#         devGradient_ind_u <- which(devGradient == sort(devGradient[devGradient >= 0])[3])
#         devGradient_ind_l <- which(devGradient == sort(devGradient[devGradient <= 0])[3])
#         
#         ## retrieve respective mapped MSMS feature with lower and higher retention time
#         ## that will use as a lower and upper bound for search space
#         upper <- MSMS_mod[intersect(which(MSMS_mod[, "gradientMSMS"] ==  mappedGradients_uni[devGradient_ind_u]),indmapped), ]
#         lower <- MSMS_mod[intersect(which(MSMS_mod[, "gradientMSMS"] ==  mappedGradients_uni[devGradient_ind_l]), indmapped), ]
#         
#         ## implement a rule for boundary values
#         if (dim(upper)[1] == 0) {
#             upper <- lower
#             upper[,"rtWOS"] <- lower[,"rtWOS"] + 10
#         }
#         if (dim(lower)[1] == 0) {
#             lower <- upper
#             lower[,"rtWOS"] <- upper[,"rtWOS"] - 10
#         }
#         
#         upperRT <- max(unique(upper[, "rtWOS"]))
#         lowerRT <- min(unique(lower[, "rtWOS"]))
#         ind_tr <- intersect(which(peaklist[, "rt"] <= upperRT), which(peaklist[, "rt"] >= lowerRT))
#         
#         peaklist_tr <- peaklist[ind_tr, ]
#         
#         ind_mapped <- which.min(abs(peaklist_tr[, "mz"] - precursorMZ[i]))
#         mapped <- peaklist_tr[ind_mapped, ]
#         
#         if (abs(mapped[, "mz"] - precursorMZ[i]) <= 0.01) {
#             
#             indWOS_minMZ <- ind_tr[ind_mapped]
#             
#             ## set mapped2 to 1
#             MSMS_mod[i, "mapped2WOS"] <- 1
#             ## write gradient of profiling to column gradientWOS
#             MSMS_mod[i, "gradientWOS"] <- gradientProfiling[indWOS_minMZ]
#             
#             ## write mz of profiling to column mzWOS
#             MSMS_mod[i, "mzWOS"] <- peaklist[indWOS_minMZ, "mz"]
#             ## write retention time of profiling to column rtWOS
#             MSMS_mod[i, "rtWOS"] <- peaklist[indWOS_minMZ, "rt"]
#             
#             ## write numbers of biological replicates where feature is present
#             MSMS_mod[i, "att72WOS"] <- peaklist[indWOS_minMZ, "att72"]
#             MSMS_mod[i, "obt72WOS"] <- peaklist[indWOS_minMZ, "obt72"]
#             MSMS_mod[i, "clev72WOS"] <- peaklist[indWOS_minMZ, "clev72"]
#             MSMS_mod[i, "quad72WOS"] <- peaklist[indWOS_minMZ, "quad72"]
#             MSMS_mod[i, "x102772WOS"] <- peaklist[indWOS_minMZ, "x1072"]
#             MSMS_mod[i, "x5712672WOS"] <- peaklist[indWOS_minMZ, "x5772"]
#             
#         }
#     }
# 
#     ## for MeJA
#     indmapped <- which(MSMS_mod[, "mapped1MJ"] == 1)
#     mappedGradients <- MSMS_mod[indmapped, "gradientMSMS"]
#     mappedGradients_uni <- unique(mappedGradients)
#     
#     if (MSMS_mod[i,"mapped1MJ"] == 0) {
#         
#         ## get feature that has gradient of rank +-3 to calculated one
#         devGradient <- MSMS_mod[i, "gradientMSMS"] - mappedGradients_uni
#         ## upper and lower feature which is +-3
#         devGradient_ind_u <- which(devGradient == sort(devGradient[devGradient >= 0])[3])
#         devGradient_ind_l <- which(devGradient == sort(devGradient[devGradient <= 0])[3])
#         
#         ## retrieve respective mapped MSMS feature with lower and higher retention time
#         ## that will use as a lower and upper bound for search space
#         upper <- MSMS_mod[intersect(which(MSMS_mod[, "gradientMSMS"] ==  mappedGradients_uni[devGradient_ind_u]),indmapped), ]
#         lower <- MSMS_mod[intersect(which(MSMS_mod[, "gradientMSMS"] ==  mappedGradients_uni[devGradient_ind_l]), indmapped), ]
#         
#         ## implement a rule for boundary values
#         if (dim(upper)[1] == 0) {
#             upper <- lower
#             upper[,"rtMJ"] <- lower[,"rtMJ"] + 10
#         }
#         if (dim(lower)[1] == 0) {
#             lower <- upper
#             lower[,"rtMJ"] <- upper[,"rtMJ"] - 10
#         }
#         
#         upperRT <- max(unique(upper[, "rtMJ"]))
#         lowerRT <- min(unique(lower[, "rtMJ"]))
#         ind_tr <- intersect(which(peaklistMJ[, "rt"] <= upperRT), which(peaklistMJ[, "rt"] >= lowerRT))
#         
#         peaklist_tr <- peaklistMJ[ind_tr, ]
#         
#         ind_mapped <- which.min(abs(peaklist_tr[, "mz"] - precursorMZ[i]))
#         mapped <- peaklist_tr[ind_mapped, ]
#         
#     
#         if (abs(mapped[, "mz"] - precursorMZ[i]) <= 0.01) {
# 
#             indMJ_minMZ <- ind_tr[ind_mapped]
# 
#             ## set mapped2 to 1
#             MSMS_mod[i, "mapped2MJ"] <- 1
#             ## write gradient of profiling to column gradientMJ
#             MSMS_mod[i, "gradientMJ"] <- gradientProfilingMJ[indMJ_minMZ]
# 
#             ## write mz of profiling to column mzMJ
#             MSMS_mod[i, "mzMJ"] <- peaklistMJ[indMJ_minMZ, "mz"]
#             ## write retention time of profiling to column rtMJ
#             MSMS_mod[i, "rtMJ"] <- peaklistMJ[indMJ_minMZ, "rt"]
# 
#             ## write numbers of biological replicates where feature is present
#             MSMS_mod[i, "att0MJ"] <- peaklistMJ[indMJ_minMZ, "att0"]
#             MSMS_mod[i, "obt0MJ"] <- peaklistMJ[indMJ_minMZ, "obt0"]
#             MSMS_mod[i, "clev0MJ"] <- peaklistMJ[indMJ_minMZ, "clev0h"]
#             MSMS_mod[i, "quad0MJ"] <- peaklistMJ[indMJ_minMZ, "quad0h"]
#             MSMS_mod[i, "x10270MJ"] <- peaklistMJ[indMJ_minMZ, "x10270"]
#             MSMS_mod[i, "x571260MJ"] <- peaklistMJ[indMJ_minMZ, "x570"]
# 
#             MSMS_mod[i, "att72MJ"] <- peaklistMJ[indMJ_minMZ, "att72"]
#             MSMS_mod[i, "obt72MJ"] <- peaklistMJ[indMJ_minMZ, "obt72"]
#             MSMS_mod[i, "clev72MJ"] <- peaklistMJ[indMJ_minMZ, "clev72"]
#             MSMS_mod[i, "quad72MJ"] <- peaklistMJ[indMJ_minMZ, "quad72"]
#             MSMS_mod[i, "x102772MJ"] <- peaklistMJ[indMJ_minMZ, "x1072"]
#             MSMS_mod[i, "x5712672MJ"] <- peaklistMJ[indMJ_minMZ, "x5772"]
# 
#         }
#     }
# }



## truncate MSMS_mod: remove entries which have sum of 0 in the mentioned columns
## i.e. remove the ones that could not be mapped
MSMS_mod <- MSMS_mod[apply(data.matrix(MSMS_mod[, c("mapped1WOS", "mapped1MJ", "mapped2WOS", "mapped2MJ")]) - 1, 1, sum) > 0,]

## change entries of biological replicates to binary values:
## set the entries with less than 6 replicates to 0, with more than 6 replicates to 1
entriesC <- MSMS_mod[,which(colnames(MSMS_mod) == "att0"):which(colnames(MSMS_mod) == "x571260")]
entriesC[entriesC < 6] <- 0
entriesC[entriesC >= 6] <- 1
## set the entries with less than 3 replicates to 0, with more than 3 replicates to 1
entries72 <- MSMS_mod[,which(colnames(MSMS_mod) == "att72WOS"):which(colnames(MSMS_mod) == "x5712672MJ")]
entries72[entries72 < 3] <- 0
entries72[entries72 >= 3] <- 1

## write entries to MSMS_mod: replace by binary matrix entries
MSMS_mod[,which(colnames(MSMS_mod) == "att0"):which(colnames(MSMS_mod) == "x571260")] <- entriesC
MSMS_mod[,which(colnames(MSMS_mod) == "att72WOS"):which(colnames(MSMS_mod) == "x5712672MJ")] <- entries72

## create MSMS_mod for WOS: create data frame with mappings for WOS only
#MSMS_WOS <- MSMS_mod[apply(data.matrix(MSMS_mod[, c("mapped1WOS", "mapped2WOS")]) -1, 1, sum) > 0,]
## remove column with MJ
#MSMS_WOS <- MSMS_WOS[, !(names(MSMS_WOS) %in% c("att0MJ", "obt0MJ","clev0MJ","quad0MJ","x10270MJ","x571260MJ"))]
#MSMS_WOS <- MSMS_WOS[, !(names(MSMS_WOS) %in% c("att72MJ", "obt72MJ", "clev72MJ", "quad72MJ", "x102772MJ", "x5712672MJ"))]
#MSMS_WOS <- MSMS_WOS[, !(names(MSMS_WOS) %in% c("mapped1MJ", "mapped2MJ"))]

uniquePrecursor <- unique(MSMS_mod[, "precursor"])
## first row entries with unique precursor
indsMSMS <- match(uniquePrecursor, MSMS_mod[, "precursor"])

## how many metabolites are found in each species?
apply(data.matrix(MSMS_mod[indsMSMS, c("att72WOS", "obt72WOS", "clev72WOS", "quad72WOS", "x102772WOS", "x5712672WOS")]), 2, sum)


## MSMS_mod for MJ0
#MSMS_MJ <- MSMS_mod[apply(data.matrix(MSMS_mod[, c("mapped1MJ", "mapped2MJ")]) -1, 1, sum) > 0,]
##MSMS_MJ <- MSMS_MJ[, !(names(MSMS_MJ0) %in% c("att0MJ", "obt0MJ","clev0MJ","quad0MJ","x10270MJ","x571260MJ"))]
#MSMS_MJ <- MSMS_MJ[, !(names(MSMS_MJ) %in% c("att72WOS", "obt72WOS", "clev72WOS", "quad72WOS", "x102772WOS", "x5712672WOS"))]
#MSMS_MJ <- MSMS_MJ[, !(names(MSMS_MJ) %in% c("mapped1WOS", "mapped2WOS"))]

#uniqueMJ <- unique(MSMS_MJ[, "precursor"])
#indsMJ <- match(uniqueMJ, MSMS_MJ[, "precursor"])

apply(data.matrix(MSMS_mod[indsMSMS, c("att0", "obt0","clev0","quad0","x10270","x571260")]), 2, sum)
apply(data.matrix(MSMS_mod[indsMSMS, c("att72MJ", "obt72MJ", "clev72MJ", "quad72MJ", "x102772MJ", "x5712672MJ")]), 2, sum)

## end



## append column which state if the compound could be mapped by precursor ion
MSMS_mod <- cbind(MSMS, mapped_precursor = rep(1, dim(MSMS)[1]))
mm <- matrix(0, ncol = 18, nrow = dim(MSMS)[1])
colnames(mm) <- c("Natt0", "Natt72WOS", "Natt72MJ", "Nobt0", "Nobt72WOS", "Nobt72MJ", 
    "Nclev0", "Nclev72WOS", "Nclev72MJ", "Nquad0", "Nquad72WOS", "Nquad72MJ",
    "Nxobt10270", "Nxobt102772WOS", "Nxobt102772MJ", 
    "Nxobt571260", "Nxobt5712672WOS", "Nxobt5712672MJ")
MSMS_mod <- cbind(MSMS_mod, mm)
MSMS_mod[MSMS_mod[,"precursor"] %in% as.character(unique(MSMS[, "precursor"]))[indNotMapped], "mapped_precursor"] <- 0
## append column which state if the fragment with highest intensity could be 
## mapped (use a decreased range, e.g. 0.05 for it)
i <- 10
inds <- which(MSMS[, "precursor"] == uniqueMZRTPC[indNotMapped[i]])

indWOS <- which(abs(gradientMSMS[indNotMapped[i]] - gradientProfiling) <= 0.03)
which(abs(peaklist[indWOS,"mz"] - MSMS[inds[which(MSMS[inds, "intensity"] == 100)], "mz"]) < 0.01)
## append columns that give information where the compound is found 
## (e.g. Natt0: 0, Natt72WOS: 4, Natt72MJ: 5)
## extract information if the metabolite is produced or not (if 3/5 biological replicates have it)

indNotMapped <- intersect(which(abs(devWOS) > deviances[k]), which(abs(devMJ) > deviances[k]))
i <- 4
peaklistMSMS2[which.min(1000 * abs( precursorMZ[indNotMapped][i] - peaklistMSMS2[,"mz"] ) + abs(precursorRT[indNotMapped][i] - peaklistMSMS2[, "rt"])),]
precursorMZ[indNotMapped][i]; precursorRT[indNotMapped][i]
minInPeaklistWOS <- peaklist[indWOS[which.min(abs(precursorMZ[i] - peaklist[indWOS, "mz"] ))],]
##[1]   1  12  38  55  57  60  64  66  69  76  77  79  90  91 113 115 116 117 118 121 126 131 141 145 147 148 149 150
##[29] 152 158 160 165 167 168 177 182 196 197 219 224 228 234 235 236 238 239 246 248 249 251 258 260 265 266 268 269
##[57] 272 275 277 279 286 288 291 292 293 295 301 302 303 304 305 306 307 308 309 310 311 315 316 320 322 327 329 330
##[85] 331 333 336 337 339 342 343 344 346 347 349 351 357 369 376 386 387 393 395 397 406 413 420 421 422 427 430 432
##[113] 437 442 443 445 446 451 452 454 455 456 458 462 464 465 466 467 469 472 473 474 475 477 478 479 480 481 482 484
##[141] 485 497 498 499 500 501 512 534 542 563 570 571 573 574 577
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[1]),] # single
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[12]),] 
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[38]),] 
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[55]),] ## no precursor
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[57]),] 
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[60]),]
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[64]),] ## no precursor
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[90]),] ## no precursor
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[91]),] 
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[92]),] ## no precursor
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[96]),] ## no precursor
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[97]),] 
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[100]),] 
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[104]),] 
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[105]),] ## no precursor
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[110]),] 
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[113]),] 
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[114]),] ## no precursor
MSMS[which(MSMS[, "precursor"] == as.character(unique(MSMS[, "precursor"]))[123]),] ## 


## check if precursor can be found in peaktable of MS1 (45min), if not check if 100% fragment ion can be found in MS1 
levelplot(mat)


minInPeaklist[,"pcgroup"]
hist(devWOS)
hist(devMJ)
precursorMZ[i]


as.character(unique(MSMS[,"pcgroup_precursorMZ"]))[intersect(which(abs(devWOS) >= 0.01), which(abs(devMJ) >= 0.01))]
load("/home/thomas/Documents/University/Master/MScArbeit/Metabolic_profiling/WOS_0h_72h//bootstrapping_WOS.RData")
inducedAtt
##abs(precursorMZ[i] - peaklist[, "mz"]) + abs(gradientMSMS[i] - gradientProfiling)
##ind_candidates <- intersect(which(peaklist[, "rt"] <= (gradient_mapped[i] + 45)), which(peaklist[,"rt"] >= (gradient_mapped[i] - 45)))
peaklist[ind_candidates, "mz"]
