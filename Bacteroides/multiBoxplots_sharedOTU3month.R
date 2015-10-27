# This code was written by Arnon Lieber, Oct 2015 as part of the subbmission to the STM.
# The code search for the number of observed otus in all samples, whether are child or mothers
# and sum them in a table. The variables names:
# StudyId: unique identifier of the mother/child dyad
# delivery: the way of birth - vaginal Vs. C-section
# uC: Number of unique OTUs in child
# uM: Number of unique OTUs in mothers
# M_Cdyad: Number of shared OTUs between mother and her child
# nonSelfMom: Number of shared OTUs between a child and other mothers
# 

# Set your home dir here:
setwd("")
library(phyloseq)

# If needed, subsets based on a taxonomy level (in this case, Genus)
################################################################################
# tb_bifido = subset_taxa(main, Genus=="Bifidobacterium")
# tb_bacteroides = subset_taxa(main, Genus=="Bacteroides")

### One must select here the taxa prefrence based on the list in "taxaName" ###
tb <- tb_bacteroides

# Lets validate that we chose the right taxa!
tax_table(tb)[1,1:7]

# List of all the StudyIDs that have been withdrown during the study 
# (table from Nora)
wdList = c(3, 13, 15, 19, 26, 29, 39, 48, 54, 28)
# Removes these samples from the dataset
tb <- prune_samples(!sample_data(tb)$StudyID %in% wdList, tb)
# Remove these samples of children who were exposed to antibiotic 
# before 6 months of life
abExpList = c(1, 17, 20, 24, 32, 41, 55)
tb <- prune_samples(!sample_data(tb)$StudyID %in% abExpList, tb)
# Extra NA removal
becauseNA <- 6
tb <- prune_samples(!sample_data(tb)$StudyID %in% becauseNA, tb)

sample_data(tb)$DOL <- as.numeric(as.character(sample_data(tb)$DOL))
sample_data(tb)$MOL <- as.numeric(as.character(sample_data(tb)$MOL))
sample_data(tb)$Month <- as.numeric(as.character(sample_data(tb)$Month))

tb_C <- subset_samples(tb, Mom_Child=="C")
tb_M <- subset_samples(tb, Mom_Child=="M")

# Pre-processing of Children datasets
#################################################################
# Remove OTUs with low horizontal number of counts (2)
tb_C <- prune_taxa(taxa_sums(tb_C)>=1, tb_C)
# Sets Month to the closest MOL
sample_data(tb_C)$Month <- round(sample_data(tb_C)$MOL, digits = 0)

# Pre-processing of Mothers datasets
#################################################################
# # Check the numbers of vaginal samples for each mother ID, before prune 
# tb_M_Vswab = subset_samples(tb_M, SampleType=="Vaginal_Swab")
# testSamples <- data.frame(sample_data(tb_M_Vswab))
# xtabs(~StudyID, data = testSamples)
# 
# Remove OTUs with low horizontal number of counts (2)
tb_M <- prune_taxa(taxa_sums(tb_M)>=1, tb_M)

# Sets Month to the closest MOL
sample_data(tb_M)$Month <- round(sample_data(tb_M)$MOL, digits = 0)

# # Remove samples with low (15) number of counts per sample
# tb_M <- prune_samples(sample_sums(tb_M)>=15, tb_M)

##### load("Mom_Child_tb_idetes")

tb_M_Rswab = subset_samples(tb_M, SampleType=="Rectal_Swab")
# Divide the rectal samples to 2:
    tb_M_R2 = subset_samples(tb_M_Rswab, (DOL>=-37 & DOL<0))
    tb_M_R1 = subset_samples(tb_M_Rswab, DOL< (-37))
tb_M_stool = subset_samples(tb_M, SampleType=="Stool_Stabilizer")
tb_M_fecal = 
    subset_samples(tb_M, (SampleType=="Rectal_Swab") |
                       (SampleType=="Stool_Stabilizer"))
tb_M_Vswab = subset_samples(tb_M, SampleType=="Vaginal_Swab")
    tb_M_Vpre = subset_samples(tb_M_Vswab, Month<=0)
    tb_M_Vpost = subset_samples(tb_M_Vswab, Month>0)

# subset the dataset here for even more sub-groups...
    
    ##### Internal functions section #############################################################
    
    # (1)
    uniqueOTUs <- function(currObj, otuTH = 1){
      totalOTUs <- taxa_sums(currObj)
      results[[1]] <- length(totalOTUs[totalOTUs>0])
      results[[2]] <- names(totalOTUs[totalOTUs>0])
    }
    
    # (2)
    uniqueOTUs_in3m <- function(currObj, time, otuTH){
      latestMonth <- 6
      results <- list()
      results[1] <- FALSE
      sampleIndx <- 1
      if ((time!=0) & (time!=latestMonth)) {
        try(currObj <- prune_samples(sample_data(currObj)$Month %in% c(time-1, time+1), currObj))
        if(inherits(currObj, "try-error")){
          results[1] <- TRUE
        }
      } else {
        if (time==latestMonth) {
          try(currObj <- prune_samples(sample_data(currObj)$Month==(time-1), currObj))
          if(inherits(currObj, "try-error")) {
            results[1] <- TRUE
          }
        } else {
          try(currObj <- prune_samples(sample_data(currObj)$Month==(time+1), currObj))
          if(inherits(currObj, "try-error")) {
            results[1] <- TRUE
          }
        }
      }
      noOfsamples <- dim(currObj)[2]
      totalOTUs <- taxa_sums(currObj)
      
      results[[2]] <- length(totalOTUs[totalOTUs>0])
      results[[3]] <- names(totalOTUs[totalOTUs>0])
      results
    }
    
#     # (3)
#     bootStrapping <- function(currDF_M, uStudyIDs, sampleSize = 10){
#       momsSample <- sample(uStudyIDs[uStudyIDs!=id], sampleSize)
#       counter <- 1
#       nsUniqueOTUsList <- list()
#       for (id in momsSample){
#         tb_M_nonSelf <- subset_samples(currDF_M, StudyID==id)
#         ns_MomOTUs <- uniqueOTUs(tb_M_nonSelf, 1)
#         nsUniqueOTUsList[counter] <- ns_MomOTUs
#         counter <- counter + 1
#       }
#     }
    
    ####----------------------------------------------------------------------------------------
    
# Calculation of shared OTU proportion
############################################################
    # Which of the subsets has less samples? use it in the loop!
    C_uStudyIDs <- sort(unique(sample_data(tb_C)$StudyID))
    M_uStudyIDs <- sort(unique(sample_data(tb_M)$StudyID))
    uStudyIDs <- C_uStudyIDs[C_uStudyIDs %in% M_uStudyIDs]
    # Loads the table with the indication what type of C-section occured
    load("cSecType_ext") # Data loaded as "CsecType"
    
# Initialize a list of data.frames to store results
sharedOTUsOutList <- list()
varNames = c("Fecal-All","Vaginal-All", "Vaginal-preNatal", "Vaginal-postNatal")
# varNames = c("Vaginal-All", "Vaginal-preNatal", "Vaginal-postNatal")
sampleTypeDFlist <- list(tb_M_fecal, tb_M_Vswab, tb_M_Vpre, tb_M_Vpost)
# sampleTypeDFlist <- list(tb_M_Vswab, tb_M_Vpre, tb_M_Vpost)
names(sampleTypeDFlist) <- varNames

##### Begining of main looping section ###############################################################
for (df in 1:length(sampleTypeDFlist)) {
  currMomSampleType <- sampleTypeDFlist[[df]]
  print(sprintf("%s", names(sampleTypeDFlist)[df]))
  # timePoints <- c(0, 1, 2, 3, 4, 5, 6)
  latestMonth <- 6
  timePoints <- c(0, 3, 6)
  timePoints <- c(0, 2, 5)
  # Month, ID, number of unique OTUs in the mother, number of unique OTUs in the child,
  # Mother->child dyad, non-self mother -> child, comparing two adjacent child samples
  colNames = c("timePoint", "StudyID", "delivery" ,"uM", "uC", "M_Cdyad", "nonSelfMom")
  sharedOTUtable <- data.frame(
    matrix(nrow = length(uStudyIDs)*length(timePoints), ncol= length(colNames)))
  names(sharedOTUtable) = colNames
  
  counter <- 1
  otuTH <- 1
  for (id in uStudyIDs) {
    print(sprintf("%d", id))
    ### place mother here
    
    # Subsets the current mother, add try to treat cases when mother does no exist 
    # for this time point and sub-caterory (like missing stool sample within a month..)
    tb_dyadM_StudyID <- try(subset_samples(currMomSampleType, StudyID==id))
    skipRecord <- FALSE
    if (inherits(tb_dyadM_StudyID, "try-error")){
      #     sharedOTUtable[counter,] <- c(time, id, delivery, rep("NA", length(colNames)-3))
      #     counter <- counter + 1
      #     next
      skipRecord <- TRUE
    } else {
      Mom_OTUs <- uniqueOTUs(tb_dyadM_StudyID)
      # Here we get a list of unique OTUs for n random mothers
      sampleSize <- 10
      momsSample <- sample(uStudyIDs[uStudyIDs!=id], sampleSize)
      j <- 1
      nonSelfSharedOTUsList <- list()
      # momsSample is a random-generated list of moms sample ids
      for (sid in momsSample){
        try(tb_M_nonSelf <- subset_samples(currMomSampleType, StudyID==sid))
        if (inherits(tb_M_nonSelf, "try-error")){
          next
        } else {
          ns_MomOTUs <- uniqueOTUs(tb_M_nonSelf, 1)
          nonSelfSharedOTUsList[[j]] <- ns_MomOTUs
          j <- j + 1
        }
        
      }
    }
    
    for (time in timePoints) {
      # print(sprintf("%d", time))
      if (skipRecord){
        sharedOTUtable[counter,] <- c(time, id, delivery, rep("NA", length(colNames)-3))
        counter <- counter + 1
        next
      } else {
        # Subsets the current child
        tb_C_StudyID <- subset_samples(tb_C, StudyID==id)
        # This function get phyloseq subset and output array of unique OTU names
        Child_OTUs <- uniqueOTUs_in3m(tb_C_StudyID, time, otuTH)
        # In case the selected sample has zero OTUs, just write 0
        if ( !(Child_OTUs[[2]]==0) & !(Mom_OTUs[[1]]==0) ){
          uniqueChild_OTUs <- Child_OTUs[[3]]
          
          # Foreach random mom, compare her unique OTUs with the ones of the current child
          nonSelfSharedOTUs <- mapply(function(currDF, currChild_OTUs){
            sum(currChild_OTUs %in% currDF)}, 
            nonSelfSharedOTUsList, uniqueChild_OTUs, SIMPLIFY = TRUE)
          
          dyadSharedOTUs <- sum(uniqueChild_OTUs %in% Mom_OTUs)
          # Stores how many OTUs the mather and child have
          uC <- length(uniqueChild_OTUs)
          uM <- length(Mom_OTUs)
          delivery <- CsecType[CsecType$Id==id, 2]
        }
        else {
          if ((Child_OTUs[[2]]==0) & !(Mom_OTUs[[1]]==0)){
            uC <- 0
            dyadSharedOTUs <- 0
            uM <- NA
          }
          else {
            uM <- 0
            dyadSharedOTUs <- 0
            uC <- NA
          }
        }
      }
      sharedOTUtable[counter,] <- c(time, id, delivery, uM, uC, dyadSharedOTUs, mean(nonSelfSharedOTUs))
      counter <- counter + 1
    }
  }
sharedOTUsOutList[[df]] <- sharedOTUtable
}
# save(sharedOTUsOutList,
#      file=paste("sharedOTUsOutList_final"))
#      