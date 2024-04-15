
#' plotAllSubclonalCN Plot the copy number of each subclone of a sample.
#'
#' @param sample Name of the sample.
#' @param pathOutput Path to the output folder containing the output of pipelineCNA.
#'
#' @return
#' @export
#'
#' @examples 
#' 
plotAllSubclonalCN <- function(sample, pathOutput = "./output/"){
  
  allFile <- list.files(pathOutput, pattern = paste0(sample,"_subclone[1-9]_CN.seg") )

  CNVtot <- lapply(allFile, function(i) read.table(paste0(pathOutput,i), sep="\t", header=TRUE, as.is=TRUE))
  
  png(paste0(pathOutput,sample,"_compareSubclonalCN.png",sep=""), height=2250, width=1350, res=100)
  
  par(mfrow=c(length(allFile),1),cex=1, cex.main = 1.5, cex.lab = 1.5,xaxs="i")
  
  for(i in 1:length(allFile)){
    
    CNV <- CNVtot[[i]]
    
    plotSegmentation(CNV) 
    title(gsub("_CN.seg","",allFile[i]))
    
  }
  dev.off()
}


#' Title plotAllClonalCN
#'
#' @param samples Vector with sample names to be plotted
#' @param name Analysis name
#'
#' @return
#' @export
#'
#' @examples
plotAllClonalCN <- function(samples, name){
  
  #CNVtot <- lapply(samples, function(i) read.table(paste0("./output/"," ",i," _  _CN.seg"), sep="\t", header=TRUE, as.is=TRUE))
  CNVtot <- lapply(samples, function(i) read.table(paste0("./output/",i,"_Clonal_CN.seg"), sep="\t", header=TRUE, as.is=TRUE))
  
  png(paste0("./output/",name,"_compareClonalCN.png",sep=""), height=1550, width=2350, res=100)
  
  if(length(samples)>3){
    par(mfrow=c(3,ceiling(length(samples)/3)), cex=1, cex.main = 1.5, cex.lab = 1.5,xaxs="i")
  }else{
    par(mfrow=c(length(samples),1),cex=1, cex.main = 1.5, cex.lab = 1.5,xaxs="i")
  }
  
  for(i in 1:length(samples)){
    
    CNV <- CNVtot[[i]]
    
    plotSegmentation(CNV) 
    title(samples[i])
    
  }
  dev.off()
}

#' multiSampleComparisonClonalCN Compare the clonal Copy Number of multiple samples.
#'
#' @param listCountMtx Named list of raw count matrix of samples
#' @param analysisName Name of the analysis (default "all")
#' @param organism Organism to be analysed (optional - "mouse" or "human" - default "human")
#' @param par_cores number of cores (default 20)
#'
#' @return
#' @export
#'
#' @examples 
#' 
multiSampleComparisonClonalCN <- function(listCountMtx, 
                                          listNormCells = NULL, 
                                          analysisName = "all", 
                                          organism = "human" , 
                                          par_cores = 20, 
                                          plotTree = TRUE,
                                          SUBCLONES = FALSE
                                          ) {

  # resList <- lapply(names(listCountMtx), function(x) {
  #   pipelineCNA(listCountMtx[[x]], 
  # norm_cell = listNormCells[[x]], 
  # sample = x, 
  # SUBCLONES = FALSE, 
  # ClonalCN = TRUE, 
  # par_cores = par_cores, 
  # organism=organism)
  # })

  dir <- "samples_scevan/"

  if (!dir.exists(dir)) {
  dir.create(dir)
  }

  for(x in names(listCountMtx)) {
  cat(paste("\n 1. Current File name ....\n",x))
  # Construct the file path with the sample name
  file_path <- paste0(dir, x, ".rds")
  
  # Check if the file already exists
  if(file.exists(file_path)) {
    # Skip the rest of this loop iteration if the file exists
    cat(paste("\n 2.", x ,"File Already exists....\n"))
    next

  }
  
  # # Execute the pipelineCNA function only if the file does not exist
  # result <- pipelineCNA(
  #   listCountMtx[[x]],
  #   sample = x,
  #   par_cores = par_cores,
  #   norm_cell = listNormCells[[x]],
  #   SUBCLONES = SUBCLONES,
  #   beta_vega = 0.5,
  #   ClonalCN = TRUE,
  #   plotTree = plotTree,
  #   AdditionalGeneSets = NULL,
  #   organism = organism
  # )

# Maximum retries till cores reduce to 1

run_cores <- par_cores  # Initialize 'run_cores' with 'par_cores' value
run_count <- 3  # Set the maximum number of attempts

while (run_cores >= 1 && run_count > 0) {
  if (file.exists(file_path)) {
    print("Result file already exists, no need to run the pipelineCNA function.")
    break
  }

  tryCatch({
    # Deliberate error for testing if run_cores is exactly 5
    if (run_cores == 5) {
      stop("Deliberate error for testing.")
    }

    # Attempt to run the pipelineCNA function with the current number of cores
    result <- pipelineCNA(
      listCountMtx[[x]],
      sample = x,
      par_cores = run_cores,
      norm_cell = listNormCells[[x]],
      SUBCLONES = SUBCLONES,
      beta_vega = 0.5,
      ClonalCN = TRUE,
      plotTree = plotTree,
      AdditionalGeneSets = NULL,
      organism = organism
    )
    print(paste("Success with", run_cores, "cores."))
    break  # Exit the loop on success
  }, error = function(e) {
    print(paste("Failed with", run_cores, "cores. Error:", e$message))
    run_count <- run_count - 1
  })

  # Decrement cores outside the tryCatch if an error occurred and not yet at minimum
  if (run_cores > 1 && run_count > 0) {
    run_cores <- max(run_cores - 5, 1)
    print(paste("Retrying with", run_cores, "cores. Remaining attempts:", run_count))
  } else if (run_cores <= 1 || run_count <= 0) {
    print("No more attempts available or running on minimum number of cores. Stopping execution.")
    break
  }
}


# Ensured 'run_cores' is being decremented correctly and messages use 'run_cores' instead of 'par_cores'.

  # Save the result as an RDS file
  saveRDS(result, file = file_path)
  
  # Remove the result object from memory
  rm(result)
  
  # Call the garbage collector to reclaim unused memory
  gc()
}

  cat("ALL samples saved in the disk ....")
  # Get list of .rds files in the specified directory
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  resList <- list()  # Initialize empty list

  for (i in 1:length(files)) {
  # Read the .rds file
  resList[[i]] <- readRDS(files[i])
  
  # Extract the base name without the extension and set it as the list element name
  names(resList)[i] <- file_path_sans_ext(basename(files[i]))
  
  # Output message about file addition
  cat(paste("\n 2.", i, "File Added to list as '", names(resList)[i], "'\n", sep=""))
}

  # names(resList) <- names(listCountMtx)
  runCompleted <- intersect(names(resList), names(listCountMtx))

  listCountMtx <- listCountMtx[runCompleted]

  cat("Completed Run Samples............")
  print(names(listCountMtx))
  
  sampleAlterList <- lapply(names(listCountMtx), function(x) {
    analyzeSegm(x, nSub = 0)
  })
  names(sampleAlterList) <- paste0(names(listCountMtx),"_subclone", 1:length(names(listCountMtx)))
  
  names(sampleAlterList) <- paste0(analysisName,"_subclone", 1:length(names(listCountMtx)))
  
  diffList <- diffSubclones(sampleAlterList, analysisName, nSub = length(names(listCountMtx)))
  
  diffList <- testSpecificAlteration(diffList, length(names(listCountMtx)), analysisName)
  
  genesMtx <- lapply(listCountMtx, function(x) rownames(x))
  
  genesMtx <- sort(unique(unlist(genesMtx)))
  
  genesMtx <- data.frame(x = rep(0, length(genesMtx)), row.names = genesMtx)
  
  annot_mtx <- annotateGenes(genesMtx)
  
  oncoHeat <- annoteBandOncoHeat(annot_mtx, diffList, length(names(listCountMtx)), organism)
  
  rownames(oncoHeat) <- names(listCountMtx)
  
  plotOncoHeatSubclones(oncoHeat, length(names(listCountMtx)), analysisName, NULL, organism)
  
  plotAllClonalCN(names(listCountMtx), name = analysisName)
  
  if(length(names(listCountMtx))>2 & plotTree) plotCloneTreeNew(names(listCountMtx), CLONAL_MULTI = TRUE, analysisName = analysisName)
  
  for(i in 1:length(names(listCountMtx))){
    names(diffList) <- gsub(paste0("subclone",i), names(listCountMtx)[i], names(diffList))
  }
  
  names(diffList) <- gsub("clone", "shared", names(diffList))

  outputAnalysis <- list(resList, diffList)
  
  save(outputAnalysis, file = paste0("./output/",analysisName,"_outputAnalysis.RData"))
  
  outputAnalysis
  
}


