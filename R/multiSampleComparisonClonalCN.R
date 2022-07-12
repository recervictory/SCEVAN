
plotAllClonalCN <- function(samples, name){
  
  #CNVtot <- lapply(samples, function(i) read.table(paste0("./output/"," ",i," _  _CN.seg"), sep="\t", header=TRUE, as.is=TRUE))
  CNVtot <- lapply(samples, function(i) read.table(paste0("./output/",i,"_Clonal_CN.seg"), sep="\t", header=TRUE, as.is=TRUE))
  
  png(paste0("./output/",name,"_compareClonalCN.png",sep=""), height=2250, width=1350, res=200)
  
  par(mfrow=c(length(samples),1),cex=1, cex.main = 1.5, cex.lab = 1.5,xaxs="i")
  
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
#' @param organism organism to be analysed (default human)
#' @param par_cores number of cores (default 20)
#'
#' @return
#' @export
#'
#' @examples 
#' 
multiSampleComparisonClonalCN <- function(listCountMtx, analysisName = "all", organism = "human" , par_cores = 20){

  resList <- lapply(names(listCountMtx), function(x) {
    pipelineCNA(listCountMtx[[x]], sample = x, SUBCLONES = FALSE, ClonalCN = TRUE, par_cores = par_cores, organism=organism)
  })
  names(resList) <- names(listCountMtx)
  
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
  
  oncoHeat <- annoteBandOncoHeat(annot_mtx, diffList, length(names(listCountMtx)))
  
  rownames(oncoHeat) <- names(listCountMtx)
  
  plotOncoHeatSubclones(oncoHeat, length(names(listCountMtx)), analysisName, NULL)
  
  plotAllClonalCN(names(listCountMtx), name = analysisName)
  
  for(i in 1:length(names(listCountMtx))){
    names(diffList) <- gsub(paste0("subclone",i), names(listCountMtx)[i], names(diffList))
  }
  
  names(diffList) <- gsub("clone", "shared", names(diffList))

  outputAnalysis <- list(resList, diffList)
  
  save(outputAnalysis, file = paste0("./output/",analysisName,"_outputAnalysis.RData"))
  
  outputAnalysis
  
}

