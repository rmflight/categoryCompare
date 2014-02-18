readGMT <- function(gmtfile){
  gmtContents <- scan(gmtfile, what="char", sep="\n", encoding="latin1")
  nSet <- length(gmtContents)
  gmtID <- as.character(seq(1, nSet))
  oldID <- character(nSet)
  gmtDesc <- character(nSet)
  
  gmtSplit <- strsplit(gmtContents, "\t", fixed=F)
  
  gmtData <- lapply(seq(1, nSet), function(inLoc){
    oldID[inLoc] <<- gmtSplit[[inLoc]][1]
    gmtDesc[inLoc] <<- paste(gmtSplit[[inLoc]][1], gmtSplit[[inLoc]][2], sep="::")
    nSplit <- length(gmtSplit[[inLoc]])
    unique(gmtSplit[[inLoc]][seq(3, nSplit)])
  })
  names(gmtData) <- gmtID
  old2new <- gmtID
  names(old2new) <- toupper(oldID)
  names(gmtDesc) <- gmtID
  return(list(id=gmtID, old2new=old2new, description=gmtDesc, data=gmtData))
}

readGSEAResults <- function(fileList, old2new, pUse="NOM p-val", pCut=0.01){
  pUse <- make.names(pUse)
  outResults <- lapply(fileList, function(x){
    tmpRes <- read.table(x, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
    sigLoc <- which(names(tmpRes) %in% pUse)
    sigEntry <- tmpRes[,sigLoc] <= pCut
    sigID <- tmpRes[sigEntry, 'NAME']
    newID <- old2new[sigID]
    names(newID) <- NULL
    newID <- newID[!is.na(newID)]
    return(new("ccSigList", sigID=newID))
  })
  return(outResults)
}