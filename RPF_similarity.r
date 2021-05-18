#set the working directory
#load data

library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

#method for calculating LF entropy
getEntropy <- function(seq, pat, patternExists)
{	
	entropy <- 0
	
	try({
	
	if(length(patternExists)>0 && patternExists == TRUE)
	{
		seqID <- seq
		patternID <- pat
		seqLenght <- resSeqLen[resSeqLen$ID == seqID,]$SEQLENGTH
		
		positions <- resPatternPositions[resPatternPositions$SEQID == seqID & resPatternPositions$ID_FRAGMENT1 == patternID,]
		locationVector <- sort(positions$LOCATION)
		positionsVector <- numeric(seqLenght)
			
		positionsVector[locationVector[]] = locationVector[] 
			
		lfVector <- double(seqLenght)
		lfVector[locationVector[1]] = 1 /(positionsVector[locationVector[1]])
		lfVector[locationVector[2: length(locationVector)]] = 1 /(positionsVector[locationVector[2: length(locationVector)]]-positionsVector[locationVector[1: (length(locationVector)-1)]])
			
		pslfVector <- cumsum(lfVector)
		z <- sum(pslfVector)
		dpVector <- double(length(pslfVector[pslfVector != 0])) 
		dpVector <-  pslfVector[pslfVector != 0]/z
		dplogVector <- vapply(dpVector,log2,double(1))
			
		entropy <- -sum(dpVector*dplogVector)
		
	}
	})
	return(entropy)	
}


#calucalete LF entropy for every sequence/repeat in the analysis

clusterExport(cl,"getEntropy")
clusterExport(cl,"resPatternFreq")
clusterExport(cl,"resPatternPositions")
clusterExport(cl,"resSeq")
clusterExport(cl,"resSeqLen")
clusterExport(cl,"resVocabular")
clusterExport(cl,"seqVectorsF_P")
clusterExport(cl,"totalSeq")
clusterExport(cl,"totalVoc")

system.time(seqVectorsF_P <- vapply(1:totalSeq, function(i) parSapply(cl, 1:totalVoc, function(j) getEntropy(resSeq[i,], resVocabular[j,]$ID, (is.na(resPatternFreq[resPatternFreq$SEQID==resSeq[i,] & resPatternFreq$ID_FRAGMENT1==resVocabular[j,]$ID,1:2]$ID_FRAGMENT1)==FALSE))), double(totalVoc)))
rownames(seqVectorsF_P) <- resVocabular[]$ID
colnames(seqVectorsF_P) <- resSeq[]$ID

stopCluster(cl)
#sequence vectors are ready for further analysis 

#__________________________________________________________________________
#similarity matrix and clustering
#install.packages("lsa")
library(lsa)
seqVectorsF_P_Host <- seqVectorsF_P
colnames(seqVectorsF_P_Host) <- resHost[]$HOST
cosSimMatrix <- cosine(data.matrix(seqVectorsF_P_Host))
cosDisimMatrix <- abs(cosSimMatrix-1)

#install.packages("d3heatmap")
library(d3heatmap)

hc <- hclust(as.dist(cosDisimMatrix))
hcd <- as.dendrogram(hc)
plot(hcd, cex = 0.6)
d3heatmap(cosDisimMatrix, Rowv = hcd, Colv = "Rowv")

#install.packages("factoextra")
#install.packages("purrr")
library(purrr)
library(factoextra)
res.hc <- hclust(as.dist(cosDisimMatrix))
grp <- cutree(res.hc, k = 2)
plot(res.hc, cex = 1, main = "direct non-complementary")
rect.hclust(res.hc, k = 2, border = 2:5)
fviz_cluster(list(data = cosDisimMatrix, cluster = grp), repel = TRUE, ellipse.type = "norm")