#set working directory that consists distance matrix files in phylip format
library(phangorn)
library(dendextend)

dm_dn <- readDist(file="nd5_dn.txt", format = "phylip")
dm_in <- readDist(file="nd5_in.txt", format = "phylip")
dm_co <- readDist(file="nd5_clustal.mat", format = "phylip")


#par(mar=c(2.1, 4.1, 4.1, 2.1))

hc_dn <- hclust(as.dist(dm_dn), method = "average")
hcd_dn <- as.dendrogram(hc_dn)
#plot(hcd_dn, horiz=TRUE, main = "R/P-F DN", edge.root = TRUE)

hc_in <- hclust(as.dist(dm_in), method = "average")
hcd_in <- as.dendrogram(hc_in)
#plot(hcd_in, horiz=TRUE, main = "R/P-F IN", edge.root = TRUE)

hc_co <- hclust(as.dist(dm_co), method = "average")
hcd_co <- as.dendrogram(hc_co)
#plot(hcd_co, horiz=TRUE, main = "CLUSTAL O", edge.root = TRUE)


#--------------------------------------------------
cor_cophenetic(hcd_dn, hcd_co)
cor_cophenetic(hcd_in, hcd_co)

#Robinson-foulds
t_dn <- as.phylo(hcd_dn)
t_in <- as.phylo(hcd_in)
t_co <- as.phylo(hcd_co)
treedist(t_dn,t_co)
treedist(t_in,t_co)





