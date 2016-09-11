library(diagram) # for plotting purposes
library(igraph)  # for graph tweaks


# Simple plotting tool for matrix
plotmatrix <- function(matt)
{
	# our graph model is uni-directional so
	mat <- matt
	mat[lower.tri(mat)] <- 0
	plotmat(mat, box.size=0.01, shadow.size=0, arr.lcol='blue', arr.width=0.001)
}



SubGraphIt <- function(matt)
{
	ig <- graph.adjacency(matt, mode="undirected", weighted=TRUE)
	
	
	
	
}


