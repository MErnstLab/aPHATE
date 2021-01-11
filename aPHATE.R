library(dplyr)
#' Transform circular coordinates to euclidean
#' 
#' 
#' @param data matrix(n_samples, n_coordinates)
#' 2 dimensional input data array with n_samples complexes and n_coordinates
#' coordinates.
#' Coordinates have 4 columns per atom, where 3 describe the angles between
#' the atom and the three axes in degrees, while the last is the distance
#' of the atom from the coordinate system origin. 
#' 
#' @return "euclidean" matrix with all atoms in euclidean coordinates
#' for all samples
#' 
#' @export
Euclidean <- function(data){
  euclidean <- matrix(data = NA,
                      nrow = nrow(data),
                      ncol = ncol(data)*3/4)
  
  
  for (i in seq(0, (ncol(data)/4))-1) {
    coordinates <- data[,seq((i*4)+1, (i*4)+4)]
    euclidean[,(i*3)+1] <- coordinates[,4] * cos(coordinates[,1]*pi/180)
    euclidean[,(i*3)+2] <- coordinates[,4] * cos(coordinates[,2]*pi/180)
    euclidean[,(i*3)+3] <- coordinates[,4] * cos(coordinates[,3]*pi/180)
  }
  
  return(euclidean)
}


#' Calculates root mean square deviation between all the samples
#' from euclidean coordinates.
#' 
#' @param data matrix'(n_samples, n_coordinates)
#' 2 dimensional input data array with n_samples complexes and n_coordinates
#' coordinates.
#' Coordinates are in the euclidean coordinate system, so the matrix
#' has n_atoms * 3 columns.
#' 
#' @return "RMSD" matrix(n_samples, n_samples)
#' 2 dimensional output data array with n_samples rows and columns.
#' Each value represents root mean square deviation between two samples.
#' 
#' @export

RMSD <- function(data){
  return(as.matrix(dist(data, method = "euclidean")) / sqrt(length(data)))
}


#' #Performs hierarchical clustering on the data, with cutoff and cluster
#' shrinkage. In the shrinkage step all the clusters with less than n_poses
#' poses per compound are discarded.
#' 
#' @param data matrix(n_samples, n_samples)
#' Matrix with distances between samples. Expected is the output from RMSD
#' function.
#' 
#' @param samples c(n_samples)
#' Factor containing for each row/column in data a string linking it to the
#' sample.
#' 
#' @param cutoff numeric
#' The number of included leaves when cutting the hierarchical tree.
#' 
#' @return list
#' A list consisting of Clusters data.frame with cluster information and
#' Data matrix, which is the same as the input data, but with samples in 
#' discarded clusters removed. 
#' 
#' @export
hierclust <- function(data, samples, cutoff=80){
  hier <- hclust(as.dist(data))
  data_frame <- data.frame(Samples = samples,
                           Clusters = factor(cutree(hier, cutoff)))
  x <- data_frame %>%
    dplyr::group_by(Samples, Clusters) %>%
    dplyr::summarise(Pose=n() > 10)
  included <- x %>%
    dplyr::group_by(Clusters) %>%
    dplyr::summarise(Inclusion = sum(Pose) > 0)
  
  #included <- data_frame %>%
  #  dplyr::group_by(Samples, Clusters) %>%
  #  dplyr::summarise(Pose=n() > 10) %>%
  #  dplyr::group_by(Clusters) %>%
  #  dplyr::summarise(Inclusion = sum(Pose) > 0)
  
  data_frame$Inclusion <- included$Inclusion[data_frame$Clusters]
  
  
  
  return(list(Clusters = data_frame$Clusters,
              Inclusion=data_frame$Inclusion))
}


#' Calculates the Earth mover's distance
#' 
#' @param pca matrix(n_samples, n_components)
#' 2 dimensional matrix with as many principal components as need to be
#' included in the distance calculation.
#' @param samples c(n_samples)
#' A factor delineating samples. 
#' @param clusters c(n_samples)
#' A factor holding information on which sample pertains to which cluster.
#' @param scoring c(n_samples)
#' A factor containing the score for each sample.
#' 
#' @return "dist" matrix(n_samples, n_samples)
#' 2 dimensional matrix holding the distances between samples.
#' 
#' @export
emdistance <- function(pca, samples, clusters, scoring){
  initial_data <- data.frame(Samples=samples,
                             Clusters=clusters,
                             Scoring=scoring)
  distances <- as.matrix(dist(aggregate(pca,
                                        list(initial_data$Clusters),
                                        mean)))
  cluster_count <- initial_data %>%
    group_by(Samples, Clusters) %>%
    summarise(SCount=sum(Scoring))
  cluster_count_cast <- reshape2::dcast(cluster_count,
                                        Samples~Clusters,
                                        value.var = "SCount")
  cluster_count_cast[is.na(cluster_count_cast)] <- 0
  cluster_count_cast[-1] <- t(apply(cluster_count_cast[-1],
                                    1,
                                    function(x){x/sum(x)}))
  cluster_weights <- as.matrix(cluster_count_cast[,-1])
  Y <- rep(0, (nrow(cluster_weights)-1)*nrow(cluster_weights)/2)
  counter <- 1
  for (i in seq(1, nrow(cluster_weights))) {
    cur_f1_weights <- cluster_weights[i,]
    for (j in seq(i+1, nrow(cluster_weights))) {
      if (j > nrow(cluster_weights)) break
      cur_f2_weights <- cluster_weights[j,]
      flow <- transport::transport(cur_f1_weights,
                                   cur_f2_weights,
                                   distances,
                                   method="primaldual")
      curdist <- 0
      for (k in seq(1, nrow(flow))) {
        cur_penalty <- distances[flow[k, 1], flow[k, 2]]
        curdist <- curdist + cur_penalty*flow[k, 3]
      }
      Y[counter] <- curdist
      counter <- counter+1
    }
  }
  Y_sq <- pracma::squareform(Y)
  colnames(Y_sq) <- cluster_count_cast$Samples
  rownames(Y_sq) <- cluster_count_cast$Samples
  return(Y_sq)
}


#' Implementation of modified PHATE pipeline from Fabjan et al., 2021
#' Originally published by Moon et al., 2019
#' This workflow starts with circular coordinates of atoms. These get used
#' to perform PCA and hierarchical clustering. The distance used in clustering
#' is root mean square deviation of atoms. Clusters and principal components
#' are then used to calculate earth mover's distance between samples.
#' @○param data list
#' List containing 3 elements:
#' Circular matrix(n_samples, n_coordinates)
#' 2 dimensional matrix containing for every sample the circular coordinates
#' of all atoms used in the calculation. For every atom the matrix contains
#' 4 columns - 3 angles formed by the atom and 3 axes and a distance between
#' the atom and the coordinate origin.
#' Samples c(n_samples)
#' A factor containing sample information.
#' Scoring c(n_samples)
#' Scoring of the data points. Is used by the earth mover's distance
#' calculation.
#' 
#' @return data list
#' List containing 12 elements:
#' The two original elements (Circular and Samples).
#' Euclidean matrix
#' A matrix of coordinates transformed from the circular to the euclidean.
#' RMSDs matrix(n_samples, n_samples)
#' Matrix of pairwise root mean square deviations between data points.
#' Clusters c(n_samples)
#' A factor containing clustering output for all data points.
#' Inclusion c(n_samples)
#' A boolean factor where the included data points have a value of TRUE
#' Clusters_selected
#' A factor containing only the clustering output for the retained data points.
#' Samples_selected
#' A factor containing only the sample information for the retained data points.
#' pca_selected
#' pca matrix containing only information for the retained data points.
#' Scoring_selected
#' A factor with only scoring for the retained data points.
#' EMDistances matrix
#' A matrix with earth mover's distances between the retained samples.
#' 
#' @export
aPHATE <- function(data){
  
  # PCA
  data$pca <- prcomp(data$Circular, scale=TRUE)
  
  
  data$Euclidean <- Euclidean(data$Circular)
  data$RMSDs <- RMSD(data$Euclidean)
  
  clustering <- hierclust(data$RMSDs, data$Samples)
  data$Clusters <- clustering$Clusters
  data$Inclusion <- clustering$Inclusion
  
  data$Clusters_selected <- data$Clusters[data$Inclusion]  
  data$Samples_selected <- data$Samples[data$Inclusion]
  data$pca_selected <- data$pca$x[data$Inclusion,]
  data$Scoring_selected <- data$Scoring[data$Inclusion]
  
  data$Clusters_selected <- as.factor(as.integer(as.factor(as.integer(data$Clusters_selected))))
  
  data$EMDistances <- emdistance(data$pca_selected,
                                 data$Samples_selected,
                                 data$Clusters_selected,
                                 data$Scoring_selected)
  
  
  return(data)
}