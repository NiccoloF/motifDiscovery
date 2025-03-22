#' @title Find Recommended Path in a Tree Structure
#'
#' @description
#' This function identifies a recommended path in a tree structure based on user-defined minimum criteria. 
#' It filters nodes based on the number of leaves and calculates scores for potential recommendations.
#'
#' @param minidend A tree structure representing the hierarchy of nodes, 
#'                 from which leaves will be derived for analysis.
#' @param window_data A matrix or data frame containing data associated with each node.
#' @param min_card An integer specifying the minimum number of leaves that must be present in a node 
#'                  for it to be considered for recommendation.
#'
#' @return A list containing:
#'         - \code{seed_path_info}: A data frame summarizing the valid seed nodes, their number of leaves, 
#'         heights, and recommended nodes.
#'         - \code{seed_path_list}: A list containing the names of elements for each node in the seed path.
#'         - \code{score_path_list}: A list of scores adjusted for each node in the seed path.
#'         - \code{recommended_node_labels}: A list of labels for the recommended nodes.
#'         - \code{recommended_node_scores}: A numeric vector containing scores for the recommended nodes.
#'
#' @details
#' The function follows these steps:
#' 1. Extracts the leaves and heights of each node in the tree.
#' 2. Filters nodes based on the minimum number of leaves specified by \code{min_card}.
#' 3. Identifies nodes that are parents of the filtered nodes, marking them for deletion.
#' 4. For each seed node, calculates crossing points with other seed nodes.
#' 5. Calculates adjusted scores for nodes in the seed path using the C++ function \code{fMSR_adj}.
#' 6. Recommends nodes based on calculated scores and returns a summary of results.
#'
#' @importFrom dplyr %>%
#' @importFrom data.table setnames
#' @importFrom Rcpp evalCpp
#' @export
find_recommended_path <- function(minidend, window_data, min_card){
  
  # get leaves for each node
  node_list <- minidend %>% partition_leaves()
  # get number of leaves in each node as a vector
  node_leaves  <- node_list %>%
    lapply(length) %>%
    unlist() %>%
    data.table::as.data.table() %>%
    data.table::setnames('node_leaves')
  
  # get branches heights
  node_heights <- minidend %>%
    get_nodes_xy() %>%
    data.table::as.data.table() %>%
    data.table::setnames(c('x','y'))
  
  node_heights <- node_heights[,'y']
  
  # bind together on a dataframe with deep-first search order
  temp <- cbind(node_leaves, node_heights) %>%
    mutate(depth_order = 1:length(node_leaves)) # deep-first search order
  
  # generate for each nodes (using function get_parents)
  node_parents <- lapply(node_list, get_parents, node_list = node_list )
  
  # minimum cardinality (minimum number of node leaves - selected by the user)
  min_card <- min_card
  colnames(temp) <- c('node_leaves','y','depth_order')
  # add column interesting according to the number of leaves and the min_card
  temp <- temp %>%
    filter(node_leaves >= min_card) %>%
    dplyr::arrange(node_leaves, y)
  
  # no node with at least min_card branches -> EXIT (return list with NULL)
  if(dim(temp)[1] == 0){
    return(NULL)
  }
  
  # find the nodes that need to be deleted because parents of minimal ones
  delete_nodes <- c()
  for(i in 1:dim(temp)[1]){
    node_index   <- temp$depth_order[i] # focus node
    delete_now   <- node_parents[[node_index]] # all parents nodes of the focus
    delete_nodes <- c(delete_nodes,
                      delete_now[delete_now != node_index])
  }
  
  # seed nodes
  temp <- temp %>%
    filter(!(depth_order %in% unique(delete_nodes)))
  
  # crossing points between seed nodes
  seed_parents <- node_parents[temp$depth_order]
  all_parents  <- seed_parents %>% unlist() %>% unique()
  
  seed_path <- list()
  for(i in 1:length(seed_parents)){
    my_parent      <- seed_parents[[i]]
    the_rest       <- seed_parents[-i] %>% unlist() %>% unique() %>% sort()
    if(!(is.null(the_rest))){
      cross_point    <- intersect(my_parent, the_rest) %>% max()
      seed_path[[i]] <- my_parent[my_parent > cross_point]
    }else{
      seed_path[[i]] <- my_parent
    }
  }
  
  # add column with parents till the new crossing
  temp$parents <- lapply(seed_path, toString) %>% unlist()
  
  # get the list of elements for each node in the seed_path
  seed_path_list <- lapply(seed_path, function(x){node_list[x]}) # names of elements
  
  
  cppFunction('double fMSR_adj(NumericMatrix& mat) {
  
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  const arma::mat matRef(mat.begin(), nrow, ncol,false,true);
  if(!nrow) // degenerate case with zero curves
    return 0.0;
  
  double score = arma::accu(arma::square(
                 (matRef.each_col() - arma::mean(matRef, 1)).each_row() - arma::mean(matRef,0)
                 + arma::mean(arma::vectorise(matRef))))/static_cast<double>(nrow * ncol);
  if (nrow > 2) 
  {
    score /= arma::as_scalar(arma::prod(arma::regspace<arma::vec>(2, 1, nrow).transform(
      [](double k){return k * k / (k * k - 1);})));
  }
  return score;
}',depends="RcppArmadillo",includes="#include <ranges>",plugins="cpp20")
  
  # get the list of h-score adjusted for each node in the seed_path
  score_path_list <- lapply(seed_path, function(x){
    lapply(x, function(x){window_data[node_list[[x]],] %>% fMSR_adj()})}
  )
  
  recommendation <- lapply(score_path_list, recommend_node) %>% unlist()
  temp$recommended <- recommendation
  
  # get information about recommended nodes (labels and scores)
  recommended_node_labels <- list()
  recommended_node_scores  <- c()
  for(i in 1:length(seed_path_list)){
    recommended_index <- recommendation[i]
    recommended_node_labels[[i]] <- seed_path_list[[i]][[recommended_index]]
    recommended_node_scores[i]   <- score_path_list[[i]][[recommended_index]]
  }
  
  
  res <- list("seed_path_info"  = temp,
              "seed_path_list"  = seed_path_list,
              "score_path_list" = score_path_list,
              "recommended_node_labels" = recommended_node_labels,
              "recommended_node_scores" = recommended_node_scores
  )
  return(res)
}