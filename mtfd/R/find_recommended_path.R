#' @title get_path_complete
#'
#' @description Run multiple times probKMA function with different K,c and initializations,
#' with the aim to find a set of candidate motifs.
#' If the folder name_KK_cc is already present and n result files are already present,
#' load them and continue with the n_init-n runs.
#'
#' @param Y0 list of N vectors, for univariate curves y_i(x), or list of N matrices with d columns,
#' for d-dimensional curves y_i(x),  with the evaluation of curves (all curves should be evaluated
#' on a uniform grid). When y_j(x)=NA in the dimension j, then y_j(x)=NA in ALL dimensions
#' @param Y1 list of N vectors, for univariate derivative curves y'_i(x), or
#' list of N matrices with d columns, for d-dimensional derivatibe curves y'_i(x),
#' with the evaluation of the curves derivatives (all curves should be evaluated on a uniform grid).
#' When y'_j(x)=NA in the dimension j, then y'_j(x)=NA in ALL dimensions.
#' Must be provided when diss='d1_L2' or diss='d0_d1_L2'.
#' @param K vector with numbers of motifs that must be tested.
#' @param c vector with minimum motifs lengths that must be tested.
#' @param n_init number of random initialization for each combination of K and c.
#' @param name name of the folders when the results are saved.
#' @param names_var vector of length d, with names of the variables in the different dimensions.
#' @param probKMA_options list with options for probKMA (see the help of probKMA).
#' @param silhouette_align True or False. If True, try all possible alignments between the curve pieces
#' when calculating the adapted silhouette index on the results of probKMA
#' @param plot if TRUE, summary plots are drawn.
#' @return A list containing: K, c, n_init and name;...
#' @return \item{times}{ list of execution times of ProbKMA for each combination of K, c, and n_init}
#' @return \item{silhouette_average_sd}{ list of the mean (silhouette_average) and standard deviation (silhouette_sd) of the silhouette indices for each execution of the ProbKMA function}
#' @author Marzia Angela Cremona & Francesca Chiaromonte

find_recommended_path <- function(minidend, window_data, min_card){
  # get leaves for each node
  node_list <- minidend %>% partition_leaves()
  # get number of leaves in each node as a vector
  node_leaves  <- node_list %>%
    lapply(length) %>%
    unlist() %>%
    as.data.table() %>%
    setnames('node_leaves')
  
  # get branches heights
  node_heights <- minidend %>%
    get_nodes_xy() %>%
    as.data.table() %>%
    setnames(c('x','y'))
  
  node_heights <- node_heights[,'y']

  # bind together on a dataframe with deep-first search order
  temp <- cbind(node_leaves, node_heights) %>%
    mutate(depth_order = 1:length(node_leaves)) # deep-first search order
  
  # generate for each nodes (using function get_parents)
  node_parents <- lapply(node_list, get_parents, node_list = node_list )
  
  # minimum cardinality (minimum number of node leaves - selected by the user)
  min_card <- min_card
  # add column interesting according to the number of leaves and the min_card
  temp <- temp %>%
    filter(node_leaves >= min_card) %>%
    arrange(node_leaves, y)
  
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
}',depends="RcppArmadillo",includes="#include <ranges>",plugin="cpp20")
  
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