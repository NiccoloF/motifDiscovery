#' @title Get Parent Nodes from a Given Node
#'
#' @description
#' This function identifies the parent nodes of a given node in a hierarchical structure. 
#' It checks which nodes in the provided list contain the specified node, thereby determining 
#' its parent nodes.
#'
#' @param node A vector representing a single node whose parents are to be found. 
#'              It is expected to contain the elements that may be present in the parent nodes.
#' @param node_list A list of vectors, where each vector represents a node in the hierarchical structure. 
#'                  The function will search through these nodes to find parents of the specified node.
#'
#' @return A numeric vector containing the indices of the parent nodes in the `node_list` that include the specified node.
#'
#' @details
#' The function works by applying a logical check across the `node_list`. It returns the indices of all nodes 
#' that contain all elements of the specified node vector. If no parent nodes are found, the function will return 
#' an empty integer vector.
#'
#' @export
get_parents <- function(node, node_list){
  which(lapply(node_list, function(x){all(node %in% x)}) %>% unlist())
}