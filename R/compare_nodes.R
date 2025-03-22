#' @title Compare Nodes
#'
#' @description This function compares two nodes to determine if all elements of the first node (`node_1`) are present within the accolites of the second node (`node_2`). Accolites are defined as a set of elements retrieved from the `get_accolites` function, which operates on `node_1` within a specified window of data. The comparison is successful if all elements in `node_1` are found in the accolites of `node_2`.
#'
#' @param node_1 A vector representing the first node whose elements will be compared against the accolites of `node_2`.
#' @param node_2 A vector representing the second node, from which accolites will be derived for comparison.
#'
#' @return A logical value. The function returns `TRUE` if all elements of `node_1` are found in the accolites of `node_2`. It returns `FALSE` otherwise.
#'
#' @export
compare_nodes <- function(node_1, node_2){
  accolites_1    <- lapply(node_1, get_accolites, window_data, 50) %>% unlist()
  common_elements <- (node_2 %in% accolites_1) %>% sum()
  if(length(node_1) == common_elements){
    TRUE # if the number of common elements is equal to node_1 cardinality
  }else{
    FALSE
  }
}