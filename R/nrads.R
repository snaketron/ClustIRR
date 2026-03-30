
get_nrads <- function(community_occupancy_matrix, B = 1000) {
    
    check_inputs <- function(community_occupancy_matrix, B) {
        
        # check community_occupancy_matrix
        if(missing(community_occupancy_matrix)) {
            stop("community_occupancy_matrix must be a matrix")
        }
        if(is.matrix(community_occupancy_matrix)==FALSE) {
            stop("community_occupancy_matrix must be a matrix")
        }
        n <- ncol(community_occupancy_matrix)
        if(n==1) {
            stop("ncol(community_occupancy_matrix) must be >1")
        }
        
        # check B
        if(missing(B)) {
            stop("B must be a number > 0")
        }
        if(length(B)!=1) {
            stop("B must be a number > 0")
        }
        if(is.numeric(B)==FALSE) {
            stop("B must be a number > 0")
        }
        if(is.finite(B)==FALSE) {
            stop("B must be a number > 0")
        }
        if(B<=0) {
            stop("B must be a number > 0")
        }
    }
    
    check_inputs(community_occupancy_matrix = community_occupancy_matrix, B = B)
    
    rads <- apply(X = community_occupancy_matrix, MARGIN = 2, 
                  FUN = function(x) {return(sort(x, decreasing = TRUE))})
    max_rank <- min(apply(X = rads, MARGIN = 2, 
                          FUN = function(x) {sum(x!=0)}))
    
    nrads <- RADnormalization_matrix(input = t(rads), max_rank = max_rank, 
                                     average_over = B, sample_in_row = TRUE)
    
    nd <- nrads$norm_matrix
    rownames(nd) <- colnames(rads)
    nd <- reshape2::melt(nd)
    colnames(nd) <- c("rownames", "rank", "norm.abundance")
    nd$sample <- nd$rownames
    return(nd)
}
