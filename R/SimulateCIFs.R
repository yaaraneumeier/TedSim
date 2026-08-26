#' Simulate Cell Identity Factor Matrix With Lineage Barcodes
#' @param ncells Number of Cells
#' @param phyla Cell state tree
#' @param cif_center Mean of CIFs, default is 1
#' @param Sigma Standard deviation of non-diff CIFs, difault is 0.5
#' @param p_a The asymmetric division rate, default is 0.8
#' @param p_edge The edge transition probability table, default is NULL so that the child edges are chosen with equal possibilities
#' @param n_CIF Total number of SIFs, both non-diff and diff combined
#' @param n_diff Number of diff SIFs
#' @param step Sampling stepsize of diff-SIF Brownian motion
#' @param p_d Dropout rate of CRISPR/Cas9 lineage barcodes
#' @param mu mutation rate of CRISPR/Cas9 lineage barcodes
#' @param N_char number of character sites for CRISPR/Cas9 lineage barcodes
#' @param N_ms number of mutated states for CRISPR/Cas9 lineage barcodes
#' @param unif_on sampling from synthetic uniform-distributed mutated states. When set FALSE, mutated states will be drawn from a real experimental dataset
#' @param SIF_res Optional input for State Identity Factors. If not, the function will generate SIF internally.
#' @param max_walk maximum walk distance of one asymmetric division on the cell state tree
#' @param lambda a num vector that indicates the max and the min value of lambda that weights the additional random walk value
#' @param T_cell optional, a cell division tree
#' @param variable_branch_lengths optional, logical, if TRUE enables variable branch length mode. Default FALSE.
#' @param branch_length_dist optional, character, distribution for branch lengths: "exponential", "lognormal", or "gamma". Default "exponential".
#' @param branch_length_params optional, list, parameters for the distribution. Default list(rate=1) for exponential.
#' @param branch_length_seed optional, integer, seed for branch length generation. If NULL, uses current random state.
#' @param scale_state_walk optional, logical, if TRUE scales state walk with branch length. Default FALSE.
#' @param state_walk_method optional, character, method for scaling: "poisson" or "repeat". Default "poisson".
#' @param walk_rate optional, numeric, expected steps per unit time for Poisson method. Default NULL (uses max_walk/2).
#' @param scale_barcode_mutations optional, logical, if TRUE scales barcode mutations with branch length. Default FALSE.
#' @param barcode_method optional, character, method for scaling: "scale_mu" or "repeat". Default "scale_mu".
#' @param ultrametric optional, logical, if TRUE makes tree ultrametric (all leaves equidistant from root). Default FALSE.
#' @param lambda_scaling optional, character, lambda scaling method: "depth" or "total_time". Default "depth".
#' @param lambda_range
#' @param ou_mode Logical. If TRUE, use Ornstein-Uhlenbeck process instead of Brownian Motion. Default FALSE.
#' @param ou_alpha Numeric. Selection strength for OU process. Higher values = stronger mean reversion. Default 1.0. optional, numeric vector of length 2, c(lambda_max, lambda_min) for total_time scaling. Default c(1, 0.1).
#' @import ape
#' @export


#' Make a tree ultrametric using node heights method
#' All leaves will be equidistant from the root
make_ultrametric <- function(tree) {
  n_tips <- length(tree$tip.label)
  n_total <- n_tips + tree$Nnode
  
  # Initialize heights (leaves = 0)
  node_heights <- rep(0, n_total)
  
  # Compute heights bottom-up
  # height[node] = max over children of (height[child] + edge_length_to_child)
  # Iterate until stable (handles any edge ordering)
  for (iter in 1:n_total) {
    for (i in 1:nrow(tree$edge)) {
      parent <- tree$edge[i, 1]
      child <- tree$edge[i, 2]
      candidate <- node_heights[child] + tree$edge.length[i]
      if (candidate > node_heights[parent]) {
        node_heights[parent] <- candidate
      }
    }
  }
  
  # Recalculate edge lengths: parent_height - child_height
  for (i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    tree$edge.length[i] <- node_heights[parent] - node_heights[child]
  }
  
  return(tree)
}

SimulateCIFs <- function(ncells, phyla, cif_center=1, Sigma=0.5, p_a=0.8, p_edge=NULL, n_CIF, n_diff, step=1, p_d=0.1, mu=0.1, N_char=9, N_ms=100, unif_on=FALSE, SIF_res=NULL, max_walk=2, lambda=0.05, T_cell=NULL, variable_branch_lengths=FALSE, branch_length_dist="exponential", branch_length_params=list(rate=1), branch_length_seed=NULL, ultrametric = FALSE, scale_state_walk=FALSE, state_walk_method="poisson", walk_rate=NULL, scale_barcode_mutations=FALSE, barcode_method="scale_mu", lambda_scaling="depth", lambda_range=c(1, 0.1), ou_mode=FALSE, ou_alpha=1.0, evolve_params=c("s")){  
if (is.null(T_cell)){
    T_cell <- stree(ncells, type = "balanced")
  }
  
    # Branch length generation
  if (variable_branch_lengths) {
    # Only generate new branch lengths if T_cell was not provided by user
    # (if user provided T_cell, keep their branch lengths)
    if (is.null(T_cell$edge.length) || length(unique(T_cell$edge.length)) == 1) {
      n_edges <- length(T_cell$edge[,1])
      
      # Handle seed for reproducibility
      if (!is.null(branch_length_seed)) {
        old_seed <- .Random.seed
        set.seed(branch_length_seed)
      }
      
      # Sample branch lengths from specified distribution
      if (branch_length_dist == "exponential") {
        T_cell$edge.length <- rexp(n_edges, rate = branch_length_params$rate)
      } else if (branch_length_dist == "lognormal") {
        T_cell$edge.length <- rlnorm(n_edges, 
                                      meanlog = branch_length_params$meanlog,
                                      sdlog = branch_length_params$sdlog)
      } else if (branch_length_dist == "gamma") {
        T_cell$edge.length <- rgamma(n_edges,
                                      shape = branch_length_params$shape,
                                      rate = branch_length_params$rate)
      } else {
        stop("Unknown branch_length_dist. Use 'exponential', 'lognormal', or 'gamma'.")
      }
      
      # Restore random state if we used a separate seed
      if (!is.null(branch_length_seed)) {
        .Random.seed <<- old_seed
      }
    }
    
    # Make ultrametric if requested (applies whether branch lengths were generated or provided)
    if (ultrametric) {
      T_cell <- make_ultrametric(T_cell)
    }
  } else {
    # Current behavior: uniform branch lengths
    T_cell$edge.length <- rep(1, length(T_cell$edge[,1]))
  }
  #  }
  #} else {
  #  # Current behavior: uniform branch lengths
  #  T_cell$edge.length <- rep(1, length(T_cell$edge[,1]))
  #}
  N_nodes <- length(T_cell$edge[,1])+1
  cell_edges <- cbind(T_cell$edge,T_cell$edge.length)
  cell_edges <- cbind(c(1:length(cell_edges[,1])),cell_edges)
  cell_connections <- table(c(cell_edges[,2],cell_edges[,3]))
  cell_root <- as.numeric(names(cell_connections)[cell_connections==2])
  Node_cell <- cell_root

  if (is.null(SIF_res)){
    returnlist <- SIFGenerate(phyla,n_diff,step = step)
  }else{
    returnlist <- SIF_res
  }

  sif_mean <- returnlist$sif_mean
  sif_mean_raw <- returnlist$sif_mean_raw
  state_tree <- returnlist$tree
  state_edges <- cbind(state_tree$edge,state_tree$edge.length)
  state_edges <- cbind(c(1:length(state_edges[,1])),state_edges)
  state_connections <- table(c(state_edges[,2],state_edges[,3]))
  state_root <- as.numeric(names(state_connections)[state_connections==2])
  Node_state <- state_root

  sif_label <- sif_mean_raw[[1]][,2]

  S <- sif_mean[[1]]
  S <- S[S[,3]==0,]
  if (length(S)>4){
    S <- t(matrix(c(S[1,1:3],cell_root),byrow = TRUE))
  } else{
    S <- t(matrix(c(S[1:3],cell_root),byrow = TRUE))
  }

  State_table <- SimulateCellStates(cell_root, cell_edges, state_edges, sif_mean = sif_mean[[1]], S = S, p_a = p_a, p_edge, max_walk = max_walk, scale_state_walk = scale_state_walk, state_walk_method = state_walk_method, walk_rate = walk_rate)
# Lambda processing
  if (lambda_scaling == "depth") {
    # Current behavior: discrete depth-based lambda
    if(length(lambda)==1){
      lambda <- rep(lambda, log2(ncells))
    }else if(length(lambda)==2){
      temp <- approx(lambda, n=log2(ncells))
      lambda <- temp$y
    }else{
      stop("Input of lambda has more than 2 entries.")
    }
    total_time_to_node <- NULL  # Not needed for depth scaling
  } else if (lambda_scaling == "total_time") {
    # New behavior: continuous total-time-based lambda
    total_time_to_node <- node.depth.edgelength(T_cell)
    lambda <- lambda_range  # Pass the range, will be interpolated in SampleEdgeNew
  } else {
    stop("Unknown lambda_scaling. Use 'depth' or 'total_time'.")
  }

  root_barcode <- rep(0,N_char)
  neutral <- Samplelineage(Node_cell, 0, cif_center, edges = cell_edges, edges_state = state_edges, sif_mean = sif_mean, S = State_table, cif = 1, p_a = p_a, p_d = p_d, mu = mu, flag = 1, barcode = root_barcode, N_ms = N_ms, unif_on = unif_on, lambda = lambda, lambda_scaling = lambda_scaling, total_time_to_node = total_time_to_node, scale_barcode_mutations = scale_barcode_mutations, barcode_method = barcode_method, ou_mode = ou_mode, ou_alpha = ou_alpha) 
  muts <- neutral[,5:length(neutral[1,])]

  param_names <- c("kon", "koff", "s")
  N_DE_cifs = c(ifelse("kon" %in% evolve_params, n_diff, 0), ifelse("koff" %in% evolve_params, n_diff, 0), ifelse("s" %in% evolve_params, n_diff, 0))
  N_ND_cifs = c(n_CIF - N_DE_cifs[1], n_CIF - N_DE_cifs[2], n_CIF - N_DE_cifs[3])


  cifs <- lapply(c(1:3),function(parami){
    nd_cif <- lapply(c(1:N_ND_cifs[parami]),function(icif){
      rnorm(N_nodes-1,cif_center,Sigma)
    })
    nd_cif <- do.call(cbind,nd_cif)
    if(N_DE_cifs[parami]!=0){
      #if there is more than 1 de_cifs for the parameter we are looking at
      de_cif <- lapply(c(1:N_DE_cifs[parami]),function(cif_i){
        Samplelineage(Node_cell, 0, cif_center, edges = cell_edges, edges_state = state_edges, sif_mean = sif_mean, S = State_table, p_a = p_a, cif = cif_i, flag = 0, lambda = lambda, lambda_scaling = lambda_scaling, total_time_to_node = total_time_to_node, scale_barcode_mutations = scale_barcode_mutations, barcode_method = barcode_method, ou_mode = ou_mode, ou_alpha = ou_alpha)      })

      de_cif <- lapply(de_cif,function(X){X[,4]})
      de_cif <- do.call(cbind,de_cif)
      cifs <- cbind(nd_cif,de_cif)
      colnames(cifs)<-c(
        paste(param_names[parami],rep('nonDE',length(nd_cif[1,])),c(1:length(nd_cif[1,])),sep='_'),
        paste(param_names[parami],rep('DE',length(de_cif[1,])),c(1:length(de_cif[1,])),sep='_'))
    }else{
      cifs <- nd_cif
      colnames(cifs)<-paste(param_names[parami],rep('nonDE',length(nd_cif[1,])),c(1:length(nd_cif[1,])),sep='_')
    }
    return(cifs)
  })
  muts[muts == Inf] <- '-'
  colnames(State_table) <- c("parent","cluster","depth","cellID")
  cif_res <- list(cifs,State_table,state_tree,T_cell,sif_mean_raw,sif_label,muts)
  return(cif_res)

}
