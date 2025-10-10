##### Functions ELEFANT #####


threshold_proba_interaction <- function(network, V, tol=0.01, echo=FALSE){
  
  # maximize Youden's J statistic
  # J = sensitivity + specificity - 1 
  
  cutoffs <- seq(tol, 1, tol)
  
  youden_statistic <- c()
  
  for (cutoff in cutoffs){
    V_net <- V
    V_net[V_net>=cutoff] <- 1
    V_net[V_net<cutoff] <- 0
    table <- table(network,V_net)
    if ("1" %in% colnames(table)) {true_positives <- table["1","1"] } else {true_positives <- 0}
    false_negatives <- table["1","0"]
    true_negatives <- table["0","0"]
    if ("1" %in% colnames(table)) {false_positives <- table["0","1"] } else {false_positives <- 0}
    
    youden_statistic <- c(youden_statistic, true_positives/(true_positives+false_negatives)+true_negatives/(true_negatives+false_positives)-1)
  }
  
  if (echo) {
    names(youden_statistic) <- cutoffs
    print(youden_statistic)
  }
  return(cutoffs[which.max(youden_statistic)][1])
  
}


sample_network_threshold <- function(V, threshold, obligate_A, obligate_B , count=FALSE){
  
  V_net <- V
  V_net[V_net>=threshold] <- 1
  V_net[V_net<threshold] <- 0
  
  sum <- 0
  
  # if obligate, take maximum value 
  if (obligate_A==TRUE){
    list_index <- which(colSums(V_net)==0)
    for (index in list_index) {
      V_net[which.max(V[,index]), index] <- 1
      sum <- sum+1
    }
  }
  
  if (obligate_B==TRUE){
    list_index <- which(rowSums(V_net)==0)
    for (index in list_index) {
      V_net[index, which.max(V[index,])] <- 1
      sum <- sum+1
    }
  }
  
  if (count) print(sum/(ncol(V_net)+nrow(V_net)))
  return(V_net)
}


##### Data augmentation on extant tree #####

rate_node_age <- function(lambda, mu, rho, t){
  r <- lambda-mu
  exp_rt <- exp(r*t)
  return(2*lambda*(1-exp_rt/(1/rho + lambda/r*(exp_rt-1))))
}

lambda_augmented <- function(lambda, mu, rho, t){
  r <- lambda-mu
  exp_rt <- exp(r*t)
  return(lambda*(1-exp_rt/(1/rho + lambda/r*(exp_rt-1))))
}

mu_augmented <- function(lambda, mu, rho, t){
  r <- lambda-mu
  exp_rt <- exp(r*t)
  return(mu/(1-exp_rt/(1/rho + lambda/r*(exp_rt-1))))
}

lambda_mu_augmented <- function(lambda, mu, rho, nb_lineages, t){
  r <- lambda-mu
  exp_rt <- exp(r*t)
  Phi_t <- 1-exp_rt/(1/rho + lambda/r*(exp_rt-1))
  return(nb_lineages*(lambda*Phi_t + mu/Phi_t))
}



sim_nhpp <- function(f_rate, min_t, max_t, tol=1e-3){
  
  rate_min_t <- f_rate(min_t)
  rate_max_t <- f_rate(max_t)
  infinite <- c()
  if (!is.finite(rate_min_t)) {
    rate_min_t <- f_rate((tol*max_t+(1-tol)*min_t))
    infinite <- (tol*max_t+(1-tol)*min_t)
  }
  if (!is.finite(rate_max_t)) rate_max_t <- f_rate(((1-tol)*max_t+tol*min_t))
  rate_max <- max(c(rate_min_t, rate_max_t))
  
  # Simulate homogeneous poisson process with constant rate 
  # Generate the number of events N during time t from Poisson distribution
  N <- rpois(1, rate_max * (max_t-min_t))
  
  # Generate N arrival times from a uniform distribution on [0, t]
  times <- sort(runif(N, min=min_t, max=max_t), decreasing = TRUE)
  
  # Thinning of the times
  
  proba <- sapply(times, function(i) f_rate(i)/rate_max  )
  
  return(c(times[which(runif(N, min=0, max=1)<proba)],infinite))
}


data_augmentation_tree <- function(tree, lambda, mu, rho, seed){
  
  set.seed(seed)
  
  age <- max(node.depth.edgelength(tree))
  ntip <- Ntip(tree)
  
  f_rate_nodes <- pryr::partial(rate_node_age, lambda=lambda, mu=mu, rho=rho)
  lambda_t <- pryr::partial(lambda_augmented, lambda = lambda, mu = mu, rho = rho)
  mu_t <- pryr::partial(mu_augmented, lambda = lambda, mu = mu, rho = rho)
  
  nb_augmented <- 0
  augmented_tree <- tree
  
  # For each branch
  for (branch in 1:nrow(tree$edge)){
    
    node_1 <- tree$edge[branch,1]
    node_2 <- tree$edge[branch,2]
    age_1 <- age-node.depth.edgelength(tree)[node_1]
    age_2 <- age-node.depth.edgelength(tree)[node_2]
    
    # Augmented nodes
    node_times <- sim_nhpp(f_rate_nodes, min_t = age_2, max_t = age_1)
    
    # Simulate birth-death process with homogeneous rates 
    if (length(node_times)>0){
      for (i in 1:length(node_times)){
        
        t_stem <- node_times[i]
        if (node_2>ntip){
          where <- getMRCA(augmented_tree, tip = extract.clade(tree, node=node_2)$tip.label)
        } else {
          where <- which(augmented_tree$tip.label==tree$tip.label[node_2])
        }
        
        f_rate_events <- pryr::partial(lambda_mu_augmented, lambda=lambda, mu=mu, rho=rho, nb_lineages=1)
        time <- sim_nhpp(f_rate_events, min_t = 0, max_t = t_stem)[1]
        
        if (is.na(time)){ # 1) One augmented lineage reaches the present (unsampled)
          time <- 0
          stem_branch <- t_stem
          nb_augmented <- nb_augmented+1
          augmented_tree <- bind.tip(augmented_tree, tip.label = paste0("DA_t",nb_augmented),
                                     edge.length = stem_branch, 
                                     where = where, 
                                     position = t_stem - age_2)
        } else {
          if (runif(1)<(1/(1+lambda_t(time)/mu_t(time)))){ # 2) Extinction (only one augmented branch)
            stem_branch <- t_stem - time
            nb_augmented <- nb_augmented+1
            augmented_tree <- bind.tip(augmented_tree, tip.label = paste0("DA_t",nb_augmented),
                                       edge.length = stem_branch, 
                                       where = where, 
                                       position = t_stem - age_2)
            
          } else { # 3) Speciation: there is an augmented tree
            stem_branch <- t_stem - time
            
            nb_nodes <- 3
            edge <- rbind(c(1,2),c(1,3))
            lineages_alive <- as.character(c(2,3))
            edge.length <- c(time, time) 
            
            while ((length(lineages_alive)>0)&(time>0)){
              
              f_rate_events <- pryr::partial(lambda_mu_augmented, lambda=lambda, mu=mu, rho=rho, nb_lineages=length(lineages_alive))
              time <- sim_nhpp(f_rate_events, min_t = 0, max_t = time)[1]
              
              if (is.na(time)){ # branch(es) alive 
                time <- 0
              }else{
                node <- as.numeric(sample(size=1,lineages_alive)) # select one lineage
                if (runif(1)<(1/(1+lambda_t(time)/mu_t(time)))){ # extinction
                  lineages_alive <- lineages_alive[lineages_alive!=node]
                  edge.length[which(edge[,2]==node)] <- edge.length[which(edge[,2]==node)] - time
                } else { # speciation
                  edge <- rbind(edge, c(node, nb_nodes+1), c(node, nb_nodes+2))
                  lineages_alive <- lineages_alive[lineages_alive!=node]
                  lineages_alive <- c(lineages_alive, as.character(nb_nodes+1), as.character(nb_nodes+2))
                  nb_nodes <- nb_nodes + 2
                  edge.length[which(edge[,2]==node)] <- edge.length[which(edge[,2]==node)] - time
                  edge.length <- c(edge.length, time, time)
                }
              }
            }
            
            # rename nodes
            n <- length(edge[which(!edge[,2] %in% edge[,1]),2])
            edge <- edge+n
            edge[which(!edge[,2] %in% edge[,1]),2] <- 1:n
            edge[, 1] <- n + dplyr::dense_rank(edge[, 1])
            edge[!edge[,2] %in% 1:n, 2] <- n + 1 + dplyr::dense_rank(edge[!edge[,2] %in% 1:n, 2])
            
            sub_tree <- list(edge = edge, edge.length = edge.length, 
                             tip.label = paste0("DA_t",nb_augmented+1:n), Nnode = n - 1)
            class(sub_tree) <- "phylo"
            nb_augmented <- nb_augmented+n
            sub_tree$root.edge <- stem_branch
            
            augmented_tree <- bind.tree(augmented_tree, sub_tree, 
                                        where = where, 
                                        position = t_stem - age_2)
          }
        }
      }
    }
  }
  return(augmented_tree)
}


is_crown_age_extant <- function(tree, age, max_tips=Inf, min_tips=2){
  if (is.null(tree)) return(FALSE) 
  list_tips <- tree$tip.label[abs(node.depth.edgelength(tree)[1:Ntip(tree)]-age)<1e-6]
  if (length(list_tips)<min_tips) return(FALSE)
  if (length(list_tips)>max_tips) return(FALSE)
  return(abs(max(node.depth.edgelength(drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% list_tips])))-age)<1e-6)
}


Ntip_extant <- function(tree, tol=1e-6){
  if (is.null(tree)) return(FALSE) 
  tip_ages <- node.depth.edgelength(tree)[1:Ntip(tree)]
  return(length(tree$tip.label[abs(tip_ages-max(tip_ages))<1e-6]))
}


##### Reconstruct past traits ######


# Build a matrix with tip and internal covariances with nodes
vcvPhyloInternal <- function(tree){
  nbtip <- Ntip(tree)
  dis <- dist.nodes(tree)
  MRCA <- mrca(tree, full = TRUE)
  M <- dis[as.character(nbtip + 1), MRCA]
  dim(M) <- rep(sqrt(length(M)), 2)
  return(M)
}


ancestral_nodes <- function(tree, traits, theta, sigma, d, stochastic_mapping=TRUE)    {
  
  V <- vcvPhyloInternal(tree) # between each pair of tip/node it gives the distance between the root and the MRCA of the pair
  nV <- nrow(V)
  rownames(V) <- colnames(V) <- 1:nV
  indice <- which(tree$tip.label %in% rownames(traits)) # extant and sampled species 
  traits <- traits[tree$tip.label[indice],] # order traits
  AY <- V[-indice,indice,drop=FALSE] # covariances between extant sampled tips and other unsampled species/nodes
  vY <- V[indice,indice,drop=FALSE] # covariances between extant sampled tips
  inv_vY <- pseudoinverse(vY)
  
  # Ancestral state at the root
  one <- t(rbind(rep(1,(nV-length(indice)))))
  one_traits <- t(rbind(rep(1,length(indice))))
  
  # Expected value of the states at the nodes (cf Clavel & Morlon 2019)
  state_nodes <- (AY%*%inv_vY%*%(traits-one_traits%*%theta))+(one%*%theta)
  colnames(state_nodes) = colnames(traits)
  rownames(state_nodes) = rownames(AY)
  
  # Covariances between nodes/unsampled tips (cf equations 12 from Martins & Hansen 1997)
  
  if (stochastic_mapping){
    vA <- V[-indice,-indice,drop=FALSE]
    cov_nodes <- vA - AY%*%inv_vY%*%t(AY)
    state_nodes_stochastic_mapping <- c()
    for (i in 1:d){
      state_nodes_stochastic_mapping <- cbind(state_nodes_stochastic_mapping, MASS::mvrnorm(n=1, mu=state_nodes[,i], Sigma = sigma[i]*cov_nodes, tol = 1e-8))
    }
    return(state_nodes_stochastic_mapping)
  } else {
    return(state_nodes)
  }
}




branches_t <- function(tree,t){ # give list branches at time t
  
  list_node1 <- c()
  list_node2 <- c()
  
  node_ages <- node.depth.edgelength(tree)
  node_ages <- max(node_ages) - node_ages # 0 corresponds to the present
  
  for (i in 1:nrow(tree$edge)){
    node1 <- tree$edge[i,1]
    node2 <- tree$edge[i,2]
    age1 <- node_ages[node1]
    age2 <- node_ages[node2]
    if ((round(t-age2,6)>=0) & (0<=round(age1-t,6))){
      list_node1 <- c(list_node1, node1)
      list_node2 <- c(list_node2, node2)
    }
  }
  
  branches_t <- data.frame(node1=list_node1, node2=list_node2)
  
  # change names
  n_tip <- Ntip(tree)
  node_label <- is.null(tree$node.label)
  if (node_label){
    branches_t$node1 <- paste0("node_",branches_t$node1)
  } else {
    branches_t$node1 <- tree$node.label[branches_t$node1-n_tip]
  }
  for (i in 1:nrow(branches_t)){
    if (as.numeric(branches_t$node2[i])>n_tip) {
      if (node_label){
        branches_t$node2[i] <- paste0("node_", branches_t$node2[i])
      }else{
        branches_t$node2[i] <- tree$node.label[as.numeric(branches_t$node2[i])-n_tip]
      }
    }else{
      branches_t$node2[i] <- tree$tip.label[as.numeric(branches_t$node2[i])]
    }
  }
  return(branches_t)
}

theta_t <- function(tree, states, n, t){  # compute theta at time t
  
  theta_t <- c() 
  list_branches_t <- branches_t(tree,t)
  node_ages <- node.depth.edgelength(tree)
  node_ages <- max(node_ages) - node_ages # 0 corresponds to the present
  
  if (is.null(tree$node.label)){
    names(node_ages) <- c(tree$tip.label, paste0("node_", (Ntip(tree) + 1:Nnode(tree)) ))
  } else {
    names(node_ages) <- c(tree$tip.label, tree$node.label)
  }
  
  for(i in 1:nrow(list_branches_t)) {
    theta_1 <- states[list_branches_t[i,1],] 
    theta_2 <- states[list_branches_t[i,2],]
    age_1 <- node_ages[list_branches_t[i,1]]
    age_2 <- node_ages[list_branches_t[i,2]]
    ratio <- (age_1-t)/(age_1-age_2)
    theta_t <- rbind(theta_t, theta_1 + ratio*(theta_2 - theta_1))
  }
  
  rownames(theta_t) <- list_branches_t[,2]
  
  return(theta_t)
}


build_past_traits <- function(tree_B, tree_A, states_B, states_A, E, t, evolution_A=TRUE) { # reconstruct latent traits 
  theta_B_t <- theta_t(tree=tree_B, states=states_B, n=Ntip(tree_B), t=t)
  if (evolution_A==TRUE) {
    theta_A_t <- theta_t(tree=tree_A, states=states_A, n=Ntip(tree_A), t=t)
  } else {
    theta_A_t <- theta_t(tree=tree_A, states=states_A, n=Ntip(tree_A), t=0)
  }
  V_past <- as.matrix(theta_B_t)%*%E%*%t(as.matrix(theta_A_t)) #new
  return(V_past)
}


tipnames_past_networks <- function(past_network, tree_A, tree_B, list_extant_A, list_extant_B){
  
  n_A <- Ntip(tree_A)
  n_B <- Ntip(tree_B)
  
  for (i in 1:ncol(past_network)){
    name <- colnames(past_network)[i]
    if (!name %in% list_extant_A){
      if (length(grep("DA_A",name))!=1){
        if (length(grep("^fossil_A",name))!=1){ # extant branch
          if (length(grep("nodeA_",name))==1){
            node <- as.numeric(gsub("nodeA_", "", name))
            list_tips <- extract.clade(tree_A, node=n_A + node)$tip.label
            list_tips_aug <- list_tips[grep("DA_A|fossil_A", list_tips)]
            if (length(list_tips)!=length(list_tips_aug)){
              list_tips <- list_tips[!list_tips %in% list_tips_aug]
            }
            list_tips_ext <- list_tips[list_tips %in% list_extant_A]
            if (length(list_tips_ext)>0){
              list_tips <- list_tips_ext
            }
            colnames(past_network)[i] <- paste0(sort(list_tips), collapse="-")
          }
        } else { # fossil branch from simulations
          
          branch <- which(tree_A$edge[,2]==which(tree_A$tip.label==name))
          branches <- which(tree_A$edge[,1]==tree_A$edge[branch,1])
          node <- tree_A$edge[branches[branches!=branch],2]
          if (node>n_A) { list_tips <- extract.clade(tree_A, node=node)$tip.label } else { list_tips <- tree_A$tip.label[node] }
          list_tips_ext <- list_tips[list_tips %in% list_extant_A]
          if (length(list_tips_ext)>0){
            list_tips <- list_tips_ext
          }
          colnames(past_network)[i] <- paste0(sort(list_tips), collapse="-")
          
        }
      }
    }
  }
  
  
  for (i in 1:nrow(past_network)){
    name <- rownames(past_network)[i]
    if (!name %in% list_extant_B){
      if (length(grep("DA_B",name))!=1){
        if (length(grep("^fossil_B",name))!=1){ # extant branch
          if (length(grep("nodeB_",name))==1){
            node <- as.numeric(gsub("nodeB_", "", name))
            list_tips <- extract.clade(tree_B, node=n_B + node)$tip.label
            list_tips_aug <- list_tips[grep("DA_B|fossil_B", list_tips)]
            if (length(list_tips)!=length(list_tips_aug)){
              list_tips <- list_tips[!list_tips %in% list_tips_aug]
            }
            list_tips_ext <- list_tips[list_tips %in% list_extant_B]
            if (length(list_tips_ext)>0){
              list_tips <- list_tips_ext
            }
            rownames(past_network)[i] <- paste0(sort(list_tips), collapse="-")
          }
          
        } else { # fossil branch from simulations
          
          branch <- which(tree_B$edge[,2]==which(tree_B$tip.label==name))
          branches <- which(tree_B$edge[,1]==tree_B$edge[branch,1])
          node <- tree_B$edge[branches[branches!=branch],2]
          if (node>n_B) { list_tips <- extract.clade(tree_B, node=node)$tip.label } else { list_tips <- tree_B$tip.label[node] }
          list_tips_ext <- list_tips[list_tips %in% list_extant_B]
          if (length(list_tips_ext)>0){
            list_tips <- list_tips_ext
          }
          rownames(past_network)[i] <- paste0(sort(list_tips), collapse="-")
          
        }
      }
    }
  }
  return(past_network)
}


get_network_metrics <- function(network, median_degree=FALSE, 
                                phylo_module=FALSE, tree_A=NULL, tree_B=NULL,
                                past=NULL, list_extant_A=NULL, list_extant_B=NULL){
  
  # only keep interacting species 
  network <- network[which(rowSums(network)>0), which(colSums(network)>0), drop=FALSE]
  
  nb_A <- ncol(network)
  nb_B <- nrow(network)
  
  connectance <- length(which(network>0))/(nrow(network)*ncol(network))
  
  nestedness <- tryCatch({bipartite::nested(network, method="NODF2")}, error = function(e) {NA})
  
  if (!phylo_module){
    modularity <- tryCatch({bipartite::computeModules(network, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, forceLPA=FALSE)@likelihood}, error = function(e) {NA})
  } else {
    
    trees <- update_tree_past_network(network, past, tree_A, tree_B, list_extant_A, list_extant_B)
    tree_A <- trees$tree_A
    tree_B <- trees$tree_B
    
    modules <- tryCatch({bipartite::computeModules(network, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, forceLPA=FALSE)}, error = function(e) {NA})
    if (isS4(modules)){
      modularity <- modules@likelihood
      list_modules <- listModuleInformation(modules)[[2]]
      MPD_module_A <- c()
      MPD_module_B <- c()
      phylo_dist_A <- cophenetic(tree_A)
      phylo_dist_B <- cophenetic(tree_B)
      for (i in 1:length(list_modules)){
        module <- list_modules[[i]]
        MPD_module_A <- c(MPD_module_A, mean(as.dist(phylo_dist_A[module[[2]],module[[2]]])))
        MPD_module_B <- c(MPD_module_B, mean(as.dist(phylo_dist_B[module[[1]],module[[1]]])))
      }
      MPD_module_A <- mean(MPD_module_A, na.rm = TRUE)
      MPD_module_B <- mean(MPD_module_B, na.rm = TRUE)
      MPD_A <- mean(as.dist(phylo_dist_A))
      MPD_B <- mean(as.dist(phylo_dist_B))
      diff_MPD_A <- MPD_module_A - MPD_A
      diff_MPD_B <- MPD_module_B - MPD_B
    }
  }
  if (phylo_module){
    if (median_degree){
      return(c(nb_A=nb_A, nb_B=nb_B, connectance=connectance, nestedness=nestedness, modularity=modularity, median_degree_A=median(colSums(network)), median_degree_B=median(rowSums(network)), diff_MPD_A=diff_MPD_A, diff_MPD_B=diff_MPD_B))
    } else {
      return(c(nb_A=nb_A, nb_B=nb_B, connectance=connectance, nestedness=nestedness, modularity=modularity, diff_MPD_A=diff_MPD_A, diff_MPD_B=diff_MPD_B))
    }
  } else {
    if (median_degree){
      return(c(nb_A=nb_A, nb_B=nb_B, connectance=connectance, nestedness=nestedness, modularity=modularity, median_degree_A=median(colSums(network)), median_degree_B=median(rowSums(network))))
    } else {
      return(c(nb_A=nb_A, nb_B=nb_B, connectance=connectance, nestedness=nestedness, modularity=modularity))
    }
  }
}


get_network_metrics_null_model <- function(null_network_r00, null_network_quasiswap, median_degree=TRUE){
  
  # only keep interacting species 
  null_network_r00 <- null_network_r00[which(rowSums(null_network_r00)>0), which(colSums(null_network_r00)>0), drop=FALSE]
  null_network_quasiswap <- null_network_quasiswap[which(rowSums(null_network_quasiswap)>0), which(colSums(null_network_quasiswap)>0), drop=FALSE]
  
  nb_A <- ncol(null_network_r00)
  nb_B <- nrow(null_network_r00)
  
  connectance <- length(which(null_network_r00>0))/(nrow(null_network_r00)*ncol(null_network_r00))
  

  nestedness <- tryCatch({bipartite::nested(null_network_quasiswap, method="NODF2")}, error = function(e) {NA})
  
  modularity <- tryCatch({bipartite::computeModules(null_network_quasiswap, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, forceLPA=FALSE)@likelihood}, error = function(e) {NA})
  
  if (median_degree){
    return(c(nb_A=nb_A, nb_B=nb_B, connectance=connectance, nestedness=nestedness, modularity=modularity, median_degree_A=median(colSums(null_network_r00)), median_degree_B=median(rowSums(null_network_r00))))
  } else {
    return(c(nb_A=nb_A, nb_B=nb_B, connectance=connectance, nestedness=nestedness, modularity=modularity))
  }
}


get_phylo_signal_old <- function(network, past, tree_A, tree_B, list_extant_A, list_extant_B){
  
  # only keep interacting species 
  network <- network[which(rowSums(network)>0), which(colSums(network)>0), drop=FALSE] # test 09/02
  
  f = file() # silence the function
  sink(file = f)
  
  age_A <- max(node.depth.edgelength(tree_A))
  age_B <- max(node.depth.edgelength(tree_B))
  
  tree_past_A <- tree_A
  tree_past_A <- drop.tip(tree_past_A, tip=tree_past_A$tip.label[grep(paste0("fossil_A_",past,"_"),tree_past_A$tip.label)])
  sliced_trees <- treeSlice(tree_past_A, slice = age_A-past, trivial = TRUE)
  
  for (i in 1:length(sliced_trees)){
    list_tips <- sliced_trees[[i]]$tip.label
    saved_tip <- list_tips[1]
    if (length(list_tips)>1){
      tree_past_A <- drop.tip(tree_past_A, tip=list_tips[2:length(list_tips)])
    }
    list_tips_aug <- list_tips[grep("DA_", list_tips)]
    if (length(list_tips)!=length(list_tips_aug)){
      list_tips <- list_tips[!list_tips %in% list_tips_aug]
    }
    list_tips_ext <- list_tips[list_tips %in% list_extant_A]
    if (length(list_tips_ext)>0){
      list_tips <- list_tips_ext
    }
    tree_past_A$tip.label[which(tree_past_A$tip.label==saved_tip)] <- paste0(sort(list_tips), collapse="-")
  }
  
  
  tree_past_B <- tree_B
  tree_past_B <- drop.tip(tree_past_B, tip=tree_past_B$tip.label[grep(paste0("fossil_B_",past,"_"),tree_past_B$tip.label)])
  sliced_trees <- treeSlice(tree_past_B, slice = age_B-past, trivial = TRUE)
  
  for (i in 1:length(sliced_trees)){
    list_tips <- sliced_trees[[i]]$tip.label
    saved_tip <- list_tips[1]
    if (length(list_tips)>1){
      tree_past_B <- drop.tip(tree_past_B, tip=list_tips[2:length(list_tips)])
    }
    list_tips_aug <- list_tips[grep("DA_", list_tips)]
    if (length(list_tips)!=length(list_tips_aug)){
      list_tips <- list_tips[!list_tips %in% list_tips_aug]
    }
    list_tips_ext <- list_tips[list_tips %in% list_extant_B]
    if (length(list_tips_ext)>0){
      list_tips <- list_tips_ext
    }
    tree_past_B$tip.label[which(tree_past_B$tip.label==saved_tip)] <- paste0(sort(list_tips), collapse="-")
  }
  
  # prepare network
  tree_past_A <- drop.tip(tree_past_A, tip=tree_past_A$tip.label[!tree_past_A$tip.label %in% colnames(network)])
  tree_past_B <- drop.tip(tree_past_B, tip=tree_past_B$tip.label[!tree_past_B$tip.label %in% rownames(network)])
  
  if ((!is.null(tree_past_A))&(!is.null(tree_past_B))){
    # correct edge length
    correction <- node.depth.edgelength(tree_past_A)-(age_A-past)
    tree_past_A$edge.length[which(tree_past_A$edge[,2] %in% which(correction>0))] <- tree_past_A$edge.length[which(tree_past_A$edge[,2] %in% which(correction>0))] - correction[which(correction>0)]
    
    correction <- node.depth.edgelength(tree_past_B)-(age_B-past)
    tree_past_B$edge.length[which(tree_past_B$edge[,2] %in% which(correction>0))] <- tree_past_B$edge.length[which(tree_past_B$edge[,2] %in% which(correction>0))] - correction[which(correction>0)]
    
    network <- network[tree_past_B$tip.label, tree_past_A$tip.label,drop=FALSE]
    
    if (min(dim(network))>1){
      mantel_test <- phylosignal_network(network, tree_A = tree_past_A, tree_B = tree_past_B, method = "Jaccard_binary")
      
      mantel_A <- mantel_test[c(3)]
      pvalue_A <- mantel_test[c(4)]
      mantel_B <- mantel_test[c(6)]
      pvalue_B <- mantel_test[c(7)]
      
      res <- c(mantel_A, pvalue_A, mantel_B, pvalue_B)
    } else {
      res <- c(NA,NA,NA,NA)
    }} else {
      res <- c(NA,NA,NA,NA)
    }

  sink()
  close(f)
  return(res)
}


update_tree_past_network <- function(network, past, tree_A, tree_B, list_extant_A, list_extant_B){
  
  age_A <- max(node.depth.edgelength(tree_A))
  age_B <- max(node.depth.edgelength(tree_B))
  
  tree_past_A <- tree_A
  tree_past_A <- drop.tip(tree_past_A, tip=tree_past_A$tip.label[grep(paste0("fossil_A_",past,"_"),tree_past_A$tip.label)])
  sliced_trees <- treeSlice(tree_past_A, slice = age_A-past, trivial = TRUE)
  
  for (i in 1:length(sliced_trees)){
    list_tips <- sliced_trees[[i]]$tip.label
    saved_tip <- list_tips[1]
    if (length(list_tips)>1){
      tree_past_A <- drop.tip(tree_past_A, tip=list_tips[2:length(list_tips)])
    }
    list_tips_aug <- list_tips[grep("DA_", list_tips)]
    if (length(list_tips)!=length(list_tips_aug)){
      list_tips <- list_tips[!list_tips %in% list_tips_aug]
    }
    list_tips_ext <- list_tips[list_tips %in% list_extant_A]
    if (length(list_tips_ext)>0){
      list_tips <- list_tips_ext
    }
    tree_past_A$tip.label[which(tree_past_A$tip.label==saved_tip)] <- paste0(sort(list_tips), collapse="-")
  }
  
  
  tree_past_B <- tree_B
  tree_past_B <- drop.tip(tree_past_B, tip=tree_past_B$tip.label[grep(paste0("fossil_B_",past,"_"),tree_past_B$tip.label)])
  sliced_trees <- treeSlice(tree_past_B, slice = age_B-past, trivial = TRUE)
  
  for (i in 1:length(sliced_trees)){
    list_tips <- sliced_trees[[i]]$tip.label
    saved_tip <- list_tips[1]
    if (length(list_tips)>1){
      tree_past_B <- drop.tip(tree_past_B, tip=list_tips[2:length(list_tips)])
    }
    list_tips_aug <- list_tips[grep("DA_", list_tips)]
    if (length(list_tips)!=length(list_tips_aug)){
      list_tips <- list_tips[!list_tips %in% list_tips_aug]
    }
    list_tips_ext <- list_tips[list_tips %in% list_extant_B]
    if (length(list_tips_ext)>0){
      list_tips <- list_tips_ext
    }
    tree_past_B$tip.label[which(tree_past_B$tip.label==saved_tip)] <- paste0(sort(list_tips), collapse="-")
  }
  
  # prepare network
  tree_past_A <- drop.tip(tree_past_A, tip=tree_past_A$tip.label[!tree_past_A$tip.label %in% colnames(network)])
  tree_past_B <- drop.tip(tree_past_B, tip=tree_past_B$tip.label[!tree_past_B$tip.label %in% rownames(network)])
  
  if ((!is.null(tree_past_A))&(!is.null(tree_past_B))){
    # correct edge length
    correction <- node.depth.edgelength(tree_past_A)-(age_A-past)
    tree_past_A$edge.length[which(tree_past_A$edge[,2] %in% which(correction>0))] <- tree_past_A$edge.length[which(tree_past_A$edge[,2] %in% which(correction>0))] - correction[which(correction>0)]
    
    correction <- node.depth.edgelength(tree_past_B)-(age_B-past)
    tree_past_B$edge.length[which(tree_past_B$edge[,2] %in% which(correction>0))] <- tree_past_B$edge.length[which(tree_past_B$edge[,2] %in% which(correction>0))] - correction[which(correction>0)]
  }
  return(list(tree_A=tree_past_A, tree_B=tree_past_B))
}

get_phylo_signal <- function(network, past, tree_A, tree_B, list_extant_A, list_extant_B){
  
  # only keep interacting species 
  network <- network[which(rowSums(network)>0), which(colSums(network)>0), drop=FALSE]
  
  f = file() # silence the function
  sink(file = f)
  
  
  trees <- update_tree_past_network(network, past, tree_A, tree_B, list_extant_A, list_extant_B)
  tree_past_A <- trees$tree_A
  tree_past_B <- trees$tree_B
  
  if ((!is.null(tree_past_A))&(!is.null(tree_past_B))){

    network <- network[tree_past_B$tip.label, tree_past_A$tip.label,drop=FALSE]
    
    if (min(dim(network))>1){
      
      res <- c()
      
      # Test number of partners
      mantel_test <- phylosignal_network(network, tree_A = tree_past_A, tree_B = tree_past_B, method = "degree", nperm = 250)
      mantel_A <- mantel_test[c(3)]
      pvalue_A <- mantel_test[c(4)]
      mantel_B <- mantel_test[c(6)]
      pvalue_B <- mantel_test[c(7)]
      res <- c(res, c(mantel_A, pvalue_A, mantel_B, pvalue_B))
      
      # Test species interactions (constraining the number of partners)
      mantel_test <- phylosignal_network(network, tree_A = tree_past_A, tree_B = tree_past_B, method = "Jaccard_binary", nperm = 250, permutation = "nbpartners")
      mantel_A <- mantel_test[c(3)]
      pvalue_A <- mantel_test[c(4)]
      mantel_B <- mantel_test[c(6)]
      pvalue_B <- mantel_test[c(7)]
      res <- c(res, c(mantel_A, pvalue_A, mantel_B, pvalue_B))
      
    } else {
      res <- c(NA,NA,NA,NA,NA,NA,NA,NA)
    }} else {
      res <- c(NA,NA,NA,NA,NA,NA,NA,NA)
    }
  
  sink()
  close(f)
  return(res)
}

rates_coextinction <- function(network, interval=5, replicates=10, seed=3, max_extinctions=50, 
                               coext_A=TRUE, coext_B=TRUE){
  
  set.seed(seed)
  
  network <- network[rowSums(network)>0,,drop=FALSE]
  network <- network[,colSums(network)>0,drop=FALSE]
  
  nb_A <- ncol(network)
  nb_B <- nrow(network)
  
  list_percentage <- seq(0,max_extinctions,interval)
  
  # Extinction in clade A - Coextinction in clade B
  
  if (coext_B){
    coextinctions <- c()
    for (perc in list_percentage){
      nb_coextinctions <- 0
      for (rep in 1:replicates){
        sub_network <- network[,sample(1:nb_A, size = round((100-perc)*nb_A/100), replace = FALSE),drop=FALSE]
        sub_network <- sub_network[rowSums(sub_network)>0,,drop=FALSE]
        nb_coextinctions <- nb_coextinctions + nb_B - nrow(sub_network)
      }
      coextinctions <- rbind(coextinctions, c(perc, nb_coextinctions/replicates/nb_B*100))
    }
    
    rate_coext_B <- as.vector(lm(coextinctions[,2]~coextinctions[,1]+0)$coefficients)
    #plot(coextinctions[,1], coextinctions[,2])
    #abline(0,rate_coext_B)
  } else {rate_coext_B <- NA}
  
  
  
  # Extinction in clade B - Coextinction in clade A
  
  if (coext_A){
    coextinctions <- c()
    for (perc in list_percentage){
      nb_coextinctions <- 0
      for (rep in 1:replicates){
        sub_network <- network[sample(1:nb_B, size = round((100-perc)*nb_B/100), replace = FALSE),,drop=FALSE]
        sub_network <- sub_network[,colSums(sub_network)>0,drop=FALSE]
        nb_coextinctions <- nb_coextinctions + nb_A - ncol(sub_network)
      }
      coextinctions <- rbind(coextinctions, c(perc, nb_coextinctions/replicates/nb_A*100))
    }
    
    rate_coext_A <- as.vector(lm(coextinctions[,2]~coextinctions[,1]+0)$coefficients)
    #plot(coextinctions[,1], coextinctions[,2])
    #abline(0,rate_coext_A)
  } else {rate_coext_A <- NA}
  
  return(c(rate_coext_A=rate_coext_A, rate_coext_B=rate_coext_B))
}


estim_past_latent_traits <- function(tree_extant, tree_augmented, traits_extant, d, stochastic_mapping=TRUE){
  
  sigma <- vector("numeric", length = d) # rates BM
  theta <- vector("numeric", length = d) # ancestral latent trait values
  for (i in 1:d) {
    mvBM_trait <- mvBM(tree_extant, traits_extant[,i], echo = FALSE, diagnostic=FALSE, param=list(constraint="diagonal"), model = "BM1")
    sigma[i] <- mvBM_trait[[5]]
    theta[i] <- mvBM_trait[[4]][1]
  }
  
  # augmented tree
  ancestral_states = ancestral_nodes(tree=tree_augmented, traits=traits_extant, theta=theta, sigma=sigma, d=d, stochastic_mapping=stochastic_mapping)
  rownames(ancestral_states) <- c(tree_augmented$tip.label, tree_augmented$node.label)[as.numeric(rownames(ancestral_states))]
  states <- rbind(traits_extant, ancestral_states) # add latent traits of tips
  states <- states[c(tree_augmented$tip.label, tree_augmented$node.label),] # order
  
  return(states)
}



cross_validation <- function(tree_extant_A, tree_extant_B, network, 
                             threshold, d, perc_cv_A, perc_cv_B, nb_cv=100, 
                             seed=1, echo=FALSE, obligate_A, obligate_B){

  tree_extant_A <- makeNodeLabel(tree_extant_A, method = "number", prefix = "nodeA_")
  tree_extant_B <- makeNodeLabel(tree_extant_B, method = "number", prefix = "nodeB_")
  
  age_A <- max(node.depth.edgelength(tree_extant_A))
  age_B <- max(node.depth.edgelength(tree_extant_B))
  
  list_prop_recover <- c()
  list_false_positives <- c()
  list_fscores <- c()
  
  set.seed(seed)
  
  for (i in 1:nb_cv){
    
    if (echo) print(i)
    
    # Subsample network
    
    # Must keep the same root after subsampling
    root <- FALSE
    
    while (root==FALSE){
      
      sub_network <- network
      
      if (perc_cv_A>0) {sub_network <- sub_network[,sort(sample(1:ncol(sub_network), size = round(ncol(sub_network)*(1-perc_cv_A))))]  }
      if (perc_cv_B>0) {sub_network <- sub_network[sort(sample(1:nrow(sub_network), size = round(nrow(sub_network)*(1-perc_cv_B)))),]  }
      
      sub_network <- sub_network[which(rowSums(sub_network)>0),]
      sub_network <- sub_network[,which(colSums(sub_network)>0)]
      
      tree_subsample_A <- drop.tip(tree_extant_A, tip=tree_extant_A$tip.label[!tree_extant_A$tip.label %in% colnames(sub_network)])
      tree_subsample_B <- drop.tip(tree_extant_B, tip=tree_extant_B$tip.label[!tree_extant_B$tip.label %in% rownames(sub_network)])
      
      if (abs(max(node.depth.edgelength(tree_subsample_A))-age_A)<1e-6){
        if (abs(max(node.depth.edgelength(tree_subsample_B))-age_B)<1e-6){
          root <- TRUE
        }
      }
    }
    
    # Get latent traits
    s <- svd(sub_network) # singular value decomposition 
    L <- s$u
    D <- s$d 
    R <- s$v
    
    E <- matrix(0, nrow = d,  ncol = d)
    diag(E) <- D[1:d]
    l <- L[,1:d]
    r <- t(R[,1:d])
    V <- l%*%E%*%r
    
    rownames(V) <- rownames(l) <- rownames(sub_network)
    colnames(V) <- colnames(r) <- colnames(sub_network)
    
    traits_extant_B <- l
    traits_extant_A <- t(r)
    
    states_B <- estim_past_latent_traits(tree_extant=tree_subsample_B, tree_augmented=tree_extant_B, 
                                         traits_extant=traits_extant_B, d=d, stochastic_mapping=TRUE)
    
    states_A <- estim_past_latent_traits(tree_extant=tree_subsample_A, tree_augmented=tree_extant_A, 
                                         traits_extant=traits_extant_A, d=d, stochastic_mapping=TRUE)
    
    
    V_cv <- build_past_traits(tree_B=tree_extant_B, tree_A=tree_extant_A, states_B=states_B, states_A=states_A, E=E, t=0.000001)
    network_cv <- sample_network_threshold(V_cv, threshold, obligate_A = obligate_A, obligate_B = obligate_B, count=FALSE)
    
    network_stats <- network
    network_stats[which(rownames(network_stats) %in% rownames(sub_network)),which(colnames(network_stats) %in% colnames(sub_network))] <- 0
    network_stats_cv <- network_cv
    network_stats_cv[which(rownames(network_stats_cv) %in% rownames(sub_network)),which(colnames(network_stats_cv) %in% colnames(sub_network))] <- 0
    
    table_cv <- table(network_stats, network_stats_cv)
    table_cv["0","0"] <- table_cv["0","0"]-dim(sub_network)[1]*dim(sub_network)[2]
    
    
    if ("1" %in% colnames(table_cv)) {
      precision <- table_cv["1","1"]/sum(table_cv[,"1"])
      recall <- table_cv["1","1"]/sum(table_cv["1",])
      fscore <- 2*(precision*recall)/(precision+recall)
      list_prop_recover <- c(list_prop_recover, recall) # % recovered
      list_false_positives <- c(list_false_positives,  1-precision) # % false positives
      list_fscores <- c(list_fscores, fscore)
    } else {
      list_prop_recover <- c(list_prop_recover, 0)
      list_false_positives <- c(list_false_positives, 0)
      list_fscores <- c(list_fscores, 0)
    }
    
  }
  
  true_positives=mean(list_prop_recover, na.rm = TRUE)
  false_positives=mean(list_false_positives, na.rm = TRUE)
  fscore=mean(list_fscores, na.rm = TRUE)
  
  return(c(true_positives=true_positives, false_positives=false_positives, fscore=fscore))
  
}


##### Run ELEFANT  #####


run_recon <- function(rep, name, list_ages, tree_extant_A, tree_extant_B, 
                      list_extant_A = NULL, list_extant_B = NULL,
                      list_full_extant_A = NULL, list_full_extant_B = NULL,
                      traits_extant_A, traits_extant_B, 
                      d, E, thresh_proba, stochastic_mapping, 
                      obligate_A, obligate_B, evolution_A, 
                      global_metrics=TRUE, null_model=TRUE, median_degree=TRUE,
                      data_augmentation_A, data_augmentation_B, treesDA_A, treesDA_B, seed=3, path_results=NULL){
  
  print(rep)
  
  set.seed(seed*10000+rep)
  list_network_metrics <- c()
  
  if (data_augmentation_A) {
    tree_augmented_A <- treesDA_A[[rep]]
    } else {tree_augmented_A <- tree_extant_A}
  
  if (data_augmentation_B) {
    tree_augmented_B <- treesDA_B[[rep]]
    } else {tree_augmented_B <- tree_extant_B}
  
  if (is.null(list_extant_A)) list_extant_A <- tree_extant_A$tip.label
  if (is.null(list_extant_B)) list_extant_B <- tree_extant_B$tip.label
  if (is.null(list_full_extant_A)) list_full_extant_A <- tree_extant_A$tip.label
  if (is.null(list_full_extant_B)) list_full_extant_B <- tree_extant_B$tip.label
  
  # Label the nodes
  tree_augmented_A$tip.label <- gsub("DA_t", "DA_A",tree_augmented_A$tip.label)
  tree_augmented_B$tip.label <- gsub("DA_t", "DA_B",tree_augmented_B$tip.label)
  tree_augmented_A <- makeNodeLabel(tree_augmented_A, method = "number", prefix = "nodeA_")
  tree_augmented_B <- makeNodeLabel(tree_augmented_B, method = "number", prefix = "nodeB_")
  
  
  ### Reconstruct past traits on data augmented trees
  
  states_B <- estim_past_latent_traits(tree_extant=tree_extant_B, tree_augmented=tree_augmented_B, 
                                       traits_extant=traits_extant_B, d=d, stochastic_mapping=stochastic_mapping)
  
  states_A <- estim_past_latent_traits(tree_extant=tree_extant_A, tree_augmented=tree_augmented_A, 
                                       traits_extant=traits_extant_A, d=d, stochastic_mapping=stochastic_mapping)
  
  
  ##### Simulate reconstructed networks
  
  past=0.00001
  V_past_0 <- build_past_traits(tree_B=tree_augmented_B, tree_A=tree_augmented_A, states_B=states_B, states_A=states_A, E=E, t=past)
  past_network_0 <- sample_network_threshold(V_past_0, thresh_proba, obligate_A = obligate_A, obligate_B = obligate_B)
  past_network_0 <- tipnames_past_networks(past_network = past_network_0, tree_A = tree_augmented_A, tree_B = tree_augmented_B, 
                                           list_extant_A = list_extant_A, list_extant_B = list_extant_B)
  
  for (past in list_ages[list_ages!=0]){
    assign(paste0("V_past_",past), build_past_traits(tree_B=tree_augmented_B, tree_A=tree_augmented_A, states_B=states_B, states_A=states_A, E=E, t=past, evolution_A=evolution_A))
    assign(paste0("past_network_",past), sample_network_threshold(get(paste0("V_past_",past)), thresh_proba, obligate_A = obligate_A, obligate_B = obligate_B))
    assign(paste0("past_network_",past), tipnames_past_networks(past_network = get(paste0("past_network_",past)), tree_A = tree_augmented_A, tree_B = tree_augmented_B, 
                                                                list_extant_A = list_extant_A, list_extant_B = list_extant_B))
  }
  
  
  ##### Compare connectance, modularity, nestedness
  
  if (global_metrics){
    for (past in list_ages){
      if (past==0){
        if (null_model) {
          null_network_r00 <- simulate(vegan::nullmodel(past_network_0, "r00"), nsim=1, seed = 1)[,,1]
          null_network_quasiswap <- simulate(vegan::nullmodel(past_network_0, "quasiswap"), nsim=1, seed = 1)[,,1]
          list_network_metrics <- rbind(list_network_metrics, c(past, rep, get_network_metrics(past_network_0, median_degree=median_degree), get_phylo_signal(past_network_0, 0.00001, tree_A=tree_augmented_A, tree_B=tree_augmented_B, list_extant_A=list_extant_A, list_extant_B=list_extant_B),rates_coextinction(past_network_0),
                                                                get_network_metrics_null_model(null_network_r00, null_network_quasiswap, median_degree=median_degree), rates_coextinction(null_network_r00)))
        } else {
          list_network_metrics <- rbind(list_network_metrics, c(past, rep, get_network_metrics(past_network_0, median_degree=median_degree), get_phylo_signal(past_network_0, 0.00001, tree_A=tree_augmented_A, tree_B=tree_augmented_B, list_extant_A=list_extant_A, list_extant_B=list_extant_B),rates_coextinction(past_network_0)))
        }
      } else {
        if (null_model) {
          null_network_r00 <- simulate(vegan::nullmodel(get(paste0("past_network_",past)), "r00"), nsim=1, seed = 1)[,,1]
          null_network_quasiswap <- simulate(vegan::nullmodel(get(paste0("past_network_",past)), "quasiswap"), nsim=1, seed = 1)[,,1]
          list_network_metrics <- rbind(list_network_metrics, c(past, rep, get_network_metrics(get(paste0("past_network_",past)), median_degree=median_degree), get_phylo_signal(get(paste0("past_network_",past)), past, tree_A=tree_augmented_A, tree_B=tree_augmented_B, list_extant_A=list_extant_A, list_extant_B=list_extant_B),rates_coextinction(get(paste0("past_network_",past))),
                                                                get_network_metrics_null_model(null_network_r00, null_network_quasiswap, median_degree=median_degree), rates_coextinction(null_network_r00)))
        } else {
          list_network_metrics <- rbind(list_network_metrics, c(past, rep, get_network_metrics(get(paste0("past_network_",past)), median_degree=median_degree), get_phylo_signal(get(paste0("past_network_",past)), past, tree_A=tree_augmented_A, tree_B=tree_augmented_B, list_extant_A=list_extant_A, list_extant_B=list_extant_B),rates_coextinction(get(paste0("past_network_",past)))))
        }
      }
    }
  }
  
  
  ##### Look at inferred interactions 
  
  for (past in list_ages){
    interactions_past <- melt(get(paste0("past_network_",past))[sort(rownames(get(paste0("past_network_",past)))),sort(colnames(get(paste0("past_network_",past))))])
    interactions_past <- interactions_past[which(interactions_past$value>0),]
    interactions_past$value <- rep
    write.table(interactions_past, file = paste0(path_results, "temporary_list_interactions_past_",name,"_",past,"_", rep, ".csv"), append = FALSE, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  return(list_network_metrics)
}


run_ELEFANT <- function(name, nb_recon, list_ages, tree_extant_A, tree_extant_B,
                        list_extant_A = NULL, list_extant_B = NULL,
                        list_full_extant_A = NULL, list_full_extant_B = NULL,
                        traits_extant_A, traits_extant_B, 
                        d, E, thresh_proba, stochastic_mapping, 
                        obligate_A, obligate_B, evolution_A, 
                        data_augmentation_A, data_augmentation_B, treesDA_A, treesDA_B, 
                        save_DA=TRUE,
                        global_metrics=TRUE, null_model=TRUE, median_degree=TRUE,
                        seed=3, nb_cores=1, path_results=NULL){
  
  
  if (is.null(list_extant_A)) list_extant_A <- tree_extant_A$tip.label
  if (is.null(list_extant_B)) list_extant_B <- tree_extant_B$tip.label
  if (is.null(list_full_extant_A)) list_full_extant_A <- tree_extant_A$tip.label
  if (is.null(list_full_extant_B)) list_full_extant_B <- tree_extant_B$tip.label
  
  
  list_network_metrics <- do.call(rbind, mclapply(1:nb_recon, run_recon, mc.preschedule=F, mc.cores = nb_cores, name=name, list_ages=list_ages, tree_extant_A=tree_extant_A, tree_extant_B=tree_extant_B, 
                                                  list_extant_A = list_extant_A, list_extant_B = list_extant_B,
                                                  list_full_extant_A = list_full_extant_A, list_full_extant_B = list_full_extant_B,
                                                  traits_extant_A=traits_extant_A, traits_extant_B=traits_extant_B, 
                                                  d=d, E=E, thresh_proba=thresh_proba, stochastic_mapping=stochastic_mapping, 
                                                  obligate_A=obligate_A, obligate_B=obligate_B, evolution_A=evolution_A,
                                                  data_augmentation_A=data_augmentation_A, data_augmentation_B=data_augmentation_B, treesDA_A=treesDA_A, treesDA_B=treesDA_B, 
                                                  global_metrics=global_metrics, null_model=null_model, median_degree=median_degree,
                                                  seed=seed, path_results=path_results))
  
  
  if (global_metrics){
    list_network_metrics <- data.frame(list_network_metrics)
    if (null_model) {
      colnames(list_network_metrics) <- c("age", "rep", "nb_A", "nb_B", "connectance", "nestedness", "modularity", "median_degree_A", "median_degree_B", "degree_mantel_cor_A", "degree_pvalue_A", "degree_mantel_cor_B", "degree_pvalue_B", "mantel_cor_A", "pvalue_A", "mantel_cor_B", "pvalue_B", "rate_coext_A", "rate_coext_B",
                                          "nb_A_NM", "nb_B_NM", "connectance_NM", "nestedness_NM", "modularity_NM", "median_degree_A_NM", "median_degree_B_NM", "rate_coext_A_NM", "rate_coext_B_NM")
    } else {
      colnames(list_network_metrics) <- c("age", "rep", "nb_A", "nb_B", "connectance", "nestedness", "modularity", "median_degree_A", "median_degree_B", "degree_mantel_cor_A", "degree_pvalue_A", "degree_mantel_cor_B", "degree_pvalue_B", "mantel_cor_A", "pvalue_A", "mantel_cor_B", "pvalue_B", "rate_coext_A", "rate_coext_B")
    }
  }


  
  # new names for saving
  new_tips_B <- paste0("b", 1:length(list_full_extant_B))
  names(new_tips_B) <- sort(list_full_extant_B)
  new_tips_A <- paste0("a", 1:length(list_full_extant_A))
  names(new_tips_A) <- sort(list_full_extant_A)
  write.table(cbind(c("species_name", names(new_tips_A), names(new_tips_B)), c("storage_name", new_tips_A, new_tips_B)), file = paste0(path_results, "list_interactions_past_storage_name_",name,".csv"), append = FALSE, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Make final file with all interactions
  for (past in list_ages){
    write.table(t(c("species_B", "species_A", "rep")), file = paste0(path_results, "list_interactions_past_",name,"_",past,".csv"), append = FALSE, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
    for (rep in 1:nb_recon){
      interactions_past <- read.table(file = paste0(path_results, "temporary_list_interactions_past_",name,"_",past,"_", rep, ".csv"), sep=";", header=FALSE)
      write.table(interactions_past, file = paste0(path_results, "list_interactions_past_",name,"_",past,".csv"), append = TRUE, sep=";", quote = FALSE, row.names = FALSE, col.names = FALSE)
      file.remove(paste0(path_results, "temporary_list_interactions_past_",name,"_",past,"_", rep, ".csv"))
    }
    
    # Reduce storage space
    list_inter <- read.table(file = paste0(path_results, "list_interactions_past_",name,"_",past,".csv"), sep=";", header=TRUE)
    for (species in sort(names(new_tips_B), decreasing = TRUE)) list_inter$species_B <- gsub(species, new_tips_B[species],  list_inter$species_B, fixed=TRUE)
    for (species in sort(names(new_tips_A), decreasing = TRUE)) list_inter$species_A <- gsub(species, new_tips_A[species],  list_inter$species_A, fixed=TRUE)
    
    # Don't save DA lineages if not needed (for consensus network reconstruction)
    if (!save_DA){
      if (length(grep("DA_B", list_inter$species_B))>0) list_inter <- list_inter[-grep("DA_B", list_inter$species_B),]
      if (length(grep("DA_A", list_inter$species_A))>0) list_inter <- list_inter[-grep("DA_A", list_inter$species_A),]
    }
    
    list_inter_merged <- list_inter %>%
      group_by(species_B, species_A) %>%
      summarize(rep = paste(unique(rep), collapse = "-"), .groups = 'drop')
    write.table(list_inter_merged, file = paste0(path_results, "list_interactions_past_",name,"_",past,".csv"), append = FALSE, sep=";", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  return(list_network_metrics)
}



ELEFANT <- function(name, nb_recon=250, list_ages, 
                    tree_A, tree_B,
                    list_full_extant_A = NULL, list_full_extant_B = NULL,
                    only_CV=FALSE,
                    stochastic_mapping=TRUE, 
                    obligate_A=TRUE, obligate_B=TRUE, evolution_A=TRUE,
                    data_augmentation_A=TRUE, data_augmentation_B=TRUE, 
                    treesDA_A=NULL, treesDA_B=NULL,
                    lambda_A=NULL, mu_A=NULL, rho_A=NULL,
                    lambda_B=NULL, mu_B=NULL, rho_B=NULL,
                    save_DA=TRUE,
                    global_metrics=TRUE,
                    null_model=TRUE, 
                    perc_cv_A=0.1,
                    perc_cv_B=0.1,
                    seed=3, nb_cores=1, path_results=NULL){
  
  #### Step 0: Prepare the data
  
  ####  Check data:
  
  if (!inherits(tree_A, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  if (!inherits(tree_B, "phylo")) {stop("object \"tree_B\" is not of class \"phylo\".")}
  
  if (!is.ultrametric(tree_A)){stop("object \"tree_A\" is not ultrametric.")}
  if (!is.ultrametric(tree_B)){stop("object \"tree_B\" is not ultrametric.")}
  
  if (!is.numeric(nb_recon) || nb_recon %% 1 != 0 || nb_recon <= 0){stop("object \"nb_recon\" must be greater than 0.")}
  
  if (!(is.numeric(list_ages) && all(list_ages >= 0))){stop("object \"list_ages\" does not contained only positive ages.")}
  
  min_MCRA <- min(c(max(node.depth.edgelength(tree_A)), max(node.depth.edgelength(tree_B))))
  if (max(list_ages)>min_MCRA){stop("object \"list_ages\" contains ages older than the youngest MRCA of clades A or B.")}
  
  
  if (data_augmentation_A) {
    if (!inherits(treesDA_A, "multiPhylo")) {
      if ((lambda_A<0)|(mu_A<0)|(rho_A<0)|(rho_A>1)) {stop("Please provide valid data augmentation parameters for \"tree_A\". If you do not want to use data augmentation, set lambda_A=1, mu_A=0, and rho_A=1.")}
    }
  }
  if (data_augmentation_B) {
    if (!inherits(treesDA_B, "multiPhylo")) {
      if ((lambda_B<0)|(mu_B<0)|(rho_B<0)|(rho_B>1)) {stop("Please provide valid data augmentation parameters for \"tree_B\". If you do not want to use data augmentation, set lambda_B=1, mu_B=0, and rho_B=1.")}
    }
  }
  
  if (is.null(path_results))  path_results <- paste0(getwd(),"/ELEFANT_",name,"/")
  dir.create(file.path(path_results), showWarnings = FALSE)
  dir.create(file.path(path_results, "figures/"), showWarnings = FALSE)
  
  
  
  # Make sure network is binary 
  network[network>0] <- 1
  
  # Only keep the species interacting in the network at present in the phylogenetic trees
  tree_A <- drop.tip(tree_A, tip=tree_A$tip.label[!tree_A$tip.label %in% colnames(network)])
  tree_B <- drop.tip(tree_B, tip=tree_B$tip.label[!tree_B$tip.label %in% rownames(network)])
  
  network <- as.matrix(network)[tree_B$tip.label, tree_A$tip.label]
  
  network <- network[,which(colSums(network)>0)]
  network <- network[which(rowSums(network)>0),]
  
  tree_A <- drop.tip(tree_A, tip=tree_A$tip.label[!tree_A$tip.label %in% colnames(network)])
  tree_B <- drop.tip(tree_B, tip=tree_B$tip.label[!tree_B$tip.label %in% rownames(network)])
  
  
  print("Step 1: Get the latent traits at present")
  
  s <- svd(network) # singular value decomposition 
  L <- s$u
  D <- s$d 
  R <- s$v
  
  # Percentage variance explained by each singular value
  explained_variance <- D^2/sum(D^2)*100
  
  d <- min(which(cumsum(explained_variance)>=90)) # at least 90% of the present-day variance
  d <- max(c(d,3)) # at least 3 traits
  
  E <- matrix(0, nrow = d,  ncol = d)
  diag(E) <- D[1:d]
  l <- L[,1:d]
  r <- t(R[,1:d])
  V <- l%*%E%*%r
  
  rownames(V) <- rownames(l) <- rownames(network)
  colnames(V) <- colnames(r) <- colnames(network)
  
  traits_extant_B <- l[tree_B$tip.label,]
  traits_extant_A <- t(r)[tree_A$tip.label,]
  
  # Thresholding approach maximizing Youden's J statistic
  thresh_proba <- threshold_proba_interaction(network, V, tol=0.001)
  
  results_summary <- c()
  nb_A <- ncol(network)
  nb_B <- nrow(network)
  
  tree_extant_A <- tree_A # extent trees of the network 
  tree_extant_B <- tree_B
  
  list_extant_A <- tree_extant_A$tip.label   # List extant species in the network
  list_extant_B <- tree_extant_B$tip.label
  
  
  write.tree(tree_extant_A, paste0(path_results,"/tree_ELEFANT_A_", name, "_.tre"))
  write.tree(tree_extant_B, paste0(path_results,"/tree_ELEFANT_B_", name, "_.tre"))
  write.table(network, paste0(path_results,"/network_ELEFANT_", name, "_.tre"), row.names=TRUE, sep=";")
  
  
  if (!only_CV){
    if (global_metrics){
      print("Compute global metrics at present")
      
      ## Measure global metrics and phylogenetic signal
      
      original_network_metrics <- get_network_metrics(network)
      mantel_test <- phylosignal_network(network, tree_A = tree_A, tree_B=tree_B, method = "Jaccard_binary")
      
      results_summary <- rbind(results_summary, c("d",  0, d, NA))
      results_summary <- rbind(results_summary, c("thresh_proba",  0, thresh_proba, NA))
      results_summary <- rbind(results_summary, c("nb_A",  "0_observed", original_network_metrics[1], NA))
      results_summary <- rbind(results_summary, c("nb_B",  "0_observed", original_network_metrics[2], NA))
      results_summary <- rbind(results_summary, c("connectance",  "0_observed", original_network_metrics[3], NA))
      results_summary <- rbind(results_summary, c("nestedness",  "0_observed", original_network_metrics[4], NA))
      results_summary <- rbind(results_summary, c("modularity",  "0_observed", original_network_metrics[5], NA))
      results_summary <- rbind(results_summary, c("phylosig_cladeA",  "0_observed", mantel_test[c(3,4)]))
      results_summary <- rbind(results_summary, c("phylosig_cladeB",  "0_observed", mantel_test[c(6,7)]))
    }
  }
  
  if (!only_CV){
    
    print("Step 2: Perform data augmentation")
    
    if (data_augmentation_A) {
      if (is.null(treesDA_A)){
        for (rep in 1:nb_recon){
          tree_augmented <- data_augmentation_tree(tree_extant_A, lambda=lambda_A, mu=mu_A, rho=rho_A, seed=seed*10000+rep)
          tree_augmented$tip.label <- gsub("DA_t", "DA_A",tree_augmented$tip.label) 
          tree_augmented <- makeNodeLabel(tree_augmented, method = "number", prefix = "nodeA_")
          treesDA_A[[rep]] <- tree_augmented
        }
      }
    }
    if (data_augmentation_B) {
      if (is.null(treesDA_B)){
        for (rep in 1:nb_recon){
          tree_augmented <- data_augmentation_tree(tree_extant_B, lambda=lambda_B, mu=mu_B, rho=rho_B, seed=seed*10000+rep)
          tree_augmented$tip.label <- gsub("DA_t", "DA_B",tree_augmented$tip.label) 
          tree_augmented <- makeNodeLabel(tree_augmented, method = "number", prefix = "nodeB_")
          treesDA_B[[rep]] <- tree_augmented
        }
      }
    }
    
    
    print("Steps 3 & 4: Reconstruct ancestral traits and networks")
    
    
    list_network_metrics <- run_ELEFANT(name=name, nb_recon=nb_recon, list_ages=list_ages, tree_extant_A=tree_extant_A, tree_extant_B=tree_extant_B,
                                        list_extant_A=list_extant_A, list_extant_B=list_extant_B,
                                        list_full_extant_A=list_full_extant_A, list_full_extant_B=list_full_extant_B,
                                        traits_extant_A=traits_extant_A, traits_extant_B=traits_extant_B, 
                                        d=d, E=E, thresh_proba=thresh_proba, stochastic_mapping=stochastic_mapping, 
                                        obligate_A=obligate_A, obligate_B=obligate_B, evolution_A=evolution_A, 
                                        data_augmentation_A=data_augmentation_A, data_augmentation_B=data_augmentation_B, treesDA_A=treesDA_A, treesDA_B=treesDA_B, 
                                        save_DA=save_DA,
                                        global_metrics, null_model,
                                        seed=seed, nb_cores=nb_cores, path_results=path_results)
    
    
    # Process outputs
    
    if (global_metrics){
      
      for (age in list_ages){
        
        metrics_age <- list_network_metrics[which(list_network_metrics$age==age),]
        
        # Global metrics 
        results_summary <- rbind(results_summary, c("nb_A", age, mean(metrics_age$nb_A, na.rm=TRUE), sd(metrics_age$nb_A, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nb_B", age, mean(metrics_age$nb_B, na.rm=TRUE), sd(metrics_age$nb_B, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("median_degree_A", age, mean(metrics_age$median_degree_A, na.rm=TRUE), sd(metrics_age$median_degree_A, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("median_degree_B", age, mean(metrics_age$median_degree_B, na.rm=TRUE), sd(metrics_age$median_degree_B, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("connectance", age, mean(metrics_age$connectance, na.rm=TRUE), sd(metrics_age$connectance, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nestedness", age, mean(metrics_age$nestedness, na.rm=TRUE), sd(metrics_age$nestedness, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("modularity", age, mean(metrics_age$modularity, na.rm=TRUE), sd(metrics_age$modularity, na.rm=TRUE)))
        
        results_summary <- rbind(results_summary, c("phylosig_mantel_cor_A", age, mean(metrics_age$mantel_cor_A, na.rm=TRUE), sd(metrics_age$mantel_cor_A, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("phylosig_pvalue_A", age, mean(metrics_age$pvalue_A, na.rm=TRUE), sd(metrics_age$pvalue_A, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("phylosig_mantel_cor_B", age, mean(metrics_age$mantel_cor_B, na.rm=TRUE), sd(metrics_age$mantel_cor_B, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("phylosig_pvalue_B", age, mean(metrics_age$pvalue_B, na.rm=TRUE), sd(metrics_age$pvalue_B, na.rm=TRUE)))
        
        results_summary <- rbind(results_summary, c("phylosig_degree_mantel_cor_A", age, mean(metrics_age$degree_mantel_cor_A, na.rm=TRUE), sd(metrics_age$degree_mantel_cor_A, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("phylosig_degree_pvalue_A", age, mean(metrics_age$degree_pvalue_A, na.rm=TRUE), sd(metrics_age$degree_pvalue_A, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("phylosig_degree_mantel_cor_B", age, mean(metrics_age$degree_mantel_cor_B, na.rm=TRUE), sd(metrics_age$degree_mantel_cor_B, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("phylosig_degree_pvalue_B", age, mean(metrics_age$degree_pvalue_B, na.rm=TRUE), sd(metrics_age$degree_pvalue_B, na.rm=TRUE)))
        
        
        ### 95% confidence interval on global metrics
        results_summary <- rbind(results_summary, c("nb_A_CI", age, quantile(metrics_age$nb_A, 0.025, na.rm=TRUE), quantile(metrics_age$nb_A, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nb_B_CI", age, quantile(metrics_age$nb_B, 0.025, na.rm=TRUE), quantile(metrics_age$nb_B, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("median_degree_A_CI", age, quantile(metrics_age$median_degree_A, 0.025, na.rm=TRUE), quantile(metrics_age$median_degree_A, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("median_degree_B_CI", age, quantile(metrics_age$median_degree_B, 0.025, na.rm=TRUE), quantile(metrics_age$median_degree_B, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("connectance_CI", age, quantile(metrics_age$connectance, 0.025, na.rm=TRUE), quantile(metrics_age$connectance, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nestedness_CI", age, quantile(metrics_age$nestedness, 0.025, na.rm=TRUE), quantile(metrics_age$nestedness, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("modularity_CI", age, quantile(metrics_age$modularity, 0.025, na.rm=TRUE), quantile(metrics_age$modularity, 0.975, na.rm=TRUE)))
        
        results_summary <- rbind(results_summary, c("phylosig_mantel_cor_A_CI", age, quantile(metrics_age$mantel_cor_A, 0.025, na.rm=TRUE), quantile(metrics_age$mantel_cor_A, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("phylosig_mantel_cor_B_CI", age, quantile(metrics_age$mantel_cor_B, 0.025, na.rm=TRUE), quantile(metrics_age$mantel_cor_B, 0.975, na.rm=TRUE)))
        
        results_summary <- rbind(results_summary, c("phylosig_degree_mantel_cor_A_CI", age, quantile(metrics_age$degree_mantel_cor_A, 0.025, na.rm=TRUE), quantile(metrics_age$degree_mantel_cor_A, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("phylosig_degree_mantel_cor_B_CI", age, quantile(metrics_age$degree_mantel_cor_B, 0.025, na.rm=TRUE), quantile(metrics_age$degree_mantel_cor_B, 0.975, na.rm=TRUE)))
        
        
        #### Null model
        results_summary <- rbind(results_summary, c("nb_A_NM", age, mean(metrics_age$nb_A_NM, na.rm=TRUE), sd(metrics_age$nb_A_NM, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nb_B_NM", age, mean(metrics_age$nb_B_NM, na.rm=TRUE), sd(metrics_age$nb_B_NM, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("median_degree_A_NM", age, mean(metrics_age$median_degree_A_NM, na.rm=TRUE), sd(metrics_age$median_degree_A_NM, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("median_degree_B_NM", age, mean(metrics_age$median_degree_B_NM, na.rm=TRUE), sd(metrics_age$median_degree_B_NM, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("connectance_NM", age, mean(metrics_age$connectance_NM, na.rm=TRUE), sd(metrics_age$connectance_NM, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nestedness_NM", age, mean(metrics_age$nestedness_NM, na.rm=TRUE), sd(metrics_age$nestedness_NM, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("modularity_NM", age, mean(metrics_age$modularity_NM, na.rm=TRUE), sd(metrics_age$modularity_NM, na.rm=TRUE)))
        
        ### 95% confidence interval on null models
        results_summary <- rbind(results_summary, c("nb_A_CI_NM", age, quantile(metrics_age$nb_A_NM, 0.025, na.rm=TRUE), quantile(metrics_age$nb_A_NM, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nb_B_CI_NM", age, quantile(metrics_age$nb_B_NM, 0.025, na.rm=TRUE), quantile(metrics_age$nb_B_NM, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("median_degree_A_CI_NM", age, quantile(metrics_age$median_degree_A_NM, 0.025, na.rm=TRUE), quantile(metrics_age$median_degree_A_NM, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("median_degree_B_CI_NM", age, quantile(metrics_age$median_degree_B_NM, 0.025, na.rm=TRUE), quantile(metrics_age$median_degree_B_NM, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("connectance_CI_NM", age, quantile(metrics_age$connectance_NM, 0.025, na.rm=TRUE), quantile(metrics_age$connectance_NM, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nestedness_CI_NM", age, quantile(metrics_age$nestedness_NM, 0.025, na.rm=TRUE), quantile(metrics_age$nestedness_NM, 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("modularity_CI_NM", age, quantile(metrics_age$modularity_NM, 0.025, na.rm=TRUE), quantile(metrics_age$modularity_NM, 0.975, na.rm=TRUE)))
        
      }
    }
    
    results_summary <- data.frame(results_summary)
    colnames(results_summary) <- c("metric", "age", "mean", "sd")
    results_summary$mean <- round(as.numeric(results_summary$mean), 5)
    results_summary$sd <- round(as.numeric(results_summary$sd), 5)
    
    save(list = ls(), file=paste0(path_results, "results_", name, ".Rdata"))
  }
  
  
  print("Step 5: Perform the cross validation at present")
  
  res_cv <- cross_validation(tree_extant_A=tree_A, tree_extant_B=tree_B, network=network, threshold=thresh_proba,
                             d=d, perc_cv_A, perc_cv_B, nb_cv=100, seed=1, echo=FALSE, obligate_A=obligate_A, obligate_B=obligate_B)
  print(res_cv)
  
  results_summary <- rbind(results_summary, c("cv_true_positive", 0, NA, res_cv[1]))
  results_summary <- rbind(results_summary, c("cv_false_positive",  0, NA, res_cv[2]))
  results_summary <- rbind(results_summary, c("cv_f1score",  0, NA, res_cv[3]))
  
  
  
  return(results_summary)
  
}




plot_networks_ELEFANT <- function(name, results, path_results=NULL,
                                  list_thresh=NULL,
                                  list_full_extant_A = NULL, list_full_extant_B = NULL){
  
  
  #### Consider an interaction as true is present in >XX% (keep % reconstructions to plot confidence)
  
  thresh_proba <- as.numeric(results$mean[which(results$metric=="thresh_proba")])
  list_ages <- unique(as.numeric(results$age))
  list_ages <- list_ages[!is.na(list_ages)]
  
  if (is.null(path_results))  path_results <- paste0(getwd(),"/ELEFANT_",name,"/")
  
  tree_A <- read.tree(paste0(path_results,"/tree_ELEFANT_A_", name, "_.tre"))
  tree_B <- read.tree(paste0(path_results,"/tree_ELEFANT_B_", name, "_.tre"))
  
  print("List of abbreviations for species names in the network:")
  print(read.table(paste0(path_results, "list_interactions_past_storage_name_",name,".csv"), header = TRUE, sep=";"))
  
  if (is.null(list_full_extant_A)) list_full_extant_A <- tree_A$tip.label
  if (is.null(list_full_extant_B)) list_full_extant_B <- tree_B$tip.label
  
  for (past in list_ages){
    
    for (threshold_prop_inter in c(list_thresh, thresh_proba)){
      
      list_interactions_past <- read.table(file = paste0(path_results, "list_interactions_past_",name,"_",past,".csv"), sep=";", header=TRUE)
      
      # Convert back with one row per replicate 
      list_interactions_past <- list_interactions_past %>%
        separate_rows(rep, sep = "-")
      
      colnames(list_interactions_past) <- c("species_B", "species_A", "rep")
      list_interactions_past$species_B <- as.character(list_interactions_past$species_B)
      list_interactions_past$species_A <- as.character(list_interactions_past$species_A)
      list_interactions_past$rep <- as.character(list_interactions_past$rep)
      
      list_interactions_past_extant <- list_interactions_past
      if (length(grep("DA_A", list_interactions_past_extant$species_A))>0) {
        list_interactions_past_extant <- list_interactions_past_extant[-grep("DA_A", list_interactions_past_extant$species_A),]
        list_interactions_past$species_A[grep("DA_A", list_interactions_past$species_A)] <- paste0(list_interactions_past$species_A[grep("DA_A", list_interactions_past$species_A)], "_", list_interactions_past$rep[grep("DA_A", list_interactions_past$species_A)])
      }
      if (length(grep("DA_B", list_interactions_past_extant$species_B))>0) {
        list_interactions_past_extant <- list_interactions_past_extant[-grep("DA_B", list_interactions_past_extant$species_B),]
        list_interactions_past$species_B[grep("DA_B", list_interactions_past$species_B)] <- paste0(list_interactions_past$species_B[grep("DA_B", list_interactions_past$species_B)], "_", list_interactions_past$rep[grep("DA_B", list_interactions_past$species_B)])
      }
      
      list_interactions_past_extant <- list_interactions_past_extant %>%
        group_by(species_B, species_A) %>% summarise(count = n(), .groups = "drop") %>% as.data.frame()
      list_interactions_past_extant$species_B <- as.character(list_interactions_past_extant$species_B)
      list_interactions_past_extant$species_A <- as.character(list_interactions_past_extant$species_A)
      list_interactions_past_extant$count <- as.numeric(as.character(list_interactions_past_extant$count))
      
      list_interactions_past_extant$count <- list_interactions_past_extant$count/nb_recon
      
      # only keep interactions > threshold_prop_inter %
      list_interactions_past_extant <- list_interactions_past_extant[which(list_interactions_past_extant$count>=threshold_prop_inter),]
      
      # Plot augmented graph:
      
      if (nrow(list_interactions_past_extant)>0){
        
        new_tips_B <- paste0("b", 1:length(list_full_extant_B))
        names(new_tips_B) <- sort(list_full_extant_B)
        new_tips_A <- paste0("a", 1:length(list_full_extant_A))
        names(new_tips_A) <- sort(list_full_extant_A)
        
        g <- graph_from_edgelist(as.matrix(list_interactions_past_extant[,1:2]), directed = FALSE)
        edge.attributes(g)$weight <- as.numeric(list_interactions_past_extant[,3])
        
        V(g)$type <- rep("Clade A", length(V(g)))
        V(g)$type[which(names(V(g)) %in% c(tree_B$tip.label, new_tips_B))] <- "Clade B"
        for (tip in c(tree_B$tip.label, new_tips_B)) { 
          V(g)$type[grep(paste0(tip,"-"), names(V(g)))] <- "Clade B"
          V(g)$type[grep(paste0("-",tip), names(V(g)))] <- "Clade B"
        }
        V(g)$type[grep("DA_B", names(V(g)))] <- "Clade B"
        V(g)$type_col <- V(g)$type
        V(g)$augmented <- "FALSE"
        V(g)$augmented[grep("DA_A|augmented_B", names(V(g)))] <- "TRUE"
        
        # define color and shape mappings.
        shape <- c("square", "circle")
        size <- c(5, 2.5)
        col <- c( "#229954", "#6e2c00")
        names(shape) <- c("Clade A", "Clade B")
        names(col) <- c("Clade A", "Clade B")
        names(size) <- c("FALSE", "TRUE")
        
        color_simul=FALSE
        color_edge <- rep("grey",nrow(list_interactions_past_extant))
        if (color_simul==TRUE){
          list_interactions_past_extant$pair <- paste0(list_interactions_past_extant$species_B," - ", list_interactions_past_extant$species_A)
          color_edge[list_interactions_past_extant$pair %in% simulate_interactions$pair] <- "#212f3d"
        }
        
        pdf(paste0(path_results, "/figures/ancestral_network_",name, "_thresh_", round(threshold_prop_inter, 2), "_t",past,".pdf"), width=7, height = 7)
        set.seed(1)
        plot(g, vertex.color = col[V(g)$type_col], vertex.shape = shape[V(g)$type], vertex.label=NA, vertex.size=size[V(g)$augmented],
             edge.width=5*E(g)$weight, edge.color=color_edge, layout= layout_with_fr(g, weights = 5*E(g)$weight,niter = 5000 ) )
        set.seed(1)
        plot(g, vertex.color = col[V(g)$type_col], vertex.shape = shape[V(g)$type],vertex.label = V(g)$name,
             vertex.label.cex = 0.8,vertex.label.color = "black", vertex.label.dist = 0.6, vertex.label.degree = pi/4, vertex.size = size[V(g)$augmented],
             edge.width = 5 * E(g)$weight, edge.color = color_edge, layout = layout_with_fr(g, weights = 5 * E(g)$weight, niter = 5000))
        plot(1,1,col="white")
        legend("bottomleft", pch=19,names(col),col=col,cex=0.8)
        dev.off()
      } else {
        pdf(paste0(path_results, "/figures/ancestral_network_",name, "_thresh_", round(threshold_prop_inter, 2), "_t",past,".pdf"), width=7, height = 7)
        dev.off()
      }
    }
    
  }
}

plot_metrics_ELEFANT <- function(name, results,
                                 clade_A="clade_A",
                                 clade_B="clade_B",
                                 min_age=NULL){
  
  results$mean <- as.numeric(results$mean)
  results$sd <- as.numeric(results$sd)
  
  list_ages <- sort(unique(results$age))
  list_ages <- list_ages[list_ages!="0_observed"]
  list_ages <- sort(as.numeric(list_ages))
  if (is.null(min_age)) min_age <- max(list_ages)
  list_ages <- list_ages[list_ages<=min_age]
  
  ### Plot trends in global metrics 
  
  data_plot <- results[which(results$metric %in% c("nb_A", "nb_A_CI")),]
  data_plot <- rbind(data_plot, results[which(results$metric %in% c("nb_B", "nb_B_CI")),])
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric=="nb_A")] <- data_plot$mean[which(data_plot$metric==paste0("nb_A_CI"))]
  data_plot$CI_sup[which(data_plot$metric=="nb_A")] <- data_plot$sd[which(data_plot$metric==paste0("nb_A_CI"))]
  data_plot$CI_inf[which(data_plot$metric=="nb_B")] <- data_plot$mean[which(data_plot$metric==paste0("nb_B_CI"))]
  data_plot$CI_sup[which(data_plot$metric=="nb_B")] <- data_plot$sd[which(data_plot$metric==paste0("nb_B_CI"))]
  data_plot <- data_plot[which(data_plot$metric %in% c("nb_A", "nb_B")),]
  data_plot$metric <- gsub("nb_", "", data_plot$metric)
  data_plot$age <- as.numeric(data_plot$age)
  max <- max(data_plot$CI_sup)
  p1 <- (ggplot(data_plot, aes(x=age, y=mean, color=metric))+
           geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2, position=position_dodge(0.1))+geom_line()+
           geom_point()+theme_bw()+xlab("Myr")+ylab("Nb. of lineages")+
           scale_color_manual(values=c("#229954", "#6e2c00"))+
           ylim(0, max)+scale_x_reverse()) + guides(color = guide_legend(title = NULL))+
    theme(legend.position = c(0.15, 0.8), legend.background = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank())+ theme(legend.position="none")
  
  metric="median_degree_A" # Z-score = (obs - mean NM) / sd NM
  data_plot <- results[grep(metric,results$metric),]
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$age <- as.numeric(data_plot$age)
  data_plot$mean[data_plot$metric==metric] <- (data_plot$mean[data_plot$metric==metric] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]
  data_plot$mean[data_plot$metric==paste0(metric, "_CI")] <- (data_plot$mean[data_plot$metric==paste0(metric, "_CI")] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]  # add error bars 18/09/25
  data_plot$sd[data_plot$metric==paste0(metric, "_CI")] <- (data_plot$sd[data_plot$metric==paste0(metric, "_CI")] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric==metric)] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_sup[which(data_plot$metric==metric)] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_inf[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  data_plot$CI_sup[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  data_plot <- data_plot[which(data_plot$metric==metric),]
  max <- max(abs(c(data_plot$CI_sup, data_plot$CI_inf)))
  data_plot$signif=FALSE
  data_plot$signif[abs(data_plot$mean)>1.96] <- TRUE
  shapes <- c()
  if (FALSE %in% data_plot$signif) shapes <- c(shapes, 4)
  if (TRUE %in% data_plot$signif) shapes <- c(shapes, 15)
  p2 <- ggplot(data_plot, aes(x=age, y=mean, group=metric, color=metric, shape=signif)) + 
    annotate("rect", xmin = 0, xmax = max(list_ages), ymin = -min(c(1.96,max)), ymax = min(c(1.96,max)),alpha = .1,fill = "orange")+
    geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2, position=position_dodge(0.1))+ # add error bars
    geom_hline(yintercept = 0) + geom_line() + geom_point()+
    theme_bw()+xlab("Myr")+ylab(paste0("Median degree\nin ", clade_A, " (Z-score)"))+ ylim(-max, max)+scale_x_reverse()+
    scale_color_manual(values=c("black","orange"))+ theme(legend.position="none")+ scale_shape_manual(values=shapes)+ theme(legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank())
  
  metric="median_degree_B" # Z-score = (obs - mean NM) / sd NM
  data_plot <- results[grep(metric,results$metric),]
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$age <- as.numeric(data_plot$age)
  data_plot$mean[data_plot$metric==metric] <- (data_plot$mean[data_plot$metric==metric] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]
  data_plot$mean[data_plot$metric==paste0(metric, "_CI")] <- (data_plot$mean[data_plot$metric==paste0(metric, "_CI")] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]  # add error bars 18/09/25
  data_plot$sd[data_plot$metric==paste0(metric, "_CI")] <- (data_plot$sd[data_plot$metric==paste0(metric, "_CI")] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric==metric)] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_sup[which(data_plot$metric==metric)] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_inf[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  data_plot$CI_sup[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  data_plot <- data_plot[which(data_plot$metric==metric),]
  max <- max(abs(c(data_plot$CI_sup, data_plot$CI_inf)))
  data_plot$signif=FALSE
  data_plot$signif[abs(data_plot$mean)>1.96] <- TRUE
  shapes <- c()
  if (FALSE %in% data_plot$signif) shapes <- c(shapes, 4)
  if (TRUE %in% data_plot$signif) shapes <- c(shapes, 15)
  p3 <- ggplot(data_plot, aes(x=age, y=mean, group=metric, color=metric, shape=signif)) + 
    annotate("rect", xmin = 0, xmax = max(list_ages), ymin = -min(c(1.96,max)), ymax = min(c(1.96,max)),alpha = .1,fill = "orange")+
    geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2, position=position_dodge(0.1))+ # add error bars
    geom_hline(yintercept = 0) + geom_line() + geom_point()+
    theme_bw()+xlab("Myr")+ylab(paste0("Median degree\nin ", clade_B, " (Z-score)"))+ ylim(-max, max)+scale_x_reverse()+
    scale_color_manual(values=c("black","orange"))+ theme(legend.position="none")+ scale_shape_manual(values=shapes)+ theme(legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank())
  
  metric="connectance"
  data_plot <- results[grep(metric,results$metric),]
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$age <- as.numeric(data_plot$age)
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric==metric)] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_sup[which(data_plot$metric==metric)] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI"))]
  #data_plot$CI_inf[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  #data_plot$CI_sup[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  #data_plot <- data_plot[which(data_plot$metric %in% c(metric, paste0(metric,"_NM"))),]
  data_plot <- data_plot[which(data_plot$metric %in% c(metric)),]
  max <- max(data_plot$CI_sup)
  data_plot <- data_plot[order(data_plot$metric, decreasing = TRUE),]
  p4 <- ggplot(data_plot, aes(x=age, y=mean, group=metric, color=metric)) + 
    geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2, position=position_dodge(0.1))+geom_line() + geom_point()+ 
    theme_bw()+xlab("Myr")+ylab("Connectance")+ ylim(0, max)+scale_x_reverse()+
    scale_color_manual(values=c("black","orange"))+ theme(legend.position="none")+ theme(legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank())
  
  metric="nestedness" # Z-score = (obs - mean NM) / sd NM
  data_plot <- results[grep(metric,results$metric),]
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$age <- as.numeric(data_plot$age)
  data_plot$mean[data_plot$metric==metric] <- (data_plot$mean[data_plot$metric==metric] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]
  data_plot$mean[data_plot$metric==paste0(metric, "_CI")] <- (data_plot$mean[data_plot$metric==paste0(metric, "_CI")] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]  # add error bars 18/09/25
  data_plot$sd[data_plot$metric==paste0(metric, "_CI")] <- (data_plot$sd[data_plot$metric==paste0(metric, "_CI")] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric==metric)] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_sup[which(data_plot$metric==metric)] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_inf[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  data_plot$CI_sup[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  data_plot <- data_plot[which(data_plot$metric==metric),]
  max <- max(abs(c(data_plot$CI_sup, data_plot$CI_inf)))
  data_plot$signif=FALSE
  data_plot$signif[abs(data_plot$mean)>1.96] <- TRUE
  shapes <- c()
  if (FALSE %in% data_plot$signif) shapes <- c(shapes, 4)
  if (TRUE %in% data_plot$signif) shapes <- c(shapes, 15)
  p5 <- ggplot(data_plot, aes(x=age, y=mean, group=metric, color=metric, shape=signif)) + 
    annotate("rect", xmin = 0, xmax = max(list_ages), ymin = -min(c(1.96,max)), ymax = min(c(1.96,max)),alpha = .1,fill = "orange")+
    geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2, position=position_dodge(0.1))+ # add error bars
    geom_hline(yintercept = 0) + geom_line() + geom_point()+
    theme_bw()+xlab("Myr")+ylab("Nestedness (Z-score)")+ ylim(-max, max)+scale_x_reverse()+
    scale_color_manual(values=c("black","orange"))+ theme(legend.position="none")+ scale_shape_manual(values=shapes)+ theme(legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank())
  
  metric="modularity" # Z-score = (obs - mean NM) / sd NM
  data_plot <- results[grep(metric,results$metric),]
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$age <- as.numeric(data_plot$age)
  data_plot$mean[data_plot$metric==metric] <- (data_plot$mean[data_plot$metric==metric] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]
  data_plot$mean[data_plot$metric==paste0(metric, "_CI")] <- (data_plot$mean[data_plot$metric==paste0(metric, "_CI")] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]  # add error bars 18/09/25
  data_plot$sd[data_plot$metric==paste0(metric, "_CI")] <- (data_plot$sd[data_plot$metric==paste0(metric, "_CI")] - data_plot$mean[data_plot$metric==paste0(metric,"_NM")]) /  data_plot$sd[data_plot$metric==paste0(metric,"_NM")]
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric==metric)] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_sup[which(data_plot$metric==metric)] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_inf[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  data_plot$CI_sup[which(data_plot$metric==paste0(metric,"_NM"))] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI_NM"))]
  data_plot <- data_plot[which(data_plot$metric==metric),]
  max <- max(abs(c(data_plot$CI_sup, data_plot$CI_inf)))
  data_plot$signif=FALSE
  data_plot$signif[abs(data_plot$mean)>1.96] <- TRUE
  shapes <- c()
  if (FALSE %in% data_plot$signif) shapes <- c(shapes, 4)
  if (TRUE %in% data_plot$signif) shapes <- c(shapes, 15)
  p6 <- ggplot(data_plot, aes(x=age, y=mean, group=metric, color=metric, shape=signif)) + 
    annotate("rect", xmin = 0, xmax = max(list_ages), ymin = -min(c(1.96,max)), ymax = min(c(1.96,max)),alpha = .1,fill = "orange")+
    geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2, position=position_dodge(0.1))+ # add error bars
    geom_hline(yintercept = 0) + geom_line() + geom_point()+
    theme_bw()+xlab("Myr")+ylab("Modularity (Z-score)")+ ylim(-max, max)+scale_x_reverse()+
    scale_color_manual(values=c("black","orange"))+ theme(legend.position="none")+ scale_shape_manual(values=shapes)+ theme(legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank())
  
  metric="phylosig_mantel_cor_A"
  data_plot <- results[grep(metric,results$metric),]
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$age <- as.numeric(data_plot$age)
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric==metric)] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_sup[which(data_plot$metric==metric)] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot <- data_plot[which(data_plot$metric==metric),]
  data_plot$signif <- results$sd[intersect(which(results$metric==metric),which(results$age %in% list_ages))]<=0.05
  data_plot$clade <- "clade_A"
  data_plot$age <- data_plot$age + diff(range(data_plot$age))/100
  data_plot_1 <- data_plot
  metric="phylosig_mantel_cor_B"
  data_plot <- results[grep(metric,results$metric),]
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$age <- as.numeric(data_plot$age)
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric==metric)] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_sup[which(data_plot$metric==metric)] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot <- data_plot[which(data_plot$metric==metric),]
  data_plot$signif <- results$sd[intersect(which(results$metric==metric),which(results$age %in% list_ages))]<=0.05
  data_plot$clade <- "clade_B"
  data_plot <- rbind(data_plot_1, data_plot)
  max <- max(data_plot$CI_sup)
  min <- min(data_plot$CI_inf)
  data_plot <- data_plot[order(data_plot$metric, decreasing = TRUE),]
  shapes <- c()
  if (FALSE %in% data_plot$signif) shapes <- c(shapes, 4)
  if (TRUE %in% data_plot$signif) shapes <- c(shapes, 15)
  p7 <- ggplot(data_plot, aes(x=age, y=mean, group=clade, color=clade, shape=signif)) + 
    geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2, position=position_dodge(0.1))+geom_line() + geom_point()+ 
    theme_bw()+xlab("Myr")+ylab(paste0("Phylo. signal \nin species interactions"))+ ylim(min, max)+scale_x_reverse()+
    scale_color_manual(values=c("#229954", "#6e2c00"))+ scale_shape_manual(values=shapes)+ theme(legend.position="none")+ theme(legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank())
  
  
  metric="phylosig_degree_mantel_cor_A"
  data_plot <- results[grep(metric,results$metric),]
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$age <- as.numeric(data_plot$age)
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric==metric)] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_sup[which(data_plot$metric==metric)] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot <- data_plot[which(data_plot$metric==metric),]
  data_plot$signif <- results$sd[intersect(which(results$metric==metric),which(results$age %in% list_ages))]<=0.05
  data_plot$clade <- "clade_A"
  data_plot$age <- data_plot$age + diff(range(data_plot$age))/100
  data_plot_1 <- data_plot
  metric="phylosig_degree_mantel_cor_B"
  data_plot <- results[grep(metric,results$metric),]
  data_plot <- data_plot[which(data_plot$age %in% list_ages),]
  data_plot$age <- as.numeric(data_plot$age)
  data_plot$CI_inf <- NA
  data_plot$CI_sup <- NA
  data_plot$CI_inf[which(data_plot$metric==metric)] <- data_plot$mean[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot$CI_sup[which(data_plot$metric==metric)] <- data_plot$sd[which(data_plot$metric==paste0(metric, "_CI"))]
  data_plot <- data_plot[which(data_plot$metric==metric),]
  data_plot$signif <- results$sd[intersect(which(results$metric==metric),which(results$age %in% list_ages))]<=0.05
  data_plot$clade <- "clade_B"
  data_plot <- rbind(data_plot_1, data_plot)
  max <- max(data_plot$CI_sup)
  min <- min(data_plot$CI_inf)
  data_plot <- data_plot[order(data_plot$metric, decreasing = TRUE),]
  shapes <- c()
  if (FALSE %in% data_plot$signif) shapes <- c(shapes, 4)
  if (TRUE %in% data_plot$signif) shapes <- c(shapes, 15)
  p8 <- ggplot(data_plot, aes(x=age, y=mean, group=clade, color=clade, shape=signif)) + 
    geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2, position=position_dodge(0.1))+geom_line() + geom_point()+ 
    theme_bw()+xlab("Myr")+ylab(paste0("Phylo. signal \nin nb. of partners"))+ ylim(min, max)+scale_x_reverse()+
    scale_color_manual(values=c("#229954", "#6e2c00"))+ scale_shape_manual(values=shapes)+ theme(legend.position="none")+ theme(legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank())
  
  
  figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, NULL,
                      labels = c("", "", "", "", "", ""),
                      ncol = 3, nrow = 3)
  
  print(figure)
  
}






##### Simulations ######


fossil_tips <- function(tree, time, age, name="A", edge_length_fossil=0){
  
  subtrees <- treeSlice(tree, slice=age-time, trivial = TRUE)
  
  for (i in 1:length(subtrees)){
    
    node_ages <- node.depth.edgelength(tree)[1:Ntip(tree)]
    names(node_ages) <- tree$tip.label
    tips <- subtrees[[i]]$tip.label
    
    if (length(tips)>1) {
      node <- getMRCA(tree, tips)
      if (abs(max(node_ages[tips])-age)<1e-6){ # if clade still alive
        position <- time - max(node.depth.edgelength(subtrees[[i]]))
      } else { # if extinct clade
        position <- time - max(node.depth.edgelength(subtrees[[i]])) - (age - max(node_ages[tips]))
      }
    } else { # one tip
      node <- which(tree$tip.label==tips)
      position <- time - (age - node_ages[node])
    }
    tree <- bind.tip(tree, where = node, edge.length = edge_length_fossil, tip.label = paste0("fossil_",name,"_", time, "_", i), position = position)
  }
  return(tree)
}



sim_host_repertoire_evolution <- function(tree_B, species_A, rate_gain, rate_loss, initial_state){
  
  nb_A <- length(species_A)
  
  host_repertoire <- matrix(0, nrow=Ntip(tree_B)+Nnode(tree_B), ncol=nb_A)
  colnames(host_repertoire) <- species_A
  rownames(host_repertoire) <- c(tree_B$tip.label, paste0("nodeB_", 1:Nnode(tree_B))) #Ntip(tree_B)+
  
  host_repertoire[Ntip(tree_B)+1,] <- initial_state
  
  for (edge in 1:nrow(tree_B$edge)){
    
    state <- host_repertoire[tree_B$edge[edge,1],]
    branch_length <- tree_B$edge.length[edge]
    
    waiting_time <- 0
    
    while(waiting_time<branch_length){
      nb_inter <- sum(state)
      
      if (nb_inter<nb_A) {rate_gain_t <- rate_gain} else {rate_gain_t <- 0}
      if (nb_inter>1) {rate_loss_t <- rate_loss} else {rate_loss_t <- 0}
      
      waiting_time <- rexp(n=1, rate = rate_loss_t + rate_gain_t)
      
      if (waiting_time < branch_length){
        if (runif(1)<(rate_loss_t/(rate_loss_t+rate_gain_t))){ # loss
          state[sample(size=1, names(which(state==1)))] <- 0
        } else { # gain
          state[sample(size=1, names(which(state==0)))] <- 1
        }
        branch_length <- branch_length - waiting_time
        waiting_time <- 0
      }
    }
    host_repertoire[tree_B$edge[edge,2],] <- state
  }
  
  return(host_repertoire)
}





sim_network_evolution <- function(tree_B, tree_A, rate_gain, rate_loss, initial_state, seed=1){
  
  set.seed(seed)
  
  network_nodes <- matrix(0, nrow=Ntip(tree_B)+Nnode(tree_B), ncol=Ntip(tree_A)+Nnode(tree_A))
  colnames(network_nodes) <- c(tree_A$tip.label, paste0("nodeA_", 1:Nnode(tree_A)))
  rownames(network_nodes) <- c(tree_B$tip.label, paste0("nodeB_", 1:Nnode(tree_B)))
  
  node_age_A <- max(node.depth.edgelength(tree_A)) - node.depth.edgelength(tree_A)
  names(node_age_A) <- c(tree_A$tip.label, paste0("nodeA_", 1:Nnode(tree_A)))
  node_age_A[node_age_A<1e-6] <- 0
  node_age_B <- max(node.depth.edgelength(tree_B)) - node.depth.edgelength(tree_B)
  names(node_age_B) <- c(tree_B$tip.label, paste0("nodeB_", 1:Nnode(tree_B)))
  node_age_B[node_age_B<1e-6] <- 0
  
  ntip_A <- Ntip(tree_A)
  
  # initiation (origin clade B)
  list_nodes <- branches_t(tree_A, max(node_age_B))[,2]
  for (node in list_nodes[grep("node", list_nodes)]){if (initial_state %in% extract.clade(tree_A, node=Ntip(tree_A)+as.numeric(gsub("nodeA_","",node)))$tip.label){
    network_nodes[Ntip(tree_B)+1,node] <- 1
  }}
  
  
  for (edge in 1:nrow(tree_B$edge)){
    
    node_B1 <- tree_B$edge[edge,1]
    node_B2 <- tree_B$edge[edge,2]
    state <- network_nodes[node_B1,]
    
    age_B1 <- node_age_B[node_B1]
    age_B2 <- node_age_B[node_B2]
    
    time <- age_B1
    
    while (time>age_B2) {
      
      nb_inter <- sum(state)
      nb_A <- nrow(branches_t(tree_A, time)) # nb_A may change
      if (nb_inter<nb_A) {rate_gain_t <- rate_gain} else {rate_gain_t <- 0}
      if (nb_inter>1) {rate_loss_t <- rate_loss} else {rate_loss_t <- 0}
      
      waiting_time <- rexp(n=1, rate = rate_loss_t + rate_gain_t)
      
      # speciation and/or loss events in clade A: update the repertoire of clade B
      list_node_A <- intersect(which(node_age_A<time), which(node_age_A>max(c(age_B2,(time-waiting_time)))))
      if (length(list_node_A)>0) {
        for (node_A in names(node_age_A)[list_node_A]){
          if (length(grep("nodeA",node_A))==1){ # speciation in clade A = keep interactions with daughter lineages
            if (state[node_A]==1){
              state[node_A] <- 0
              state[tree_A$edge[which(tree_A$edge[,1]==(ntip_A+as.numeric(gsub("nodeA_","",node_A)))),2]] <- 1
            }
          } else { # extinction in clade A = loss interaction 
            state[node_A] <- 0
            if (sum(state)==0){ # if loss of all interactions, sample of interaction at random
              state[sample(size=1, branches_t(tree_A, max(c(age_B2,(time-waiting_time))))[,2])] <- 1
            }
            
          }
        }
      }
      
      time <- time - waiting_time
      
      if (time > age_B2){ # loss or acquisition in the repertoire of clade B
        
        # update rates
        nb_inter <- sum(state)
        nb_A <- nrow(branches_t(tree_A, time))
        if (nb_inter<nb_A) {rate_gain_t <- rate_gain} else {rate_gain_t <- 0}
        if (nb_inter>1) {rate_loss_t <- rate_loss} else {rate_loss_t <- 0}
        
        if (runif(1)<(rate_loss_t/(rate_loss_t+rate_gain_t))){ # loss
          state[sample(size=1, names(which(state==1)))] <- 0
        } else { # gain
          list_nodes_A <- branches_t(tree_A, time)[,2]
          list_nodes_A <- names(which(state[list_nodes_A]==0))
          if (length(list_nodes_A)>0) state[sample(size=1, list_nodes_A)] <- 1
        }
      }
    }
    network_nodes[node_B2,] <- state
  }
  
  return(network_nodes)
}



sim_clado_shift <- function(tree_B, tree_A, rate_gain, rate_loss, initial_state, seed=1){
  
  set.seed(seed)
  
  network_nodes <- matrix(0, nrow=Ntip(tree_B)+Nnode(tree_B), ncol=Ntip(tree_A)+Nnode(tree_A))
  colnames(network_nodes) <- c(tree_A$tip.label, paste0("nodeA_", 1:Nnode(tree_A)))
  rownames(network_nodes) <- c(tree_B$tip.label, paste0("nodeB_", 1:Nnode(tree_B)))
  
  node_age_A <- max(node.depth.edgelength(tree_A)) - node.depth.edgelength(tree_A)
  names(node_age_A) <- c(tree_A$tip.label, paste0("nodeA_", 1:Nnode(tree_A)))
  node_age_A[node_age_A<1e-6] <- 0
  node_age_B <- max(node.depth.edgelength(tree_B)) - node.depth.edgelength(tree_B)
  names(node_age_B) <- c(tree_B$tip.label, paste0("nodeB_", 1:Nnode(tree_B)))
  node_age_B[node_age_B<1e-6] <- 0
  
  ntip_A <- Ntip(tree_A)
  ntip_B <- Ntip(tree_B)
  
  # initiation (origin clade B)
  list_nodes <- branches_t(tree_A, max(node_age_B))[,2]
  for (node in list_nodes[grep("node", list_nodes)]){if (initial_state %in% extract.clade(tree_A, node=Ntip(tree_A)+as.numeric(gsub("nodeA_","",node)))$tip.label){
    network_nodes[Ntip(tree_B)+1,node] <- 1
  }}
  
  list_fossil_A <- tree_A$tip.label[grep("fossil_A_",tree_A$tip.label)]
  list_fossil_nodes_B <- c()
  for (fossil in tree_B$tip.label[grep("fossil_B_",tree_B$tip.label)]){
    list_fossil_nodes_B <- c(list_fossil_nodes_B, tree_B$node.label[tree_B$edge[which(tree_B$edge[,2]==which(tree_B$tip.label==fossil)),1]-ntip_B])
  }
  
  list_node_with_shift <- c() 
  
  for (edge in 1:nrow(tree_B$edge)){
    
    node_B1 <- tree_B$edge[edge,1]
    node_B2 <- tree_B$edge[edge,2]
    state <- network_nodes[node_B1,]
    
    name_B1 <- tree_B$node.label[node_B1-ntip_B]
    
    age_B1 <- node_age_B[node_B1]
    age_B2 <- node_age_B[node_B2]
    
    # Speciation by ecological interactions: Shift in 50% of the cases 
    
    if (!name_B1 %in% list_node_with_shift){
      list_node_with_shift <- c(list_node_with_shift, name_B1)
      if (!name_B1 %in% list_fossil_nodes_B){
        inter <- names(which(state==1))
        lineages_A <- branches_t(tree_A, t=age_B1)[,2]
        lineages_A <- lineages_A[!lineages_A %in% list_fossil_A] 
        lineages_A <- lineages_A[!lineages_A %in% inter]
        if (length(lineages_A)>0){
          state[state==1] <- 0
          state[sample(size=1, lineages_A)] <- 1 # one interaction at speciation
        }
      }
    }
    
    time <- age_B1
    
    while (time>age_B2) {
      
      nb_inter <- sum(state)
      nb_A <- nrow(branches_t(tree_A, time)) # nb_A may change
      if (nb_inter<nb_A) {rate_gain_t <- rate_gain} else {rate_gain_t <- 0}
      if (nb_inter>1) {rate_loss_t <- rate_loss} else {rate_loss_t <- 0}
      
      waiting_time <- rexp(n=1, rate = rate_loss_t + rate_gain_t)
      
      # speciation and/or loss events in clade A: update the repertoire of clade B
      list_node_A <- intersect(which(node_age_A<time), which(node_age_A>max(c(age_B2,(time-waiting_time)))))
      if (length(list_node_A)>0) {
        for (node_A in names(node_age_A)[list_node_A]){
          if (length(grep("nodeA",node_A))==1){ # speciation in clade A = keep interactions with daughter lineages
            if (state[node_A]==1){
              state[node_A] <- 0
              state[tree_A$edge[which(tree_A$edge[,1]==(ntip_A+as.numeric(gsub("nodeA_","",node_A)))),2]] <- 1
            }
          } else { # extinction in clade A = loss interaction 
            state[node_A] <- 0
            if (sum(state)==0){ # if loss of all interactions, sample of interaction at random
              state[sample(size=1, branches_t(tree_A, max(c(age_B2,(time-waiting_time))))[,2])] <- 1
            }
            
          }
        }
      }
      
      time <- time - waiting_time
      
      if (time > age_B2){ # loss or acquisition in the repertoire of clade B
        
        # update rates
        nb_inter <- sum(state)
        nb_A <- nrow(branches_t(tree_A, time))
        if (nb_inter<nb_A) {rate_gain_t <- rate_gain} else {rate_gain_t <- 0}
        if (nb_inter>1) {rate_loss_t <- rate_loss} else {rate_loss_t <- 0}
        
        if (runif(1)<(rate_loss_t/(rate_loss_t+rate_gain_t))){ # loss
          state[sample(size=1, names(which(state==1)))] <- 0
        } else { # gain
          list_nodes_A <- branches_t(tree_A, time)[,2]
          list_nodes_A <- names(which(state[list_nodes_A]==0))
          if (length(list_nodes_A)>0) state[sample(size=1, list_nodes_A)] <- 1
        }
      }
    }
    network_nodes[node_B2,] <- state
  }
  
  return(network_nodes)
}




sim_host_repertoire_diversification <- function(time_tot, species_A, rate_gain, rate_loss, initial_state,
                                                f_lambda, f_mu, f_rate_gain, f_rate_loss, seed){
  
  set.seed(seed)
  
  nb_A <- length(species_A)
  
  host_repertoire_t <- matrix(0, nrow=2, ncol=nb_A)
  colnames(host_repertoire_t) <- species_A
  rownames(host_repertoire_t) <- paste0("B", 1:2)
  host_repertoire_t[,initial_state] <- 1
  host_repertoire <- c() # final host_repertoire
  time=0.000001
  
  tree_B <- phytools::pbtree(n=2)
  tree_B$edge.length <- c(time,time)
  list_extant_B <- tree_B$tip.label <- c("B1", "B2")
  nb_B <- 2
  nb_alive <- 2
  
  list_events <- c("extinction", "speciation", "gain", "loss")
  
  while ((time < time_tot) & (nb_alive > 0)){
    
    rate_extinction <- f_mu(host_repertoire_t)
    rates_speciation <- f_lambda(host_repertoire_t) 
    rates_gain <- sapply(rowSums(host_repertoire_t), function(n) f_rate_gain(n, rate_gain, nb_A))
    rates_loss <- sapply(rowSums(host_repertoire_t), function(n) f_rate_loss(n, rate_loss))
    
    waiting_time <- rexp(n=1, rate = sum(rate_extinction) + sum(rates_speciation) + sum(rates_gain) + sum(rates_loss))
    
    if (time + waiting_time < time_tot){
      
      event <- sample(list_events, size = 1, prob = c(sum(rate_extinction), sum(rates_speciation), sum(rates_gain), sum(rates_loss)))
      
      if (event=="extinction"){
        # extend branch lengths
        tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] <- tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] + waiting_time
        # select species
        species <- sample(rownames(host_repertoire_t), size = 1, prob = rate_extinction) 
        host_repertoire <- rbind(host_repertoire, host_repertoire_t[which(rownames(host_repertoire_t)==species), , drop=FALSE])
        host_repertoire_t <- host_repertoire_t[which(rownames(host_repertoire_t)!=species), , drop=FALSE]
        # drop extinct species
        list_extant_B <- list_extant_B[list_extant_B!=species]
        nb_alive <- nb_alive - 1
      }
      if (event=="speciation"){
        # extend branch lengths
        tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] <- tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] + waiting_time
        # select species
        species <- sample(rownames(host_repertoire_t), size = 1, prob = rates_speciation)
        nb_B <- nb_B+1
        host_repertoire_t <- rbind(host_repertoire_t, host_repertoire_t[species,])
        
        rownames(host_repertoire_t)[nrow(host_repertoire_t)] <- paste0("B",nb_B)
        # add species
        tree_B <- phytools::bind.tip(tree_B, paste0("B",nb_B), edge.length = 0, position = 0, where = which(tree_B$tip.label==species))
        list_extant_B <- c(list_extant_B, paste0("B",nb_B))
        nb_alive <- nb_alive + 1
      }
      if (event=="gain"){
        # extend branch lengths
        tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] <- tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] + waiting_time
        # select species
        species <- sample(rownames(host_repertoire_t), size = 1, prob = rates_gain)
        host_repertoire_t[species,sample(size=1, names(which(host_repertoire_t[species,]==0)))] <- 1
      }
      if (event=="loss"){
        # extend branch lengths
        tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] <- tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] + waiting_time
        # select species
        species <- sample(rownames(host_repertoire_t), size = 1, prob = rates_loss)
        host_repertoire_t[species,sample(size=1, names(which(host_repertoire_t[species,]==1)))] <- 0
      }
      
    } else { # no events
      # extend branch lengths
      tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] <- tree_B$edge.length[which(tree_B$edge[,2] %in% which(tree_B$tip.label %in% list_extant_B))] + time_tot - time
    }
    
    time <- time+waiting_time
    
  }
  
  host_repertoire <- rbind(host_repertoire, host_repertoire_t)
  
  if (nb_alive==0) {alive=FALSE} else {alive=TRUE} 
  
  return(list(host_repertoire=host_repertoire, tree_B=tree_B, alive=alive))
}



