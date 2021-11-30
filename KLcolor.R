KLcolor <- function(
  object, 
  group.by = NULL,
  reduction = "umap", 
  min.alpha = 0.2, 
  seed = 4213,
  color.initial.set = NULL, 
  col.distance.metric = "manhattan",
  return.imtermediate = FALSE
) {
  
  if (!is.null(group.by)) {
    Idents(object) <- group.by
  }
  
  smooth.matrix <- as.sparse(
    x = fastDummies::dummy_cols(
      object@active.ident
    )[, -1]
  ) 
  
  colnames(smooth.matrix) <- gsub( pattern =".data_",
                                   replacement = "",
                                   x = colnames(smooth.matrix)
  )
  
  cluster.freq <- colSums(smooth.matrix)
  smooth.matrix <- sweep(x = smooth.matrix,
                         MARGIN = 2, 
                         STATS = cluster.freq,
                         FUN = "/"
  )
  
  cluster.embeddings <- as.matrix(t(smooth.matrix) %*% object[[reduction ]]@cell.embeddings)
  cluster.embeddings <- apply(X = cluster.embeddings, MARGIN = 2, function(x) {
    x <- (x - min(x))/(max(x) - min(x))
    return(x)
  })
  
  cluster.dist <- rdist::cdist(X = cluster.embeddings, Y = cluster.embeddings)   
  rownames(cluster.dist) <- colnames(cluster.dist) <- rownames( cluster.embeddings)
  cluster.kernel <- exp((-10)*cluster.dist**2)
  if (!is.null(color.initial.set)) {
    extra.color <-  length(cluster.freq) - length(color.initial.set)
    if (extra.color > 0) {
      warning("Input color set is smaller than the number of clusters. \n",
              "The rest ", extra.color, 
              " colors will be sampled from the input color set")
      
    }
    color.set.0 <- t(col2rgb(color.initial.set, alpha = F))
    f.order <- function(x) {
      x.int <- round(x)
      x.int <- MinMax(data = x.int, min = 1, max = nrow(color.set.0))
      W.f <- color.set.0[x.int,]/255
      W.dist <-  rdist::pdist(W.f, metric = col.distance.metric) + 0.01
      diag(W.dist) <- 1
      cost <- cluster.kernel*log(cluster.kernel *  (1/W.dist))
      penalty <- 1e2*(sum((x.int - x) **2) )
      cost <- sum(cost) + penalty
      return(cost)
    }
    ## discrete optimization
    set.seed(seed)
    nextfun <- function(x) sample(1:nrow(color.set.0), ncol(cluster.kernel), replace=TRUE)
    random.order <-  sample(1:nrow(color.set.0), size =  ncol(cluster.kernel), replace = T)
    r <- optim(par = random.order , 
               fn = f.order,
               gr = nextfun,
               method="SANN", 
               control =  list(maxit = 100, fnscale=1)
    )
    order.r <- MinMax(data = round(r$par), min = 1, max = nrow(color.set.0))
    w <- color.set.0[order.r,]
    rownames(w) <- rownames(cluster.embeddings)
    
  } else {
    f.W <- function(x) {
      W.f <- matrix(data = x, ncol = 3)
      W.dist <-  rdist::pdist(X = W.f, metric = col.distance.metric)
      W.dist <- W.dist/255
      diag(W.dist) <- 1
      cost <- cluster.kernel*log(cluster.kernel *  (1/W.dist))
      cost <- sum(cost) + sum(1e5*ReLu(W.f - 1)) + sum(1e5*ReLu(W.f*(-1))) + sum(1e3*W.f)
      return(cost)
    }
    set.seed(seed)
    W.init <- matrix(data = sample(x = 1:255,
                                   size = 3*ncol(cluster.kernel),
                                   replace = T),
                     nrow = nrow(cluster.kernel), 
                     ncol = 3
    )
    r <- optim(as.vector(W.init), f.W)
    w <- matrix(data = r$par, ncol = 3)
  }
  w <- MinMax(data = w/255, min = 0, max = 1)
  # w <- (w - min(w))/(max(w) - min(w))
  cluster.alpha <- 1/sqrt(cluster.freq)
  cluster.alpha <- min.alpha + 
    (1 - min.alpha)*
    (cluster.alpha - min(cluster.alpha))/
    (max(cluster.alpha) - min(cluster.alpha))
  
  cols.set <- sapply(X = 1:nrow(w), FUN = function(x) {
    col.x <-  rgb(red = w[x,1], 
                  green =  w[x,2],
                  blue =  w[x,3],
                  alpha = cluster.alpha[x] , 
                  maxColorValue = 1
    )
    return(col.x)
  })
  names(cols.set) <- rownames(cluster.embeddings)
  
  if (return.imtermediate) {
    W.dist <-  rdist::pdist(X = w, metric = col.distance.metric)
    output.list <- list(color = cols.set,
                        W_col =  w,
                        W_dist = W.dist,
                        cluster.kernel = cluster.kernel)
    return(output.list)
  } else {
    return (cols.set)
  }
  
}
