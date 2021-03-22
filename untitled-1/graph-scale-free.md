# Graph scale free

The following code generate a scale free graph for any type of k.

```r
nb_init_nodes <- function(nb_init_edges)
{
  n_cal <- (1+sqrt(1+8*nb_init_edges))/2
  if( n_cal == floor(n_cal) )
  {
    return(n_cal)
  }else{
    n_cal <- floor(n_cal)+1
  }
  n_cal
}

scale_free_init <- function(k)
{
  nb_node <- nb_init_nodes(k)
  nb_tri <- nb_node*(nb_node-1)/2
  mat <- matrix( nrow = nb_node, ncol = nb_node)
  diag(mat) <- 0
  vect_init <- c( rep(1, k) , rep(0, (nb_tri-k) ))
  mat[lower.tri(mat)] <- vect_init
  mat <- t(mat)
  mat[lower.tri(mat)] <- vect_init  
  igraph <- graph_from_adjacency_matrix(mat,mode = c("undirected"))
  igraph
}

scale_free_degree_range <- function(igraph,random)
{
  random <- runif(1)
  vect_deg <- degree(igraph)
  nbedges <- sum(vect_deg)/2
  sum <- 0
    for(i in 1:length(vect_deg))
    {
        sum <- sum + vect_deg[i]/(nbedges*2)
        if( sum > random)
        {
          return(i)
        }
    }
}

scale_free <- function(n,k,q)
{
  igraph <- scale_free_init(k)
  add <- 0
  nodenb <- k
  while( nodenb < n)
  {
    igraph <- add_vertices(igraph,1)
    for(i in 1:q)
    {
      igraph <- add_edges(igraph,c( (nodenb-1+q) ,scale_free_degree_range(igraph)))
    }
    nodenb <- nodenb + 1
  }
  igraph
}
```

