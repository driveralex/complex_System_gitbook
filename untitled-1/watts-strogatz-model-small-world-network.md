# Watts-Strogatz model : Small world network

The Watts–Strogatz model produces graphs with small-world properties.

During the implementation of this algorithm the first version created completed the task efficiently but was very slow \(multiple minutes for n=1000\) due to the cascade of for loops.

After rewriting the r code here is implementation that perform well and can produce the required the task quasi instantly.

```r
Watts_Strogatz_opt <- function(n,p,m)
{
  igraph <- make_ring(n, directed = FALSE)
  ## End step1
  for(j in 1:m)
  {
    for( i in 1: n)
    {
      igraph <- add_edges(igraph, c(i,((i+j)%%n+1)  ))
    }
  }
  ## End step2
  nb_tri <- n*(n-1)/2
  vect_rand <- runif(nb_tri , min = 0, max = 1) - p
  for(i in 1:nb_tri)
  {
    if(vect_rand[i]>0)
    {
      vect_rand[i] <- 1 
    }
    else
    {
      vect_rand[i] <- 0 
    }
  }
  sp <- matrix( nrow = n, ncol = n)
  sp[lower.tri(sp)] <- vect_rand
  sp <- t(sp)
  sp[lower.tri(sp)] <- vect_rand
  diag(sp) <- 0

  adj_mat <-  as_adjacency_matrix(igraph,type = c("both"))
  res_mat <- adj_mat*sp
  igraph <- graph_from_adjacency_matrix(res_mat,mode = c("undirected"))

  ## End Step 3 (remouving edges)

  act_size <- gsize(igraph)
  aim_size <- n*(2+2*m)/2
  nb_edges_to_add <- aim_size - act_size
  for( i in 1: nb_edges_to_add)
  {
    edge_added <- FALSE
    while(edge_added == FALSE)
    {
      proposal <- floor(runif(2,min = 1, max = n))
      if( (are_adjacent(igraph, proposal[1], proposal[2] ) == FALSE ) && ( proposal[1] != proposal[2] ) )
      {
        igraph <- add_edges(igraph, c( proposal[1], proposal[2] )  )
        edge_added <- TRUE
      }
    }
  }
  igraph
}
```

