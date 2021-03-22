# Erdos-Renyi model

The Graph Erdos-Renyi $$G(n, p) = (V,E)$$ constructed from a set $$V$$ of $$n$$ vertices. The edge between 2 vertices $$i$$ and $$j$$ exists with probability $$p$$.

Below is a proposition of algorithm :

```r
Erdos_Renyi <- function(n,p)
{ 
  igraph  <- make_empty_graph(directed = FALSE)
  igraph <- add_vertices(igraph,n, color = "red")
  for(i in 1 : n)
  {
    for(j in i : n)
    {
      if( runif(1) < p )
      {
        igraph <- add_edges(igraph, c(i,j))
      }
    }
  }
  igraph
}
```

This implementation gives the right output. However due to the 2 for loop performance is quite bad and it can take several minutes with large graph.

A faster implementation is proposed below with approach closer to adjacency matrix.

Note the transpose trick used to make the matrix symmetric.

```r
Erdos_Renyi_optimized <- function(n,p)
{
  nb_tri <- n*(n-1)/2
  mat <- diag(0,nrow = n, ncol = n)
  vect_rand <- runif(nb_tri , min = 0, max = 1) - p
  for(i in 1:nb_tri)
  {
    if(vect_rand[i]>0)
    {
      vect_rand[i] <- 0 
    }
    else
    {
      vect_rand[i] <- 1 
    }
  }
  mat[lower.tri(mat)] <- vect_rand
  mat <- t(mat)
  mat[lower.tri(mat)] <- vect_rand
  igraph <- graph_from_adjacency_matrix(mat,mode = c("undirected"))
  igraph
}
```



