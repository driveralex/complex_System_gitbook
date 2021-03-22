# Representation

The model described is a simple agent-based on an undirected graph. Each node has only one attribute that represente his vote.

First we initialise the network with a Bernoulli's law of parameter $$0.5$$ :

```r
initVoteBernoulli <- function(igraph)
{
  for( i in 1: vcount(igraph))
  {
    if( runif(1) < 1/2 )
    {
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=1)
    }
    else
    {
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=0)
    }
  }
  igraph
}
```

Then each node is influence by his neigbors. Here is fast implementation that work close to adjency matrix :

```r
Nb_neighbors_vect <- function(igraph)
{
  mat <- as_adjacency_matrix(igraph)
  len <- nrow(mat)
  N_vect <- vector(mode = "numeric",len)
  for(i in 1 : len )
  {
    res <- sum(mat[i,])
    if(res > 0)
    {
      N_vect[i] <- sum(mat[i,])
    }else if(res == 0)
    {
      N_vect[i] <- 1
    }
    else{
      print("BUG - Nb_neighbors_vect")
    }
  }
  N_vect
}

probability_vect <- function(igraph,N_vect)
{
  adj_mat <-  as_adjacency_matrix(igraph, sparse = FALSE)
  vect_vote <- get.vertex.attribute(igraph, "vote")
  prob_vect <- adj_mat%*%vect_vote
  out <- prob_vect/N_vect
  out
}

getOneVote <- function(igraph)
{
  vc <- vcount(igraph)
  vect <- vector(mode = "numeric",vc)
  for( i in 1 : vc)
  {
    vect[i] <- get.vertex.attribute(igraph, "vote",i) 
  }
  vect
}

vote_probability <- function(igraph,noise,N_vect)
{
  p_v <- probability_vect(igraph,N_vect)
  len <- length(p_v)
  f_p_opti <- vector(mode = "numeric",len)
  f_p_opti <- (1-2*noise)*p_v + noise
  f_p_opti
}

vote <- function(igraph,noise,N_vect)
{

  f_p <- vote_probability(igraph,noise,N_vect)
  vc <- vcount(igraph)
  for( i in 1: vc)
  {
    if( runif(1) < f_p[i] )
    {
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=0)
    }
    else
    {
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=1)
    }
  }
  igraph
}
declare_winner <- function(mat)
{
  end_vote_vect <- mat[,ncol(mat)]

  if( sum(end_vote_vect)/nrow(mat) == 0.5  )
  {
    cat("Deuce\n")
  }else if( sum(end_vote_vect)/nrow(mat) < 0.5 )
  {
    cat("Jerry win\n")
  }else if( sum(end_vote_vect)/nrow(mat) > 0.5 )
  {
    cat("Tom win\n")
  }
}

simulation <- function(igraph,noise,time)
{
  N_vect <- Nb_neighbors_vect(ig)

  igraph <- initVoteBernoulli(igraph)
  mat <- vector(mode = "numeric",(vcount(igraph)*time))
  dim(mat) <- c(vcount(igraph),time)
  for(i in 1: time)
  {
    mat[,i] <- getOneVote(igraph)
    igraph <- vote(igraph,noise,N_vect)
    cat(i," over ",time,"\n")
  }
  declare_winner(mat)
  mat
}
```

Then you just have to give the network of your choice to the simulation and you have completed the matrix that represent the evolution of votes.

Before the development of this implementation I have built another model, which is way less time-efficient. Here is the code:

```r
probability_node <- function(igraph, node)
{
  adj_lst_node <- as_adj_list(igraph)[[node]]
  if( length( adj_lst_node ) == 0 )
  {
    return(0)  
  }
  sum <- 0
  for(i in 1 : length( adj_lst_node ) )
  {
    sum <- sum + get.vertex.attribute(igraph, "vote",    adj_lst_node[i]    ) 
  }
  out <- sum/ length( adj_lst_node )
  out
}

probability_vect <- function(igraph)
{
  vc <- vcount(igraph)
  vect <- vector(mode = "numeric",vc)
  sum <- 0
  for( i in 1 : vc )
  {
    sum <- sum + probability_node(igraph,i)
    vect[i] <- probability_node(igraph,i)
  }
  vect
}

probability_tot <- function(igraph)
{
  vect <- probability_vect(igraph)
  out <- sum(vect) / vcount(igraph)
}

getOneVote <- function(igraph)
{
  vc <- vcount(igraph)
  vect <- vector(mode = "numeric",vc)
  for( i in 1 : vc)
  {
    vect[i] <- get.vertex.attribute(igraph, "vote",i) 
  }
  vect
}

vote_probability <- function(igraph,noise)
{
  p_v <- probability_vect(igraph)
  len <- length(p_v)
  f_p <- vector(mode = "numeric",len)
  for(i in 1 : len  )
  {
    f_p[i] <- (1-2*noise)*p_v[i] + noise
  }
  f_p
}

vote <- function(igraph,noise)
{
  start_time <- Sys.time()
  end_time <- Sys.time()
  f_p <- vote_probability(igraph,noise)
  end_time <- Sys.time()
  exec_time <- (end_time-start_time)
  print(exec_time)
  vc <- vcount(igraph)
  for( i in 1: vc)
  {
    if( runif(1) < f_p[i] )
    {
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=0)
    }
    else
    {
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=1)
    }
  }
  igraph
}

simulation <- function(igraph,noise,time)
{
  igraph <- initVoteBernoulli(igraph)
  mat <- vector(mode = "numeric",(vcount(igraph)*time))
  dim(mat) <- c(vcount(igraph),time)
  for(i in 1: time)
  {
    mat[,i] <- getOneVote(igraph)
    igraph <- vote(igraph,noise)
    cat(i," over ",time,"\n")
  }
  mat
}
```

Notice the difference introduce by working with matrix in comparison with for loops.

The evolution of the voting rate can be represented via the output matrix; however, this data still needs to be simplified to be more human readable. The role of the declare winner function is to inform us about the winner if the vote is done at the end of the simulation.

