# Influenceing scenarios

With a scale free graph of 501 nodes \(to avoid Deuce\), $$k = 3$$, $$m = 2$$, and for simulation of $$0.01$$ noise and $$time = 3000$$, we will try to influence the vote. We will work in Jerry's Team.

### **Scenario A**

We have the possibility convince 10 people to vote for Jerry. We will use the metrics to identify them.

We want to have:

* Node connected to the giant component.
* Node with a small length to all the others nodes.
* Node with high closeness centrality.

Functions n\_length, reacheable\_node and how\_closness are used to identify nodes with a high influence.

On the other side we must take in consideration the limitation on the sum of degree of 100. First, we use the how\_bad\_is\_degree function. Then the function high\_cut also play a role to eliminate the nodes with highest degree with a nonlinear function.

All those parameters can be tuned by hand with coefficient to improve the ranking.

Then to ensure that we do respect the degree condition we cycle down the top result until the condition is respected.

Here is the code to do the selection of node to influence:

```r
selection <- function(igraph)
{
  vect_len        <- n_length(igraph)
  vect_reach      <- reacheable_node(igraph)
  vect_close      <- how_closness(igraph)
  vect_bad_degree <- how_bad_is_degree(igraph)

  coef_len   <- 3
  coef_reach <- 5
  coef_close <- 1.5

  coef_bad_degree <- 1
  alpha <- 0.013

  vect_rank_positif <- vect_len*coef_len + vect_reach*coef_reach + vect_close*coef_close

  vect_rank_negatif <- vect_bad_degree*coef_bad_degree
  vect_rank_negatif <- high_cut(vect_rank_negatif,alpha)
  print(vect_rank_negatif)
  print(degree(igraph))

  vect_rank_tot <- vect_rank_positif - vect_rank_negatif

  nb_elect <- 10
  pool <- vector(mode = "numeric",nb_elect)
  ref <- 1
  condition_is_ok <- FALSE
  while(condition_is_ok == FALSE)
  {
    sorted <- sort.int(vect_rank_tot,decreasing = TRUE, index.return = TRUE)
    for(i in ref : (ref+nb_elect) )
    {
      pool[(i-ref)] <- sorted$ix[i] 
    }
    vect_deg <- degree(igraph)
    somm <- 0
    #print( vect_deg[pool] )
    for(i in 1:length(pool))
    {

      somm <- somm + vect_deg[pool[i]]
    }
    if(100 > somm)
    {
      cat("Number of rank down to respect condition=",ref," (SUM of degree=",somm,")\n")
      condition_is_ok <- TRUE
    }
    ref <- ref + 1
  }
  pool
}

high_cut <- function(inputvect,alpha)
{
  outputvect <- inputvect + alpha*inputvect*inputvect
}

how_closness <- function(igraph)
{
  tryCatch( vect_close <- closeness(igraph,normalized = TRUE) , warning=function() print("-") )
  if( min(vect_close) == max(vect_close) )
  {
    print("May have pb w: how_closness")
    return( c(rep(0,vcount(igraph)) ))
  }
  n_comp_min <- min(vect_close)
  vect_close <- vect_close-n_comp_min
  n_comp_max <- max(vect_close)
  vect_close <- vect_close*100/n_comp_max
  vect_close
}

how_bad_is_degree <- function(igraph)
{
  vect_deg <- degree(igraph) 
  if( min(vect_deg) == max(vect_deg) )
  {
    print("May have pb w: how_bad_is_degree")
    return( c(rep(0,vcount(igraph)) ))
  }
  n_comp_min <- min(vect_deg)
  vect_deg <- vect_deg-n_comp_min
  n_comp_max <- max(vect_deg)
  vect_deg <- (vect_deg*100/n_comp_max)+1
  vect_deg
}



reacheable_node <- function(igraph)
{
  Comp = clusters(igraph)
  vect_reach <- Comp$csize[Comp$membership]
  if( min(vect_reach) == max(vect_reach) )
  {
    #print("One giant Component")
    return( c(rep(0,vcount(igraph)) ))
  }
  n_comp_min <- min(vect_reach)
  vect_reach <- vect_reach-n_comp_min
  n_comp_max <- max(vect_reach)
  vect_reach <- vect_reach*100/n_comp_max
  vect_reach
}


n_length <- function(igraph)
{
  nb_node <- vcount(igraph)
  n_length_vect <- vector(mode = "numeric",nb_node)
  for( i in 1:  nb_node)
  {
    bfs_i <- bfs(igraph,i,
                 unreachable = FALSE,
                 dist = TRUE);
    #print(bfs_i$dist)
    n_length_vect[i] <- sum(bfs_i$dist,na.rm = TRUE)/nb_node  
  }
  n_length_max <- max(n_length_vect)

  for( i in 1:  nb_node)
  {
    if( n_length_vect[i] == 0 )
    {
      n_length_vect[i] <- n_length_max
    }
  }
  n_length_min <- min(n_length_vect)


  n_length_vect <- (n_length_vect) - n_length_min
  n_length_max <- max(n_length_vect)
  n_length_vect <- 100 - (n_length_vect*100/n_length_max)
}
```

Then the set node chosen is tested in simulation, with they vote force to 1. To ensure that that our methods is ok we repeat simulation multiple time and observe the result. Below is the code to test the model.

```r
initVoteBernoulli <- function(igraph,sway_vector)
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
  nb_sway <- length(sway_vector)
  for( i in 1 : nb_sway )
  {
    igraph <- set.vertex.attribute(igraph,"vote", sway_vector[i] ,value=0)
  }
  igraph
}

vote <- function(igraph,noise,N_vect,sway_vector)
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
  nb_sway <- length(sway_vector)
  for( i in 1:nb_sway )
  {
    igraph <- set.vertex.attribute(igraph,"vote", sway_vector[i] ,value=0)
  }
  igraph
}
declare_winner <- function(mat)
{
  end_vote_vect <- mat[,ncol(mat)]
  out <- vector(mode = "numeric",2)

  if( sum(end_vote_vect)/nrow(mat) == 0.5  )
  {
    cat("Deuce\n")
  }else if( sum(end_vote_vect)/nrow(mat) < 0.5 )
  {
    cat("Jerry win with ",(sum(end_vote_vect)/nrow(mat)),"\n")
    out[1] <- 0
  }else if( sum(end_vote_vect)/nrow(mat) > 0.5 )
  {
    cat("Tom win with ",(sum(end_vote_vect)/nrow(mat)),"\n")
    out[1] <- 1
  }
  out[2] <- (sum(end_vote_vect)/nrow(mat)) 
  out
}


simulation <- function(igraph,noise,time,sway_vector)
{
  N_vect <- Nb_neighbors_vect(igraph)

  igraph <- initVoteBernoulli(igraph,sway_vector)
  mat <- vector(mode = "numeric",(vcount(igraph)*time))
  dim(mat) <- c(vcount(igraph),time)
  for(i in 1: time)
  {
    mat[,i] <- getOneVote(igraph)
    igraph <- vote(igraph,noise,N_vect,sway_vector)
    #cat(i," over ",time,"\n")
  }
  win <- declare_winner(mat)
}

n <- 501
k <- 3
q <- 2

igraph <- scale_free(n,k,q)
sway_vector <- selection(igraph)

noise <- 0.01
time <- 3000

out <-  0
nb_simu <- 100
winner <- 0
score <- 0
for(i in 1: nb_simu )
{
  out <- simulation(igraph,noise,time,sway_vector)
  winner <- winner + out[1]
  score <- score + out[2]
  cat("Progress:",i*100/nb_simu,"%\n")
}
cat("TOM win in ",winner/nb_simu,"%\n")
cat("Avreage score is",score/nb_simu,"%\n")
```

Unfortunaly the results did not seem to be influenced by our method.

### Introducting Zealots

Zealots are people that can't change of opinion. To add them into our model we first define on witch node they are:

```r
zelot <- function(igraph)
{
  nb_node <- floor(0.4*vcount(igraph))
  nb_node
  vect_zelot <- vector(mode = "numeric",vcount(igraph))
  for(i in 1: vcount(igraph) )
  {
    if( runif(1)>0.8 )
    {
      vect_zelot[i] <- 1
    }else if( runif(1) < 0.2 )
    {
      vect_zelot[i] <- 0
    }else
    {
      vect_zelot[i] <- -1
    }
  }
  vect_zelot
}

initZelot <- function(igraph,zelot_vector)
{
  for( i in 1: vcount(igraph))
  {
    if( zelot_vector[i] == 1 )
    {
      igraph <- set.vertex.attribute(igraph,"zelot", i ,value=1)
    }
    else if (zelot_vector[i] == 0 )
    {
      igraph <- set.vertex.attribute(igraph,"zelot", i ,value=0)
    }else
    {
      igraph <- set.vertex.attribute(igraph,"zelot", i ,value=-1)
    }
  }
  igraph
}
```

The voting function is adapted:

```r
initVoteBernoulli <- function(igraph,sway_vector,zelot_vector)
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
  for( i in 1: vcount(igraph))
  {
    if( zelot_vector[i] == 1 )
    {
      igraph <- set.vertex.attribute(igraph,"zelot", i ,value=1)
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=1)
    }
    else if (zelot_vector[i] == 0 )
    {
      igraph <- set.vertex.attribute(igraph,"zelot", i ,value=0)
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=0)
    }else
    {
      igraph <- set.vertex.attribute(igraph,"zelot", i ,value=-1)
    }
  }
  nb_sway <- length(sway_vector)
  for( i in 1 : nb_sway )
  {
    igraph <- set.vertex.attribute(igraph,"vote", sway_vector[i] ,value=1)
  }
  igraph
}

vote <- function(igraph,noise,N_vect,sway_vector,zelot_vector)
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
  for( i in 1: vcount(igraph))
  {
    if( zelot_vector[i] == 1 )
    {
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=1)
    }
    else if (zelot_vector[i] == 0 )
    {
      igraph <- set.vertex.attribute(igraph,"vote", i ,value=0)
    }
  }

  nb_sway <- length(sway_vector)
  for( i in 1:nb_sway )
  {
    igraph <- set.vertex.attribute(igraph,"vote", sway_vector[i] ,value=1)
  }
  igraph
}
```

Then as for the others parameters we remove them from the ranking with the following trick:

```r
is_zelot <- function(igraph)
{
  nb <- vcount(igraph)
  vect_zelot <- vector(mode = "numeric",nb)
  for( i  in 1: nb)
  {
    if( get.vertex.attribute(igraph, "zelot",i) == 1 || get.vertex.attribute(igraph, "zelot",i) == 0 )
    {
      vect_zelot[i] <- 99999999
    }else if(get.vertex.attribute(igraph, "zelot",i) == -1 )
    {
      vect_zelot[i] <- 0
    }
  }
  vect_zelot
}

selection <- function(igraph)
{
  vect_len        <- n_length(igraph)
  vect_reach      <- reacheable_node(igraph)
  vect_close      <- how_closness(igraph)
  vect_bad_degree <- how_bad_is_degree(igraph)
  vect_zelot      <- is_zelot(igraph)

  coef_len   <- 3
  coef_reach <- 5
  coef_close <- 1.5

  coef_bad_degree <- 1
  alpha <- 0.013

  vect_rank_positif <- vect_len*coef_len + vect_reach*coef_reach + vect_close*coef_close

  vect_rank_negatif <- vect_bad_degree*coef_bad_degree + vect_zelot
  vect_rank_negatif <- high_cut(vect_rank_negatif,alpha)
  vect_rank_tot <- vect_rank_positif - vect_rank_negatif

  print(vect_rank_tot)

  nb_elect <- 10
  pool <- vector(mode = "numeric",nb_elect)
  ref <- 1
  condition_is_ok <- FALSE
  while(condition_is_ok == FALSE)
  {
    sorted <- sort.int(vect_rank_tot,decreasing = TRUE, index.return = TRUE)
    for(i in ref : (ref+nb_elect) )
    {
      pool[(i-ref)] <- sorted$ix[i] 
    }
    vect_deg <- degree(igraph)
    somm <- 0
    #print( vect_deg[pool] )
    for(i in 1:length(pool))
    {

      somm <- somm + vect_deg[pool[i]]
    }
    if(100 > somm)
    {
      cat("Number of rank down to respect condition=",ref," (SUM of degree=",somm,")\n")
      condition_is_ok <- TRUE
    }
    ref <- ref + 1
  }
  pool
}
```

Then to run the simulation do the following call in the next order:

```r
n <- 501
k <- 3
q <- 2

igraph <- scale_free(n,k,q)
zelot_vector <- zelot(igraph)
igraph <- initZelot(igraph,zelot_vector)
sway_vector <- selection(igraph)


noise <- 0.001
time <- 3000

out <-  0
nb_simu <- 3
winner <- 0
score <- 0
for(i in 1: nb_simu )
{
  out <- simulation(igraph,noise,time,sway_vector,zelot_vector)
  winner <- winner + out[1]
  score <- score + out[2]
  cat("Progress:",i*100/nb_simu,"%\n")
}
cat("TOM win in ",winner/nb_simu,"%\n")
cat("Avreage score is",score/nb_simu,"%\n")
```

With the zealot the effect of the 10 influenced nodes seem increassed.

### Effect of remouving node

As for the nodes, we can influence the network by changing its topology. However, finding the edges to remove might be more complicated that for nodes in the voter model. Because my proposition of node influencing did not work well, and limitation in time to execute this function: I will not implement it.

