# Go further: Real data input

For the project about the voter model we will use some indicators based on real data. We will also use a new faster voter model based on a scale free graph.

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

vote <- function(igraph,noise,N_vect)
{
  nb_node <- vcount(igraph)
  vect_influence <- vector(mode = "numeric",nb_node)
  vect_nb_neighb <- vector(mode = "numeric",nb_node)
  f_p <-  vector(mode = "numeric",nb_node)

  adj_mat <- as_adjacency_matrix(igraph)
  vect_vote <- get.vertex.attribute(igraph,"vote")
  for(i in 1: nb_node)
  {
    vect_influence[i] <-  vect_vote%*%adj_mat[i,]
    vect_nb_neighb[i] <- sum(adj_mat[i,])

    if(vect_nb_neighb[i] == 0)
    {
      f_p[i] <- vect_vote[i]
    }else{
      f_p[i] <- (1-2*noise)*(vect_influence[i]/vect_nb_neighb[i])+noise
    }

    if( f_p[i] > 0.5 )
    {
      igraph <- set.vertex.attribute(igraph,"vote",i,value = 1)
    }else{
      igraph <- set.vertex.attribute(igraph,"vote",i,value = 0)
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
  igraph <- initVoteBernoulli(igraph)
  vote_plot(igraph,"Init")
  mat <- vector(mode = "numeric",(vcount(igraph)*time))
  dim(mat) <- c(vcount(igraph),time)
  for(i in 1: time)
  {
    igraph <- vote(igraph,noise)
    cat(i," over ",time,"\n")
  }
  declare_winner(mat)
  mat
}

n <- 20
p <- 0.1

vote_plot <- function(igraph,mainstr)
{
  nb_node <- vcount(igraph)
  for(i in 1: nb_node)
  {
    node_status <- get.vertex.attribute(igraph,"vote", i)
    if( node_status == 1 )
    {
      igraph <- set.vertex.attribute(igraph, 'color', i, rgb(1,0,0) )
    }
    else if (node_status == 0 )
    {
      igraph <- set.vertex.attribute(igraph, 'color', i, rgb(0.9,0.9,0.7) )
    }
    else{
      print("BUG PLOT attribute")
    }
  }
  plot(igraph, main = mainstr )
}

igraph <- sample_gnp(n, p, directed = FALSE, loops = FALSE)


noise <- 0.1
time <- 3

simulation(igraph,noise,time)
```

Now on top of this model we can add some real data. First, we will load some date about age repartition of the US population. Here is the source of data : [https://www.census.gov/data/tables/time-series/demo/popest/2010s-national-detail.html\#par\_textimage\_1537638156](https://www.census.gov/data/tables/time-series/demo/popest/2010s-national-detail.html#par_textimage_1537638156) This how you can load .cvs files.

```r
path <- file.path("yourpathtofile.csv")

pop_data <- read.table(path,
            header = TRUE,
            sep = "\t",
            na.strings = "n/a",
            stringsAsFactors = FALSE)


nb_node <- 1000
p <- 0.01
igraph <- sample_gnp(nb_node, p, directed = FALSE, loops = FALSE)

vect_age <- pop_data[,5]
vect_pop_age <- pop_data[,6]

poptot <- sum(vect_pop_age)
popbynode <- poptot/nb_node
vect_node_pop_age <- floor(vect_pop_age/popbynode)

age_it <- 1
node_it <- 1
while(sum(vect_node_pop_age) > 0)
{
  while(vect_node_pop_age[age_it] > 0)
  {
    igraph <- set.vertex.attribute(igraph,"age",node_it, value = age_it)
    vect_node_pop_age[age_it] <- vect_node_pop_age[age_it] - 1
    node_it <- node_it + 1
  }
  age_it <- age_it + 1
}
```

Notice that modeling a population of 300 million inhabitants by only 1000 nodes creates uncertainties on our model.

This is the kind of problem that must be considered during development on large and complex systems.

You can now include data into your simulation and try to predict the future of political opinion, spread of epidemic or population growth.

### Conclusion:

 Graph network simulation is great way to modelized complex phenomena. Creating a network and perform simulation on it is simple and fun. Moreover, working with nodes attributes make the developpement flexible and allow you to improve the model incrementally and to infinity. The only limits you will met is the time you have to create the model, the computing power you have \( It is therefore necessary to optimize the algorithms or the code\) and the quality of the data you will use as inputs.



A faire : Comparé l’effet des confinements avec la limite de percolation \(immunité collective, surtout avec les multiples confinements\)

Mettre plus en avant les resultats..

