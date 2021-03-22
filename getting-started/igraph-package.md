# Igraph package

The R Igraph package is dedicated for creation, manipulation and analysis of networks.

First you need to download the igraph library via the R packet manager:

```r
install.packages("igraph")
```

Then load the package at the begging of the script

```r
library(igraph)
```

### 1\) From adjacency matrix to Igraph

The function `undirectedIgraphFromAdjancyMatrix` defines an undirected graph in Igraph from an adjacency matrix.

```r
undirectedIgraphFromAdjancyMatrix <- function(input_mat)
{
  nbnode <- nrow(input_mat)
  graph  <- make_empty_graph(directed = FALSE)
  igraph <- add_vertices(graph,nbnode, color = "red")
  for(i in 1: nbnode)
  {
    for(j in i:nbnode)
    {
      if(input_mat[i,j] == 1)
      {
        igraph <- add_edges(igraph, c(i,j))
      }
    }
  }
  igraph
}
```

The function `directedIgraphFromAdjancyMatrix` do the same job but with a directed graph.

```r
directedIgraphFromAdjancyMatrix <- function(input_mat)
{
  nbnode <- nrow(input_mat)
  graph  <- make_empty_graph(directed = TRUE)
  igraph <- add_vertices(graph,nbnode, color = "red")
  for(i in 1: nbnode)
  {
    for(j in 1:nbnode)
    {
      if(input_mat[i,j] == 1)
      {
        igraph <- add_edges(igraph, c(i,j))
      }
    }
  }
  igraph
}
```

### 2\) Usual functions

#### Degrees

Using igraph you can have **degree** with

```r
igraphdegrees <- function(igraph)
{
  degree(igraph)
}
```

#### **Global clustering coefficient**

$$
C := \frac{number\ of\ triangles * 3}{number\ of\ connected\ triples}
$$

```r
global_clustering_coefficients <- function(igraph)    
{
  transitivity(igraph)
}
```

#### **Local clustering coefficient** for a given node

$$
C_{i} := \frac{number\ of \ pairs \ of\ neighbors\ of\ i\ that\ are\ connected}{number\ of\ pairs\ of\ neighbors\ of\ i}
$$

```r
local_clustering_coefficient <- function(igraph,node)
{
  transitivity( igraph, type = "local")[node] 
}
```

#### **Normalized closeness centrality**

```r
igraph_Closeness_centrality_node <- function(mattestbis,node)
{
  closeness(ig, mode="in",normalized = FALSE)[node]*nrow(mattestbis)
}
```

#### **Betweenness centrality**

```r
igraph_betweenness_centrality_node <- function(input_mat,node)
{
  igraph <- undirectedIgraphFromAdjancyMatrix(input_mat)
  betweenness(igraph,normalized = FALSE)[node]
}
```

### 3\) Plotting ability

It's possible the change the color of the node accordingly to their characteristic.

#### **Color**

To determine the color of a node we use a linear repartition of attribute with rgb color.

You can use a color of reference like the one used in the example: Esisar's purple color.

```r
inputcolor.R <- 146
inputcolor.G <- 39
inputcolor.B <- 143

findcolor <- function(input_val,low_val, high_val,inputcolor )
{
  range <- high_val - low_val;
  step.R <- floor(inputcolor.R * (input_val-low_val) / range) 
  step.G <- floor(inputcolor.G * (input_val-low_val) / range)
  step.B <- floor(inputcolor.B * (input_val-low_val) / range)
  step <- c(step.R/256,step.G/256,step.B/256)
}
```

The we just have to modify the vertex attribute.

```r
colorgraph <- function(igraph,inputcolor,inputtype)
{
  if(inputtype == "degree" )
  {
    degreevect <- igraphdegrees(igraph)
    for(i in 1:length(degreevect))
    {
      vcolor <- findcolor( degree(igraph)[i] , min(degree(igraph)), max(degree(igraph)) , inputcolor )
      igraph <- set.vertex.attribute(igraph, 'color', i, rgb(vcolor[1],vcolor[2],vcolor[3]))
    }
  }else if( inputtype == "closeness" )
  {
    closnessvect <- closeness(igraph,normalized = FALSE)
    for(i in 1 : length(closnessvect) )
    {
      vcolor <- findcolor( closnessvect[i] , min(closnessvect), max(closnessvect) , inputcolor )
      igraph <- set.vertex.attribute(igraph, 'color', i, rgb(vcolor[1],vcolor[2],vcolor[3]))
    }
  }else if( inputtype == "clustering" )
  {
    local_cluster <- transitivity(igraph, type = "local")
    for(i in 1 : length(local_cluster) )
    {
      if( is.nan(local_cluster[i]) == FALSE )
      {
        # /!\ Points were local clustering can't be etablish are not colored here
        vcolor <- findcolor( local_cluster[i] , min(local_cluster, na.rm = TRUE), max(local_cluster, na.rm = TRUE) , inputcolor )
        igraph <- set.vertex.attribute(igraph, 'color', i, rgb(vcolor[1],vcolor[2],vcolor[3]))
      }
    }
  }else if( inputtype == "betweenness" )
  {
    betw <- betweenness(igraph,normalized = FALSE)
    for(i in 1 : length(betw) )
    {
      if( is.nan(betw[i]) == FALSE )
      {
        vcolor <- findcolor( betw[i] , min(betw, na.rm = TRUE), max(betw, na.rm = TRUE) , inputcolor )
        igraph <- set.vertex.attribute(igraph, 'color', i, rgb(vcolor[1],vcolor[2],vcolor[3]))
      }
    }
  }
  igraph
}
```

#### **Size**

The process to modify the node size is almost the same as for color. We have replace the input color by a size coefficient to change the scale of plotting.

```r
findsize <- function(input_val,low_val, high_val, sizeCoef )
{
  range <- high_val - low_val;
  step  <-  ((input_val-low_val) / range)*sizeCoef
}
```

```r
sizegraph <- function(igraph,inputtype,sizeCoef)
{
  if(inputtype == "degree" )
  {
    degreevect <- igraphdegrees(igraph)
    for(i in 1:length(degreevect))
    {
      wsize <- findsize( degree(igraph)[i] , min(degree(igraph)), max(degree(igraph)) , sizeCoef )
      igraph <- set.vertex.attribute(igraph, 'size', i, wsize)
    }
  }else if( inputtype == "closeness" )
  {
    closnessvect <- closeness(igraph,normalized = FALSE)
    for(i in 1 : length(closnessvect) )
    {
      wsize <- findsize( closnessvect[i] , min(closnessvect), max(closnessvect) , sizeCoef )
      igraph <- set.vertex.attribute(igraph, 'size', i, wsize)
    }
  }else if( inputtype == "clustering" )
  {
    local_cluster <- transitivity(igraph, type = "local")
    for(i in 1 : length(local_cluster) )
    {
      if( is.nan(local_cluster[i]) == FALSE )
      {

        wsize <- findsize( local_cluster[i] , min(local_cluster, na.rm = TRUE), max(local_cluster, na.rm = TRUE) , sizeCoef )
        igraph <- set.vertex.attribute(igraph, 'size', i, wsize)
      }else
      {
        # /!\ Points were local clustering can't be etablish we set default size to 1
        igraph <- set.vertex.attribute(igraph, 'size', i, 1)
      }
    }
  }else if( inputtype == "betweenness" )
  {
    betw <- betweenness(igraph,normalized = FALSE)
    for(i in 1 : length(betw) )
    {
      if( is.nan(betw[i]) == FALSE )
      {
        wsize <- findsize( betw[i] , min(betw, na.rm = TRUE), max(betw, na.rm = TRUE) , sizeCoef )
        igraph <- set.vertex.attribute(igraph, 'size', i, wsize)
      }
    }
  }
  igraph
}
```

Those functions can be performed more efficiently.

