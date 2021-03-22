# Breadth-first search algorithm and applications

### 1\) Breadth-first search algorithms

Breadth-first search algorithm give the distance between a given node to all the other node belonging to the same component.

#### 1. _Naive_ implementation

The first implementation we use is _"naive"_ and use the following principles:

Let $$n$$ = number of nodes.

Let $$m$$ = number of edges.

 Let $$i$$ = studied node.

 Let $$s$$ = destination node.

 Let $$d$$ = number of rounds.

**1/** Creation of an array $$ D_{s}$$ of $$n$$ _****_ ****integers to store the distances$$ d(i,s)$$.

**2/** Initialization of __$$ D_{s}$$ : $$ \forall \ i \neq s \  ,  \ D_{s}(i)=-1 \$$ else $$ D_{s}(i) = 0$$.

**3/** Find all nodes with distance d. If there is no, then stop. 

**4/** Find the neighbors of these nodes, assign those neighbors which don’t have a distance yet, with the distance $$d + 1$$. 

**5/** Set $$d = d + 1$$ and go to **3/**.

Below is a R implementation:

```r
breadth_first <- function(input_mat,node)
{
  k <- nrow(input_mat)
  outvect <- vector(mode = "numeric",k)
  for( i in 1:k )
  {
    if( i == node )
    {
      outvect[i] = 0
    }
    else{
      outvect[i] = -1
    }
  }
  end <- FALSE
  d <- 0
  while( end == FALSE  )
  {
    positivecondition <- list()
    for(i in 1:k)
    {
      if( outvect[i]== d)
      {
        positivecondition[length(positivecondition)+1] = i
      }
    }
    if(length(positivecondition) == 0 )
    {
      end <- TRUE
    }else{
      for(i in 1:length(positivecondition))
      {
        for(j in 1:k)
        {
          if( (areneighbours(input_mat,positivecondition[[i]],j) == TRUE) && (outvect[j] == -1) )
          {
            outvect[j] <- d + 1
          }
        }
      }
    }
    d <- d + 1
  }
  outvect
}
```

With the details of the `areneighbours`function bellow.

```r
areneighbours <- function(input,i,j)
{
  if( is.vector(input) == TRUE )
  {
    matinput <- fromListToMatrice(input)
    if(matinput[i,j] == 1 )
    {
      output <- TRUE
    }else
    {
      output <- FALSE
    }
  }
  else if( is.matrix(input) == TRUE )
  {
    if(input[i,j] == 1 )
    {
      output <- TRUE
    }else
    {
      output <- FALSE
    }
  }
  output
}
```

The way we implement this could be better executed \(notably the with the `areneighbours`\) but with this algorithm, for a typical network, complexity is $$ O(m + n *log(n))$$.

#### 2. _**Stack**_ implementation

A better implementation could be done with using a queue. Here is a implementation with the stack algorithm:

```r
breadth_first_stack <- function(input_mat,node)
{
  aplist <- fromMatriceToList(input_mat)
  stack <- vector(mode = "numeric")
  stack[1] <- node
  read <- 1
  write <- 2
  k <- nrow(input_mat)
  ds <- vector(mode = "numeric",k)
  for( i in 1:k )
  {
    if( i == node )
    {
      ds[i] = 0
    }
    else{
      ds[i] = -1
    }
  }
  while(read != write)
  {
    for(i in 1:length(aplist[[stack[read]]]))
    {
      if(ds[ aplist[[stack[read]]][i] ] == -1 )
      {
        ds[ aplist[[stack[read]]][i] ] <- ds[ stack[read]] + 1
        stack[write] <- aplist[[stack[read]]][i]
        write <- write +1
      }
    }
    read <- read + 1 
  }
  ds
}
```

With the stack approach the complexity is smaller: $$ O(m + n)$$ 

For this tow approach we have the result of a unique node. To have the complete matrices of distance as required in Q1, we repeat the operation for all nodes and aggregate the outputs. This is the purpose of the `mat_D_breadth_first_stack` function using the stack algorithm.

```r
mat_D_breadth_first_stack <- function(input_mat)
{
  k <- nrow(input_mat)
  outputmat <- vector(mode = "numeric")

  for(i in 1:k)
  {
    outputmat <- c(outputmat,breadth_first_stack(input_mat,i));
  }
  dim(outputmat) <- c(k,k)
  outputmat
}
```

### 2\) Diameter 

Diameter is a record of the largest distance observed in a component. Once the Breadth-first search executed, finding the diameter for a given node is straightforward.

```r
diammeter <- function(input_mat,node)
{
  k <- nrow(input_mat)
  bfs_node <- breadth_first_stack(input_mat,node);
  diameter <- max(bfs_node)
  lst_compenent_memeber <- list()
  for(i in 1:k)
  {
    if(bfs_node[i] > 0 )
    {
      lst_compenent_memeber[[length(lst_compenent_memeber)+1]] <- i
    }
  }
  l <- length(lst_compenent_memeber)
  for(j in 1:l)
  {
    if( max(breadth_first_stack(input_mat,lst_compenent_memeber[[j]])) > diameter )
    {
      diameter <- max(breadth_first_stack(input_mat,lst_compenent_memeber[[j]]))
    }
  }
  diameter
}
```

### 3\) Closeness centrality 

First, we have to underline that if the network is not composed of a single unique component, the result of the closeness centrality must be taken with care. It is possible to handle each component individually, but this can bias the values. Indeed, nodes in smaller component may have higher value. This is what the following example do:

```r
closeness_centrality_node <- function( input_mat , node ) 
{
  dist <- breadth_first_stack_list(input_mat,node)
  n <- length(dist)
  sum <- 0
  for(i in 1:n)
  {
    if( dist[[i]] > 0 ) #Here nodes that are not in the component are not take in account.faire la remarque page 47 (A modfifier ?)
    {
      sum <- sum + dist[[i]]
    }
  }
  output <- length(dist)/(sum)
  output
}
```

We use the very inelegant \(but functional\) breadth\_first\_stack\_list function described below.

```r
breadth_first_stack_list <- function(input_mat,node)
{
  output_vect <- breadth_first_stack(input_mat,node)
  k <- length(output_vect)
  output_list <- list()
  for(i in 1:k)
  {
    if(output_vect[i] >= 0)
    {
      output_list[[length(output_list)+1]] <-  output_vect[i]
    }
  }
  output_list
}
```

Note that the igraph package throw a warning when the components are not all linked together.

### 4\) Betweennes centrality

Implementing the Betweennes centrality was the hardest task to complete. Finding all the shortest path is the blocking point. I tried to do this job with a recusive function knowing the lenght of every shortest path with a previous breadfirst search. Unfortunately I did not maneged to make this function work. ~~The following programm visit all the nodes, but returning the target point and aggregate the coresting path is not so easy.~~

