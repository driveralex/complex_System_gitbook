# Representations of graph and first metrics

### **1\) Matrix** $$\leftrightarrow$$ **List**

A graph may be represented by its adjacency matrix or its adjacency list. Below is a function that gets the adjacency list from the adjacency matrix :

```r
fromMatriceToList <- function(mat)
{
  k <- nrow(mat)
  lst <- list()
  for(i in 1:k)
  {
    lst[[i]] <- numeric()
  }
  for (i in 1:k)
  {
    for(j in i:k)
    {
      if( mat[i,j] == 1)
      {
        lst[[i]] <- c(lst[[i]],j)
        lst[[j]] <- c(lst[[j]],i)
      }
    }
  }
  lst
}
```

And now the opposite: a function that gets the adjacency matrix form the adjacency list :

```r
fromListToMatrice <- function(lst)
{
k <- length(lst)
mat <- vector(mode = "numeric",k*k)
dim(mat) <- c(k,k)
for( i in 1:k )
{
    l <- length(lst[[i]])
    for(j in 1:l)
    {
    mat[i,lst[[i]][j]] <- 1
    }
}
mat
}
```

### **2\) List and distribution of degrees**

Here is a function which return the list of degrees from the adjacency matrix or adjacency list :

```r
listOfDegrees <- function(input)
{
  if( is.vector(input) == TRUE )
  {
    k <- length(input)
    lst <- list()
    for( i in 1:k)
    {
      lst[[i]] <- length(input[[i]])
    }
  }
  else if( is.matrix(input) == TRUE )
  {
    templst <- fromMatriceToList(input)
    k <- length(templst)
    lst <- list()
    for( i in 1:k)
    {
      lst[[i]] <- length(templst[[i]])
    }
  } else
  {
    print("Wrong input of function listOfDegrees")
  }
  lst
}
```

Here is a function which return the degrees distribution from the adjacency list.

```r
degreeDistribution <- function(input)
{
  lstdeg <- listOfDegrees(input)
  print(lstdeg[[which.max(lstdeg)]])

  lstdistrib <- list()
  lstdistrib <- replicate(lstdeg[[which.max(lstdeg)]],0)
  k <- length(lstdeg)
  for(i in 1:k)
  {
    lstdistrib[[lstdeg[[i]]]] <- lstdistrib[[lstdeg[[i]]]] + 1
  }
  print(lstdistrib)
  print("Endfunction")
}
```

### **3\) Clustering coefficient**

Here is a function that return the list of cluster coefficients from the adjacency matrix.

```r
clustering_coef <- function(input)
{
  len <- nrow(input)
  lst <- list()
  for(i in 1:len)
  {
    output <- 0
    for(k in 1:len)
    {
      for(p in 1 :len )
      {
        output <- output + input[i,k]*input[k,p]*input[p,i]
      }
    }
    lst[[i]] <- output
  }
  lst
}
```

