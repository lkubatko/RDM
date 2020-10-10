## This is a code for estimating splits of a tree-topology with zero-valued RD differences

check.equivalence<-function(array_zero_ID){

P<-dim(array_zero_ID)[1]
increment <- 2*P
max_split <- increment
split_set_indicator <- matrix(nrow=max_split,ncol=P)
n_split_set <- 0
# This will be incremented as we find new splits
for (i in 1:P)
 for (j in 1:P)
  for (k in 1:P)
{

   if ( ((i != j) && (i != k)) && ((j < k) && (array_zero_ID[i,j,k] == 0)))   
 {

i_split <- 0
repeat
  {

## i_split is counting the number of already-identified splits compared 
i_split <- i_split + 1

if (i_split > n_split_set) 
   {
    n_split_set <- i_split

    if (n_split_set > max_split)
     {
       max_split <- max_split + increment
       new_array <- matrix(nrow=max_split,ncol=P)
       new_array[(1:(max_split-increment)),] <- split_set_indicator
       split_set_indicator <- new_array
       rm(new_array)
     }

    this_split_set_indicator <- integer(P)
    this_split_set_indicator[i] <- 1
    this_split_set_indicator[j] <- 2
    this_split_set_indicator[k] <- 2   

    split_set_indicator[n_split_set,] <- this_split_set_indicator

    break
   }
## that is, if no other equivalent split is found, define (i,j,k) as a new split; then get out  

  if ( (split_set_indicator[i_split,i] > 0) && (split_set_indicator[i_split,k]==0) && (( (split_set_indicator[i_split,j] > 0)  &&  (split_set_indicator[i_split,j] !=  split_set_indicator[i_split,i]) ) ) )
   {
    split_set_indicator[i_split,k] <- split_set_indicator[i_split,j]

    break
   }
  else if ( (split_set_indicator[i_split,i] > 0) && (split_set_indicator[i_split,j]==0) && ( ( (split_set_indicator[i_split,k] > 0)  &&  (split_set_indicator[i_split,k] !=  split_set_indicator[i_split,i]) ) ) )
   {
    split_set_indicator[i_split,j] <- split_set_indicator[i_split,k]

    break
   }

  } 
 }

}

split_set_indicator <- matrix(ncol=P,split_set_indicator[(1:n_split_set),])

return(split_set_indicator)

}


## next we need to eliminate the repetitions

eliminate.rep<-function(split_set_indicator){
  
  P<-ncol(split_set_indicator)
  n_split_set<-nrow(split_set_indicator)
  increment <- P-1
  max_split <- increment
  next_split_set_indicator <- matrix(nrow=max_split,ncol=P)
  next_n_split_set <- 0
  merge_record <- NULL
  
  for (i_split in 1:n_split_set)
  {
    
    j_split <- 0
    
    repeat
    {
      j_split <- j_split + 1
      
      if (j_split > next_n_split_set)
      {
        
        
        next_n_split_set <- next_n_split_set + 1
        
        
        if (next_n_split_set > max_split)
        {
          max_split <- max_split + increment
          next_array <- matrix(nrow=max_split,ncol=P)
          next_array[(1:(max_split-increment)),] <- next_split_set_indicator
          next_split_set_indicator <- next_array
          rm(next_array)
        }
        
        
        
        next_split_set_indicator[next_n_split_set,] <- split_set_indicator[i_split,]
        
        break
      }
      
      if (  ((length( intersect( which(next_split_set_indicator[j_split,] == 1),which(split_set_indicator[i_split,] == 1) ) ) > 0) && (length( intersect( which(next_split_set_indicator[j_split,] == 2),which(split_set_indicator[i_split,] == 2) ) ) > 0) )
            &&    ((length( intersect( which(next_split_set_indicator[j_split,] == 1),which(split_set_indicator[i_split,] == 2) ) ) == 0) && (length( intersect( which(next_split_set_indicator[j_split,] == 2),which(split_set_indicator[i_split,] == 1) ) ) == 0) )  )
      {
        
        
        next_split_set_indicator[j_split,(union( which(next_split_set_indicator[j_split,] == 1),which(split_set_indicator[i_split,] == 1) ))] <- 1
        next_split_set_indicator[j_split,(union( which(next_split_set_indicator[j_split,] == 2),which(split_set_indicator[i_split,] == 2) ))] <- 2
        
        
        break
      }
      else if (  ((length( intersect( which(next_split_set_indicator[j_split,] == 1),which(split_set_indicator[i_split,] == 2) ) ) > 0) && (length( intersect( which(next_split_set_indicator[j_split,] == 2),which(split_set_indicator[i_split,] == 1) ) ) > 0) )
                 &&    ((length( intersect( which(next_split_set_indicator[j_split,] == 1),which(split_set_indicator[i_split,] == 1) ) ) == 0) && (length( intersect( which(next_split_set_indicator[j_split,] == 2),which(split_set_indicator[i_split,] == 2) ) ) == 0) )  )
      {
        
        
        next_split_set_indicator[j_split,(union( which(next_split_set_indicator[j_split,] == 1),which(split_set_indicator[i_split,] == 2) ))] <- 1
        next_split_set_indicator[j_split,(union( which(next_split_set_indicator[j_split,] == 2),which(split_set_indicator[i_split,] == 1) ))] <- 2
        
        
        break
      }
    }
  }
  
  new_split_set_indicator <- matrix(nrow=max_split,ncol=P)
  
  
  new_n_split_set <- 0
  
  
  for (i_split in 1:next_n_split_set)
  {
    
    j_split <- 0
    
    repeat
    {
      j_split <- j_split + 1
      
      if (j_split > new_n_split_set)
      {
        
        
        new_n_split_set <- new_n_split_set + 1
        
        
        if (new_n_split_set > max_split)
        {
          max_split <- max_split + increment
          new_array <- matrix(nrow=max_split,ncol=P)
          new_array[(1:(max_split-increment)),] <- new_split_set_indicator
          new_split_set_indicator <- new_array
          rm(new_array)
        }
        
        
        
        new_split_set_indicator[new_n_split_set,] <- next_split_set_indicator[i_split,]
        
        break
      }
      
      if (setequal(which(new_split_set_indicator[j_split,] == 1),which(next_split_set_indicator[i_split,] == 1)))
      {
        new_split_set_indicator[j_split,(union( which(new_split_set_indicator[j_split,] == 2),which(next_split_set_indicator[i_split,] == 2) ))] <- 2
        
        break
      }
      else if (setequal(which(new_split_set_indicator[j_split,] == 2),which(next_split_set_indicator[i_split,] == 2)))
      {
        new_split_set_indicator[j_split,(union( which(new_split_set_indicator[j_split,] == 1),which(next_split_set_indicator[i_split,] == 1) ))] <- 1
        
        break
      }
      else if (setequal(which(new_split_set_indicator[j_split,] == 1),which(next_split_set_indicator[i_split,] == 2)))
      {
        new_split_set_indicator[j_split,(union( which(new_split_set_indicator[j_split,] == 2),which(next_split_set_indicator[i_split,] == 1) ))] <- 2
        
        break
      }
      else if (setequal(which(new_split_set_indicator[j_split,] == 2),which(next_split_set_indicator[i_split,] == 1)))
      {
        new_split_set_indicator[j_split,(union( which(new_split_set_indicator[j_split,] == 1),which(next_split_set_indicator[i_split,] == 2) ))] <- 1
        
        break
      }
      
    }
  }
  
  return(new_split_set_indicator)
}
## done with elimination


zeroRDD_to_splits<-function(array_zero_ID){
  split_set_indicator<-check.equivalence(array_zero_ID)
  new_split_set_indicator<-eliminate.rep(split_set_indicator)
  return(new_split_set_indicator)
}