## This is a code for computing Brownina Motion allele frequencies zero-valued RD differences of a tree-topology

## trans_mat_allele_freq is a (P+1) X L matrix, containing transformed allele frequencies of P taxa and an outgroup for L loci

BMaf_to_zeroRDD<-function(trans_mat_allele_freq,
                          use=c("complete.obs","pairwise.complete.obs","everything","all.obs","na.or.complete")){
  
  P <- nrow(trans_mat_allele_freq) - 1
  
  H <- P*(P-1)*(P-2)/6
  ## H is the smallest number of  RD differences that have to be zero in order for it
  ## to agree to a tree-topology
  
  use<-match.arg(use)
  mat_normalized <- t(trans_mat_allele_freq)
  mat_sigma <- cov(mat_normalized,use = use)
  for (x in 1:(P+1))
    for (y in 1:(P+1))
      if (mat_sigma[x,y] < 0)
        mat_sigma[x,y] <- 0
  
  #######################################################################
  
  raw_array_ID <- array(dim=c(P,P,P),0)
  
  n_list_ID <- P*(P-1)*(P-2)/2
  
  list_ID <- matrix(nrow=n_list_ID,ncol=4)
  ## list_ID lists the relevant absolute IDs from array_ID
  
  array_zero_ID <- array(dim=dim(raw_array_ID),1)
  
  i_list_ID <- 0
  
  
  for (i in 1:P)
    for (j in 1:P)
      for (k in 1:P)
      {
        raw_array_ID[i,j,k] <- abs(mat_sigma[i,j] - mat_sigma[i,k])
        
        if ( ((i != j) && (i != k)) && (j < k))
        {
          i_list_ID <- i_list_ID + 1
          
          list_ID[i_list_ID,] <- c(raw_array_ID[i,j,k],i,j,k)      
        }
        
        
        if ( ((i < j) && (i < k)) && (j < k))
        {
          
          maxrd <- which.max(c(mat_sigma[i,j],mat_sigma[i,k],mat_sigma[j,k]))
          
          if (maxrd == 1)
          {
            array_zero_ID[k,i,j] <- 0
            array_zero_ID[k,j,i] <- 0
            
          }
          if (maxrd == 2)
          {
            array_zero_ID[j,i,k] <- 0
            array_zero_ID[j,k,i] <- 0
            
          }
          if (maxrd == 3)
          {
            array_zero_ID[i,j,k] <- 0
            array_zero_ID[i,k,j] <- 0
            
          }
          
        }
        
      }
  
  ## array_zero_ID
  ## This is an array where zero valued RD differences have a value zero, and the rest have a value of 1
  ## Note that the array is P X P X P, and for each i, array_zero_ID[i,,] is symmetric
  
  return(array_zero_ID)
  
  }





