### Define libraries, parameters and functions

library("MASS")
library("ape")
library("abind")
library("stats")



P <- 40
## 2000 sites 
compute_L_by_multiplying_to_P_square <- 1.25
## 4000 sites
#compute_L_by_multiplying_to_P_square <- 2.5
## 6000 sites
#compute_L_by_multiplying_to_P_square <- 3.75
## 8000 sites
#compute_L_by_multiplying_to_P_square <- 5
L <- ceiling(P*P*compute_L_by_multiplying_to_P_square)
outgroup_root_br_length_scale_to_max <- 1
outgroup_extr_br_length_scale_to_max <- 10
Mu_0_element <- -1.49
#fixed_br_length <- 0.05
fixed_br_length <- 0.02
#fixed_br_length <- 0.01


N_Iterations <- 100
version_number <- 10.1
### SimID is a unique identifier for each simulation and has the format P40_N25_V2
SimID <- paste("P",P,"_N",N_Iterations,"_V",version_number,sep="")
Foldername <- paste("Folder_",SimID,sep="")
dir.create(Foldername)


### Done defining

### Do simulations N_Iterations times

for (i_it in 1:N_Iterations)
{
  ### Generate the tree
  
  simtree <- rtree(P,br=fixed_br_length)
  #simtree <-  rtree(n=P,br=runif,min=0.01,max=0.1)
  
  new_label<-gsub("t","",simtree$tip.label)
  new1<-as.character(as.numeric(new_label)+1)
  simtree$tip.label<- paste("species",new1,sep='')
  
  # max_edge_length <- max(simtree$edge.length)
  max_edge_length <- 0.02
  
  outgroup_root_br_length <- outgroup_root_br_length_scale_to_max*max_edge_length
  outgroup_extr_br_length <- outgroup_extr_br_length_scale_to_max*max_edge_length
  outgroup_info <- c(outgroup_root_br_length,outgroup_extr_br_length)
  
  
  Mu_0 <- rep(Mu_0_element,(P+1))
  Sigma_0=matrix(nrow=(P+1),ncol=(P+1),0)
  dist_nodes_simtree <- dist.nodes(simtree)
  mrca_simtree <- mrca(simtree)
  
  for (i in 1:P){
    for (j in 1:P){
      Sigma_0[i,j] <- dist_nodes_simtree[mrca_simtree[i,j],(P+1)]
    }
  }
  
  # row (P+1) represents the outgroup
  Sigma_0[(1:P),(1:P)] <- Sigma_0[(1:P),(1:P)] + outgroup_root_br_length
  Sigma_0[(P+1),(P+1)] <- outgroup_extr_br_length  
  
  
  ### Simulate info for L loci for the tree
  total_sampled_per_taxon<- 100   ##   # of individuals sampled for each population
  
  logit_mat_allele_freq <- (t(mvrnorm(n=L,mu=Mu_0,Sigma=Sigma_0)))
  mat_allele_count <- round(total_sampled_per_taxon*(exp(logit_mat_allele_freq)/(exp(logit_mat_allele_freq)+1)))
  mat_allele_freq <- mat_allele_count/total_sampled_per_taxon
  label<- c(simtree$tip.label,"species1")
  row.names(mat_allele_freq)<- label
  
  Simr_filename <- paste(Foldername,"/","Sim_",SimID,"_Iteration_",i_it,".txt",sep="")
  Simr_tree <- paste(Foldername,"/","Sim_",SimID,"_tree_",i_it,".phy",sep="")
  
  write.tree(simtree,file=Simr_tree) 
  write.table(mat_allele_freq, file=Simr_filename)
  
}

##################################### Analysis

dir.create(paste(Foldername,"/RDM",sep = ""))
for (iter in 1:N_Iterations){
  
  Simr_filename <- paste(Foldername,"/","Sim_",SimID,"_Iteration_",iter,".txt",sep="")
  mat_allele_freq<- read.table(file=Simr_filename)
  mat_allele_freq[mat_allele_freq==1]<-0.99
  mat_allele_freq[mat_allele_freq==0]<-0.01
  trans_mat_allele_freq<-log(mat_allele_freq/(1-mat_allele_freq))
  
  rd_tre<-RDM(trans_mat_allele_freq,outgroup = 'species1')
  name_tre<-paste(Foldername,"/RDM/RDouttree",iter,sep='')
  write.tree(rd_tre,file=name_tre)
}




    
 

    
    



