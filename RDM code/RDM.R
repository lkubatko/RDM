RDM<-function(trans_mat_allele_freq,outgroup,
              use=c("complete.obs","pairwise.complete.obs","everything","all.obs","na.or.complete")){
  
  use<-match.arg(use)
  names<-row.names(trans_mat_allele_freq)
  if (is.character(names)==0){
    print("Please use the population names as the row names of your transformed allele frequency matrix")
  }else if (is.character(outgroup)){
    index<-which(names==outgroup)
  }else {
    index<-outgroup
  }
  trans_mat_allele_freq<-rbind(trans_mat_allele_freq[-index, ],trans_mat_allele_freq[index,])
  label<-row.names(trans_mat_allele_freq)
  array_zero_ID<-BMaf_to_zeroRDD(trans_mat_allele_freq,use=use)
  base_tree<-zeroRDD_to_splits(array_zero_ID)
  rd_tre<-topology_to_newick(base_tree,label=label)
  return(rd_tre)
}
