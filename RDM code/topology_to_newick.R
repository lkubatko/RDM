## please check if the last row in the freq matrix is the outgroup

multiphylo<-function(base_tree,label){

require(Matrix)
require(phangorn)
  
outgroup<-label[length(label)]

for (i in 1:dim(base_tree)[1]){
  ind3<-which(base_tree[i,]==0)
  ind4<-which(base_tree[i,]==1)
  ind5<-which(base_tree[i,]==2)
  if (length(ind3)==0){
    grp1<-label[ind4]
    if (length(ind4)>1){
      grp1<-paste("(",paste(grp1,collapse = ","),")",sep = "")
    }
    grp2<-label[ind5]
    if (length(ind5)>1){
      grp2<-paste("(",paste(grp2,collapse = ","),")",sep = "")
    }
    grp<-paste("(",outgroup,",",grp1,",",grp2,");",sep="")
  }else{
    grp1<-label[ind4]
    if (length(ind4)>1){
      grp1<-paste("(",paste(grp1,collapse = ","),")",sep = "")
    }
    grp2<-label[ind5]
    if (length(ind5)>1){
      grp2<-paste("(",paste(grp2,collapse = ","),")",sep = "")
    }
    grp3<-label[ind3]
    if (length(ind3)>1){
      grp3<-paste("(",paste(grp3,collapse = ","),")",sep = "")
    }
    grp<-paste("(",outgroup,",",grp1,",",grp2,",",grp3,");",sep="")
  }
  tree<-read.tree(text = grp)
  mp<-list(tree)
  class(mp)<- "multiPhylo"
  if (i==1) mp2<-mp else mp2<-c(mp2,mp)
}
return(mp2)
}

topology_to_newick<-function(new_split_set_indicator,label){
  mp2<-multiphylo(new_split_set_indicator,label=label)
  tre1<-superTree(mp2, rooted = TRUE)
  return(tre1)
}
  

  



