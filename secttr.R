DIM<-dim(var_sect_sect_time_1)
for (i in 1:DIM[1]){
  for (j in 1:DIM[2]){
    for (k in 1:DIM[3]){
      for (l in 1:DIM[4]){
        V<-var_sect_sect_time_1[i,j,k,l]
        var_secttr_time_1[i,paste(j,k,sep="."),l]<-V
      }
    }
  }
}
