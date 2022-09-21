cluster.counts.long.proportion.1E6 <- function(cluster.counts.directory="./data_results/cluster_counts/",write.csv.to.file=T){

  cluster.count.files <- list.files(cluster.counts.directory,full.names = T,pattern = "cluster_counts_*.*csv")
  cell.types <- stringr::str_extract(cluster.count.files,"[A-Za-z]+[0-9]+\\+|[A-Za-z]+\\+|PBMC")
  if(any(table(cell.types)>1)){
    stop(paste("Cluster count files conflict: more than one instance of",names(which(table(cell.types)>1))))
  }
  cluster.counts <- sapply(cluster.count.files,read.csv,check.names = F)
  ##
  cluster.counts.list <- lapply(cluster.counts,function(cluster.counts.frame){
    cluster.cols <- grep("[0-9]+",colnames(cluster.counts.frame))
    dat.prop <- dat.1E6 <- cluster.counts.frame
    ##
    dat.1E6[,cluster.cols] <-dat.1E6[,cluster.cols]*(1E6/rowSums(dat.1E6[,cluster.cols]))
    ##
    dat.prop[,cluster.cols] <- prop.table(as.matrix(dat.prop[,cluster.cols]),1)
    ##
    dat.melt <- list(count = cluster.counts.frame,
                     proportion = dat.prop,
                     `1E6`= dat.1E6)
    dat.melt <- lapply(names(dat.melt),function(i){
      reshape2::melt(dat.melt[[i]],
                     measure.vars = grep("[0-9]+",colnames(dat.melt[[i]]),value = T),
                     value.name = i,
                     variable.name = 'cluster'
      )
    })
    return(suppressMessages(plyr::join_all(dat.melt)))
  })
  ##
  cluster.counts.long <- do.call(rbind,unname(cluster.counts.list))
  ##temp fix for missing batch date in older files; CD4+ and CD8+
  batch.dates <- stats::setNames(unique(cluster.counts.long$batch.date[cluster.counts.long$cell.type=="PBMC"]),
                                 nm=unique(cluster.counts.long$batch[cluster.counts.long$cell.type=="PBMC"]))
  batch.vec.cd4 <- cluster.counts.long$batch[cluster.counts.long$cell.type=="CD4+"]
  cluster.counts.long$batch.date[cluster.counts.long$cell.type=="CD4+"] <- batch.dates[as.character(batch.vec.cd4)]
  batch.vec.cd8 <- cluster.counts.long$batch[cluster.counts.long$cell.type=="CD8+"]
  cluster.counts.long$batch.date[cluster.counts.long$cell.type=="CD8+"] <- batch.dates[as.character(batch.vec.cd8)]
  ##
  if(write.csv.to.file){
    out.name <- paste("ECHO",sprintf("%03d",min(cluster.counts.long$batch)),"thru",sprintf("%03d",max(cluster.counts.long$batch)),sep = "_")
    out.name <- paste0(paste(out.name,"cluster_counts_long",Sys.Date(),sep = "_"),".csv")
    write.csv(cluster.counts.long,
              file = file.path(cluster.counts.directory,out.name),
              row.names = F)
  }else{
    return(cluster.counts.long)
  }
}
