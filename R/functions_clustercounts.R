cluster.counts.long.from.path <- function(cluster.counts.path){
  cluster.counts.frame <- utils::read.csv(cluster.counts.path,check.names = F,stringsAsFactors = T)
  if('batch'%in%colnames(cluster.counts.frame)){
    cluster.counts.frame$batch <- factor(cluster.counts.frame$batch)
  }
  cluster.cols <- grep("[0-9]+",colnames(cluster.counts.frame))
  ##
  dat.prop <- dat.per1million <- cluster.counts.frame
  ##
  dat.per1million[,cluster.cols] <-dat.per1million[,cluster.cols]*(1E6/rowSums(dat.per1million[,cluster.cols]))
  ##
  dat.prop[,cluster.cols] <- prop.table(as.matrix(dat.prop[,cluster.cols]),1)*100
  ##
  dat.melt <- list(count = cluster.counts.frame,
                   proportion = dat.prop,
                   per1million= dat.1E6)
  dat.melt <- lapply(names(dat.melt),function(i){
    reshape2::melt(dat.melt[[i]],
                   measure.vars = grep("[0-9]+",colnames(dat.melt[[i]]),value = T),
                   value.name = i,
                   variable.name = 'cluster'
    )
  })
  return(suppressMessages(plyr::join_all(dat.melt)))
}

cluster.counts.long.proportion.1E6 <- function(cluster.counts.directory="./data_results/cluster_counts/",write.csv.to.file=T){

  cluster.count.files <- list.files("./data_results/cluster_counts/",full.names = T,pattern = "[A-Za-z]+[0-9]+\\+|[A-Za-z]+\\+|PBMC*.*csv")
  cell.types <- stringr::str_extract(cluster.count.files,"[A-Za-z]+[0-9]+\\+|[A-Za-z]+\\+|PBMC")
  if(any(table(cell.types)>1)){
    stop(paste("Cluster count files conflict: more than one instance of",names(which(table(cell.types)>1))))
  }
  cluster.counts <- sapply(cluster.count.files,utils::read.csv,check.names = F)
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
    utils::write.csv(cluster.counts.long,
                     file = file.path(cluster.counts.directory,out.name),
                     row.names = F)
  }else{
    return(cluster.counts.long)
  }
}

cluster.counts.long.generate.frames <- function(cluster.counts.long.filepath){
  #value column: count = raw cluster counts per-sample/per cell.type
  #value column: proportion = count converted to proportion (sum to 1) of cluster per-sample/per cell.type
  #value column: 1E6 = raw cluster counts per-sample/per cell.type normalized to count per-million
  clusters.long <- utils::read.csv(cluster.counts.long.filepath,check.names = F,
                                   colClasses = list(
                                     "subject"="factor",
                                     "condition"="factor",
                                     "cell.type"="factor",
                                     "batch.date"="factor"
                                   )
  )
  clusters.long$visit <- factor(clusters.long$visit,levels = c(paste0("V",c(4,6,7)),"Adult"))
  clusters.long$batch <- factor(clusters.long$batch,levels=seq(max(clusters.long$batch)))
  clusters.long$cluster <- factor(clusters.long$cluster,levels=seq(max(clusters.long$cluster)))
  ##
  clusters.long.split <- split(clusters.long,clusters.long$cell.type)
  clusters.long.split <- lapply(clusters.long.split,droplevels)
  ##
  clusters.wide.split <- lapply(clusters.long.split,function(i){
    left.side <- names(which(!sapply(i,is.numeric)));left.side <- left.side[!left.side %in% 'cluster']
    left.side <- paste(left.side,collapse = "+")

    value.vars <- names(which(sapply(i,is.numeric)))

    right.side <- 'cluster'

    dcast.formula <- stats::as.formula(paste(left.side,"~",right.side))

    sapply(value.vars,function(v){reshape2::dcast(i,formula = dcast.formula,value.var = v)},simplify = F)

    #dat <- reshape2::dcast(i,formula = dcast.formula,value.var = "1E6")
  })
  return(list(split = clusters.long.split,
              wide = clusters.wide.split)
  )
}
