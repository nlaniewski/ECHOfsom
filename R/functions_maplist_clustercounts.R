
counts.from.maplist <- function(maplist.rds.path,count.type=c('cluster','node')){
  #read in map list
  map.list <- readRDS(maplist.rds.path)
  ##list contains both node and cluster indices; retrieve based on 'count.type' argument
  count.list <- lapply(map.list,'[[',count.type)
  ##table/counts
  count.table <- lapply(count.list,table)
  ##rbind;data.frame
  counts <- as.data.frame(do.call(rbind,count.table))
  ##
  return(counts)
}

counts.add.mdat <- function(counts.frame){
  if(is.null(rownames(counts.frame))){
    stop("Need rownames as filepaths")
  }
  ##cluster names
  c.names <- colnames(counts.frame)
  ##sample/unique id column;for merging
  counts.frame <- cbind(sample = gsub("_CONCATENATED|_MadMedian.fcs|_HistCut|.fcs", "",
                                      basename(rownames(counts.frame))),
                        counts.frame)
  ##add a conditional in order to auto-generate a 'PBMC' label
  if(!any(grepl("+",counts.frame$sample,fixed = T))){counts.frame$sample <- paste0(counts.frame$sample,"_PBMC")}
  ##generate meta-data from row names (file paths); merge; drop rownames (file paths)
  mdat <- ECHOfcs::mdat.frame.from.paths.echo(rownames(counts.frame))
  counts.frame <- plyr::join(counts.frame, mdat, by = "sample")
  rownames(counts.frame) <- NULL
  ##logical test for cluster total/.fcs event total match
  fcs.total.cluster.total.mismatch <- which(rowSums(counts.frame[c.names])!=counts.frame$total.events)
  sample.mismatch <- counts.frame$sample[fcs.total.cluster.total.mismatch]
  message(paste("The following samples have a mismatch between cluster totals and .fcs '$TOT' value:",
                paste(sample.mismatch, collapse = " "),
                sep = "\n")
  )
  ##fix totals; cluster total takes priority over .fcs total (most likely due to 'name-fix')
  counts.frame[fcs.total.cluster.total.mismatch,'total.events'] <- rowSums(counts.frame[fcs.total.cluster.total.mismatch,c.names])
  ##logical test
  if(all(rowSums(counts.frame[c.names])==counts.frame$total.events)){
    return(counts.frame)
  }else{
    message("cluster total (rowSums) and 'total.events' mismatch")
  }
}

counts.colnames.ordering <- function(counts.frame){
  count.names <- colnames(counts.frame)[!is.na(suppressWarnings(as.numeric(colnames(counts.frame))))]
  c.names <- c(
    "global.id",
    "sample",
    "subject",
    "visit",
    "condition",
    "cell.type",
    "batch",
    "batch.date"
  )
  c.names <- c(c.names,count.names)
  counts.frame <- counts.frame[,colnames(counts.frame) %in% c.names]
  counts.frame <- counts.frame[,c.names]
  return(counts.frame)
}

counts.out.filename <- function(counts.frame,count.type){
  out.filename <- paste0(paste("ECHO",
                               min(levels(counts.frame$batch)),
                               "thru",
                               max(levels(counts.frame$batch)),
                               unique(counts.frame$cell.type),
                               paste(count.type,"counts",sep = "_"),
                               Sys.Date(),
                               sep = "_"),
                         ".csv")
  return(out.filename)
}

maplist.to.counts <- function(maplist.rds.path,count.type=c('cluster','node'),write.counts=TRUE){
  counts <- counts.from.maplist(maplist.rds.path,count.type = count.type)
  counts <- counts.add.mdat(counts)
  counts <- counts.colnames.ordering(counts)
  if(write.counts){
    out.filename <- counts.out.filename(counts,count.type = count.type)
    f.path <- file.path("./data_results/cluster_counts/")
    if(count.type=='node'){
      f.path <- file.path(f.path,"node_counts")
    }
    message("Writing counts .csv to './data_results/cluster_counts/'")
    utils::write.csv(counts,
                     file.path(f.path,out.filename),
                     row.names = F)
  }else{
    return(counts)
  }
}
