#' Build an fsom object using a list of flowFrames
#'
#' @param fcs.list a list of flowFrames (named); returned from 'ECHOfcs::read.fcs.selected.markers'
#'
#' @return fsom object; data columns renamed and transformed; meta-data added
#' @export
#'
fsom.initialize <- function(fcs.list){
  if(length(unique(lapply(fcs.list,class)))==1&unique(lapply(fcs.list,class))=="flowFrame"){
    #extract markernames
    fcs.list.markernames <- sapply(strsplit(ECHOfcs::get.metal.markers(fcs.list),"_"),'[',2)
    #extract expression matrix from each flowFrame
    fcs.list <- lapply(fcs.list, flowCore::exprs)
    #initialize fsom list
    fsom <- list()
    class(fsom) <- "FlowSOM"
    #row bind the list of matrices; data
    fsom$data <- do.call(rbind,fcs.list)
    #transform the data
    message("transforming data...")
    fsom$data <- apply(fsom$data,2,function(x) asinh(x)/5)
    #
    gc()
    #store markernames
    fsom$markers <- fcs.list.markernames
    #test for name match; rename data column names
    if(all(colnames(fsom$data)==names(fsom$markers))){
      colnames(fsom$data) <- fsom$markers
    }else{
      message("data column names and names(markernames) mismatch")
    }
    #store seed value
    fsom$seed <- 20040501
    #store mdat
    fsom$mdat <- ECHOfcs::mdat.frame.from.paths.echo(names(fcs.list))
    if(!all(as.numeric(sapply(fcs.list,nrow))==fsom$mdat$total.events)){
      message("Updating 'total.events' in 'mdat'")
      fsom$mdat$total.events <- as.numeric(sapply(fcs.list,nrow))
    }
    #return the fsom
    return(fsom)
  }else{
    stop("Need a list of flowFrames")
  }
}

#' FlowSOM::metaClustering_consensus modification; generates cluster factor
#'
#' @param dat.mat numeric matrix; fsom$map$codes
#' @param max.clusters final number of clusters; factored
#'
#' @return numeric factor; length = nrow/# of fsom$map$codes, levels = seq(max.clusters)
#' @export
#'
generate.clusters <- function(dat.mat, max.clusters){
  ccp.res <- ConsensusClusterPlus::ConsensusClusterPlus(d = t(dat.mat),
                                                        maxK = max.clusters,
                                                        reps = 200,
                                                        pItem = 1,
                                                        pFeature = 1,
                                                        title = tempdir(),
                                                        plot = "png",
                                                        clusterAlg = "hc",
                                                        distance = "euclidean",
                                                        seed = 4242)
  ccp.res <- as.factor(ccp.res[[max.clusters]]$consensusClass)
  return(ccp.res)
}

generate.cluster.medians <- function(fsom){
  cluster.medians <- stats::aggregate(. ~ cluster,
                                      data = data.frame(cluster = as.numeric(fsom$metaclustering),fsom$map$codes,check.names = F),
                                      FUN = stats::median
  )
  return(cluster.medians)
}

echo.fsom.colnames <- function(c.names){
  ##prepare column names
  c.names <- unname(c.names)
  c.names <- c.names[order(c.names)]
  c.names[grep("CD",c.names)] <- grep("CD",c.names,value = T)[order(as.numeric(stringr::str_extract(grep("CD",c.names,value = T),"[0-9]+")))]
  c.names[grep("IL",c.names)] <- grep("IL",c.names,value = T)[order(as.numeric(stringr::str_extract(grep("IL",c.names,value = T),"[0-9]+")))]
  c.names <- c(c.names,'cluster','node')
  return(c.names)
}

fsom.cluster.counts <- function(fsom=fsom){
  if(is.null(fsom$mdat)&is.null(fsom$map$mapping)&is.null(fsom$metaclustering)){
    stop("Missing list features")
  }
  dat <- data.frame(sample = rep(fsom$mdat$sample,fsom$mdat$total.events),
                    cluster = fsom$metaclustering[fsom$map$mapping[,1]])
  message("Generating per-sample cluster counts...")
  dat <- data.frame(do.call(rbind,lapply(split(dat,dat$sample),function(i) table(i$cluster))),
                    check.names = F)
  dat$sample <- rownames(dat); rownames(dat) <- NULL
  return(dat)
}

fsom.node.counts <- function(fsom=fsom){
  if(is.null(fsom$mdat)&is.null(fsom$map$mapping)){
    stop("Missing list features")
  }
  dat <- data.frame(sample = rep(fsom$mdat$sample,fsom$mdat$total.events),
                    node = factor(fsom$map$mapping[,1]))
  message("Generating per-sample node counts...")
  dat <- data.frame(do.call(rbind,lapply(split(dat,dat$sample),function(i) table(i$node))),
                    check.names = F)
  dat$sample <- rownames(dat); rownames(dat) <- NULL
  return(dat)
}

#' Generate per-sample counts (cluster/node)
#'
#' @param fsom fsom (fsom.object in environment)
#'
#' @return a list of length 2; counts$cluster and counts$node
#' @export
#'
generate.counts <- function(fsom=fsom){
  counts <- list(cluster=fsom.cluster.counts(fsom),
                 node=fsom.node.counts(fsom))
  return(counts)
}

fsom.get.nodes.from.cluster <- function(metaclusters,cluster.number){
  if(!is.factor(metaclusters)){
    stop("Need factored metaclusters")
  }
  which(metaclusters %in% cluster.number)
}

fsom.merge.nodes <- function(metaclusters,nodes){
  if(!is.factor(metaclusters)){
    stop("Need factored metaclusters")
  }
  m <- metaclusters
  levels(m) <- c(levels(m),length(levels(m))+1)
  m[nodes] <- length(levels(m))
  m <- factor(m)
  levels(m) <- c(1:length(levels(m)))
  return(m)
}

echo.melt <- function(fsom,counts=c('cluster','node'),drop.factor=NULL){
  if(is.null(fsom$counts)){
    stop("Need cluster/node counts...")
  }
  v.name <- match.arg(counts)
  counts <- switch(v.name,
                   cluster = fsom$counts$cluster,
                   node = fsom$counts$node)
  c.names <- colnames(counts)[sapply(counts,is.numeric)]
  counts.frame <- suppressMessages(plyr::join(fsom$mdat,counts))
  props.frame <- counts.frame
  props.frame[,c.names] <- prop.table(as.matrix(props.frame[,c.names]),1)*100
  if(length(c.names)){
    counts.melt <- reshape2::melt(counts.frame,
                                  measure.vars = c.names,
                                  variable.name = v.name,
                                  value.name = paste("count of",unique(fsom$mdat$cell.type))
    )
    prop.melt <- reshape2::melt(props.frame,
                                measure.vars = c.names,
                                variable.name = v.name,
                                value.name = paste("% of",unique(fsom$mdat$cell.type))
    )
    dat.melt <- suppressMessages(plyr::join(counts.melt,prop.melt))
  }
  if(!is.null(drop.factor)){
    dat.melt <- droplevels(dat.melt[dat.melt$condition!=drop.factor,])
  }
  return(dat.melt)
}

fsom.somnambulation <- function(fsom.rds.path,c.names=NULL){
  ##load fsom.rds
  if(!exists('fsom')){
    message(paste("loading fsom.rds:",fsom.rds.path))
    fsom <- readRDS(fsom.rds.path)
  }else{
    message("fsom is already defined (in environment)")
  }
  ##total.events
  if(!is.null(fsom$mdat)){
    if(length(fsom$map$mapping[,1])==sum(fsom$mdat$total.events)){
      total.events <- sum(fsom$mdat$total.events)
    }else{
      total.events <- sum(fsom$mdat$total.events)
      message("length(fsom$map$mapping[,1]) DOES NOT EQUAL sum(fsom$mdat$total.events)")
    }
  }else{
    total.events <- length(fsom$map$mapping[,1])
  }
  ##generate mdats
  if(!is.null(fsom$mdat)){
    mdats = list(cluster = echo.melt(fsom,counts = 'cluster'),
                 node = echo.melt(fsom,counts = 'node')
    )
  }else{
    mdats <- NULL
  }

  if(is.null(c.names)){
    c.names <- colnames(fsom$data)
    c.names <- c(c.names, "cluster", "node")
  }else if(c.names=="ECHO"){
    ##prepare column names
    message(paste("Column name style/ordering:",c.names))
    c.names <- echo.fsom.colnames(c.names = colnames(fsom$data))
  }

  ##prepare data.frame for plotting all events; sub-sample if too many rows to speed up plotting
  message("dat.all")
  if(nrow(fsom$data)>1E6){
    message("All events > 1E6; sub-sampling...")
    set.seed(20040501)
    sample.val <- sample(1:nrow(fsom$data),1E6)
    dat.all <- cbind(fsom$data[sample.val,],
                     cluster = as.numeric(fsom$metaclustering[fsom$map$mapping[,1]][sample.val]),
                     node = fsom$map$mapping[,1][sample.val])[,c.names]
  }else{
    dat.all <- cbind(fsom$data,
                     cluster = as.numeric(fsom$metaclustering[fsom$map$mapping[,1]]),
                     node = fsom$map$mapping[,1])[,c.names]
  }

  ##prepare cluster matrices
  message("cluster.mats")
  cluster.mats <- sapply(levels(fsom$metaclustering),function(i){
    c.index <- which(fsom$metaclustering[fsom$map$mapping[,1]]==i)
    if(length(c.index)>2E5){
      set.seed(20040501)
      sample.val <- sample(c.index,2E5)
      dat <- cbind(fsom$data[sample.val,],cluster = as.numeric(i),node = fsom$map$mapping[,1][sample.val])[,c.names]
    }else{
      dat <- cbind(fsom$data[c.index,],cluster = as.numeric(i),node = fsom$map$mapping[,1][c.index])[,c.names]
    }
    return(dat)
  },simplify = F)

  ##pre-calculate per-column x,y limits
  ##decimal value ceiling for setting x,y upper limits and lower limits
  ##can probably simplify this chunk of code but works fine...
  message("xy.limits")
  col.maxes.all <- apply(dat.all,2,max)
  col.mins.all <- apply(dat.all,2,min)
  col.maxes.mats <- apply(sapply(cluster.mats,function(i) apply(i,2,max)),1,max)
  col.mins.mats <- apply(sapply(cluster.mats,function(i) apply(i,2,min)),1,min)
  col.maxes <- apply(rbind(col.maxes.all,col.maxes.mats),2,max)
  col.mins <- apply(rbind(col.mins.all,col.mins.mats),2,min)
  #
  ceiling_dec <- function(x, level=1){round(x + 5*10^(-level-1), level)}
  #
  xy.upper.lim <- c(sapply(col.maxes[!names(col.maxes) %in% c('cluster','node')],ceiling_dec)+0.05,
                    col.maxes[c('cluster','node')]+1
  )*1
  xy.lower.lim <- c(sapply(abs(col.mins[!names(col.mins) %in% c('cluster','node')]),ceiling_dec)+0.05,
                    col.mins[c('cluster','node')]-1
  )*-1
  xy.lims <- rbind(xy.lower.lim,xy.upper.lim)

  message("making list")
  gc()
  fsom.somnambulated <- list(total.events = total.events,
                             dat.all = dat.all,
                             cluster.mats = cluster.mats,
                             xy.lims = xy.lims,
                             mdats = mdats,
                             nc.vals = col.maxes[c('cluster','node')],
                             heatmaps = list(p.heat = pheat(dat.mat = fsom$cluster.medians,color.type = 'sequential',border_color = NA,silent=T),
                                             pl.heat = plheat(fsom$cluster.medians),
                                             cluster.corr.heat = pheatmap::pheatmap(stats::cor(log2(fsom$counts$cluster[grep("HD",fsom$counts$cluster$sample,invert=T),sapply(fsom$counts$cluster,is.numeric)]+1)),
                                                                                    border_color=NA,silent=T
                                             )
                             ),
                             metaclustering = fsom$metaclustering,
                             cluster.titles = sapply(levels(fsom$metaclustering),function(i){
                               nodes.in.cluster <- which(fsom$metaclustering==as.numeric(i))
                               if(length(nodes.in.cluster)>10){
                                 title.sub <- paste0(paste0("(",paste(nodes.in.cluster[1:10],collapse = " "),"...)"),"(",length(nodes.in.cluster)," total",")")
                               }else{
                                 title.sub <- paste0("(",paste(nodes.in.cluster,collapse = " "),")")
                               }
                             }),
                             node.titles = stats::setNames(
                               sapply(as.numeric(fsom$metaclustering),function(i){
                                 title.sub <- paste0("(",i,")")
                               }),
                               nm = seq(length(fsom$metaclustering))
                             )
  )
  class(fsom.somnambulated) <- "FlowSOM Somnambulated"

  return(fsom.somnambulated)
}

#' Map samples to existing FlowSOM codes/clusters
#'
#' @param fcs.file.paths character vector of .fcs file paths; returned from 'ECHOfcs::get.fcs.paths.echo()'
#' @param fsom.codes.path existing codes/clusters list saved as .rds
#'
#' @return a per-sample list of mapped nodes/clusters; essentially a per-sample,per-event/row node/cluster assignment/index
#' @export
#'
map.samples.to.fsom.codes <- function(fcs.file.paths,fsom.codes.path){
  fsom.codes <- readRDS(fsom.codes.path)
  if(!is.list(fsom.codes)){
    stop("Expect a list")
  }
  ##
  message("Parallel reading/mapping of .fcs files...")
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  map.list <- parallel::parSapply(cl,fcs.file.paths,function(i,dims=unname(colnames(fsom.codes$codes)),
                                                             codes=fsom.codes$codes,
                                                             clusters=fsom.codes$clusters,
                                                             f.levels=seq(fsom.codes$clusters)){
    #
    fcs.tmp <- ECHOfcs::read.fcs.selected.markers(fcs.path = i,selected.markers = dims)
    if(!all(match(sapply(strsplit(flowCore::markernames(fcs.tmp),"_"),'[[',2),dims)==seq(dims))){
      stop('column/marker name error')
    }
    flowCore::colnames(fcs.tmp) <- sapply(strsplit(flowCore::markernames(fcs.tmp),"_"),'[[',2)
    fcs.tmp@exprs <- asinh(fcs.tmp@exprs)/5
    map <- EmbedSOM::MapDataToCodes(data = fcs.tmp@exprs,
                                    codes = codes)
    list(node = factor(map[,1],levels = f.levels),
         cluster= clusters[map[,1]])
  },simplify = F)
  parallel::stopCluster(cl)
  return(map.list)
}


# echo.stats.pval <- function(mdat = fsom$mdat,cluster.counts=fsom$cluster.counts,factor.group='condition'){
#   dat.melt <- echo.melt()
#   y.name <- grep("% of", colnames(dat.melt),value = T)
#
#   p.val <- sapply(levels(dat.melt$cluster),function(i){
#     round(t.test(dat.melt[dat.melt$cluster==i,y.name] ~ dat.melt[dat.melt$cluster==i,factor.group])$p.value,3)
#   })
#
#   p.vals.frame <- data.frame(cluster=names(p.val),p.val)
#   p.vals.frame$p.val.adjust <-p.adjust(p.vals.frame$p.val,
#                                        method = "bonferroni"
#   )
#   p.vals.frame$mark <- ""
#   p.vals.frame$mark[p.vals.frame$p.val.adjust<0.05] <- "*"
#   return(p.vals.frame)
# geom_text(data = p.vals.frame, aes(label=mark),size=10,
#           x = Inf, y = Inf, hjust = 1.25, vjust = 1.25,
#           inherit.aes = FALSE)
# }
