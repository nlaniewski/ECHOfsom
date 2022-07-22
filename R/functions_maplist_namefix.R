##this chunk of code is to handle reassignment of "089-2_V4_SEB"/"HD0191_Adult_SEB_014" (barcode issue);
##and drop low-viability samples mistakenly plated
map.list.namefix.echo <- function(map.list){
  i <- names(map.list)
  ##store index
  index.hd0191_adult_seb_014 <- grep("089-2_V4_SEB",i)
  f.name.089_v4_seb <- names(map.list)[index.hd0191_adult_seb_014]
  index.089_v4_seb <- grep("HD0191_Adult_SEB_014",i)
  f.name.hd0191_adult_seb_014 <- names(map.list)[index.089_v4_seb]
  ##reassign names
  names(map.list)[index.hd0191_adult_seb_014] <- f.name.hd0191_adult_seb_014
  names(map.list)[index.089_v4_seb] <- f.name.089_v4_seb
  ##mistakenly plated;low viability following thaw/overnight rest
  low.viability.samples <- c("022-2_V6","038-2_V6","058-2_V4")
  map.list <- map.list[grep(paste0(low.viability.samples,collapse = "|"),i,invert = T)]
  ##
  return(map.list)
}

##this chunk of code is to handle reassignment of "089-2_V4_SEB"/"HD0191_Adult_SEB_014" (barcode issue);
##and drop low-viability samples mistakenly plated
cluster.namefix.echo <- function(cluster.counts.frame){
  i <- cluster.counts.frame
  ##store counts
  cluster.counts.hd0191_adult_seb_014 <- i[grep("089-2_V4_SEB",rownames(i)),]
  cluster.counts.089_v4_seb <- i[grep("HD0191_Adult_SEB_014",rownames(i)),]
  ##reassign
  i[grep("089-2_V4_SEB",rownames(i)),] <- cluster.counts.089_v4_seb
  i[grep("HD0191_Adult_SEB_014",rownames(i)),] <- cluster.counts.hd0191_adult_seb_014
  ##
  i <- i[grep("ECHO_020",rownames(i),invert = T),]
  ##
  low.viability.samples <- c("022-2_V6","038-2_V6","058-2_V4")#mistakenly plated;low viability following thaw/overnight rest
  i <- i[grep(paste0(low.viability.samples,collapse = "|"),rownames(i),invert = T),]
  ##
  return(i)
}
