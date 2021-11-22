LengthCheck <- function(values, cutoff=0){
  return(vapply(X = values, FUN = function(x) {
    return(length(x = x) > cutoff)
  }, FUN.VALUE = logical(1)))
}

testModuleScore <- function (expMatrix, genes.list = NULL, genes.pool = NULL, n.bin = 25,
                             seed.use = 1, ctrl.size = 100, enrich.name = "Cluster", uniqueOnly=T)
{
  genes.old <- genes.list
  if (0) {
    1==1
  }
  else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.list <- lapply(X = genes.list, FUN = function(x) {
      return(intersect(x = x, y = rownames(x = expMatrix)))
    })
    cluster.length <- length(x = genes.list)
  }
  if (!all(LengthCheck(values = genes.list))) {
    warning(paste("Could not find enough genes in the object from the following gene lists:",
                  paste(names(x = which(x = !LengthCheck(values = genes.list)))),
                  "Attempting to match case..."))
    genes.list <- lapply(X = genes.old, FUN = CaseMatch,
                         match = rownames(x = expMatrix))
  }
  if (!all(LengthCheck(values = genes.list))) {
    stop(paste("The following gene lists do not have enough genes present in the object:",
               paste(names(x = which(x = !LengthCheck(values = genes.list)))),
               "exiting..."))
  }
  if (is.null(x = genes.pool)) {
    genes.pool = rownames(x = expMatrix)
  }
  data.avg <- apply(X = expMatrix[genes.pool, ], MARGIN = 1,
                    FUN = mean)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg)/n.bin)))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector("list", cluster.length)
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(x = genes.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut ==
                                                                              data.cut[genes.use[j]])], size = ctrl.size, replace = FALSE)))
    }
  }
  if (uniqueOnly){
    ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  }

  ctrl.scores <- c()
  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores <- rbind(ctrl.scores, apply(X = expMatrix[genes.use,
    ], MARGIN = 2, FUN = mean))
  }
  genes.scores <- c()
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    genes.scores <- rbind(genes.scores, apply(X = expMatrix[genes.use,
    ], MARGIN = 2, FUN = mean))
  }
  genes.scores.use <- genes.scores - ctrl.scores
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  genes.scores.use <- t(x = as.data.frame(x = genes.scores.use))
  return(genes.scores.use)
}
