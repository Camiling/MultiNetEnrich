#' Perform functional network enrichment for multiple data sets
#'
#' This function performs functional network enrichment for multiple data sets
#'
#' @param X list of the \eqn{K} observed \eqn{n_k} by \eqn{p} data matrices
#' @param database gene set database(s) to use from MSigDb. Provided as one or more strings containing the URL(s), from \url{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/}.
#'                 If not specified, Wikipathways, PID and Hallmark are used.
#' @param nperm the number of random permuations to use. Default \eqn{1000}.
#' @param maxitr maximum number of iterations. Default \eqn{1000}.
#' @param scale should variables be scaled? Default \code{TRUE}
#'
#' @return a fitted \code{MultiNetEnrich} object
#'
#' @importFrom foreach %dopar%
#'
#'
#' @export
#'
MultiNetEnrich <- function(X, database = NULL,nperm=1000, maxitr = 1e3, scale=TRUE){

  # If only one data set is provided, single-network version is used instead.
  if(!is.list(X)){
    stop('X must be a list of data matrices.')
  if(! length(unique(unlist(lapply(X,ncol))))==1){
    stop('All data matrices must have the same number of variables.')
  }
  }
  if(is.null(database)){
    database = c('https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/c2.cp.wikipathways.v7.5.1.symbols.gmt',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/h.all.v7.5.1.symbols.gmt',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/c2.cp.pid.v7.5.1.symbols.gmt')
  }
  else{
    cat('Checking validity of URLs provided...')
    if(!all(sapply(database ,valid_url))){
      stop('\n Please provide valid URLs.')
    }
    cat('\n URLs valid. Commencing network analysis. \n')
  }

  p <- dim(X[[1]])[2]
  K <- length(X)
  # Perform joint network analysis
  res.joint = jointGHS::jointGHS(X,epsilon=1e-3,AIC_eps = 0.1,maxitr=maxitr,scale=scale,verbose=F)
  cat('Network analysis done.')

  # Write ranked gene lists to files
  dir.create("MultiNetEnrich_tmp")
  ranked.lists = list()
  for(k in 1:K){
    theta.est = stats::cov2cor(res.joint$theta[[k]])
    theta.est[which(abs(theta.est) < 1e-5, arr.ind = T)] = 0
    g.tmp = igraph::graph.adjacency(theta.est!=0, mode='undirected', diag=F)
    df.deg = data.frame(gene = colnames(X[[k]]),
                        degree=igraph::degree(g.tmp)+1)
    df.deg.ordered = df.deg[rev(order(df.deg$degree)),]
    ranked.lists[[k]] = df.deg.ordered
    df.deg.unique = df.deg.ordered[which(!duplicated(df.deg.ordered$gene)),]
    utils::write.table(df.deg.unique, file = paste0('MultiNetEnrich_tmp/df_degree_', k, '.rnk'), sep = "\t", quote=FALSE,
                row.names = F, col.names = F)
  }
  # Perform GSEA
  cat('Commencing functional network enrichment analysis. \n')
  enriched.list = list()
  for(k in 1:K){
    for(i in 1:length(database)){
      dir.create("MultiNetEnrich_tmp_gsea")
      tmp.gsea = GSEA::GSEA(input.ds=paste0('MultiNetEnrich_tmp/df_degree_', k, '.rnk'), input.cls = NULL,gs.db = database[i],
                      reshuffling.type = "gene.labels", output.directory = 'MultiNetEnrich_tmp_gsea', gs.size.threshold.min = 2,
                      gsea.type='preranked', nperm=nperm)
      unlink("MultiNetEnrich_tmp_gsea", recursive = TRUE)
      enriched.tmp = tmp.gsea$report1[tmp.gsea$report1$`FDR q-val`<0.25 & tmp.gsea$report1$`NOM p-val`<0.05,]
      if(nrow(enriched.tmp)>0){
        if(length(enriched.list) < k){
          enriched.list[[k]] = enriched.tmp
        }
        else{
          enriched.list[[k]] = rbind(enriched.list[[k]], enriched.tmp)
        }
      }
    }
  }
  unlink("MultiNetEnrich_tmp", recursive = TRUE)
  out=list()
  out$enrichment.list = enriched.list
  out$networks = lapply(res.joint$theta, FUN = function(s) abs(s)>1e-5)
  out$ranked.lists = ranked.lists
  class(out) = "MultiNetEnrich"
  return(out)
}

