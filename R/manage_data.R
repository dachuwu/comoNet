
set.dzlevel <- function(df, 
                       var_dz=c("ICD_1","ICD_2","ICD_3","ICD_4","ICD_5"), 
                       dz.lv=NULL){
  if(is.null(dz.lv)){ 
    stop("Must specify dz.lv!")
  } else if(class(dz.lv) != "character"){
    stop("Incorrect class for dz.lv; must be a character vector.")
  }
  if(!all(var_dz%in%colnames(df))) stop("Invalid var_dz.")
  
  for(v in var_dz){
    df[,v] <- as.integer(factor(df[,v], levels = dz.lv))
  }
  ind <- apply(df[,var_dz], 1, function(x){
    any(!is.na(x))})
  cat("\nremove ",NROW(df)-length(ind)," records with NA dz's.\n")
  df <- df[ind,]
  row.names(df)<-NULL
  return(df)
} 

if(F){
  sub_ipseq <- function(x, df, 
                     var_id ="ID",
                     var_t = "date", 
                     var_dz = c("ICD_1","ICD_2","ICD_3","ICD_4","ICD_5")){
  sdf <- df[df[,var_id] == x, c(var_t,var_dz)]
  
  out <- lapply(var_dz, function(nv){
    ind <- !is.na(sdf[,nv])
    data.frame(
      t = sdf$date[ind], 
      dz = sdf[ind,nv], stringsAsFactors = F)
  }) 
  out <- do.call(rbind, out)
  out <- out[order(out$t),]
  out <- out[!duplicated(out$dz),]
  row.names(out)<-NULL
  return(out)
}
}


get_ipseq <- function(df, var_id ="ID",
                      var_t = "date", 
                      var_dz = c("ICD_1","ICD_2","ICD_3","ICD_4","ICD_5"), 
                      Par= F){
  
  df <- df[,c(var_id,var_t,var_dz)]
  sub_ipseq <- function(sdf, svar_t = var_t, svar_dz = var_dz){
    
    out <- lapply(svar_dz, function(nv){
      ind <- which(!is.na(sdf[,nv]))
      if(length(ind)>0){
        data.frame(
          t = sdf[,svar_t][ind], 
          dz = sdf[ind, nv], stringsAsFactors = F)}
    }) 
    out <- do.call(rbind, out)
    out <- out[order(out$t),]
    out <- out[!duplicated(out$dz),]
    row.names(out)<-NULL
    return(out)
  }
  
  df <- split(df, as.factor(df[,var_id]))
  # ipseq <- by(df,  as.factor(df[,var_id]), function(k) sub_ipseq(sdf = k) )
    
  if(Par==F){
    ipseq <- lapply(df, function(k) sub_ipseq(sdf = k) )
  
  } else if (Par==T){

    cl <- parallel::makeCluster(getOption("cl.cores", parallel::detectCores()))
    parallel::clusterExport(cl, envir =environment(), 
                            varlist = c("sub_ipseq","var_t","var_dz") )
    ipseq <- parallel::parLapply(cl, df,  function(k) sub_ipseq(k))
    parallel::stopCluster(cl)
    
  } else if (is.integer(Par)){
    
    cl <- parallel::makeCluster(getOption("cl.cores", min(parallel::detectCores(),Par)))
    parallel::clusterExport(cl, envir =environment(), 
                            varlist = c("sub_ipseq","var_t","var_dz"))
    ipseq <- parallel::parLapply(cl, df,  function(k) sub_ipseq(k))
    parallel::stopCluster(cl)
    
  } else {
    stop("incorrect specification of Par (must be True, False or integer)")
  }
  return(ipseq)

}


