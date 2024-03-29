if(F){
comoNet_from_prog <- function(ipseq, dz.lv, 
                              bidir=F, conditional=F){
  
  
  if(class(ipseq)!="list") stop("ipseq must be a list containing individual progression data.frame.")
  
  dti <- F
  for(i in 1:length(ipseq)){
    
    if(!all(colnames(ipseq[[i]])[1]=="t", colnames(ipseq[[i]])[2]=="dz")){
      stop("incorrect variable names of individual progression data.frame [",i,"].")
    }
    if(class(ipseq[[i]]$t)%in%c("Date")){
      ipseq[[i]]$t <- as.integer(ipseq[[i]]$t)
      dti <- T
    }else if(!class(ipseq[[i]]$t)%in%c("integer","numeric")){
      stop("incorrect progression time class [",i,"]; should be integer, numeric or Date class.")
    }
    
    if(class(ipseq[[i]]$dz)%in%c("factor", "character")){
      ipseq[[i]]$dz <- as.integer(factor(as.character(ipseq[[i]]$dz), levels = dz.lv))
      warning("coerce disease (factor, character) to integer based in dz.lv.")
    }else if(!class(ipseq[[i]]$dz)%in%c("integer")){
      stop("incorrect disease class [",i,"]; should be integer, factor or character.")
    }
    
  }
  if(dti) warning("coerce progression time (Date) to integer")
  cat("\n ipseq structure checked. \n")
  
  
  if(conditional){
    build_CondRaw(ipseq, dz.lv, bidir)
  } else {
    build_Raw(ipseq, dz.lv, bidir)
  }
  
  
}


}

require(Matrix)
comoNet_from_prog_sp <- function(ipseq, dz.lv, 
                              bidir=F, conditional=F){
  
  
  if(class(ipseq)!="list") stop("ipseq must be a list containing individual progression data.frame.")
  
  dti <- F
  for(i in 1:length(ipseq)){
    
    if(!all(colnames(ipseq[[i]])[1]=="t", colnames(ipseq[[i]])[2]=="dz")){
      stop("incorrect variable names of individual progression data.frame [",i,"].")
    }
    if(class(ipseq[[i]]$t)%in%c("Date")){
      ipseq[[i]]$t <- as.integer(ipseq[[i]]$t)
      dti <- T
    }else if(!class(ipseq[[i]]$t)%in%c("integer","numeric")){
      stop("incorrect progression time class [",i,"]; should be integer, numeric or Date class.")
    }
    
    if(class(ipseq[[i]]$dz)%in%c("factor", "character")){
      ipseq[[i]]$dz <- as.integer(factor(as.character(ipseq[[i]]$dz), levels = dz.lv))
      warning("coerce disease (factor, character) to integer based in dz.lv.")
    }else if(!class(ipseq[[i]]$dz)%in%c("integer")){
      stop("incorrect disease class [",i,"]; should be integer, factor or character.")
    }
    
  }
  if(dti) warning("coerce progression time (Date) to integer")
  cat("\n ipseq structure checked. \n")
  
  
  if(conditional){
    build_CondRaw_sp(ipseq, dz.lv, bidir)
  } else {
    build_Raw_sp(ipseq, dz.lv, bidir)
  }
  
  
}
