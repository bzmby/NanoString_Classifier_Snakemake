###########################################################################
###########################################################################
#               Normalization - mRNA
###########################################################################
###########################################################################

#### function to plot norm factors vs. deviance from typical HK profile:
# HKnormalized the normalized data from the HKs
# norm.factors: the normalization factors from each sample
# analyte.type: "gene" "protein" etc
normalization.summary.plot = function(HKnormalized,norm.factors,analyte.type,col="black",main="",cex=cex,pch=pch)
{
  # for each sample, calculate the MSE of its HKs from the average profile:
  avg.hk.profile = apply(HKnormalized,2,function(x) mean(x,na.rm = T))
  HK.MSE = rowMeans((sweep(HKnormalized,2,avg.hk.profile))^2)
  
  # identify high HK.MSE outliers:
  outlier.HKMSE = HK.MSE > mean(HK.MSE,na.rm=T)+2*sd(HK.MSE,na.rm=T)

  # identify hkgeomean outliers:
  outlier.normfactor = abs(norm.factors-mean(norm.factors,na.rm=T))>2*sd(norm.factors,na.rm=T)
  
  # plot the result:
  par(las = 1)
  plot(norm.factors,HK.MSE,xlab = "Normalization factor",ylab = paste0("MSE of reference ",analyte.type,"s from mean profile"),col=col,pch=pch,main=main,cex=cex)
  highlight = outlier.HKMSE|outlier.normfactor
  if(sum(highlight,na.rm = T)>0)
  {
    text(norm.factors[highlight],HK.MSE[highlight],rownames(HKnormalized)[highlight])
  }
  
  mse.df <- data.frame("normalization.factor" = norm.factors, "HK.MSE" = HK.MSE)
  colnames(mse.df)[1] <- paste(analyte.type, ".normalization.factors", sep = "")
  write.csv(mse.df, file = paste(path.to.normalization.results,"//", analyte.type, "_normalization_summary.csv",sep=""))
  
}



#### function to normalize the raw data given a list of HKs:
normalize.given.HKs = function(raw,HKdata,HKs)
{   
  # use only the selected HKs:
  HKdata.selected = HKdata[,is.element(dimnames(HKdata)[[2]],HKs)]
  # get mean of HKs
  HKmeans = apply(HKdata.selected,1,mean)
  # apply norm factors:
  norm.factors = HKmeans-mean(HKmeans,na.rm = T)
  normalized = raw-norm.factors
  HKnormalized = HKdata-norm.factors
  out = list(normalized=normalized,HKnormalized=HKnormalized,norm.factors=norm.factors)
  return(out)
}

normalizerawdata = function(raw,
                            processed,
                            data.type,
                            method=c("choose HKs","use selected","none"),
                            choose.method = c("geNorm"),
                            n.HKs = 10,
                            auto.HKs,
                            codeset.HKs,
                            path.to.normalization.results,
                            log,
                            plottypearg,
                            path.results,
                            path.inc,
                            prb.annots)
{
  print("Starting normalization")
  cat("LOG:Starting normalization",file=log,sep='\n\n',append=TRUE)
  #cat("document.write('<p>Starting normalization</p>');", file=paste(path.inc,"//status.js",sep=""),append=TRUE)
  
  codeset.genes = dimnames(raw)[[2]]
  
  ### basic processing: threshold and log-transform:
  raw = replace(raw,raw<1,1)
  raw = log2(raw)
  
  #processed <- replace(raw,raw<1,1)
  processed <- replace(processed,processed<1,1)
  processed <- log2(processed)
  
  ### split data into HKs and informative genes:

  suppressWarnings(write.table(table(codeset.HKs==1),file=log,sep='\t',append=TRUE))
  suppressWarnings(write.table(table(codeset.HKs),file=log,sep='\t',append=TRUE))
  cat("dim(raw):",file=log,sep='\t',append=TRUE)
  cat(dim(raw),file=log,sep='\t',append=TRUE)
  HKdata = raw[,codeset.HKs,drop=F]
  raw = raw[,setdiff(dimnames(raw)[[2]],codeset.HKs)]

  # if it's already normalized data, do nothing:
  #---------------------------------------------
  if(method=="none")
  {
    warning("Warning: No mRNA normalization will be performed - assuming the input data is already normalized")
    cat("LOG: Warning: No mRNA normalization will be performed - assuming the input data is already normalized",file=log,sep='\n\n',append=TRUE)
    cat("document.write('<p>Warning: No mRNA normalization will be performed - assuming the input data is already normalized</p>');", file=paste(path.inc,"//status.js",sep=""),append=TRUE)
    
    normalized <- processed
    raw <- cbind(raw,HKdata)
    raw <- raw[,colnames(processed)[colnames(processed) %in% colnames(raw)],drop=F]
    if(!isTRUE(all.equal(raw,processed,check.attributes=F)) & data.type == "raw")
      normalized <- raw
    HKnormalized = NULL#HKdata
    HKs = NULL
    normalizedall = list(); normalizedall$norm.factors = rep(0,dim(raw)[1])
  }

  # use selected user selected HK genes
  #-------------------------------
  if(method == "use selected"){
    message("Message: All user selected HK will be used. No HK selection will be performed")
    cat("LOG: Message: All user selected HK will be used. No HK selection will be performed",file=log,sep='\n\n',append=TRUE)
    cat("document.write('<p>Message:All user selected HK will be used. No HK selection will be performed</p>');", file=paste(path.inc,"//status.js",sep=""),append=TRUE)
    
    HKs <- codeset.HKs
    # now normalize:
    normalizedall = normalize.given.HKs(raw,HKdata,HKs) 
    normalized = normalizedall$normalized
    HKnormalized = normalizedall$HKnormalized
    
    # plot the normalization results:
    for(r in 1:length(plottypearg)){
      plottype=plottypearg[r];
      tempfilename = drawplot(filename=paste(path.to.normalization.results,"//gx normalization results",sep=""),plottype)
      tempfilename=gsub(path.results,"results",tempfilename)
      normalization.summary.plot(HKnormalized=HKnormalized,norm.factors=normalizedall$norm.factors,analyte.type="mRNA",
                                 col="slateblue2",main = "Normalization summary: mRNA",cex=2,pch=16)
      dev.off()}
    
  }
  
  
  # select normalization genes in a data-driven manner:
  #-------------------------------------------------------
  if(method == "choose HKs")
  {
    # make sure there's an allowed method chosen:
    if(!is.element(choose.method,c("geNorm")))
    {
      warning("Warning: No method for choosing housekeepers selected - using the geNorm algorithm")
      cat("LOG: Warning: No method for choosing housekeepers selected - using the geNorm algorithm",file=log,sep='\n\n',append=TRUE)
      cat("document.write('<p>Warning: No method for choosing housekeepers selected - using the geNorm algorithm</p>');", file=paste(path.inc,"//status.js",sep=""),append=TRUE)
      choose.method="geNorm"
    }
    
    if(choose.method=="geNorm")
    {
      minNrHKs = 2
      # run geNorm:
      temp = selectHKs(HKdata,method="geNorm",minNrHKs=minNrHKs,log=TRUE,Symbols = dimnames(HKdata)[[2]],trace=FALSE)
      ## select top number of genes:
      # if you decided beforehand how many to take:
      if(!auto.HKs)
      {
        HKs = as.vector(temp$ranking)[1:n.HKs]
      }
      # if you want to choose the top housekeepers dynamically:
      nselect=n.HKs
      if(auto.HKs)
      {
        nselect = dim(HKdata)[2]-(1:length(temp$variation))[temp$variation==min(temp$variation)]+1
        # require at least 6 HKs:
        nselect = max(nselect,6)
        HKs = as.vector(temp$ranking)[1:nselect] 
      }
      
      
      # plot the results as a sanity check:
      for(r in 1:length(plottypearg)){
        plottype=plottypearg[r];
        tempfilename = drawplot(filename=paste(path.to.normalization.results,"//HK selection details - pairwise variance",sep=""),plottype)
        tempfilename=gsub(path.results,"results",tempfilename)
        
        par(las = 1)
        plot(temp$variation,
             col=c(rep("grey",dim(HKdata)[2]-nselect),rep("black",nselect)),
             pch=16,
             main="Genes selected using geNorm",
             ylab="Pairwise variation during stepwise selection",xlab="Order removed" #,ylim = c(.95,1.2)*range(temp$variation)
             )
        #text(1:length(temp$variation),temp$variation,labels = names(temp$variation),pos = 3)
        
        abline(v=dim(HKdata)[2]-nselect+.5)
        legend("topleft",col=c("black","grey"),pch=c(16,1),legend=c("selected","unselected"))
        dev.off()}
      
    }
    
    # now normalize:
    normalizedall = normalize.given.HKs(raw,HKdata,HKs) 
    normalized = normalizedall$normalized
    HKnormalized = normalizedall$HKnormalized
    
   
    # plot the normalization results:
    for(r in 1:length(plottypearg)){
      plottype=plottypearg[r];
      tempfilename = drawplot(filename=paste(path.to.normalization.results,"//gx normalization results",sep=""),plottype)
      tempfilename=gsub(path.results,"results",tempfilename)
      normalization.summary.plot(HKnormalized=HKnormalized,norm.factors=normalizedall$norm.factors,analyte.type="mRNA",
                                 col="slateblue2",main = "Normalization summary: mRNA",cex=2,pch=16)
      dev.off()}
    
  }
  
  
  # write csv of selected HKs:
  if(method!="none")
  {
    print("Normalized mRNA data")
    cat("LOG:Normalized mRNA data",file=log,sep='\n\n',append=TRUE)
    cat("document.write('<p>Normalized mRNA data</p>');", file=paste(path.inc,"//status.js",sep=""),append=TRUE)
    cat("document.write('  	  	<ul>');\n",file=paste(path.inc,"//panel2_5.js",sep=""),append=TRUE)
    
    HKstowrite = data.frame(HKs,stringsAsFactors = FALSE)
    names(HKstowrite) = "HKs in the order by which they were selected"
    #write.csv(HKstowrite,file=paste(path.to.normalization.results,"//selected housekeepers.csv",sep=""),row.names=FALSE)
    # write unselected HKs:
    unselected.HKs = setdiff(codeset.HKs,HKs)
    #write.csv(unselected.HKs,file=paste(path.to.normalization.results,"//discarded housekeepers.csv",sep=""),row.names=FALSE)
    #     hk.selection.summary = rbind(cbind(HKstowrite,1:nrow(HKstowrite)),
    #                                  cbind(unselected.HKs),rep("discarded",length(unselected.HKs)))
    hk.selection.summary <- data.frame(HK = c(HKstowrite[,1],unselected.HKs),
                                       selection = c(1:nrow(HKstowrite),rep("discarded",length(unselected.HKs))))
    # report variance of each gene after normalization:
    normgenesd = signif(apply(HKnormalized[,as.character(hk.selection.summary[,1])],2,sd),3)
    tmp.lab <- paste(prb.annots[as.character(hk.selection.summary[,1]),"Probe.Label"],prb.annots[as.character(hk.selection.summary[,1]),"Analyte.Type"],sep="-")
    hk.selection.summary = cbind(tmp.lab,hk.selection.summary,normgenesd)
    colnames(hk.selection.summary)=c("Gene Name","Gene","Order selected by geNorm","SD after normalization")
    
    write.csv(hk.selection.summary,file=paste(path.to.normalization.results,"//selected housekeepers.csv",sep=""),row.names=FALSE)
    
    file=paste("results//Normalization","//selected housekeepers.csv",sep="")
    strTemp=paste("document.write('	       <li><a href=\"",file,"\">View Selected HK Genes</a></li>')\n",sep="")
    cat(strTemp,file=paste(path.inc,"//panel2_5.js",sep=""),append=TRUE)
    
  }
  
  # threshold normalized data at 1: <----------turn of the thresholding
  # normalized = replace(normalized,normalized<0,0)
  # HKnormalized = replace(HKnormalized,HKnormalized<0,0)
  
  out = list(normalized.data=normalized,HKnormalized=HKnormalized,HKs=HKs,norm.factors=normalizedall$norm.factors)
  # write csv of normalized data:
  colnames(normalized) <- paste(prb.annots[colnames(normalized),"Probe.Label"],prb.annots[colnames(normalized),"Analyte.Type"],sep="-")
  #write.csv(normalized,file=paste(path.to.normalization.results,"//mRNA_normalized_log2_data.csv",sep=""))
  #write.csv(2^normalized,file=paste(path.to.normalization.results,"//mRNA_normalized_linear_data.csv",sep=""))

  log.temp.df <- data.frame("log" = c("mRNA - normalized log2 count data"))
  write.table(log.temp.df,file=paste(path.to.normalization.results,"//mRNA_normalized_data.csv",sep=""), col.names = F, row.names = F)
  write.table(normalized,file=paste(path.to.normalization.results,"//mRNA_normalized_data.csv",sep=""), append = T, sep = ",", col.names = NA, row.names = T)

  linear.temp.df <- data.frame("linear" = c("", "mRNA - normalized linear count data"))
  write.table(linear.temp.df,file=paste(path.to.normalization.results,"//mRNA_normalized_data.csv",sep=""), append = T, sep = ",", col.names = F, row.names = F)
  write.table(2^normalized,file=paste(path.to.normalization.results,"//mRNA_normalized_data.csv",sep=""), append = T, sep = ",", col.names = NA, row.names = T)

  print("Created mRNA normalized data files")
  cat("LOG:Created mRNA normalized data files",file=log,sep='\n\n',append=TRUE)
  cat("document.write('<p>Created mRNA normalized data files</p>');", file=paste(path.inc,"//status.js",sep=""),append=TRUE)
  return(out)
}


###########################################################################
###########################################################################
#               Normalization - protein
###########################################################################
###########################################################################




# @raw.dat: nxp numeric matrix or data.frame containing n samples of p protein measurments
# @neg.colnames: 2 characters indicating the name of the negative proteins
# @norm.ref.prot.names: optional name of proteins to be used in normalization
# @smp.annot: for now it is a place holder for sample annotation (to indicate repeates, reference sample, etc. 
#   that could inform normalization in the future versions)

# values:
# @n.prot: nxp normalized data
# @above.background: nxp boolean call on each data point being below or above background

preprocess.prot <- function(raw.dat,
                            neg.colnames, 
                            norm.ref.prot.names = NULL,
                            norm.method = c("geomean.stable1","geomean.all","geomean.stable2"),
                            normmodule.nchoos,
                            smp.annot = NULL, 
                            prb.annots,
                            path.to.normalization.results,
                            log, 
                            plottypearg,
                            path.results,
                            path.inc){
  
  get.orth.proj <- function(y,x,mod){
    er <- (y - (cbind(1,x) %*% coef(mod)) )/as.numeric(sqrt(t( c(1,-coef(mod)[2])  ) %*% c(1,-coef(mod)[2])))
    ort.proj <- cbind(x,y) - er %*% t(c(-coef(mod)[2],1)/sqrt(c(1,-coef(mod)[2])%*% c(1,-coef(mod)[2])))
    return(ort.proj)
  }
  
  # validation of input
  #--------------------
  if(is.data.frame(raw.dat)){
    numeric.cols <- apply(raw.dat,2,function(x) is.numeric(x))
    if(!all(numeric.cols)){
      warning("non-numeric protein raw.dat columns were dropped")
      raw.dat <- raw.dat[,numeric.cols]
    }
    raw.dat <- as.matrix(raw.dat)
  }
    
  
  if(!is.matrix(raw.dat))
    stop("raw.dat protein needs to be matrix")
  
  if(!all(neg.colnames %in% colnames(raw.dat)))
    stop("neg.colnames for protein data need to be in colnames(raw.dat)")
  
  if(length(neg.colnames) < 2)
    stop("for protein data, at least two neg.colnames are needed to estimate varaibility in lane background counts")
  
  if(!is.null(norm.ref.prot.names)){
    if(any(!norm.ref.prot.names %in% colnames(raw.dat)))
      stop("norm.ref.prot.names for protein data need to be in colnames(raw.dat)")
  }
  
  
  norm.method <- match.arg(norm.method)
  #   if(is.null(norm.ref.prot.names))
          
  #log2 tranform
  prot <- log2(raw.dat + 1)
  
  
  #=============================================
  # Find background (via deming reg model)
  #=============================================
  
  # Here we use the deming package with stdpat(1,0,1,0), which is the classical deming, 
  # note if the ratio of SDs are to be changed via stdpat, pay attention to scaling of sigma accordingly
  # when making CI for the IgGs
  visualize <- F
  if(length(neg.colnames)==2){
    require(deming)
    dmod <- deming(prot[,neg.colnames[2]]~prot[,neg.colnames[1]],stdpat=c(1,0,1,0))
    sigma <- dmod$sigma
    
    alpha <- atan(coef(dmod)[2])
    var.y.on.line <- sigma^2 * sin(alpha)
    var.x.on.line <- sigma^2 * coef(dmod)[2]
    var.on.line <- var.x.on.line + var.y.on.line 
    
    background.stats <- data.frame(prot[,neg.colnames])
    background.stats[[make.names(paste(neg.colnames[2],"hat",sep="."))]] <- cbind(1,prot[,neg.colnames[1]]) %*% dmod$coefficients
    background.stats[["UL"]] <- background.stats[[make.names(paste(neg.colnames[2],"hat",sep="."))]] +
      qnorm(.975)* sqrt(var.on.line * sin(alpha))
    background.stats[["LL"]] <- background.stats[[make.names(paste(neg.colnames[2],"hat",sep="."))]] -
      qnorm(.975)* sqrt(var.on.line * sin(alpha))
    background.stats[["out"]] <- prot[,neg.colnames[2]]<background.stats[["LL"]]| prot[,neg.colnames[2]]>background.stats[["UL"]]
    
    
    # Note: currently the estimate of the variability in background, sigma, 
    #       is taken wihtout excluding the outliers. The plot below shows the
    #       difference in trend after excluding the outliers (for eval purposes)
    
    if(visualize){
      par(las = 1)
      plot(prot[,neg.colnames[1]],prot[,neg.colnames[2]],col = 1+background.stats[["out"]],pch = c(1,19)[1+background.stats[["out"]]],
           main = "Deming regression of negatives", xlab =neg.colnames[1], ylab =neg.colnames[2])
      abline(coef = dmod$coefficients,lty=1,col=1)
      abline(coef = deming(prot[!background.stats[["out"]],neg.colnames[2]]~prot[!background.stats[["out"]],neg.colnames[1]],stdpat=c(1,0,1,0))$coefficients,lty=2,col=2)
      abline(coef = coef(lm(prot[,neg.colnames[2]]~prot[,neg.colnames[1]])),lty=3,col=3)
      abline(coef = c(mean(prot[,neg.colnames[2]]-prot[,neg.colnames[1]]),1),lty=4,col=4)
      legend("bottomright",legend = c("Deming","Deming w/o outliers","OLS","Best fit slope=1"),col = 1:4,lty = 1:4)
      
    }
    
    orth.proj <- get.orth.proj(y = prot[,neg.colnames[2]],x = prot[,neg.colnames[1]],mod = dmod)
    
    # Here we get a representation of the "mean" of the two negatives along with the CI around this est.
    background.stats[["orth.proj.bar"]] <- apply(orth.proj,1,mean) 
    background.stats[["opbar.UL"]] <- background.stats[["orth.proj.bar"]] + qnorm(.975)* sqrt( (var.on.line * sin(alpha)+ var.on.line * cos(alpha))/2 )
    background.stats[["opbar.LL"]] <- background.stats[["orth.proj.bar"]] + qnorm(.025)* sqrt( (var.on.line * sin(alpha)+ var.on.line * cos(alpha))/2 )
    
    # The upper limit of the orthogonal projection is chosen as the threshold for background
    background.stats[["background.threshold"]] <- background.stats[["opbar.UL"]]
    
  } else{
    
    background.stats <- data.frame(prot[,neg.colnames])
    background.stats[["mean"]] <- apply(prot[,neg.colnames],MARGIN = 1,mean)
    background.stats[["sd"]] <- apply(prot[,neg.colnames],MARGIN = 1,sd)
    background.stats[["background.threshold"]] <- background.stats[["mean"]] + qnorm(.975) * background.stats[["sd"]]
    
  }

  above.background <- sweep(prot,MARGIN = 1,STATS = background.stats[rownames(prot),"background.threshold"],FUN = ">")
  
  
  #=============================================
  # Normalize data
  #=============================================
  nchoos <- normmodule.nchoos$protein
  if(isTRUE(normmodule.nchoos$protein == "Dynamic"))
    nchoos <- 15
  
  # Assign candidate normalizing prots (15 for now)
  # Note: we don't exclude frequently background proteins
  #       from consideration here
  #------------------------------------------------
  ref.lane <- apply(prot,MARGIN = 2,FUN = function(x) median(x,na.rm=T))
  prot.RLN <- sweep(prot,MARGIN = 2,STATS = ref.lane,FUN = "-")
  
  #evaluate variability
  prot.RLN.MAD <- apply(abs(sweep(prot.RLN,MARGIN = 1,STATS = apply(prot.RLN,1,mean),FUN = "-" )),MARGIN = 2,FUN = function(x) mean(x,na.rm=T))
  sorted.prbs.by.MAD <- names(sort(prot.RLN.MAD))
  
  cand.names <- sorted.prbs.by.MAD[1:min(nchoos,ncol(prot.RLN.MAD))]
  if(!is.null(norm.ref.prot.names)){
    cand.names <- sorted.prbs.by.MAD[sorted.prbs.by.MAD %in% norm.ref.prot.names]
    cand.names <- cand.names[1:min(nchoos,length(cand.names))]
  }
    

  prot.cand <- prot[,cand.names]
  ref.cand.profile <- colSums(prot.cand,na.rm=T)/nrow(prot.cand)
  Mdel <- sweep(prot.cand,MARGIN = 2,STATS = ref.cand.profile,FUN = "-")
  Mresid<- t(scale(t(Mdel),center = T,scale = F))
  
  
  norm.factor <- data.frame(matrix(nrow = nrow(prot.cand),ncol=4,dimnames = list(rownames(prot.cand),c("user.defined","geomean.all","geomean.stable1","geomean.stable2") )))
  
  # if (!is.null(norm.ref.prot.names)){
  #   normalizers <- norm.ref.prot.names
  #   norm.factor$user.define <- apply(prot[,normalizers],1,mean)
  #   n.prot <- sweep(prot,MARGIN = 1,STATS = norm.factor$user.defined,FUN = "-")
  # }
    
  
  if(norm.method == "geomean.all"){
    normalizers <- colnames(prot)
    norm.factor$geomean.all <- apply(prot[,normalizers],1,mean)
    n.prot <- sweep(prot,MARGIN = 1,STATS =  norm.factor$geomean.all,FUN = "-")
  }
    
  
  if(norm.method == "geomean.stable1"){
    normalizers <- colnames(prot.cand)
    norm.factor$geomean.stable1 <- attr(Mresid,"scaled:center")
    n.prot <- sweep(prot,MARGIN = 1,STATS =  norm.factor$geomean.stable1,FUN = "-")
  }
  
  #reformat normalizers
  #--------------------
  by.lane.normalizers <- matrix(TRUE,nrow=nrow(prot),ncol=length(normalizers),dimnames = list(rownames(prot),normalizers))
  

  
  # Identify cases where the fit might be improved:
  # Note: right now it just uses Mclust to find the cases where we end up with a group of at least 6 member
  #       we could use at other metrics (e.g. R2) to evaluate fit improvement
  
  if(norm.method == "geomean.stable2"){
    require(cluster)
    require(mclust)
    Mclst <- apply(Mresid,MARGIN = 1,FUN = function(x) Mclust(x))
    Mtb <- unlist(lapply(Mclst,function(x) max(table(x$classification))> 5 & (x$G)>1))
    by.lane.normalizers <- matrix(TRUE,nrow=nrow(prot),ncol=ncol(Mresid),dimnames = list(rownames(prot),colnames(Mresid)))
    
    #Initialize the norm.factor
    norm.factor$geomean.stable2 <- apply(prot[,colnames(Mresid)],1,mean)
    
    #Look at the individual clusters and for appropriate lanes get the adjust
    #normalization factor based on a smaller set of noramlizers.
    for(lane.i in names(Mtb)[Mtb]){
      
      del.i <- Mdel[lane.i,]
      clst.i <- Mclst[[lane.i]]
      tb.i <- table(clst.i$classification)
      grps.i <- which(tb.i==max(tb.i))
      normalizers.i <- names(which(clst.i$classification == grps.i))
      
      
      norm.factor[lane.i,"geomean.stable2"] <- mean(del.i[normalizers.i])
      by.lane.normalizers[ lane.i,!colnames(by.lane.normalizers) %in% normalizers.i] <- FALSE
      
      
      if(visualize){
        par(las = 1)
        plot(ref.cand.profile, prot.cand[lane.i,],col=clst.i$classification, main = lane.i, ylim = range(prot.cand))
        abline(mean(del.i),1,lty=1)
        abline(norm.factor[lane.i,"geomean.stable2"],1,col = which(tb.i==max(tb.i)),lty = 2)
      }
      
    }
    
    n.prot <- sweep(prot,MARGIN = 1,STATS = norm.factor$geomean.stable2,FUN = "-")
  }

  out.prot <- list(n.prot = n.prot,
                   above.background = above.background,
                   MAD = sort(prot.RLN.MAD),
                   background.stats= background.stats,
                   norm.factors = norm.factor[,sapply(norm.factor,function(x) any(!is.na(x)))],
                   by.lane.normalizers = by.lane.normalizers)
  

  if(visualize){
    used.in.normalization <- names(which(apply(by.lane.normalizers,2,any)))
    par(las = 1)
    plot(out.prot$MAD, main = "Mean Absolute Deviance by protein",ylab = "MAD",
         pch = 19,
         col = c("gray","black")[1+names(out.prot$MAD) %in% used.in.normalization])
    legend("topleft",legend = c("selected","unselected"),pch=19,col=c("black","gray"))
    
  }
  
  


    
    
  if(TRUE)   #<------------- implement when we've got the right values for the ??'s below
  {
    # Write protein normalization summary
    prot.norm.summary <- data.frame(protein = prb.annots[colnames(by.lane.normalizers),"Probe.Label"],
                                    rank = match(colnames(by.lane.normalizers),names(out.prot$MAD)),
                                    "SD after normalization"= signif(apply(out.prot$n.prot[,colnames(by.lane.normalizers)],2,sd),digits = 3),check.names = F)
    
    if(any(!by.lane.normalizers))
      prot.norm.summary <- data.frame(protein = prb.annots[colnames(by.lane.normalizers),"Probe.Label"],
                                      "Used in normalizing samples" = unlist(lapply(as.data.frame(by.lane.normalizers),function(x) paste(rownames(by.lane.normalizers)[x],collapse = ";"))),
                                      rank = match(colnames(by.lane.normalizers),names(out.prot$MAD)),
                                      "SD after normalization"= signif(apply(out.prot$n.prot[,colnames(by.lane.normalizers)],2,sd),digits = 3),check.names = F)
    
    write.csv(prot.norm.summary,file = paste(path.to.normalization.results,"//proteins used in normalization.csv",sep=""))
    
    
    
    # plot the normalization results:
    used.in.normalization <- names(which(apply(by.lane.normalizers,2,any)))
    tmp <- 1*out.prot$above.background * out.prot$n.prot
    colnames(tmp) <- prb.annots[colnames(tmp), "Probe.Label"]
    protein.measured.g <- visualize.matrix(x = tmp,xlabel = "proteins",ylabel = "Samples",low.col = "white",high.col = "orange")
    protein.measured.g <- protein.measured.g + ggtitle("Protein Expression\n Background Thresholded")
   
    for(r in 1:length(plottypearg)){
    plottype=plottypearg[r]
    tempfilename = drawplot(filename=paste(path.to.normalization.results,"//protein normalization results",sep=""),plottype)
    tempfilename=gsub(path.results,"results",tempfilename)
    normalization.summary.plot(HKnormalized=out.prot$n.prot[,cand.names],
                               norm.factors=out.prot$norm.factors,
                               analyte.type="protein",col="darkgoldenrod2",main = "Normalization summary: Protein",cex=2,pch=17)
    dev.off()
    
    tempfilename = drawplot(filename=paste(path.to.normalization.results,"//protein normalizer selection",sep=""),plottype)
    tempfilename=gsub(path.results,"results",tempfilename)
    
    par(las = 1)
    plot(out.prot$MAD, main = "Mean Absolute Deviance by protein",ylab = "MAD",
         pch = 19,
         col = c("gray","black")[1+names(out.prot$MAD) %in% used.in.normalization])
    legend("topleft",legend = c("selected","unselected"),pch=19,col=c("black","gray"))
    dev.off()
    
    # plot the antibodies above background
    filename<-paste(path.to.normalization.results,"//protein expression thresholded",sep="")
    ggsave(filename = paste(filename,plottype,sep="."), plot = protein.measured.g,width = 8,height = 6,scale = 1,dpi = 250,units = "in")
    
    }
  }
    
	print("Normalized protein data")
	cat("LOG:Normalized protein data",file=log,sep='\n\n',append=TRUE)
	cat("document.write('<p>Normalized protein data</p>');", file=paste(path.inc,"//status.js",sep=""),append=TRUE)
	cat("document.write('  	  	<ul>');\n",file=paste(path.inc,"//panel2_5.js",sep=""),append=TRUE)

	# threshold normalized data at 1:
	# normalized = replace(normalized,normalized<0,0)
	# HKnormalized = replace(HKnormalized,HKnormalized<0,0)
	# out = list(normalized=normalized,HKnormalized=HKnormalized,HKs=HKs,norm.factors=normalizedall$norm.factors)

	# write csv of normalized protein data:
	relabeled.n.prot <- n.prot
	colnames(relabeled.n.prot) <- paste(prb.annots[colnames(n.prot),"Probe.Label"],prb.annots[colnames(n.prot),"Analyte.Type"],sep="-")
	#write.csv(relabeled.n.prot,file=paste(path.to.normalization.results,"//Protein_normalized_log2_data.csv",sep=""))
        #write.csv(2^relabeled.n.prot,file=paste(path.to.normalization.results,"//Protein_normalized_linear_data.csv",sep=""))
  
        log.temp.df2 <- data.frame("log" = c("Protein - normalized log2 count data"))
        write.table(log.temp.df2,file=paste(path.to.normalization.results,"//Protein_normalized_data.csv",sep=""), col.names = F, row.names = F)
        write.table(relabeled.n.prot,file=paste(path.to.normalization.results,"//Protein_normalized_data.csv",sep=""), append = T, sep = ",", col.names = NA, row.names = T)

        linear.temp.df2 <- data.frame("linear" = c("", "Protein - normalized linear count data"))
        write.table(linear.temp.df2,file=paste(path.to.normalization.results,"//Protein_normalized_data.csv",sep=""), append = T, sep = ",", col.names = F, row.names = F)
        write.table(2^relabeled.n.prot,file=paste(path.to.normalization.results,"//Protein_normalized_data.csv",sep=""), append = T, sep = ",", col.names = NA, row.names = T)

	print("Created normalized protein data files")
	cat("LOG:Created normalized protein data files",file=log,sep='\n\n',append=TRUE)
	cat("document.write('<p>Created normalized protein data files</p>');", file=paste(path.inc,"//status.js",sep=""),append=TRUE)


	names(out.prot)[which(names(out.prot)=="n.prot")] <- "normalized.data"
	return(out.prot)

}




# ==========================================
#         normalization util
#===========================================
# This function gets a vector of of specified probes for normalization
# and returns the probes for the analyte specified. If analyte type doesn't exist
# or the vector of specified probes is empty (NULL), it returns NULL
# @probes: character vector specifying probes ids to be used in normalization
# @probe.annotation.id.column: the column name in probe annotation specifying the unique IDs
# @pannot: prbe annotation data.frame
# @analyte.type: character specifying the analyte type of intrest (e.g. "mRNA" or "protein")
# NOTE: pannot is expected to have a column "Analyte.Type" try to fix this to be less of a magic word!

get.probes.for.analyte <- function(probes,
                                   probe.annotation.id.column ="UniqueID",
                                   pannot,
                                   analyte.type){
  #selected probes for normalization ==> redefine what's HK and how many to choose
  probes.filtered.by.analyte <- NULL
  if(!is.null(probes)){
    selected.prbs.indices <- match(probes,as.character(pannot[,probe.annotation.id.column]))
    selected.prbs.indices <- selected.prbs.indices[pannot$Analyte.Type[selected.prbs.indices] %in% analyte.type]
    if(length(selected.prbs.indices)>0)
      probes.filtered.by.analyte <- pannot[selected.prbs.indices,probe.annotation.id.column]
  }
  return(probes.filtered.by.analyte)
}


#' This function gets the normalization options based on the input provided by the wizard
#' If mRNA normalization probes are provided, it uses them to normalize. If not, it tries to use HK
#' in the annotations, it tries to find the best candidate HK. Subsequently based on the choice
#' of normmodule.nchoos (i.e. numeric value or "Dynamic"), it sets the auto.HKs option
#' @param   normmodule.norm.probes: character vector specifying probes ids to be used in normalization
#' @param normalize T/F 
#' @param normmodule.nchoos: Either charcter "Dynamic" or a numeric value specifying how selection of normalization candidate are to be done
#' @param probe.annotation.id.column: the column name in probe annotation specifying the unique IDs
#' @param pannot: prbe annotation data.frame
#' NOTE: here we are forcing ERCCs to be excluded from consideration as normalizers per team's decision
#' It returns the following values which are normalization options. see normalizerawdata function above:
#' @return  codeset.HKs: a character vector of the probe Ids to be used in normalization
#' @return  auto.HKs: TRUE/FALSE, when TRUE uses geNorm ranking to find the optimal number housekeepers.
#' @return  n.HKs: a number to determine the number of normalizers to be picked. Ignored when auto.HKs is TRUE

get.mRNA.normalization.opt <- function(raw.dat,
                                       normalize,
                                       normmodule.norm.probes,
                                       normmodule.nchoos,
                                       probe.annotation.id.column,
                                       pannot){

  codeset.HKs  <- get.probes.for.analyte (probes = normmodule.norm.probes,
                                          probe.annotation.id.column = probe.annotation.id.column,
                                          pannot = pannot,
                                          analyte.type = "mRNA")
  
  
  # Auto option: If for whatever reason codeset.HKs is empty we will try to use HK for normalization 
  # (the reason are: the user chose auto method or failed to correctly specify)
  #--------------------------------------------------------------------------------------------------
  if(is.null(codeset.HKs)){
    codeset.HKs <- as.character(pannot[which(pannot$Control.Type  == "Housekeeping"),probe.annotation.id.column])
    
    # If HK are not found in the annotations find candidate HKs using MAD
    if(length(codeset.HKs)==0 ){
      
      #probes not to be considered as normalizers
      exclude.as.normalizers <- pannot[pannot$Control.Type %in% c("Negative","Positive"),probe.annotation.id.column]
      
      tmp <- raw.dat[,!colnames(raw.dat) %in% exclude.as.normalizers,drop=F]
      if(2>=ncol(tmp))
        stop("Not including control mRNA probes, there are no more probes left to be used to select good candidate normalizer probes from")
      
      # get probes 90% of the times above 100 counts
      not.lowexpression <- apply(tmp,2,function(x) mean(x>100) >.9)
      
      tmp <- log2(tmp+1)
      tmp <- sweep(tmp,2,STATS = apply(tmp,2,median),FUN = "-")
      tmp <- sweep(tmp,1,STATS = apply(tmp,1,mean),FUN = "-")
      MAD <- apply(abs(tmp),2,mean)
      #     plot(MAD,log="y",col=1+names(MAD)%in% rownames(pannot)[pannot$Control.Type == "Housekeeping"], pch = 1+ not.lowexpression)
      not.lowexpression <- not.lowexpression[order(MAD,decreasing = F)]
      MAD <- MAD[order(MAD,decreasing = F)]
      
      
      num.candidates.to.select <- min(20, max(10,round(.2*length(MAD))))
      len.to.cover <- which(cumsum(not.lowexpression)<=num.candidates.to.select & num.candidates.to.select<= cumsum(not.lowexpression))[1]
      
      # if couldn't find the num.candidate.to.select non lowexpression probes iteratively reduce the num.candidates.to.select to at least 10
      while(10<=num.candidates.to.select & is.na(len.to.cover)){
        
        num.candidates.to.select <- num.candidates.to.select - 1
        len.to.cover <- which(cumsum(not.lowexpression)<=num.candidates.to.select & num.candidates.to.select<= cumsum(not.lowexpression))[1]
        
        # if couldn't find at least 10 non lowexpression probes just go with the best 10 by MAD criteria
        if(num.candidates.to.select<=10){
          len.to.cover <- min(10,length(not.lowexpression))
          not.lowexpression <- not.lowexpression == not.lowexpression
        }
        
      }
      
      codeset.HKs <- pannot[names(not.lowexpression[0:len.to.cover]),probe.annotation.id.column]
    }
    
  }

  # Manual option: User selected probes for normalization 
  # ==> redefine what's HK and how many to choose
  #----------------------------------------------------------------------------------------------
  exclude.as.normalizers <- pannot[pannot$Control.Type %in% c("Negative","Positive"),probe.annotation.id.column]
  non.HK <- setdiff(colnames(raw.dat)[!colnames(raw.dat) %in% exclude.as.normalizers],codeset.HKs)

  # Default
  n.HKs <- NA
  refine.selected.HKs <- FALSE
  HK.selection.method <-  "none"
  if(!normalize | length(non.HK)==0)
    return(list(codeset.HKs=codeset.HKs,
                auto.HKs= refine.selected.HKs,
                n.HKs = n.HKs,
                HK.selection.method = HK.selection.method))
    
  
  HK.selection.method <- "choose HKs"
  refine.selected.HKs <- length(codeset.HKs)>5
  
  if(is.numeric(normmodule.nchoos$mRNA)){
    n.HKs <- min(normmodule.nchoos$mRNA, length(codeset.HKs))
    if(n.HKs == length(codeset.HKs))
      refine.selected.HKs <- FALSE
  }
  
  
  if(!refine.selected.HKs){
    HK.selection.method <- "use selected"
    n.HKs <- NA
  }

    

  
  return(list(codeset.HKs=codeset.HKs,
              auto.HKs= refine.selected.HKs & isTRUE(normmodule.nchoos$mRNA == "Dynamic"),
              n.HKs = n.HKs,
              HK.selection.method = HK.selection.method))
}
