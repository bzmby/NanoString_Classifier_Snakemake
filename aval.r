############3#to get the  merged data###############

PrepData = function(infile){
	data <- read.table(infile, header=TRUE, as.is = TRUE, sep = ",")
 	data = data[, 2:ncol(data)]
		 merged.data = data
 	colnames(merged.data[1]) < "CodeCount"
	rownames(merged.data) <- merged.data$Name
 	x <- list(counts = merged.data, raw.counts = merged.data)
 return(x)
}

##################################
geo.mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#############################
#' Normalization part
Normalization_NanoS = function(nsdata, method){
  # Get counts and number of samples
 		 counts <- nsdata$counts
 		 no.cols <- length(colnames(counts))
 			 lane.norm.factors <- c()
 			 endo.counts <- counts[grepl("Endogenous", counts$CodeClass), 4:no.cols]
	  names <- counts$Name[grepl("Endogenous", counts$CodeClass)]

  # Housekeeping normalization
  if (method == "housekeeping"){
   	 hk.counts <- counts[counts$CodeClass == "Housekeeping",4:no.cols]
   	 mean.hk <- apply(
   		   X = hk.counts,
   		   MARGIN = 2,
    	  FUN = geo.mean
    )
   	 global.mean.hk <- mean(mean.hk)
    lane.norm.factors <- global.mean.hk / mean.hk
  }

  # multiply with factors
  x <- t(apply(
    X = counts[,4:no.cols],
  	  MARGIN = 1,
  		  FUN = '*',
  			  lane.norm.factors
  ))
  	  		  counts <- cbind(counts[,1:3],x)
  	  nsdata$counts <- counts
  	  nsdata$norm.factors <- lane.norm.factors
  return(nsdata)
}



#' Positive Control Normalization
#' Normalize using the positive controls (spike-ins)
Normalization_positive_NanoS <- function(nsdata, pcm = "geo.mean"){
      counts <- nsdata$counts
      no.cols <- length(colnames(counts))
      pos.counts <- counts[counts$CodeClass == "Positive",4:no.cols]
  		pos.means.samples <- apply(
    		X = pos.counts,
      MARGIN = 2,
      FUN = pcm
  )
     pos.mean <- mean(pos.means.samples)
     pos.factors <- pos.mean / pos.means.samples
  # multiply with factors
  x <- t(apply(
        X = counts[,4:no.cols],
        MARGIN = 1,
            FUN = '*',
            pos.factors
  ))
            counts <- cbind(counts[,1:3],x)
        nsdata$counts <- counts
        nsdata$pos.factors <- pos.factors
  return(nsdata)
}

#' Background Correction
# Currently: use the mean + 2sd as background and subtract that from every count value
CorrectBackground_NanoS = function(nsdata, bm = "mean", sd.factor = 2){
      counts <- nsdata$counts
      negative.counts <- counts[counts$CodeClass == "Negative", ]
      no.cols <- length(colnames(counts))
      bg_values <- c()
          for (i in 4:no.cols){
              lane <- negative.counts[,i]
          neg.sd <- sd(lane)
          neg.mean <- mean(lane)
          neg.value <- neg.mean + (sd.factor*neg.sd)
          bg_values <- c(bg_values, neg.value)
  }

     x <- t(apply(
         X = counts[,4:no.cols],
         FUN = '-',
         MARGIN = 1,
     bg_values
  ))
      x[x < 0] <- 1
      nsdata$raw.counts <- nsdata$counts
          nsdata$counts <- cbind(counts[,1:3],x)
          nsdata$bg.corr.counts <- nsdata$counts
      nsdata$background <- bg_values
  return(nsdata)
}



data <- PrepData('result15.csv')
data2 <- CorrectBackground_NanoS(data)
data3 <- Normalization_positive_NanoS(data2)
data4 <- Normalization_NanoS(data3, method="housekeeping")
write.table(data4$counts, 'Normalized_data.txt', sep = '\t')
