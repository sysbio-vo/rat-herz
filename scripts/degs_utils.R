getDEGS <- function(meta.vars, pheno.data, exprs, noannotate=FALSE) {
  # Subset groups for comparison
  ind <- c()
  for (i in meta.vars) {
    ind <- c(ind, which(pheno.data$Condition==i))
  }
  pdata.short <- pheno.data[ind, ]
  exprs.short <- exprs[,which(colnames(exprs) %in% pdata.short$Sample_ID)]
  exprs.short <- exprs.short[,order(match(colnames(exprs.short), pdata.short$Sample_ID))]
  
  # Create design matrix
  design = model.matrix(~factor(pdata.short$Condition, levels=meta.vars), data=pdata.short)
  colnames(design) <- c(meta.vars[1], paste(meta.vars[1], "vs", meta.vars[2], sep=""))
  # Fit with linear models
  fit <- lmFit(exprs.short, design)
  fit <- eBayes(fit)
  # Get all the genes with logFC, p-values, no filtering
  degs <- topTable(fit, coef=paste(meta.vars[1], "vs", meta.vars[2], sep=""), adjust.method="fdr", number=nrow(fit))

  # Merge degs with expression matrix
  exprs.degs <- merge(degs, exprs.short, by="row.names")
  if (noannotate) {
    rownames(exprs.degs) <- exprs.degs[, 1]
    exprs.degs <- exprs.degs[, -1]
  } else {
    colnames(exprs.degs)[1] <- "ENTREZID"
    # Add information about gene names
    EntrezID_Symbol<-select(org.Rn.eg.db, exprs.degs$ENTREZID, c("SYMBOL", "GENENAME"))
    exprs.degs <- cbind(EntrezID_Symbol, exprs.degs)
    exprs.degs <- exprs.degs[,-4]
  }
  
  return(exprs.degs)
}

filterDEGS <- function(degs, pval, fc, adj) {
  # Filter by p-values
  if (missing(adj)) {
    degs <- degs[degs$adj.P.Val < pval,]
  } else if (adj==FALSE) {
    degs <- degs[degs$P.Value < pval,]
  } else if (adj==TRUE) {
    degs <- degs[degs$adj.P.Val < pval,]
  }
  
  # Sort by logFC
  degs <- degs[order(abs(degs$logFC), decreasing = TRUE),]
  # Filter by logFC
  degs <- degs[abs(degs$logFC) > fc,]  
  return (degs)
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

geneBarPlot <- function(exprs, pdata, symbol) {
  pdata <- pdata[, c("Trimester", "Condition")]
  pdata <- pdata[which(pdata$Condition!="High risk"),]
  gName = select(org.Hs.eg.db, symbol, c("ENTREZID"), keytype = "SYMBOL")
  exprs <- t(exprs[which(rownames(exprs) == gName$ENTREZID),])
  
  gene <- merge(exprs, pdata, by="row.names", all=FALSE)
  rownames(gene) <- gene[, 1]
  gene <- gene[,-1]
  
  #smr <- describeBy(gene$`3952`, group = list(gene$Condition, gene$Trimester), digits=3, mat=TRUE)
  smr <- summarySE(gene, measurevar=gName$ENTREZID, groupvars=c("Trimester","Condition"))
  colnames(smr)[4] <- "meanLogFC"
  smr$Trimester <- factor(smr$Trimester)
  smr$Condition <- factor(smr$Condition, levels=c("Low risk", "Control", "Preeclampsia"))
  
  pl <- ggplot(smr, aes(x=Condition, y=meanLogFC, fill=Trimester)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=meanLogFC-se, ymax=meanLogFC+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    labs(title = symbol)
  return(pl)
}