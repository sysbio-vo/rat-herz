library(limma)
source("degs_utils.R")
library(rta10rnentrezg.db)


pd.plain <- read.table("../pdata/pdata.tsv", header = TRUE, sep = "\t")
exprs <- read.table("../exprs/exprs.tsv", header=TRUE, check.names=FALSE)
exl = which(pd.plain$Outlier==TRUE)
pd.plain <- pd.plain[-exl,]
exprs <- exprs[,-exl]

probesetsID <- rownames(exprs)
probesetsID_EntrezID <- select(rta10rnentrezg.db, probesetsID, "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(!is.na(probesetsID_EntrezID$ENTREZID)),]


exprs <- exprs[which(rownames(exprs) %in% probesetsID_EntrezID$PROBEID),]

rownames(probesetsID_EntrezID) <- probesetsID_EntrezID$PROBEID
  
exprs <- merge(exprs, probesetsID_EntrezID, by="row.names")
rownames(exprs) <- exprs$ENTREZID
exprs <- exprs[, -c(1, 17, 18)]

degs <- getDEGS(c("SD_d21", "PE_d21"), pd.plain, exprs)
degs <- filterDEGS(degs, 0.05, 0)
write.table(degs, "../degs/SDd21_PEd21.tsv", sep="\t", row.names = FALSE, quote=FALSE)

degs <- getDEGS(c("PE_NP", "PE_d21"), pd.plain, exprs)
degs <- filterDEGS(degs, 0.05, 0)
write.table(degs, "../degs/PENP_PEd21.tsv", sep="\t", row.names = FALSE, quote=FALSE)

degs <- getDEGS(c("SD_NP", "SD_d21"), pd.plain, exprs)
degs <- filterDEGS(degs, 0.05, 0)

degs <- getDEGS(c("SD_NP", "PE_NP"), pd.plain, exprs)
degs <- filterDEGS(degs, 0.05, 0)
