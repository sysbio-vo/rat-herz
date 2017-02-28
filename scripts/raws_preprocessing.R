library(affy)
library(arrayQualityMetrics)
#install.brainarray("rta10", orgnsm = "rn", force.reinstall = FALSE)
library(rta10rnentrezgcdf)
library(rta10rnentrezg.db)
library(rta10rnentrezgprobe)
source("plots_utils.R")

pd <- read.AnnotatedDataFrame("../pdata/pdata.tsv")
pd.plain <- read.table("../pdata/pdata.tsv", header = TRUE, sep = "\t")
affyData = ReadAffy(phenoData=pd, sampleNames=pd$Sample_ID, filenames=as.character(rownames(pd)),
                    celfile.path="../raws/")

affyData@cdfName <- "rta10rnentrezgcdf"
eset = rma(affyData)

arrayQualityMetrics(expressionset = eset,
                    outdir = "../plots/AQM/AQM_report",
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Condition"))

write.table(exprs(eset), "../exprs/exprs.tsv", sep="\t", quote=FALSE)

exprs <- exprs(eset)
exl = which(pd.plain$Outlier==TRUE)
pd.plain <- pd.plain[-exl,]
exprs <- exprs[,-exl]

pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pd.plain, c("Pregnancy", "Phenotype", "Date", "Condition"), ncol=2)
save_plot("../plots/PCA/pca.pdf", base_height=6, base_aspect_ratio = pl[[2]]/2, pl[[1]])
save_plot("../plots/PCA/pca.svg", base_height=6, base_aspect_ratio = pl[[2]]/2, pl[[1]])