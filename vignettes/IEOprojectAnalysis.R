## ----setup, echo=FALSE, cache=FALSE-------------------------------------------
library(knitr) ## kable()
library(kableExtra) ## kable_styling(), save_kable()
library(usethis) ## use_directory(), proj_path()

knitr::opts_chunk$set(
  collapse=TRUE,
  comment="",
  fig.align="center",
  cache=FALSE
)

## this option avoid use_directory() being verbose later on
options(usethis.quiet=TRUE)

## ----message=FALSE------------------------------------------------------------
library(SummarizedExperiment)
se <- readRDS(file.path(system.file("extdata", package="IEOproject"), "GSE45669.rds"))
#se <- readRDS("../inst/extdata/GSE45669.rds")

se


## -----------------------------------------------------------------------------
head(rowData(se))

## -----------------------------------------------------------------------------
dim(colData(se))
head(colData(se), n=3)

## ---- message=FALSE-----------------------------------------------------------
library(IEOproject)
library(edgeR)

dge <- DGEList(counts=assays(se)$counts, genes=rowData(se))
dim(dge)

## -----------------------------------------------------------------------------
assays(se)$logCPM <- cpm(dge, log=TRUE)
assays(se)$logCPM[1:5, 1:5]

## -----------------------------------------------------------------------------

# loop over the columns and count unique values
for (col in names(colData(se))) {
  unique_vals <- unique(colData(se)[[col]])
  if (any(length(unique_vals) %in% c(1, 16))) {
    next # skip to the next column
  } else {
    counts <- table(colData(se)[[col]])
    cat(paste0("Column '", col, "' has ", length(unique_vals), " unique values:\n"))
    print(counts)
  }
}


## -----------------------------------------------------------------------------
# simpler names
names(colData(se))[37] <- "id"
names(colData(se))[36] <- "lacStage"
names(colData(se))[35] <- "protocol"

# creating factors
colData(se)$idFac <- factor(colData(se)$id)
# colData(se)$lacStageFac <- factor(colData(se)$lacStage)
colData(se)$lacStageFac <- factor(colData(se)$lacStage,
                                  levels = c("Colostrum", "Transitional", "Mature"),
                                  labels = c(1, 2, 3))

#colData(se)$protocolFac <- factor(colData(se)$protocol)
colData(se)$protocolFac <- factor(colData(se)$protocol,
                                  levels = c("Soft spin, Unwashed", "Hard spin, Unwashed", "Hard spin, Washed once", "Hard spin, Washed twice"),
                                  labels = c(1, 2, 3, 4 ))


## ----pheno, echo=FALSE, message=FALSE-----------------------------------------
tmpdf <- data.frame("Patient"=colData(se)$id,
                    "lacStage"=colData(se)$lacStage,
                    "protocol"=colData(se)$protocol,
                    check.names=FALSE)
ktab <- kable(tmpdf, caption="Phenotypic variables.")
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
par(mar=c(7, 5, 2, 2))
ord <- order(dge$sample$lib.size/1e6)
ordmreads <- dge$sample$lib.size[ord]/1e6
names(ordmreads) <- colnames(se)[ord]
bp <- barplot(ordmreads, las=1, ylab="Millions of reads",
              xlab="", col=c("blue", "green", "red")[colData(se)$lacStageFac[ord]], las=2, ylim = c(0, 50))
legend("topleft", c("Colostrum", "Transitional", "Mature"), fill=c("blue", "green", "red"), inset=0.01)

## -----------------------------------------------------------------------------
par(mar=c(7, 5, 2, 2))
ord <- order(dge$sample$lib.size/1e6)
ordmreads <- dge$sample$lib.size[ord]/1e6
names(ordmreads) <- colnames(se)[ord]
bp <- barplot(ordmreads, las=1, ylab="Millions of reads",
              xlab="", col=c("blue", "green", "red", "orange")[colData(se)$protocolFac[ord]], las=2, ylim = c(0, 60))
legend("topleft", c("Soft spin, Unwashed", "Hard spin, Unwashed", "Hard spin, Washed once", "Hard spin, Washed twice"), 
       fill=c("blue", "green", "red", "orange"), inset=0.01)

## ----distRawExp, echo=FALSE, fig.height=5, fig.width=5, out.width="600px", fig.cap="Non-parametric density distribution of expression profiles per sample.", message=FALSE----
library(geneplotter)
par(mar=c(4, 5, 1, 1))
lst <- as.list(as.data.frame(assays(se)$logCPM))
multidensity(lst, xlab="log 2 CPM", legend=NULL,
             main="", las=1)

## -----------------------------------------------------------------------------
par(mar=c(7, 5, 2, 2))
boxplot(assays(se)$logCPM, col="gray", ylab=expression(log[2] * "CPM"),
cex.axis=1.2, cex.lab=1.5, las=2)


## -----------------------------------------------------------------------------
# Create logical mask excluding column containing "SRR801705" in the name
mask <- !grepl("SRR801705", colnames(se))

# Subset the data frame using the mask
se_sample_filtered <- se[, mask]

# The same for the samples in the dge object
mask <- rownames(dge$samples) != "SRR801705"
dge_sample_filtered <- dge
dge$samples_masked <- dge$samples[mask, ]
dge_sample_filtered$samples <- dge$samples_masked

## ----exprdist, echo=FALSE, out.width="600px", fig.cap="Distribution of average expression level per gene."----
avgexp <- rowMeans(assays(se_sample_filtered)$logCPM)
hist(avgexp, xlab="log2 CPM", main="", las=1, ylim=c(0,7000))

