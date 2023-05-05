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
# colData(se)$lacStageFac <- factor(colData(se)$lacStage)
se$lacStageFac <- factor(colData(se)$lacStage,
                                  levels = c("Colostrum", "Transitional", "Mature"),
                                  labels = c(1, 2, 3))

#colData(se)$protocolFac <- factor(colData(se)$protocol)
se$protocolFac <- factor(colData(se)$protocol,
                                  levels = c("Soft spin, Unwashed", "Hard spin, Unwashed", "Hard spin, Washed once", "Hard spin, Washed twice"),
                                  labels = c(1, 2, 3, 4 ))


# Create logical mask excluding column containing "SRR801705" in the name
mask <- !grepl("SRR801705", colnames(se))

# Subset the data frame using the mask
se_sample_filtered <- se[, mask]

# The same for the samples in the dge object
# in the counts table



# in the samples table
mask <- rownames(dge$samples) != "SRR801705"
dge_sample_filtered <- dge
dge$samples_masked <- dge$samples[mask, ]
dge_sample_filtered$samples <- dge$samples_masked


mask <- !grepl("SRR801705", colnames(dge$counts))
dge_sample_filtered$counts <- dge_sample_filtered$counts[, mask]


## ----exprdist, echo=FALSE, out.width="600px", fig.cap="Distribution of average expression level per gene."----
avgexp <- rowMeans(assays(se_sample_filtered)$logCPM)
hist(avgexp, xlab="log2 CPM", main="", las=1, ylim=c(0,7000))

## -----------------------------------------------------------------------------

mask <- filterByExpr(dge, group=colData(se_sample_filtered)$protocolFac)
se.filt <- se_sample_filtered[mask, ]
dim(se.filt)
dge.filt <- dge_sample_filtered[mask, ]
dim(dge.filt)

## -----------------------------------------------------------------------------
par(mar=c(4, 5, 1, 1))

h <- hist(avgexp, xlab=expression("Expression level (" * log[2] * "CPM)"), 
          main="", las=1, col="grey", cex.axis=1.2, cex.lab=1.5)

x <- cut(rowMeans(assays(se.filt)$logCPM), breaks=h$breaks)

lines(h$mids, table(x), type="h", lwd=10, lend=1, col="darkred")

legend("topright", c("All genes", "Filtered genes"), fill=c("grey", "darkred"))


## -----------------------------------------------------------------------------
dge.filt <- calcNormFactors(dge.filt)

## -----------------------------------------------------------------------------
assays(se.filt)$logCPM <- cpm(dge.filt, log=TRUE,
                              normalized.lib.sizes=TRUE)

## ----maPlots, fig.height=18, fig.width=10, dpi=100, echo=FALSE, fig.cap="MA-plots of filtered and normalized expression values."----
par(mfrow=c(5, 4), mar=c(4, 5, 3, 1))
for (i in 1:ncol(se.filt)) {
  A <- rowMeans(assays(se.filt)$logCPM)
  M <- assays(se.filt)$logCPM[, i] - A
  smoothScatter(A, M, main=colnames(se.filt)[i], las=1)
  abline(h=0, col="blue", lwd=2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col="red", lwd=2)
}

## -----------------------------------------------------------------------------
table(se.filt$idFac, se.filt$protocolFac)

## -----------------------------------------------------------------------------
table(se.filt$lacStageFac, se.filt$protocolFac)

## ----sampleClustering, fig.height=5, fig.width=8, dpi=100, echo=FALSE, fig.cap="Figure S6: Hierarchical clustering of the samples. Labels correspond to treatment and sample identifer, while colors indicate sample group."----
par(mar=c(8, 5, 1, 1))
logCPM <- cpm(dge.filt, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(se.filt$protocolFac)
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(se.filt)
outcome <- paste(se.filt$lacStageFac, colnames(se), sep="\n")
names(outcome) <- colnames(se.filt)

sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)

plot(sampleDendrogram, main="Hierarchical clustering of samples",
     cex=0.7)

legend("topright", levels(se.filt$protocolFac),
       fill=seq_len(nlevels(se.filt$protocolFac)))

## ----mdsPlot, fig.height=5, fig.width=8, dpi=100, echo=FALSE, fig.cap="Figure S7: Multidimensional scaling plot of the samples. Labels correspond to treatment and colors indicate sample group."----
outcome <- se.filt$lacStage
names(outcome) <- colnames(se.filt)
plotMDS(dge.filt, labels=outcome, col=batch)
legend("bottomleft", levels(se.filt$protocolFac),
       fill=seq_len(nlevels(se.filt$protocolFac)), inset=0.05)

