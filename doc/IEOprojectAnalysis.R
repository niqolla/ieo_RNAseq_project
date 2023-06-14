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
assays(se)$logCPM <- cpm(dge, log=TRUE,  prior.count=0.25)
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

## ----libsizes, echo=FALSE-----------------------------------------------------
par(mar=c(7, 5, 2, 2))
ord <- order(dge$sample$lib.size/1e6)
ordmreads <- dge$sample$lib.size[ord]/1e6
names(ordmreads) <- colnames(se)[ord]
bp <- barplot(ordmreads, las=1, ylab="Millions of reads",
              xlab="", col=c("blue", "green", "red")[colData(se)$lacStageFac[ord]], las=2, ylim = c(0, 50))
legend("topleft", c("Colostrum", "Transitional", "Mature"), fill=c("blue", "green", "red"), inset=0.01, cex=0.85)

## ----libsizes_2, echo=FALSE---------------------------------------------------
par(mar=c(7, 5, 2, 2))
ord <- order(dge$sample$lib.size/1e6)
ordmreads <- dge$sample$lib.size[ord]/1e6
names(ordmreads) <- colnames(se)[ord]
bp <- barplot(ordmreads, las=1, ylab="Millions of reads",
              xlab="", col=c("blue", "green", "red", "orange")[colData(se)$protocolFac[ord]], las=2, ylim = c(0, 60))
legend("topleft", c("Soft spin, Unwashed", "Hard spin, Unwashed", "Hard spin, Washed once", "Hard spin, Washed twice"), 
       fill=c("blue", "green", "red", "orange"), inset=0.01, cex=0.85)

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


## ---- echo=FALSE--------------------------------------------------------------
# Create factors for lacStage and protocol
se$lacStageFac <- factor(colData(se)$lacStage,
                                  levels = c("Colostrum", "Transitional", "Mature"),
                                  labels = c(1, 2, 3))

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
hist(avgexp, xlab="log2 CPM", main="", las=1, ylim=c(0,8000))

## -----------------------------------------------------------------------------

#cpmcutoff <- round(10/min(dge_sample_filtered$sample$lib.size/1e6), digits=1)
cpmcutoff <- 0.1
print(cpmcutoff)

#nsamplescutoff <- min(table(se_sample_filtered$protocolFac))
nsamplescutoff<- 15 
print(nsamplescutoff)

mask <- rowSums(cpm(dge_sample_filtered) > cpmcutoff) >= nsamplescutoff

se.filt <- se_sample_filtered[mask, ]
dim(se.filt)
dge.filt <- dge_sample_filtered[mask, ]
dim(dge.filt)



par(mar=c(4, 5, 1, 1))

h <- hist(avgexp, xlab=expression("Expression level (" * log[2] * "CPM)"), 
          main="", las=1, col="grey", cex.axis=0.95, cex.lab=1.2)

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
  smoothScatter(A, M, main=colnames(se.filt)[i], las=1,cex.main=2.5)
  abline(h=0, col="blue", lwd=2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col="red", lwd=2)
}

## -----------------------------------------------------------------------------
table(se.filt$lacStageFac, se.filt$protocolFac)

## ----sampleClusteringTechnique, fig.height=5, fig.width=8, dpi=100, echo=FALSE, fig.cap="Hierarchical clustering of the samples by technical factor - protocol used."----
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

plot(sampleDendrogram, main="Hierarchical clustering by protocol used.",
     cex=2.0, cex.axis=1.9, cex.main=1.9)

legend("right", levels(se.filt$protocolFac),
       fill=seq_len(nlevels(se.filt$protocolFac)), 
       legend=c("Soft spin, Unwashed", "Hard spin, Unwashed", "Hard spin, Washed once", "Hard spin, Washed twice"))

## ----sampleClusteringLacStage, fig.height=5, fig.width=8, dpi=100, echo=FALSE, fig.cap="Hierarchical clustering of the samples by biological factor - lactation stage."----
par(mar=c(8, 5, 1, 1))
logCPM <- cpm(dge.filt, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(se.filt$lacStageFac)
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(se.filt)
outcome <- paste(se.filt$protocolFac, colnames(se), sep="\n")
names(outcome) <- colnames(se.filt)

sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)

plot(sampleDendrogram, main="Hierarchical clustering by lactation stage.")

legend("topright", levels(se.filt$lacStageFac),
       fill=seq_len(nlevels(se.filt$lacStageFac)),
       legend = c("Colostrum", "Transitional", "Mature"))

## -----------------------------------------------------------------------------
outcome <- se.filt$lacStage
batch <- as.integer(se.filt$protocolFac)
names(outcome) <- colnames(se.filt)
plotMDS(dge.filt, labels=outcome, col=batch)
legend("topright", levels(se.filt$protocolFac),
       fill=seq_len(nlevels(se.filt$protocolFac)), inset=0.05, 
       legend=c("Soft spin, Unwashed", "Hard spin, Unwashed", "Hard spin, Washed once", "Hard spin, Washed twice"))

## -----------------------------------------------------------------------------
Col_exp <- rowMeans(logCPM[, se_sample_filtered$lacStage=="Colostrum"])
Tra_exp <- rowMeans(logCPM[, se_sample_filtered$lacStage=="Transitional"])
Mat_exp <- rowMeans(logCPM[, se_sample_filtered$lacStage=="Mature"])

## ----TraColfig, fig.cap="A plot analysing the difference in expression between Colostral and Transitional stage."----
plot((Tra_exp+Col_exp)/2, Tra_exp-Col_exp, pch=".", cex=4, las=1)

## ----MatColfig, fig.cap="A plot analysing the difference in expression between Colostral and Mature stage."----
plot((Mat_exp+Col_exp)/2, Mat_exp-Col_exp, pch=".", cex=4, las=1)

## ----TraMatfig, fig.cap="A plot analysing the difference in expression between Transitional and Mature stage."----
plot((Mat_exp+Tra_exp)/2, Mat_exp-Tra_exp, pch=".", cex=4, las=1)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(sva)
library(limma)

## ----CThist, fig.cap="P-value diagnostic plot for differental expression between Colostrum and Transitional."----
se.filt.all <- se.filt[,se.filt$lacStageFac!=3]
se.filt.all$stage <- droplevels(se.filt.all$lacStageFac)

mod <- model.matrix(~ se.filt.all$stage, colData(se.filt.all))
mod0 <- model.matrix(~ 1, colData(se.filt.all))

pv <- f.pvalue(assays(se.filt.all)$logCPM, mod, mod0)
hist(pv, main="", las=1)

## -----------------------------------------------------------------------------
mask <- p.adjust(pv, method="fdr") < 0.1
DEgenesEGs <- names(pv)[mask]
DEgenesSyms <- mcols(se.filt)[DEgenesEGs, "symbol"]
DEgenesPvalue <- pv[mask]
DEgenesDesc <- mcols(se.filt)[DEgenesEGs, "description"]
DEgenesDesc <- sub(" \\[.+\\]", "", DEgenesDesc)
DEgenesTab <- data.frame(EntrezID=DEgenesEGs,
                         Symbol=DEgenesSyms,
                         Description=DEgenesDesc,
                         "P value"=DEgenesPvalue,
                         stringsAsFactors=FALSE, check.names=FALSE)
DEgenesTab <- DEgenesTab[order(DEgenesTab[["P value"]]), ] 

## ----CTvolc, fig.cap="P-value diagnostic plot for differental expression between Colostrum and Transitional."----
sv <- sva(assays(se.filt.all)$logCPM, mod=mod, mod0=mod0)
mod <- cbind(mod, sv$sv)
colnames(mod) <- c(colnames(mod)[1:2], paste0("SV", 1:sv$n))
fit3 <- lmFit(assays(se.filt.all)$logCPM, mod)
fit3 <- eBayes(fit3)
tt3 <- topTable(fit3, coef=2, n=Inf)

DEgenes_no3 <- rownames(tt3)[tt3$adj.P.Val < 0.1]
sort(table(tt3[DEgenes_no3, "chr"]), decreasing=TRUE)
res3 <- decideTests(fit3, p.value=0.1)

genesmd <- data.frame(chr=as.character(seqnames(rowRanges(se.filt.all))), 
                      symbol=rowData(se.filt.all)[, 5], stringsAsFactors=FALSE)
fit3$genes <- genesmd

volcanoplot(fit3, coef=2, highlight=7, names=fit3$genes$symbol, 
            main="Volcano plot: Colostrum and Transitional (known + unknown)", las=1)


## -----------------------------------------------------------------------------
top7_no3 <- order(fit3$lods[,2], decreasing=TRUE)
top7_no3 <- fit3$genes$symbol[top7_no3][1:7]

## -----------------------------------------------------------------------------
print(paste("The top 7 genes are", paste(top7_no3, collapse = ", "), sep = ": "))

## -----------------------------------------------------------------------------
mask <- DEgenesTab$EntrezID %in% DEgenes_no3
DEgenesTab <- DEgenesTab[mask,]

## ----CTtab, echo=FALSE, warning=FALSE-----------------------------------------
## generate full table in a CSV file and store it in the 'doc' directory
## twice, once in 'doc' to enable quickly look up during vignette editing
## and building with 'devtools::build_vignettes()' and a second time in
## 'inst/doc' to make these files available at install.
use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))
fnameCSV <- "DEgenes_no3.csv"
fpathCSV <- proj_path(file.path("doc", fnameCSV))
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)
fpathCSV <- proj_path(file.path("inst", "doc", fnameCSV))
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)

DEgenesEGs_no3 <- DEgenesEGs[mask]

## generate full table in HTML and store it into the 'doc' directory
## twice, just as we did with the CSV file. note that because the
## table caption is not translated from Markdown, but directly copied
## into HTML, we need to avoid using the '<' symbol, as in FDR < 10%,
## and put its HTML code instead (&lt;)
ktab <- kable(DEgenesTab, "html", escape=FALSE, row.names=TRUE,
              caption=sprintf("Differentially expressed genes. Differentially expressed genes between between Colostrum and Transitional lactation stages. FDR &lt; 10%% (CSV <a href=\"%s\" download>file</a>).",
                              fnameCSV))
ktab <- kable_styling(ktab,
                      bootstrap_options=c("stripped", "hover", "responsive"),
                      fixed_thead=TRUE)
fnameHTML <- "DEgenes_no3.html"
fpathHTML <- proj_path(file.path("doc", fnameHTML))
save_kable(ktab, file=fpathHTML, self_contained=TRUE)
fpathHTML <- proj_path(file.path("inst", "doc", fnameHTML))
save_kable(ktab, file=fpathHTML, self_contained=TRUE)


ktab <- kable(DEgenesTab[1:10, ], "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Differentially expressed genes. Top-10 differentially expressed genes with lowest p-value between Colostrum and Transitional lactation stages.",
                              fnameHTML, fnameCSV))
kable_styling(ktab, position="center")

## ----TMhist, fig.cap="P-value diagnostic plot for differental expression between Transitional and Mature."----
se.filt.all <- se.filt[,se.filt$lacStageFac!=1]
se.filt.all$stage <- droplevels(se.filt.all$lacStageFac)

mod <- model.matrix(~ se.filt.all$stage, colData(se.filt.all))
mod0 <- model.matrix(~ 1, colData(se.filt.all))


pv <- f.pvalue(assays(se.filt.all)$logCPM, mod, mod0)
hist(pv, main="", las=1)

## -----------------------------------------------------------------------------
mask <- p.adjust(pv, method="fdr") < 0.1
DEgenesEGs <- names(pv)[mask]
DEgenesSyms <- mcols(se.filt)[DEgenesEGs, "symbol"]
DEgenesPvalue <- pv[mask]
DEgenesDesc <- mcols(se.filt)[DEgenesEGs, "description"]
DEgenesDesc <- sub(" \\[.+\\]", "", DEgenesDesc)
DEgenesTab <- data.frame(EntrezID=DEgenesEGs,
                         Symbol=DEgenesSyms,
                         Description=DEgenesDesc,
                         "P value"=DEgenesPvalue,
                         stringsAsFactors=FALSE, check.names=FALSE)
DEgenesTab <- DEgenesTab[order(DEgenesTab[["P value"]]), ] 
rownames(DEgenesTab) <- 1:nrow(DEgenesTab)

## ----TMvolc, fig.cap="P-value diagnostic plot for differental expression between Transitional and Mature"----
sv <- sva(assays(se.filt.all)$logCPM, mod=mod, mod0=mod0)
mod <- cbind(mod, sv$sv)
colnames(mod) <- c(colnames(mod)[1:2], paste0("SV", 1:sv$n))
fit3 <- lmFit(assays(se.filt.all)$logCPM, mod)
fit3 <- eBayes(fit3)

tt3 <- topTable(fit3, coef=2, n=Inf)
DEgenes_no1 <- rownames(tt3)[tt3$adj.P.Val < 0.1]
sort(table(tt3[DEgenes_no1, "chr"]), decreasing=TRUE)
res3 <- decideTests(fit3, p.value=0.1)

genesmd <- data.frame(chr=as.character(seqnames(rowRanges(se.filt.all))), symbol=rowData(se.filt.all)[, 5], stringsAsFactors=FALSE)

fit3$genes <- genesmd

volcanoplot(fit3, coef=2, highlight=7, names=fit3$genes$symbol, 
            main="Volcano plot: Transitional and Mature (known + unknown)", las=1)

## -----------------------------------------------------------------------------
top7_no1 <- order(fit3$lods[,2], decreasing=TRUE)
top7_no1 <- fit3$genes$symbol[top7_no1][1:7]

## -----------------------------------------------------------------------------
print(paste("The top 7 genes are", paste(top7_no1, collapse = ", "), sep = ": "))

## -----------------------------------------------------------------------------
mask <- DEgenesTab$EntrezID %in% DEgenes_no3
DEgenesTab <- DEgenesTab[mask,]

## ----TMtab, echo=FALSE, warning=FALSE-----------------------------------------
## generate full table in a CSV file and store it in the 'doc' directory
## twice, once in 'doc' to enable quickly look up during vignette editing
## and building with 'devtools::build_vignettes()' and a second time in
## 'inst/doc' to make these files available at install.
use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))
fnameCSV <- "DEgenes_no1.csv"
fpathCSV <- proj_path(file.path("doc", fnameCSV))
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)
fpathCSV <- proj_path(file.path("inst", "doc", fnameCSV))
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)

DEgenesEGs_no1 <- DEgenesEGs[mask]

## generate full table in HTML and store it into the 'doc' directory
## twice, just as we did with the CSV file. note that because the
## table caption is not translated from Markdown, but directly copied
## into HTML, we need to avoid using the '<' symbol, as in FDR < 10%,
## and put its HTML code instead (&lt;)
ktab <- kable(DEgenesTab, "html", escape=FALSE, row.names=TRUE,
              caption=sprintf("Differentially expressed genes. Differentially expressed genes between between Transitional and Mature lactaction stages FDR &lt; 10%% (CSV <a href=\"%s\" download>file</a>).",
                              fnameCSV))
ktab <- kable_styling(ktab,
                      bootstrap_options=c("stripped", "hover", "responsive"),
                      fixed_thead=TRUE)
fnameHTML <- "DEgenes_no1.html"
fpathHTML <- proj_path(file.path("doc", fnameHTML))
save_kable(ktab, file=fpathHTML, self_contained=TRUE)
fpathHTML <- proj_path(file.path("inst", "doc", fnameHTML))
save_kable(ktab, file=fpathHTML, self_contained=TRUE)


ktab <- kable(DEgenesTab[1:10, ], "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Differentially expressed genes. Top-10 differentially expressed genes with lowest p-value between Transitional and Mature lactation stages",
                              fnameHTML, fnameCSV))
kable_styling(ktab, position="center")

## ----CMhist, fig.cap="P-value diagnostic plot for differental expression between Colostrum and Mature."----
se.filt.all <- se.filt[,se.filt$lacStageFac!=2]
se.filt.all$stage <- droplevels(se.filt.all$lacStageFac)

mod <- model.matrix(~ se.filt.all$stage,
                    colData(se.filt.all))
mod0 <- model.matrix(~ 1, colData(se.filt.all))


pv <- f.pvalue(assays(se.filt.all)$logCPM, mod, mod0)
#sum(p.adjust(pv, method="fdr") < 0.05)
#sum(p.adjust(pv, method="fdr") < 0.1)
hist(pv, main="", las=1)

## ----warning=FALSE------------------------------------------------------------
mask <- p.adjust(pv, method="fdr") < 0.1
DEgenesEGs <- names(pv)[mask]
DEgenesSyms <- mcols(se.filt)[DEgenesEGs, "symbol"]
DEgenesPvalue <- pv[mask]
DEgenesDesc <- mcols(se.filt)[DEgenesEGs, "description"]
DEgenesDesc <- sub(" \\[.+\\]", "", DEgenesDesc)
DEgenesTab <- data.frame(EntrezID=DEgenesEGs,
                         Symbol=DEgenesSyms,
                         Description=DEgenesDesc,
                         "P value"=DEgenesPvalue,
                         stringsAsFactors=FALSE, check.names=FALSE)
DEgenesTab <- DEgenesTab[order(DEgenesTab[["P value"]]), ] ## order by p-value
rownames(DEgenesTab) <- 1:nrow(DEgenesTab)

## ----CMvolc, fig.cap="P-value diagnostic plot for differental expression between Colostrum and Mature"----
sv <- sva(assays(se.filt.all)$logCPM, mod=mod, mod0=mod0)
mod <- cbind(mod, sv$sv)
colnames(mod) <- c(colnames(mod)[1:2], paste0("SV", 1:sv$n))
fit3 <- lmFit(assays(se.filt.all)$logCPM, mod)
fit3 <- eBayes(fit3)
tt3 <- topTable(fit3, coef=2, n=Inf)

DEgenes_no2 <- rownames(tt3)[tt3$adj.P.Val < 0.1]
sort(table(tt3[DEgenes_no2, "chr"]), decreasing=TRUE)
res3 <- decideTests(fit3, p.value=0.1)

genesmd <- data.frame(chr=as.character(seqnames(rowRanges(se.filt.all))), symbol=rowData(se.filt.all)[, 5], stringsAsFactors=FALSE)
fit3$genes <- genesmd

volcanoplot(fit3, coef=2, highlight=7, names=fit3$genes$symbol, main="Volcano plot: Colostrum and Mature (known + unknown)", las=1)

## -----------------------------------------------------------------------------
top7_no2 <- order(fit3$lods[,2], decreasing=TRUE)
top7_no2 <- fit3$genes$symbol[top7_no2][1:7]

## -----------------------------------------------------------------------------
print(paste("The top 7 genes are", paste(top7_no2, collapse = ", "), sep = ": "))

## -----------------------------------------------------------------------------
mask <- DEgenesTab$EntrezID %in% DEgenes_no2
DEgenesTab <- DEgenesTab[mask,]

## ----CMtab, echo=FALSE, warning=FALSE-----------------------------------------
## generate full table in a CSV file and store it in the 'doc' directory
## twice, once in 'doc' to enable quickly look up during vignette editing
## and building with 'devtools::build_vignettes()' and a second time in
## 'inst/doc' to make these files available at install.
use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))
fnameCSV <- "DEgenes_no2.csv"
fpathCSV <- proj_path(file.path("doc", fnameCSV))
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)
fpathCSV <- proj_path(file.path("inst", "doc", fnameCSV))
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)

DEgenesEGs_no2 <- DEgenesEGs[mask]

## generate full table in HTML and store it into the 'doc' directory
## twice, just as we did with the CSV file. note that because the
## table caption is not translated from Markdown, but directly copied
## into HTML, we need to avoid using the '<' symbol, as in FDR < 10%,
## and put its HTML code instead (&lt;)
ktab <- kable(DEgenesTab, "html", escape=FALSE, row.names=TRUE,
              caption=sprintf("Differentially expressed genes. Differentially expressed genes between between Colostrum and Mature lactaction stages FDR &lt; 10%% (CSV <a href=\"%s\" download>file</a>).",
                              fnameCSV))
ktab <- kable_styling(ktab,
                      bootstrap_options=c("stripped", "hover", "responsive"),
                      fixed_thead=TRUE)
fnameHTML <- "DEgenes_no2.html"
fpathHTML <- proj_path(file.path("doc", fnameHTML))
save_kable(ktab, file=fpathHTML, self_contained=TRUE)
fpathHTML <- proj_path(file.path("inst", "doc", fnameHTML))
save_kable(ktab, file=fpathHTML, self_contained=TRUE)


ktab <- kable(DEgenesTab[1:10, ], "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Differentially expressed genes. Top-10 differentially expressed genes with lowest p-value between Colostrum and Mature lactation stages",
                              fnameHTML, fnameCSV))
kable_styling(ktab, position="center")

## ----echo=FALSE---------------------------------------------------------------
x <- c(length(DEgenesEGs_no3), length(DEgenesEGs_no1), length(DEgenesEGs_no2))
y <- c(8015, 1948, 8527)
correlation_value_with_paper <- cor(x, y)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(org.Hs.eg.db)
library(GOstats)
library(KEGGREST)

geneUniverse <- rownames(se)

## -----------------------------------------------------------------------------
params <- new("GOHyperGParams", geneIds=DEgenesEGs_no3, 
              universeGeneIds=geneUniverse, annotation="org.Hs.eg.db", 
              ontology="BP", pvalueCutoff=0.05, testDirection="over")

conditional(params) <- TRUE
hgOverCond <- hyperGTest(params)
goresults <- summary(hgOverCond)

## -----------------------------------------------------------------------------
mask <- goresults$OddsRatio != "Inf"
goresults <- goresults[mask, ]
goresults <- goresults[order(goresults$OddsRatio, decreasing=TRUE), ]
goresults <- goresults[goresults$Size >= 3 & goresults$Size <= 100 & goresults$Count >=25, ]

## ----message=FALSE, warning=FALSE---------------------------------------------
geneIDs <- geneIdsByCategory(hgOverCond)[goresults$GOBPID]
geneSYMs <- sapply(geneIDs, function(id) select(org.Hs.eg.db, columns="SYMBOL", key=id, keytype="ENTREZID")$SYMBOL)
geneSYMs <- sapply(geneSYMs, paste, collapse=", ")
goresults <- cbind(goresults, Genes=geneSYMs)
rownames(goresults) <- 1:nrow(goresults)

## ----CTGOtab, echo=FALSE------------------------------------------------------
ktab <- kable(goresults, "html", caption="GO results between Colostrum and Transitional")
ktab <- kable_styling(ktab, bootstrap_options=c("stripped", "hover", "responsive"), fixed_thead=TRUE)
save_kable(ktab, file="../doc/goresults_no3.html", self_contained=TRUE)

ktab <- kable(goresults[1:10, 2:7], "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("GO Enrichment Between Colostrum and Transitional"))
kable_styling(ktab, position="center")

use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))
fnameCSV <- "goresults_no3.csv"
fpathCSV <- proj_path(file.path("doc", fnameCSV))
write.csv(goresults, fpathCSV, row.names=FALSE)
fpathCSV <- proj_path(file.path("inst", "doc", fnameCSV))
write.csv(goresults, fpathCSV, row.names=FALSE)

## -----------------------------------------------------------------------------
KEGGparams <- new("KEGGHyperGParams", geneIds=DEgenesEGs_no3, universeGeneIds=geneUniverse, annotation="org.Hs.eg.db", pvalueCutoff=0.05, testDirection="over")
KEGGhgOver <- hyperGTest(KEGGparams)
KEGGresults <- summary(KEGGhgOver)

## -----------------------------------------------------------------------------
mask <- KEGGresults$OddsRatio != "Inf"
KEGGresults <- KEGGresults[mask,]
KEGGresults <- KEGGresults[order(KEGGresults$OddsRatio, decreasing = TRUE), ]

## -----------------------------------------------------------------------------
source("../R/custom_functions.R")

names_names <- sapply(KEGGresults$KEGGID, getKEGGName)
KEGGresults$Term <- as.data.frame(names_names)$names_names

class_names <- sapply(KEGGresults$KEGGID, getKEGGClass)
KEGGresults$KEGGClassNames <- as.data.frame(class_names)$class_names

## -----------------------------------------------------------------------------
KEGGresults <- KEGGresults[KEGGresults$Size < 150, ]

## ----message=FALSE------------------------------------------------------------
KEGGgeneIDs <- geneIdsByCategory(KEGGhgOver)[KEGGresults$KEGGID]
KEGGgeneSYMs <- sapply(KEGGgeneIDs, function(id) select(org.Hs.eg.db, columns="SYMBOL", key=id, keytype="ENTREZID")$SYMBOL)
KEGGgeneSYMs <- sapply(KEGGgeneSYMs, paste, collapse=", ")
KEGGresults_with_genes <- cbind(KEGGresults, Genes=KEGGgeneSYMs)


## ----CTKEGGtab, echo=FALSE----------------------------------------------------
ktab <- kable(KEGGresults_with_genes, "html", caption="KEGG results between Colostrum and Transitional")
ktab <- kable_styling(ktab, bootstrap_options=c("stripped", "hover", "responsive"), fixed_thead=TRUE)
save_kable(ktab, file="../doc/KEGGresults_no3.html", self_contained=TRUE)

ktab <- kable(KEGGresults[1:10, 2:8], "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("KEGG Enrichment Between Colostrum and Transitional"))
kable_styling(ktab, position="center")

use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))
fnameCSV <- "KEGGresults_with_genes_no3.csv"
fpathCSV <- proj_path(file.path("doc", fnameCSV))
write.csv(KEGGresults_with_genes, fpathCSV, row.names=FALSE)
fpathCSV <- proj_path(file.path("inst", "doc", fnameCSV))
write.csv(KEGGresults_with_genes, fpathCSV, row.names=FALSE)


## -----------------------------------------------------------------------------
new_table <- data.frame(Gene = character(), Term = character(), stringsAsFactors = FALSE)

for (i in 1:nrow(goresults)) {
  genes <- strsplit(goresults$Genes[i], ", ")[[1]]
  term <- rep(goresults$Term[i], length(genes))
  new_rows <- data.frame(Gene = genes, Term = term, stringsAsFactors = FALSE)
  new_table <- rbind(new_table, new_rows)
}

mask <-new_table$Gene %in% top7_no3 
goGenesInTop7_no3 <- new_table[mask,]

ktab <- kable(goGenesInTop7_no3, "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Top 7 fold-change genes in GO enriched pathways between Colostrum and Transitional."))
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
KEGG_new_table <- data.frame(Gene = character(), Term = character(), stringsAsFactors = FALSE)

for (i in 1:nrow(KEGGresults_with_genes)) {
  genes <- strsplit(KEGGresults_with_genes$Genes[i], ", ")[[1]]
  term <- rep(KEGGresults_with_genes$Term[i], length(genes))
  KEGG_new_rows <- data.frame(Gene = genes, Term = term, stringsAsFactors = FALSE)
  KEGG_new_table <- rbind(KEGG_new_table, KEGG_new_rows)
}

mask <-KEGG_new_table$Gene %in% top7_no2 
KEGGGenesInTop7_no3 <- KEGG_new_table[mask,]

ktab <- kable(KEGGGenesInTop7_no3, "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Top 7 fold-change genes in GO enriched pathways between Colostrum and Transitional."))
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
params <- new("GOHyperGParams", geneIds=DEgenesEGs_no1, universeGeneIds=geneUniverse, annotation="org.Hs.eg.db", ontology="BP", pvalueCutoff=0.05, testDirection="over")

conditional(params) <- TRUE
hgOverCond <- hyperGTest(params)
goresults <- summary(hgOverCond)

## -----------------------------------------------------------------------------
mask <- goresults$OddsRatio != "Inf"
goresults <- goresults[mask, ]
goresults <- goresults[order(goresults$OddsRatio, decreasing=TRUE), ]

goresults <- goresults[goresults$Size >= 3 & goresults$Size <= 100 & goresults$Count >=15, ]

## ----message=FALSE, warning=FALSE---------------------------------------------
geneIDs <- geneIdsByCategory(hgOverCond)[goresults$GOBPID]
geneSYMs <- sapply(geneIDs, function(id) select(org.Hs.eg.db, columns="SYMBOL", key=id, keytype="ENTREZID")$SYMBOL)
geneSYMs <- sapply(geneSYMs, paste, collapse=", ")
goresults <- cbind(goresults, Genes=geneSYMs)
rownames(goresults) <- 1:nrow(goresults)

## ----TMGOtab, echo=FALSE------------------------------------------------------
ktab <- kable(goresults, "html", caption="GO results between Transitional and Mature")
ktab <- kable_styling(ktab, bootstrap_options=c("stripped", "hover", "responsive"), fixed_thead=TRUE)
save_kable(ktab, file="../doc/goresults_no1.html", self_contained=TRUE)

ktab <- kable(goresults[1:10, 2:7], "html", escape=FALSE, row.names=TRUE, caption=sprintf("GO Enrichment Between Transitional and Mature"))
kable_styling(ktab, position="center")

use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))
fnameCSV <- "goresults_no1.csv"
fpathCSV <- proj_path(file.path("doc", fnameCSV))
write.csv(goresults, fpathCSV, row.names=FALSE)
fpathCSV <- proj_path(file.path("inst", "doc", fnameCSV))
write.csv(goresults, fpathCSV, row.names=FALSE)

## -----------------------------------------------------------------------------
KEGGparams <- new("KEGGHyperGParams", geneIds=DEgenesEGs_no1, universeGeneIds=geneUniverse, annotation="org.Hs.eg.db", pvalueCutoff=0.05, testDirection="over")
KEGGhgOver <- hyperGTest(KEGGparams)
KEGGresults <- summary(KEGGhgOver)

## -----------------------------------------------------------------------------
mask <- KEGGresults$OddsRatio != "Inf"
KEGGresults <- KEGGresults[mask,]
KEGGresults <- KEGGresults[order(KEGGresults$OddsRatio, decreasing = TRUE), ]

## -----------------------------------------------------------------------------
source("../R/custom_functions.R")

names_names <- sapply(KEGGresults$KEGGID, getKEGGName)
KEGGresults$Term <- as.data.frame(names_names)$names_names

class_names <- sapply(KEGGresults$KEGGID, getKEGGClass)
KEGGresults$KEGGClassNames <- as.data.frame(class_names)$class_names

## -----------------------------------------------------------------------------
KEGGresults <- KEGGresults[KEGGresults$Size < 130, ]

## ----message=FALSE------------------------------------------------------------
KEGGgeneIDs <- geneIdsByCategory(KEGGhgOver)[KEGGresults$KEGGID]
KEGGgeneSYMs <- sapply(KEGGgeneIDs, function(id) select(org.Hs.eg.db, columns="SYMBOL", key=id, keytype="ENTREZID")$SYMBOL)
KEGGgeneSYMs <- sapply(KEGGgeneSYMs, paste, collapse=", ")
KEGGresults_with_genes <- cbind(KEGGresults, Genes=KEGGgeneSYMs)

## ----TMKEGGtab, echo=FALSE----------------------------------------------------
ktab <- kable(KEGGresults_with_genes, "html", caption="KEGG results between Transitional and Mature")
ktab <- kable_styling(ktab, bootstrap_options=c("stripped", "hover", "responsive"), fixed_thead=TRUE)
save_kable(ktab, file="../doc/KEGGresults_no1.html", self_contained=TRUE)

ktab <- kable(KEGGresults[1:10, 2:8], "html", escape=FALSE, row.names=TRUE, caption=sprintf("KEGG Enrichment Between Transitional and Mature"))
kable_styling(ktab, position="center")

use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))
fnameCSV <- "KEGGresults_with_genes_no1.csv"
fpathCSV <- proj_path(file.path("doc", fnameCSV))
write.csv(KEGGresults_with_genes, fpathCSV, row.names=FALSE)
fpathCSV <- proj_path(file.path("inst", "doc", fnameCSV))
write.csv(KEGGresults_with_genes, fpathCSV, row.names=FALSE)

## -----------------------------------------------------------------------------
new_table <- data.frame(Gene = character(), Term = character(), stringsAsFactors = FALSE)

for (i in 1:nrow(goresults)) {
  genes <- strsplit(goresults$Genes[i], ", ")[[1]]
  term <- rep(goresults$Term[i], length(genes))
  new_rows <- data.frame(Gene = genes, Term = term, stringsAsFactors = FALSE)
  new_table <- rbind(new_table, new_rows)
}

mask <-new_table$Gene %in% top7_no1 
goGenesInTop7_no1 <- new_table[mask,]

ktab <- kable(goGenesInTop7_no1, "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Top 7 fold-change genes in GO enriched pathways between Transitional and Mature"))
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
KEGG_new_table <- data.frame(Gene = character(), Term = character(), stringsAsFactors = FALSE)

for (i in 1:nrow(KEGGresults_with_genes)) {
  genes <- strsplit(KEGGresults_with_genes$Genes[i], ", ")[[1]]
  term <- rep(KEGGresults_with_genes$Term[i], length(genes))
  KEGG_new_rows <- data.frame(Gene = genes, Term = term, stringsAsFactors = FALSE)
  KEGG_new_table <- rbind(KEGG_new_table, KEGG_new_rows)
}

mask <-KEGG_new_table$Gene %in% top7_no1 
KEGGGenesInTop7_no1 <- KEGG_new_table[mask,]

ktab <- kable(KEGGGenesInTop7_no1, "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Top 7 fold-change genes in GO enriched pathways between Transitional and Mature"))
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
params <- new("GOHyperGParams", geneIds=DEgenesEGs_no2,universeGeneIds=geneUniverse, annotation="org.Hs.eg.db", ontology="BP", pvalueCutoff=0.05, testDirection="over")

conditional(params) <- TRUE
hgOverCond <- hyperGTest(params)
goresults <- summary(hgOverCond)

## -----------------------------------------------------------------------------
mask <- goresults$OddsRatio != "Inf"
goresults <- goresults[mask, ]
goresults <- goresults[order(goresults$OddsRatio, decreasing=TRUE), ]

goresults <- goresults[goresults$Size >= 3 & goresults$Size <= 100 & goresults$Count >=15, ]

## ----message=FALSE, warning=FALSE---------------------------------------------
geneIDs <- geneIdsByCategory(hgOverCond)[goresults$GOBPID]
geneSYMs <- sapply(geneIDs, function(id) select(org.Hs.eg.db, columns="SYMBOL", key=id, keytype="ENTREZID")$SYMBOL)
geneSYMs <- sapply(geneSYMs, paste, collapse=", ")
goresults <- cbind(goresults, Genes=geneSYMs)
rownames(goresults) <- 1:nrow(goresults)

## ----CMGOtab, echo=FALSE------------------------------------------------------
ktab <- kable(goresults, "html", caption="GO results between Colostrum and Mature")
ktab <- kable_styling(ktab, bootstrap_options=c("stripped", "hover", "responsive"), fixed_thead=TRUE)
save_kable(ktab, file="../doc/goresults_no2.html", self_contained=TRUE)

ktab <- kable(goresults[1:10, 2:7], "html", escape=FALSE, row.names=TRUE, caption=sprintf("GO Enrichment Between Colostrum and Mature"))
kable_styling(ktab, position="center")

use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))
fnameCSV <- "goresults_no2.csv"
fpathCSV <- proj_path(file.path("doc", fnameCSV))
write.csv(goresults, fpathCSV, row.names=FALSE)
fpathCSV <- proj_path(file.path("inst", "doc", fnameCSV))
write.csv(goresults, fpathCSV, row.names=FALSE)

## -----------------------------------------------------------------------------
KEGGparams <- new("KEGGHyperGParams", geneIds=DEgenesEGs_no2, universeGeneIds=geneUniverse, annotation="org.Hs.eg.db", pvalueCutoff=0.05, testDirection="over")
KEGGhgOver <- hyperGTest(KEGGparams)
KEGGresults <- summary(KEGGhgOver)

## -----------------------------------------------------------------------------
mask <- KEGGresults$OddsRatio != "Inf"
KEGGresults <- KEGGresults[mask,]
KEGGresults <- KEGGresults[order(KEGGresults$OddsRatio, decreasing = TRUE), ]

## -----------------------------------------------------------------------------
source("../R/custom_functions.R")

names_names <- sapply(KEGGresults$KEGGID, getKEGGName)
KEGGresults$Term <- as.data.frame(names_names)$names_names

class_names <- sapply(KEGGresults$KEGGID, getKEGGClass)
KEGGresults$KEGGClassNames <- as.data.frame(class_names)$class_names

## -----------------------------------------------------------------------------
KEGGresults <- KEGGresults[KEGGresults$Size < 100, ]

## ----message=FALSE------------------------------------------------------------
KEGGgeneIDs <- geneIdsByCategory(KEGGhgOver)[KEGGresults$KEGGID]
KEGGgeneSYMs <- sapply(KEGGgeneIDs, function(id) select(org.Hs.eg.db, columns="SYMBOL", key=id, keytype="ENTREZID")$SYMBOL)
KEGGgeneSYMs <- sapply(KEGGgeneSYMs, paste, collapse=", ")
KEGGresults_with_genes <- cbind(KEGGresults, Genes=KEGGgeneSYMs)

## ----CMGKEGGtab, echo=FALSE---------------------------------------------------
ktab <- kable(KEGGresults_with_genes, "html", caption="KEGG results between Colostrum and Mature")
ktab <- kable_styling(ktab, bootstrap_options=c("stripped", "hover", "responsive"), fixed_thead=TRUE)
save_kable(ktab, file="../doc/KEGGresults_no2.html", self_contained=TRUE)

ktab <- kable(KEGGresults[1:10, 2:8], "html", escape=FALSE, row.names=TRUE, caption=sprintf("KEGG Enrichment Between Colostrum and Mature"))
kable_styling(ktab, position="center")

use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))
fnameCSV <- "KEGGresults_with_genes_no2.csv"
fpathCSV <- proj_path(file.path("doc", fnameCSV))
write.csv(KEGGresults_with_genes, fpathCSV, row.names=FALSE)
fpathCSV <- proj_path(file.path("inst", "doc", fnameCSV))
write.csv(KEGGresults_with_genes, fpathCSV, row.names=FALSE)

## -----------------------------------------------------------------------------
new_table <- data.frame(Gene = character(), Term = character(), stringsAsFactors = FALSE)

for (i in 1:nrow(goresults)) {
  genes <- strsplit(goresults$Genes[i], ", ")[[1]]
  term <- rep(goresults$Term[i], length(genes))
  new_rows <- data.frame(Gene = genes, Term = term, stringsAsFactors = FALSE)
  new_table <- rbind(new_table, new_rows)
}

mask <-new_table$Gene %in% top7_no2 
goGenesInTop7_no2 <- new_table[mask,]

ktab <- kable(goGenesInTop7_no2, "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Top 7 fold-change genes in GO enriched pathways between Colostrum and Mature"))
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
KEGG_new_table <- data.frame(Gene = character(), Term = character(), stringsAsFactors = FALSE)

for (i in 1:nrow(KEGGresults_with_genes)) {
  genes <- strsplit(KEGGresults_with_genes$Genes[i], ", ")[[1]]
  term <- rep(KEGGresults_with_genes$Term[i], length(genes))
  KEGG_new_rows <- data.frame(Gene = genes, Term = term, stringsAsFactors = FALSE)
  KEGG_new_table <- rbind(KEGG_new_table, KEGG_new_rows)
}

mask <-KEGG_new_table$Gene %in% top7_no2 
KEGGGenesInTop7_no2 <- KEGG_new_table[mask,]

ktab <- kable(KEGGGenesInTop7_no2, "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Top 7 fold-change genes in GO enriched pathways between Colostrum and Mature"))
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
sessionInfo()

