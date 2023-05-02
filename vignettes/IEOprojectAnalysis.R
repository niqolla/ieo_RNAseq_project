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

## -----------------------------------------------------------------------------
# extract the colData from the SingleCellExperiment object
coldata <- colData(se)

# loop over the columns and count unique values
for (col in names(coldata)) {
  unique_vals <- unique(coldata[[col]])
  if (any(length(unique_vals) %in% c(1, 16))) {
    next # skip to the next column
  } else {
    counts <- table(coldata[[col]])
    cat(paste0("Column '", col, "' has ", length(unique_vals), " unique values:\n"))
    print(counts)
  }
}


## ---- message=FALSE-----------------------------------------------------------
library(IEOproject)
library(edgeR)

dge <- DGEList(counts=assays(se)$counts, genes=rowData(se))
dim(dge)

## -----------------------------------------------------------------------------
assays(se)$logCPM <- cpm(dge, log=TRUE)
assays(se)$logCPM[1:5, 1:5]

## -----------------------------------------------------------------------------
table(se$characteristics_ch1)
table(se$characteristics_ch1.2)

## -----------------------------------------------------------------------------
se$cell_line <- se$characteristics_ch1
levels(se$cell_line) <- c("A549", "NHBE")
se$treatment <- se$characteristics_ch1.2
tmplevels <- gsub(" treatment", "", gsub("treatment: ", "", levels(se$treatment)))
tmplevels <- gsub(" infected ", "", tmplevels)
tmplevels <- gsub("\\(MOI ", "MOI", gsub(")", "", gsub("-", "", tmplevels)))
levels(se$treatment) <- tmplevels

## -----------------------------------------------------------------------------
table(se$extract_protocol_ch1)

## -----------------------------------------------------------------------------
se$description

