library(KEGGREST)

getKEGGName <- function(input) {
  # Construct the KEGG ID
  kegg_id <- paste0("hsa", input)

  tryCatch(
    {
      # Retrieve KEGG information
      rezKEGG <- keggGet(kegg_id)
      name <- rezKEGG[[1]]$NAME
      name <- gsub(" - Homo sapiens \\(human\\)", "", name)
      return(name)

    },
    error = function(err) {
      # Code to execute if an error occurs
      name <- "no info"
      return(name)
    }
  )
}

getKEGGClass <- function(input) {

  # Construct the KEGG ID
  kegg_id <- paste0("hsa", input)

  tryCatch(
    {
      # Retrieve KEGG information
      rezKEGG <- keggGet(kegg_id)
      class_name <- rezKEGG[[1]]$CLASS
      if (is.null(class_name)) {
        class_name <- "no info"
      }
      return(class_name)

    },
    error = function(err) {
      # Code to execute if an error occurs
      class_name <- "no info"
      return(class_name)
    }
  )
}
