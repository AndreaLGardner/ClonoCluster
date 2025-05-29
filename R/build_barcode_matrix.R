## ---- build_barcode_matrix
#' Build a sparse matrix representing the barcode graph.
#'
#' @param bt Data.table. Columns are "rn" (cell ids) and "Barcode".
#' @param value Numeric. A value to be used for edge weights within barcodes. Default NULL will use the reciprocal of barcode size.
#'
#' @return A sparse matrix of class `dgCMatrix` suitable for `clonocluster_model` or `Seurat::FindClusters`.
#'
#' @export build_barcode_matrix
#' @md
#' @importFrom dplyr filter mutate ungroup n 
#' @importFrom tidyr pivot_wider unnest
build_barcode_matrix <- function(bt, value = NULL){

  #bc <- bt[, .SD %>% unique, .SDcols = c("rn", "Barcode")] %>%
  #dcast(rn ~ Barcode,
  #  fun.aggregate = function(x) ifelse(length(x) > 0, 1, 0))
  ## Want to be able to include non-barcoded cells in transcritpomic clustering, but not have them all be assigned the same label for lineage clustering
  ## Also want to be able to include information saying that fusion cells contain more than 1 barcode
  ## To do both of these tasks, I will first assign all non-barcoded cells a dummy name 'na1, na2, na3, ..., naN'
  ## Then create the barcode matrix by first spiltting multi-barcoded cells into multiple rows where each row represents a unique barcode in a cell
  ## Then creating the one-hot encoded table by making wider.
  
  #### if there are non-bcarded (Barcode == NA_ cells) then remap NA cells to dummy barcodes
  
  if (any(is.na(bt$Barcode))) {
    #### set aside barcoded cells    
    bt1 <- bt[!is.na(bt$Barcode), ]
    #### find non-barcoded (Barcode == NA) cells and assign a unique dummy name
    bt2 <- bt[is.na(bt$Barcode), ]
    bt2$Barcode <- paste0("na", seq_len(nrow(bt2)))
    #### combine back together with new dummy barcode names for NA cells
    bt <- rbind(bt1, bt2)
  }

 #### split multi-barcoded cells on underscore such that each row is a cell and a unique barcode (there will be duplicated rows for multi-barcoded cells)

  split_barcodes <- strsplit(bt$Barcode, "_")

  bt_expanded <- bt[rep(seq_len(nrow(bt)), lengths(split_barcodes)), ]

  bt_expanded$Barcode <- unlist(split_barcodes)

  #### one-hot table for barcode identity
  bc <- bt_expanded %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = Barcode, values_from = value, values_fill = list(value = 0)) %>%
    data.table() 

  rn <- bc[, rn]

  bc <- bc[, .SD, .SDcols = 2:ncol(bc)] %>% as.matrix

  rownames(bc) <- rn
  # matrix to normalize size
  nm <- matrix(0, nrow = nrow(bc), ncol = nrow(bc),
    dimnames = list(rownames(bc), rownames(bc))) %>% as("dgCMatrix")

  l <- apply(bc, 2, function(x) which(x == 1))

  message("Making lineage matrix")

  pb <- utils::txtProgressBar(min = 0, max = length(l), style = 3)

  for (i in 1:length(l)){

    utils::setTxtProgressBar(pb, i)

    ll <- data.table::CJ(names(l[[i]]), names(l[[i]]))

    nm[ll[, V1], ll[, V2]] <- 1/length(l[[i]])

  }

  close(pb)

  if (!is.null(value)){

    nm@x[nm@x > 0] <- value

    nm %<>% as("dgCMatrix")

  }

  return(nm)

}
