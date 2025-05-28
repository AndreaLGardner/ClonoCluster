## ---- build_barcode_matrix_fast
#' Build a sparse matrix representing the barcode graph.
#'
#' @param bt Data.table. Columns are "rn" (cell ids) and "Barcode".
#' @param value Numeric. A value to be used for edge weights within barcodes. Default NULL will use the reciprocal of barcode size.
#'
#' @return A sparse matrix of class `dgCMatrix` suitable for `clonocluster_model` or `Seurat::FindClusters`.
#'
#' @export build_barcode_matrix_fast
#' @md
build_barcode_matrix_fast <- function(bt, value = NULL){

  #bc <- bt[, .SD %>% unique, .SDcols = c("rn", "Barcode")] %>%
  #dcast(rn ~ Barcode,
  # fun.aggregate = function(x) ifelse(length(x) > 0, 1, 0))

  ## Want to be able to include non-barcoded cells in transcritpomic clustering, but not have them all be assigned the same label for lineage clustering
  ## Also want to be able to include information saying that fusion cells contain more than 1 barcode
  ## To do both of these tasks, I will first assign all non-barcoded cells a dummy name 'na1, na2, na3, ..., naN'
  ## Then create the barcode matrix by first spiltting multi-barcoded cells into multiple rows where each row represents a unique barcode in a cell
  ## Then creating the one-hot encoded table by making wider.

  #### set aside barcoded cells
  bt1 <- bt[is.na(bt$Barcode)==FALSE,]
  
  #### find non-barcoded (Barcode == NA) cells and assign a unique dummy name to each so that they don't cluster together
  bt2 <- bt[is.na(bt$Barcode)==TRUE,]
  bt2 <- bt2 %>% mutate(Barcode = paste0('na',seq(1,nrow(bt2))))
  
  #### combine back together with new dummy barcode names for NA cells
  bt <- rbind(bt1, bt2)
  
  #### split multi-barcoded cells on underscore such that each row is a cell and a unique barcode (there will be duplicated rows for multi-barcoded cells)
  bc <- df %>%
    mutate(Barcode = str_split(Barcode, "_")) %>%
    unnest(Barcode)

  rn <- bc[, rn]

  bc <- bc[, .SD, .SDcols = 2:ncol(bc)] %>% as.matrix

  rownames(bc) <- rn

  # matrix to normalize size
  nm <- list()

  l <- apply(bc, 2, function(x) which(x == 1))

  message("Making lineage matrix")

  pb <- utils::txtProgressBar(min = 0, max = length(l), style = 3)

  for (i in 1:length(l)){

    utils::setTxtProgressBar(pb, i)

    ll <- matrix(1/length(l[[i]]), nrow = length(l[[i]]), ncol = length(l[[i]]),
                dimnames = list(names(l[[i]]), names(l[[i]])))

    r <- matrix(0, nrow = nrow(ll), ncol = nrow(bc) - ncol(ll),
                dimnames = list(rownames(ll),
                rownames(bc)[!rownames(bc) %chin% colnames(ll)])
                )

    ll <- cbind(ll, r)

    d <- matrix(0, nrow = nrow(bc) - nrow(ll), ncol = ncol(ll),
                dimnames = list(rownames(bc)[!rownames(bc) %chin% rownames(ll)],
                colnames(ll))
              )

    ll <- rbind(ll, d)

    ll <- ll[rownames(bc), rownames(bc)]

    ll %<>% as("dgCMatrix")

    nm[[length(nm) + 1]] <- ll

  }

  close(pb)

  nm %<>% purrr::reduce(`+`)

  if (!is.null(value)){

    nm@x[nm@x > 0] <- value

    nm %<>% as("dgCMatrix")

  }

  return(nm)

}
