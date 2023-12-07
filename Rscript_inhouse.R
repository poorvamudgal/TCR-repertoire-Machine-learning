gene_usage_internal <- function (.data, .genes = HUMAN_TRBV_MITCR, .quant = c(NA, "read.count", 
                                                       "umi.count", "read.prop", "umi.prop"), .norm = F, .ambig = F) 
{
  message("\n========================================\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!! The tcR package WILL SOON BE ORPHANED \n!! AND REMOVED FROM CRAN.\n!!\n!! A new package is available that is \n!! designed to replace tcR: \n!! immunarch  --  https://immunarch.com/\n!!\n!! We will be happy to help you to move\n!! to the new package. Feel free to contact us:\n!! http://github.com/immunomind/immunarch\n!!\n!! Sincerely, \n!!  immunarch dev team and \n!!  Vadim I. Nazarov, lead developer of tcR\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n=======================================")
  .process.df <- function(.df, .quant, .cols) {
    cast.fun <- dcast
    if (length(.cols) == 2) {
      cast.fun <- acast
      len <- 2
    }
    for (i in 1:length(.cols)) {
      .df[which(!(.df[[.cols[i]]] %in% .genes[[i]])), .cols[i]] <- "Ambiguous"
    }
    if (length(.cols) == 1) {
      .cols <- c(.cols, ".")
      len <- 1
    }
    .df = .df %>% select(na.exclude(c(.quant, .cols[1:len])))
    .df = .df %>% grouped_df(.cols[1:len])
    if (!is.na(.quant)) {
      .df = .df %>% summarise(Freq = sum(Read.count))
    }
    else {
      .df = .df %>% summarise(Freq = n())
    }
    .df
  }
  .fix.ambig <- function(.res, .ambig) {
    if (length(.genes) == 2) {
      .res <- lapply(.res, function(x) {
        x[row.names(x) != "Ambiguous", ][, colnames(x) != 
                                           "Ambiguous"]
      })
      if (length(.data) == 1) {
        .res <- .res[[1]]
      }
      .res
    }
    else {
      .res[.res[, 1] != "Ambiguous", ]
    }
  }
  quant <- NA
  if (!is.na(.quant[1])) {
    quant <- "Read.count"
  }
  if (has.class(.data, "data.frame")) {
    .data <- list(Sample = .data)
  }
  if (has.class(.genes, "list")) {
    genecols <- c(paste0(substr(.genes[[1]][1], 4, 4), ".gene"), 
                  paste0(substr(.genes[[2]][1], 4, 4), ".gene"))
  }
  else {
    genecols <- paste0(substr(.genes[1], 4, 4), ".gene")
    .genes <- list(.genes)
  }
  tbls <- lapply(.data, .process.df, .quant = quant, .cols = genecols)
  if (length(.genes) == 2) {
    tbls <- lapply(tbls, function(x) {
      genrows <- .genes[[1]][is.na(match(.genes[[1]], row.names(x)))]
      gencols <- .genes[[2]][is.na(match(.genes[[2]], colnames(x)))]
      if (length(genrows) > 0) {
        x = rbind(x, matrix(0, ncol = ncol(x), nrow = length(genrows)))
        row.names(x)[(nrow(x) - length(genrows) + 1):nrow(x)] <- genrows
      }
      if (length(gencols) > 0) {
        x = cbind(x, matrix(0, nrow = nrow(x), ncol = length(gencols)))
        colnames(x)[(ncol(x) - length(gencols) + 1):ncol(x)] <- gencols
      }
      x[is.na(x)] <- 0
      x
    })
    if (.norm) {
      tbls <- lapply(tbls, function(x) x/sum(x))
    }
    return(.fix.ambig(tbls, .ambig))
  }
  res <- tbls[[1]]
  colnames(res) <- c("Gene", names(.data)[1])
  if (length(.data) > 1) {
    for (i in 2:length(.data)) {
      colnames(tbls[[i]]) <- c("Gene", names(.data)[i])
      res <- merge(res, tbls[[i]], by = "Gene", all = T)
    }
  }
  res <- merge(res, data.frame(Gene = .genes[[1]], Something = 0, 
                               stringsAsFactors = F), by = "Gene", all = T)
  res <- res[, -ncol(res)]
  res[is.na(res)] <- 0
  if (!.ambig) {
    res <- .fix.ambig(res, .ambig)
  }
  if (.norm) {
    if (length(.genes) == 1) {
      res[, -1] <- apply(as.matrix(res[, -1]), 2, function(col) col/sum(col))
    }
    else {
      res <- res/sum(res)
    }
  }
  res
}


entropy.seq.internal <- function (.data, .genes = HUMAN_TRBV, .frame = c("all", "in", 
                                                 "out"), .quant = c(NA, "read.count", "umi.count", "read.prop", 
                                                                    "umi.prop"), .ambig = F) 
{
  if (class(.data) == "list") {
    return(sapply(.data, entropy.seq.internal, .quant = .quant, .frame = .frame, 
                  .genes = .genes, .ambig = .ambig))
  }
  .data <- get.frames(.data, .frame)
  if (has.class(.genes, "list") && length(.genes) == 2) {
    entropy(gene_usage_internal(.data, .genes = .genes, .quant = .quant, 
                      .ambig = .ambig))
  }
  else {
    entropy(as.matrix(gene_usage_internal(.data, .genes = .genes, .quant = .quant, 
                                .ambig = .ambig)[, -1]))
  }
}

