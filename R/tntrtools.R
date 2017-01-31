# Export a continous character matrix to TNT format.
# Species names should be given as row names. Character names (optional) are given by column names.
# `dec` is used to round the values.
writecont.tnt <- function(mat, file, dec = 3) {
  out_file <- file(file, "w")
  n_char <- ncol(mat)
  n_tip <- nrow(mat)
  writeLines("NSTATES CONT ;", out_file)
  writeLines("XREAD", out_file)
  tnt_comment <- paste0("\'TNT file generated with the R function `writecont.tnt` on ", date(), ".\'")
  writeLines(tnt_comment, out_file)
  writeLines(paste(n_char, n_tip), out_file)
  writeLines("&[continuous]", out_file)
  tip_names <- rownames(mat)
  tip_names <- gsub(" ", "_", tip_names)
  rownames(mat) <- tip_names
  cnames <- colnames(mat)
  write.table(round(mat, dec), file = out_file, append = TRUE, col.names = FALSE, quote = FALSE)
  if (!is.null(cnames)) {
    writeLines(";\n\nCNAMES", out_file)
    for (i in 0:(length(cnames) - 1)) {
      writeLines(paste0("{  ", i, " ", cnames[i + 1], ";"), out_file)
    }
  }
  writeLines(";\nPROC /;\n", out_file)
  close.connection(out_file)
}

# Write 2D or 3D landmark data to a file in TNT format. `A` can be either an
# array with landmarks from a single configuration or a list of arrays of
# different configurations. Some species may not be included in all the arrays.
# If that happens, `allow_missing = TRUE` will drop those species from the
# entire dataset. `allow_missing = TRUE` will keep all the species and fill in
# with question marks as necessary. A file name to write out a TNT log can be
# given with the `log` parameter. writeland.tnt will also add an ECHO command.
# `dec` is used to round the values
writeland.tnt <- function (A, file, dec = 3, allow_missing = FALSE, log = NULL) {
  .writearray.tnt <- function (A, n, global_names, file) {
  # Write an array of landmark data in TNT format to an open file connection.
  A <- round(A, dec)
  k <- dim(A)[2] # Number of dimensions
  p <- dim(A)[1] # Number of landmarks
  dimnames(A)[[3]] <- gsub(" ", "_", dimnames(A)[[3]]) -> item_names

  if (k == 2) {
    writeLines("& [landmark 2D]", file)
  } else {
    writeLines("& [landmark 3D]", file)
    }
  for (i in 1:n) { # Iterate over terminal items
    this_line <- global_names[i]
    if (global_names[i] %in% item_names) {
      for (j in 1:p) { # Iterate over landmarks
        coords <- sapply(A[j, , global_names[i]], function (x) ifelse ((x >= 0), paste0("+", x), x))
        coords <- paste(coords, collapse = ",")
        coords <- gsub("NA,NA,NA", "?", coords) # This is for missing landmarks
        this_line <- paste(this_line, coords, collapse = "  ")
      }
    } else {
        this_line <- paste(this_line, paste(rep("?", p), collapse = " "))
      }
    writeLines(this_line, file)
  }
  }

  file <- file(file, "w")
  if (! is.null(log)) writeLines(paste0("LOG ", log, ";\nECHO =;"), file)
  writeLines("NSTATES CONT;", file)
  writeLines("XREAD", file)
  tnt_comment <- paste0("\'TNT file generated with the R function `writeland.tnt` on ", date(), ".\'")
  writeLines(tnt_comment, file)
  if (class(A) == "list") {
    item_names <- vector(mode = "character")
    if (allow_missing) {
      for (a in A) {
        item_names <- c(item_names, setdiff(dimnames(a)[[3]], item_names))
      }
    } else {
      item_names <- Reduce(intersect, lapply(A, function (x) dimnames(x)[[3]]))
    }
    item_names <- gsub(" ", "_", item_names)
    n <- length(item_names) # Number of terminal items (sp, ssp, individuals, &c)
    cnames <- names(A)
    writeLines(paste(length(A), n), file)
    for (a in A) {
      .writearray.tnt(a, n, item_names, file)
    }
  } else {
    item_names <- dimnames(A)[[3]]
    item_names <- gsub(" ", "_", item_names)
    n <- length(item_names)
    cnames <- NULL
    writeLines(paste("1", n), file)
    .writearray.tnt(A, n, item_names, file)
  }
  if (! is.null(cnames)) {
    writeLines(";\n\nCNAMES", file)
    for (i in 0:(length(cnames) - 1)) {
      writeLines(paste0("{  ", i, " ", cnames[i + 1], ";"), file)
    }
  }
  writeLines(";\n", file)
  if (! is.null(log)) writeLines("LOG /;", file)
  writeLines("PROC /;\n", file)
  close.connection(file)
}

# Write a phylo or multiPhylo object in the TNT parenthetical tree format. TNT
# cannot read trees without a matrix in memory, so this function offers the
# option to create a dummy matrix. Trees can have diffent tip labels.
write.tree.tnt <- function(tr, out_file, dummy = TRUE) {
  require(ape)
  out_file <- file(out_file, "w")
  if (class(tr) == "multiPhylo") {
    n_tr <- length(tr)
    tip_labels <- Reduce(union, lapply(tr, function(x)
      x$tip.label))
  } else {
    n_tr <- 1
    tip_labels <- tr$tip.label
    tr <- list(tr)
  }

  if (dummy) {
    writeLines("NSTATES 8 ;", out_file)
    writeLines("XREAD", out_file)
    writeLines(paste("1", length(tip_labels)), out_file)
    mat <- cbind(tip_labels, 0)
    write.table(
      mat,
      file = out_file,
      col.names = FALSE,
      row.names = FALSE,
      append = TRUE,
      quote = FALSE
    )
    writeLines(";", out_file)
  }

  tr <-
    lapply(tr, function (x) {
      x$edge.length <- NULL
      x$node.label <- NULL
      x
    })

  writeLines("TREAD", out_file)
  for (i in 1:n_tr) {
    nwkstr <- write.tree(tr[[i]])
    tntstr <- gsub(",", " ", nwkstr)
    tntstr <- gsub(";", "*", tntstr)
    writeLines(tntstr, out_file)
  }
  writeLines(";\nPROC/;\n", out_file)
  close.connection(out_file)
}
