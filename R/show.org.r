beauty.org <- function(x, beauti = c("e", "m", "s")) {
  x[is.na(x)] <- "NA"
  if (beauti == "s") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("\\*.*\\*", x)+1)/2) # bold seulement si != de "" et si pas de bold
    if (length(x[!y]) != 0) x[!y] <- sub("(^ *)([:alpha]*)", "\\1\\*\\2", sub("([:alpha:]*)( *$)", "\\1\\*\\2", x[!y]))
    if (length(x[y]) != 0) x[y] <- sub("(^ *$)", "\\1 ", x[y]) # rajouter suffisamment d'espaces lorsque la case est vide pour l'alignement globale
  }
  if (beauti == "e") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("/.*/", x)+1)/2) # it seulement si != de "" et si pas de it
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1/\\2", sub("([:alpha:]*)( *$)", "\\1/\\2", x[!y]))
    if (length(x[y]) != 0) x[y] <- sub("(^ *$)", "\\1 ", x[y]) # rajouter suffisamment d'espaces lorsque la case est vide pour l'alignement globale
  }
  if (beauti == "m") {
    y <- as.logical((regexpr("^ *$", x)+1)/2) | as.logical((regexpr("=.*=", x)+1)/2) # it seulement si != de "" et si pas de mono
    if (length(x[!y]) != 0) x[!y] <-sub("(^ *)([:alpha]*)", "\\1=\\2", sub("([:alpha:]*)( *$)", "\\1=\\2", x[!y]))
    if (length(x[y]) != 0) x[y] <- sub("(^ *$)", "\\1 ", x[y]) # rajouter suffisamment d'espaces lorsque la case est vide pour l'alignement globale
  }
  return(x)
}

header.org <- function(caption = NULL, caption.level = NULL) {
  res <- ""
  if (is.null(caption.level))
    caption.level <- ""
  if (!is.null(caption)) {
    if (is.numeric(caption.level) & caption.level > 0) {
      res <- paste(paste(rep("*", caption.level), collapse = ""), caption, "\n")
    } else if (is.character(caption.level) & caption.level %in% c("s", "e", "m")) {
      if (caption.level == "s")
        res <- paste(beauty.org(caption, "s"), "\n", sep = "")
      else if (caption.level == "e")
        res <- paste(beauty.org(caption, "e"), "\n", sep = "")
      else if (caption.level == "m")
        res <- paste(beauty.org(caption, "m"), "\n", sep = "")
    } else if (caption.level == "none")
      res <- paste(caption, "\n", sep = "")
    else
      res <- paste("#+CAPTION: ", caption, "\n", sep = "")
    }
  return(res)
}

show.org.table <- function(x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = FALSE, footer = FALSE, align = NULL, col.width = 1, style = NULL, lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", ...) {

  x <- tocharac(x, include.rownames, include.colnames, rownames, colnames, format, digits, decimal.mark, na.print)
  nrowx <- nrow(x)
  ncolx <- ncol(x)
  
  if (!is.null(style)) {
    style <- expand(style, nrowx, ncolx)
    style[!(style %in% c("s", "e", "m"))] <- ""
    style[style == "s"] <- "*"
    style[style == "e"] <- "/"
    style[style == "m"] <- "="
  } else {
    style <- ""
  }
  before_cell_content <- after_cell_content <- style
  before_cell_content <- paste.matrix(" ", before_cell_content, sep = "")
  after_cell_content <- paste.matrix(after_cell_content, " ", sep = "")

  line_separator <- FALSE
  line_separator_pos <- NULL
  if (is.logical(header) & header)
    header <- 1
  if (header > 0) {
    line_separator_pos <- min(c(header, nrowx))
    line_separator <- TRUE
  }

  csep <- matrix("+", nrowx+1, ncolx+1)
  csep[, 1] <- "|"
  csep[, ncol(csep)] <- "|"
  results <- print.character.matrix(x, line_separator = line_separator, line_separator_pos = line_separator_pos, hsep = "-", vsep = "|", csep = csep, before_cell_content = before_cell_content, after_cell_content = after_cell_content, print = FALSE)

  cat(header.org(caption = caption, caption.level = caption.level))
  cat(results, sep = "\n")
}

show.org.list <- function(x, caption = NULL, caption.level = NULL, list.type = "bullet", ...) {
  indent.mark <- "  "
  if (list.type == "bullet") mark <- rep("-", length(x))
  if (list.type == "number") mark <- paste(seq(1, length(x), 1), ".", sep = "")
  if (list.type == "none")  { mark <- rep("", length(x)); indent.mark = ""}
  if (list.type == "label") {
    if (is.null(names(x))) {
      namesx <- paste("[ [ ", 1:length(x), " ] ]", sep = "")
    } else {
      namesx <- names(x)
    }
    mark <- paste("- ", namesx, " ::", sep = "")
    indent.mark = "  "
  }
  
  charac.x <- vector("character", length(x))
  for (i in 1:length(x)) {
    tmp <- x[[i]]
    tmp <- gsub('\t|(*COMMIT)(*FAIL)', indent.mark, tmp, perl = TRUE)
    charac.x[i] <- sub("(^ *)", paste("\\1", mark[i], " ", sep = ""), tmp)
  }
  cat(header.org(caption = caption, caption.level = caption.level))
  cat(charac.x, sep = "\n")
}
