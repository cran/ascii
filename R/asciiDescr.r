##' @export
##' @method ascii meanscomp
ascii.meanscomp <- function (x, header = TRUE, caption = NULL, include.rownames = TRUE, include.colnames = TRUE, ...) {
    rlab <- ifelse(is.null(x$row.label), x$row, x$row.label)
    clab <- ifelse(is.null(x$column.label), x$column, x$column.label)
    msg1 <- gettext("Mean value of", domain = "R-descr")
    msg2 <- gettext("according to", domain = "R-descr")
    msg <- switch(1 + (caption == ""), caption, paste(msg1, " \"", clab, "\" ", msg2, " \"", rlab, "\"", sep = ""))
    ascii(x$table, caption = msg, include.rownames = include.rownames, include.colnames = include.colnames, header = header, ...)
}

##' @export
##' @method ascii CrossTable
ascii.CrossTable <- function (x, ...) {
  t <- x$t
  CPR <- x$prop.row
  CPC <- x$prop.col
  CPT <- x$prop.tbl
  GT <- x$gt
  RS <- x$rs
  CS <- x$cs
  TotalN <- x$total.n
  CST <- x$chisq
  CSTc <- x$chisq.corr
  FTt <- x$fisher.ts
  FTl <- x$fisher.lt
  FTg <- x$fisher.gt
  McN <- x$mcnemar
  McNc <- x$mcnemar.corr
  ASR <- x$asr
  RowData <- x$RowData
  ColData <- x$ColData
  digits <- x$digits
  max.width <- x$max.width
  vector.x <- x$vector.x
  expected <- x$expected
  prop.r <- (is.na(CPR[1]) == FALSE)
  prop.c <- (is.na(CPC[1]) == FALSE)
  prop.t <- (is.na(CPT[1]) == FALSE)
  chisq <- (is.na(CST[1]) == FALSE)
  prop.chisq <- x$prop.chisq
  fisher <- (class(FTt) == "htest")
  resid <- x$resid
  sresid <- x$sresid
  asresid <- x$asresid
  mcnemar <- x$print.mcnemar
  missing.include <- x$missing.include
  format <- x$format
  outDec <- getOption("OutDec")
  nsep <- "  | "
  if (format == "SAS") {
    resid <- sresid <- asresid <- FALSE
    hdd <- 1
    psep <- "  | "
  } else {
    if (format == "SPSS") {
      hdd <- 100
      psep <- "% | "
    } else {
      stop("unknown format")
    }
  }
  if (vector.x)
    expected <- prop.chisq <- prop.c <- prop.t <- resid <- sresid <- asresid <- FALSE
  ColTotal <- gettext("Total", domain = "R-descr")
  RowTotal <- ColTotal
  CWidth <- max(digits + 2, c(nchar(t), nchar(dimnames(t)[[2]]),
                              nchar(RS), nchar(CS), nchar(RowTotal)))
  RWidth <- max(c(nchar(dimnames(t)[[1]]), nchar(ColTotal)))
  if (is.na(RowData) == FALSE)
    RWidth <- max(RWidth, nchar(RowData))
  RowSep <- paste(rep("-", CWidth + 2), collapse = "")
  RowSep1 <- paste(rep("-", RWidth + 1), collapse = "")
  SpaceSep1 <- paste(rep(" ", RWidth), collapse = "")
  SpaceSep2 <- paste(rep(" ", CWidth), collapse = "")
  FirstCol <- formatC(dimnames(t)[[1]], width = RWidth, format = "s")
  ColTotal <- formatC(ColTotal, width = RWidth, format = "s")
  RowTotal <- formatC(RowTotal, width = CWidth, format = "s")
  caption <- gettext("Cell Contents", domain = "R-descr")
  if (format == "SAS") {
    content <- c("N",
                 ifelse(expected, "Expected N", NA),
                 ifelse(prop.chisq, "Chi-square contribution", NA),
                 ifelse(prop.r, "N / Row Total", NA),
                 ifelse(prop.c, "N / Col Total", NA),
                 ifelse(prop.t, "N / Table Total", NA)
                 )
    content <- as.list(content[!is.na(content)])
  } else if (format == "SPSS") {
    content <- c("Count",
                 ifelse(expected, "Expected Values", NA),
                 ifelse(prop.chisq, "Chi-square contribution", NA),
                 ifelse(prop.r, "Row Percent", NA),
                 ifelse(prop.c, "Column Percent", NA),
                 ifelse(prop.t, "Total Percent", NA),
                 ifelse(resid, "Residual", NA),
                 ifelse(sresid, "Std Residual", NA),
                 ifelse(asresid, "Adj Std Resid", NA)
                 )
    content <- as.list(content[!is.na(content)])
  }
  if (vector.x) {
    stop("For single vector description, use freq() function in package descr.\n")
  }
  nelements <- 1 + expected + prop.chisq + prop.r + prop.c +
    prop.t + resid + sresid + asresid
  nr <- nrow(t) * nelements + 1 + prop.c
  nc <- ncol(t)
  m <- matrix(nrow = nr, ncol = nc + 1)
  rnames <- rownames(t)
  n.rnames <- vector(mode = "numeric", length = nrow(t) + 1)
  k <- 1
  for (i in 1:nrow(t)) {
    for (l in 1:nc) m[k, l] <- formatC(t[i, l], format = "d")
    m[k, nc + 1] <- RS[i]
    n.rnames[i] <- 1
    k <- k + 1
    if (expected) {
      for (l in 1:nc) m[k, l] <- formatC(CST$expected[i,
                                                      l], digits = 1, format = "f", decimal.mark = outDec)
      m[k, nc + 1] <- " "
      n.rnames[i] <- n.rnames[i] + 1
      k <- k + 1
    }
    if (prop.chisq) {
      for (l in 1:nc) m[k, l] <- formatC((((CST$expected[i,
                                                         l] - t[i, l])^2)/CST$expected[i, l]), digits = digits,
                                         format = "f", decimal.mark = outDec)
      m[k, nc + 1] <- " "
      n.rnames[i] <- n.rnames[i] + 1
      k <- k + 1
    }
    if (prop.r) {
      for (l in 1:nc) m[k, l] <- formatC(CPR[i, l] * hdd,
                                         digits = digits, format = "f", decimal.mark = outDec)
      m[k, nc + 1] <- formatC(hdd * RS[i]/GT, digits = digits,
                              format = "f", decimal.mark = outDec)
      n.rnames[i] <- n.rnames[i] + 1
      k <- k + 1
    }
    if (prop.c) {
      for (l in 1:nc) m[k, l] <- formatC(CPC[i, l] * hdd,
                                         digits = digits, format = "f", decimal.mark = outDec)
      m[k, nc + 1] <- " "
      n.rnames[i] <- n.rnames[i] + 1
      k <- k + 1
    }
    if (prop.t) {
      for (l in 1:nc) m[k, l] <- formatC(CPT[i, l] * hdd,
                                         digits = digits, format = "f", decimal.mark = outDec)
      m[k, nc + 1] <- " "
      n.rnames[i] <- n.rnames[i] + 1
      k <- k + 1
    }
    if (resid) {
      for (l in 1:nc) m[k, l] <- formatC(CST$observed[i,
                                                      l] - CST$expected[i, l], digits = digits, format = "f",
                                         decimal.mark = outDec)
      m[k, nc + 1] <- " "
      n.rnames[i] <- n.rnames[i] + 1
      k <- k + 1
    }
    if (sresid) {
      for (l in 1:nc) m[k, l] <- formatC(CST$residual[i,
                                                      l], digits = digits, format = "f", decimal.mark = outDec)
      m[k, nc + 1] <- " "
      n.rnames[i] <- n.rnames[i] + 1
      k <- k + 1
    }
    if (asresid) {
      for (l in 1:nc) m[k, l] <- formatC(ASR[i, l], digits = digits,
                                         format = "f", decimal.mark = outDec)
      m[k, nc + 1] <- " "
      n.rnames[i] <- n.rnames[i] + 1
      k <- k + 1
    }
  }
  ColTotal <- gettext("Total", domain = "R-descr")
  RowTotal <- ColTotal
  rnames[length(rnames) + 1] <- RowTotal
  for (l in 1:nc) m[k, l] <- formatC(c(CS[l]), format = "d")
  m[k, nc + 1] <- formatC(GT, format = "d")
  n.rnames[length(n.rnames)] <- 1
  if (prop.c) {
    k <- k + 1
    for (l in 1:nc) m[k, l] <- formatC(hdd * CS[l]/GT, digits = digits,
                                       format = "f", decimal.mark = outDec)
    n.rnames[length(n.rnames)] <- 2
  }
  colnames(m) <- c(colnames(t), ColTotal)
  nc <- nc + 1
  mcolnames <- colnames(m)
  res.m <- ascii(m, include.colnames = T, header = T, lgroup = rnames, n.lgroup = n.rnames, lstyle = "s")
  res.t <- NULL
  if (chisq) {
    res.t <- list("Pearson's Chi-squared test" = paste(gettext("Chi^2 = ", domain = "R-descr"), format(CST$statistic),
                    ", ", gettext("d.f. = ", domain = "R-descr"), format(CST$parameter),
                    ", ", gettext("p = ", domain = "R-descr"), format(CST$p.value), sep = ""))
    if (all(dim(t) == 2)) {
      res.t <- c(res.t, "Pearson's Chi-squared test with Yates' continuity correction" = paste(gettext("Chi^2 = ", domain = "R-descr"),
                          format(CSTc$statistic),
                          ", ", gettext("d.f. = ", domain = "R-descr"),
                          format(CSTc$parameter), ", ", gettext("p = ", domain = "R-descr"),
                          format(CSTc$p.value), sep = ""))
    }
  }
  if (is.na(McN[1]) == FALSE) {
    res.t <- c(res.t, "McNemar's Chi-squared test" = paste(gettext("Chi^2 = ", domain = "R-descr"), format(McN$statistic),
                        ", ", gettext("d.f. = ", domain = "R-descr"), format(McN$parameter),
                        ", ", gettext("p = ", domain = "R-descr"), format(McN$p.value), sep = ""))
    if (is.na(McNc[1]) == FALSE) {
      res.t <- c(res.t, "McNemar's Chi-squared test with continuity correction" = paste(gettext("Chi^2 = ", domain = "R-descr"), format(McNc$statistic),
                          ", ", gettext("d.f. = ", domain = "R-descr"),
                          format(McNc$parameter), ", ", gettext("p = ", domain = "R-descr"),
                          format(McNc$p.value), sep = ""))
    }
  }
  if (fisher) {
    res.t <- c(res.t, "Fisher's Exact Test for Count Data" = paste(gettext("Alternative hypothesis: two.sided",
                        domain = "R-descr"), gettext(", p = ", domain = "R-descr"), format(FTt$p.value), sep = ""))
  }
  if (format == "SPSS") {
    if (any(dim(t) >= 2) & any(chisq, mcnemar, fisher)) {
      MinExpF = min(CST$expected)
      res.t <- c(res.t, list("Minimum expected frequency" = format(MinExpF)))
      NMinExpF = length(CST$expected[which(CST$expected <
        5)])
      if (NMinExpF > 0) {
        NCells = length(CST$expected)
        res.t <- c(res.t, "Cells with Expected Frequency < 5" = paste(format(NMinExpF), " ", "of ",
                            NCells, " (", format(100 * NMinExpF/NCells),
                            "%)", sep = ""))
      }
      ## cat("\n")
    }
  }
  content <- ascii(content, caption = "Cell Contents", caption.level = "s")
  if (!is.null(res.t))
    res.t <- ascii(res.t, caption = "Statistics for All Table Factors", caption.level = "s", list.type = "label")
  res <- asciiMixed$new(args = list(content, res.m, res.t))
  return(res)
}

##' @export
##' @method ascii freqtable
ascii.freqtable <- function (x, header = TRUE, footer = TRUE, digits = c(0, 2, 2), format = "f", na.print = "", include.rownames = TRUE, include.colnames = TRUE, caption = attr(x, "xlab"), ...) {
  class(x) <- "matrix"
  ascii(x, header = header, footer = footer, include.rownames = include.rownames, include.colnames = include.colnames, caption = caption, digits = digits, format = format, na.print = na.print, ...)
}

##' Ascii formatting for a microbenchmark
##'
##' The default implementation returns an asciiMixed object with the
##' units for the first element.
##' @export
##' @method ascii microbenchmark
##' @param x an object of class 'microbenchmark'
##' @param unit What unit to print the timings in. Default value taken
##'     from the option 'microbenchmark.unit'
##' @param order If present, order results according to this column of
##'     the output.
##' @param signif If present, limit the limit of significant digits
##'     shown.
##' @param row.names Argument passed to ascii
##' @param caption logical; if not NULL, then add caption with units
##'     specified; otherwise, add units as part of an asciiMixed
##'     object.
##' @param ... Other parameters to pass to ascii for the summary table
##' @return ascii object
ascii.microbenchmark <- function (x, unit, order, signif, row.names=FALSE,
                                  caption=NULL, ...) {
    s <- summary(x, unit = unit)
    timing_cols <- c("min", "lq", "median", "uq", "max", "mean")
    if (!missing(signif)) {
        s[timing_cols] <- lapply(s[timing_cols], base::signif,
                                 signif)
    }
    if (!missing(order)) {
        if (order %in% colnames(s)) {
            s <- s[order(s[[order]]), ]
        }
        else {
            warning("Cannot order results by", order, ".")
        }
    }
    unit.caption <- sprintf("Unit: %s", attr(s, "unit"))
    if (!is.null(caption)) {
        ascii(s, ..., row.names = row.names, caption = paste0(caption, " (", unit.caption, ")"))
    } else {
        asciiMixed$new(args=list(ascii(unit.caption),
                                 ascii(s, ..., row.names = row.names)))
    }
}

##' Translation of the printCoefmat function for ascii
##'
##' Compared with printCoefmat, this drops the quote and right
##' arguments, and adds include.rownames, include.colnames and header
##' default arguments.
##' @export
##' @param x coefficient summary table that is suitable for
##'     printCoefmat
##' @param digits minimum number of significant digits to be used for
##'     most numbers.
##' @param signif.stars locial; if 'TRUE', P-values are additionally
##'     encoded visually as 'significance stars' in order to help
##'     scanning of long coefficient tables.  It defaults to the
##'     'show.signif.stars' slot of 'options'.
##' @param signif.legend logical; if 'TRUE', a legend for the
##'     'significance stars' is printed provided 'signif.stars =
##'     TRUE'.
##' @param dig.tst minimum number of significant digits for the test
##'     statistics, see 'tst.ind'.
##' @param cs.ind indices (integer) of column numbers which are (like)
##'     *c*oefficients and *s*tandard errors to be formatted together.
##' @param tst.ind indices (integer) of column numbers for test
##'     statistics.
##' @param zap.ind indices (integer) of column numbers which should be
##'     formatted by zapsmall, i.e., by 'zapping' values close to 0.
##' @param P.values logical or 'NULL'; if 'TRUE', the last column of
##'     'x' is formatted by format.pval as P values.  If 'P.values =
##'     NULL', the default, it is set to 'TRUE' only if
##'     'options("show.coef.Pvalue")' is 'TRUE' _and_ 'x' has at least
##'     4 columns _and_ the last column name of 'x' starts with
##'     '"Pr("'.
##' @param has.Pvalue logical; if 'TRUE', the last column of 'x'
##'     contains P values; in that case, it is printed if and only if
##'     'P.values' (above) is true.
##' @param eps.Pvalue lower threshold for reporting p-values.
##' @param na.print a character string to code NA values in printed
##'     output.
##' @param include.rownames argument passed to ascii
##' @param include.colnames argument passed to ascii
##' @param header argument passed to ascii
##' @param ... other argments passed to ascii
##' @importFrom stats symnum
##' @return ascii object. This is character, rather than numeric.
asciiCoefmat <- 
function (x, digits = max(3L, getOption("digits") - 2L), signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, dig.tst = max(1L, min(5L, digits - 
        1L)), cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(), 
    P.values = NULL, has.Pvalue = nc >= 4L && length(cn <- colnames(x)) && 
        substr(cn[nc], 1L, 3L) %in% c("Pr(", "p-v"), eps.Pvalue = .Machine$double.eps, 
    na.print = "NA", 
    include.rownames=TRUE, include.colnames=TRUE, header=TRUE, ...) 
{
    if (is.null(d <- dim(x)) || length(d) != 2L) 
        stop("'x' must be coefficient matrix/data frame")
    nc <- d[2L]
    if (is.null(P.values)) {
        scp <- getOption("show.coef.Pvalues")
        if (!is.logical(scp) || is.na(scp)) {
            warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
            scp <- TRUE
        }
        P.values <- has.Pvalue && scp
    }
    else if (P.values && !has.Pvalue) 
        stop("'P.values' is TRUE, but 'has.Pvalue' is not")
    if (has.Pvalue && !P.values) {
        d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
        nc <- nc - 1
        has.Pvalue <- FALSE
    }
    else xm <- data.matrix(x)
    k <- nc - has.Pvalue - (if (missing(tst.ind)) 
        1
    else length(tst.ind))
    if (!missing(cs.ind) && length(cs.ind) > k) 
        stop("wrong k / cs.ind")
    Cf <- array("", dim = d, dimnames = dimnames(xm))
    ok <- !(ina <- is.na(xm))
    for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
    if (length(cs.ind)) {
        acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
        if (any(ia <- is.finite(acs))) {
            digmin <- 1 + if (length(acs <- acs[ia & acs != 0])) 
                floor(log10(range(acs[acs != 0], finite = TRUE)))
            else 0
            Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - 
                digmin)), digits = digits)
        }
    }
    if (length(tst.ind)) 
        Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst), 
            digits = digits)
    if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc)))) 
        for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
    ok[, tst.ind] <- FALSE
    okP <- if (has.Pvalue) 
        ok[, -nc]
    else ok
    x1 <- Cf[okP]
    dec <- getOption("OutDec")
    if (dec != ".") 
        x1 <- chartr(dec, ".", x1)
    x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
    if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
        Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L, 
            digits - 1L))
    }
    if (any(ina)) 
        Cf[ina] <- na.print
    if (any(inan <- is.nan(xm))) 
        Cf[inan] <- "NaN"
    if (P.values) {
        if (!is.logical(signif.stars) || is.na(signif.stars)) {
            warning("option \"show.signif.stars\" is invalid: assuming TRUE")
            signif.stars <- TRUE
        }
        if (any(okP <- ok[, nc])) {
            pv <- as.vector(xm[, nc])
            Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst, 
                eps = eps.Pvalue)
            signif.stars <- signif.stars && any(pv[okP] < 0.1)
            if (signif.stars) {
                Signif <- symnum(pv, corr = FALSE, na = FALSE, 
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  symbols = c("***", "**", "*", ".", " "))
                Cf <- cbind(Cf, format(Signif))
            }
        }
        else signif.stars <- FALSE
    }
    else signif.stars <- FALSE
    out <- ascii(Cf, na.print = na.print,
                 include.rownames=include.rownames, include.colnames=include.colnames,
                 header=header,
                 ...)
    if (signif.stars && signif.legend) {
        sleg <- attr(Signif, "legend")
        out <- asciiMixed$new(args=list(out,
                                        ascii(paste("Signif. codes: ", sleg))))
    }
    out
}
