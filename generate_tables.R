# =====================================================================
#  generate_tables.R
#
#  Builds the manuscript tables straight from the raw simulation CSVs.
#  For every table family it writes one LaTeX file (tables.tex) with
#  one table per page and, if a LaTeX installation is on the PATH, the
#  matching tables.pdf.
#
#  Four families (one folder each):
#   1. PLS_linear_H1_rpgrp_tables      Table 1 and  Table 25-31     
#                                                     
#   2. PLS_logistic_H2_grp_tables      Table 2 and Table S18     
#                                                     
#   3. PLS_linear_H1_bandwidth_tables  Table S1-S24  
#   4. PLS_logistic_H2_bandwidth_tables Table S32-S37  
#
#  Usage
#     Rscript generate_tables.R
#     PLS_CSV_DIR=/path/to/csv PLS_OUT_DIR=/path/to/out Rscript generate_tables.R
#

## ---------------------------------------------------------------------
##   input / output location
## ---------------------------------------------------------------------
csv_dir   <- Sys.getenv("PLS_CSV_DIR",
                        "E:/桌面/代码原始输出结果/代码原始输出结果")
out_root  <- Sys.getenv("PLS_OUT_DIR", ".")
build_pdf <- TRUE
latex_cmd <- Sys.getenv("PLS_LATEX", "pdflatex")

csv_file <- list(
  h1_rpgrp = file.path(csv_dir, "PLS_linear_H1_rpgrp_result.csv"),
  h2_grp   = file.path(csv_dir, "PLS_logistic_H2_grp_result.csv"),
  h1_bw    = file.path(csv_dir, "PLS_linear_H1_result.csv"),
  h2_bw    = file.path(csv_dir, "PLS_logistic_H2_result.csv")
)
out_dir <- list(
  h1_rpgrp = file.path(out_root, "PLS_linear_H1_rpgrp_tables"),
  h2_grp   = file.path(out_root, "PLS_logistic_H2_grp_tables"),
  h1_bw    = file.path(out_root, "PLS_linear_H1_bandwidth_tables"),
  h2_bw    = file.path(out_root, "PLS_logistic_H2_bandwidth_tables")
)

## ---------------------------------------------------------------------
## 1.  helpers
## ---------------------------------------------------------------------
PCOLS    <- c("p=50", "p=100", "p=300", "p=600", "p=900", "p=1200")
PVALS    <- c(50, 100, 300, 600, 900, 1200)
EST_NAME <- c("Lasso", "SCAD", "MCP")            # estimate_type = 1, 2, 3
EST_TEX  <- paste0("post-", EST_NAME)
ERR_TEX  <- c("N(0,1)", "t_3")                   # error_type   = 1, 2
ERR_TAG  <- c("N", "t3")

## read one simulation file; the first (unnamed) column holds the statistic
read_sim <- function(path) {
  if (!file.exists(path)) stop("CSV not found: ", path, call. = FALSE)
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- "label"
  d
}

## map a unique key (built from `keycols`) to the row of six p-values
make_index <- function(dat, keycols) {
  k <- do.call(paste, c(lapply(keycols, function(cc) as.character(dat[[cc]])),
                        list(sep = "|")))
  if (anyDuplicated(k))
    warning("duplicated keys for ", paste(keycols, collapse = ", "), call. = FALSE)
  idx <- seq_len(nrow(dat))
  names(idx) <- k
  list(dat = dat, idx = idx)
}

lookup <- function(tb, ...) {
  key <- paste(vapply(list(...), as.character, character(1)), collapse = "|")
  i <- tb$idx[[key]]
  if (is.null(i)) stop("no row for key: ", key, call. = FALSE)
  unlist(tb$dat[i, PCOLS], use.names = FALSE)
}

## first estimate_type whose row is complete (used for RP_n / GRP_n)
lookup_first <- function(tb, label, ...) {
  for (est in 1:3) {
    v <- lookup(tb, label, ..., est)
    if (!anyNA(v)) return(v)
  }
  stop("no non-missing values for ", label, call. = FALSE)
}

f3      <- function(x) sprintf("%.3f", x)
num_row <- function(v) paste(vapply(v, f3, character(1)), collapse = "&")

## Sigma of the reference tables
pho_tex <- function(ph) {
  if (ph == 0) return("I_p")
  paste0("(", format(ph, nsmall = 1), "^{|i-j|})_{p\\times p}")
}

## ---------------------------------------------------------------------
## 2.  preamble, header blocks and one block of rows
## ---------------------------------------------------------------------
doc_preamble <- function(thetable = NULL) {
  paste0(
    "\n\\documentclass[10pt]{article}\n",
    "\\usepackage[a4paper,margin=19mm]{geometry}\n",
    "\\usepackage[T1]{fontenc}\n",
    "\\usepackage{amsmath,amssymb}\n",
    "\\usepackage{color}\n",
    "\\usepackage[font=small,labelfont=bf,labelsep=colon,justification=raggedright,\n",
    "            singlelinecheck=false,skip=6pt]{caption}\n",
    "\\usepackage{newtxtext,newtxmath}\n",
    "\\newsavebox{\\TAB}\n",
    if (!is.null(thetable)) paste0("\\renewcommand{\\thetable}{", thetable, "}\n"),
    "\\begin{document}\n")
}

doc_assemble <- function(chunks, thetable = NULL) {
  paste(c(doc_preamble(thetable), chunks, "\n\\end{document}\n"), collapse = "\n")
}

## header of family 1: "a" column, n and p rows with multicolumn cells
head_h1_rpgrp <- function(n) {
  nc <- rep(paste0("\\multicolumn{1}{c}{n=", n, "}"), 6)
  pc <- paste0("\\multicolumn{1}{c}{p=", PVALS, "}")
  paste0(
    "\\begin{table}[ht!]\n\\color{blue}\n\\savebox{\\TAB}{\\scriptsize\n",
    "\\renewcommand{\\arraystretch}{0.3}\\tabcolsep 0.3cm\n",
    "\\begin{tabular}{*{8}{c}}\n\\hline\n",
    paste(c("", "\\multicolumn{1}{c}{a}", nc[1:4]), collapse = "&"), "\n",
    paste(c("", nc[5:6]), collapse = "&"), "\\\\\n",
    paste(c("", "", pc[1:4]), collapse = "&"), "\n",
    paste(c("", pc[5:6]), collapse = "&"), " \\\\\n\\hline")
}

## header of family 2 (same style, n on a single line)
head_h2_grp <- function(n) {
  nc <- rep(paste0("\\multicolumn{1}{c}{n=", n, "}"), 6)
  pc <- paste0("\\multicolumn{1}{c}{p=", PVALS, "}")
  paste0(
    "\\begin{table}[ht!]\n\\color{blue}\n\\savebox{\\TAB}{\\scriptsize\n",
    "\\renewcommand{\\arraystretch}{0.3}\\tabcolsep 0.3cm\n",
    "\\begin{tabular}{*{8}{c}}\n\\hline\n",
    paste(c("", "\\multicolumn{1}{c}{a}", nc), collapse = "&"), "\\\\\n",
    paste(c("", "", pc[1:4]), collapse = "&"), "\n",
    paste(c("", pc[5:6]), collapse = "&"), " \\\\\n\\hline")
}

## header of family 3 (c_h column, small font)
head_h1_band <- function(n) {
  nc <- rep(paste0("\\multicolumn{1}{c}{n=", n, "}"), 6)
  pc <- paste0("\\multicolumn{1}{c}{p=", PVALS, "}")
  paste0(
    "\\begin{table}[ht!]\n\\color{blue}\n\\savebox{\\TAB}{\\small\n",
    "\\renewcommand{\\arraystretch}{0.86}\\tabcolsep 0.15cm\n",
    "\\begin{tabular}{*{8}{c}}\n\\hline\n",
    paste(c("", "$c_h$", nc), collapse = "&"), "\\\\\n",
    paste(c("", "", pc[1:4]), collapse = "&"), "\n",
    paste(c("", pc[5:6]), collapse = "&"), " \\\\\n\\hline")
}

## header of family 4 (c_h column, footnote size)
head_h2_band <- function(n) {
  paste0(
    "\\begin{table}[!htbp]\n\\color{blue}\n\\savebox{\\TAB}{\\footnotesize\n",
    "\\renewcommand{\\arraystretch}{1.0} \\tabcolsep 0.42cm\n",
    "\\begin{tabular}{@{}cccccccc@{}}\n\\hline\n",
    paste(c("", "$c_h$", rep(paste0("n=", n), 6)), collapse = "&"), "\\\\\n",
    paste(c("", "", paste0("p=", PVALS)), collapse = "&"), "\\\\\n\\hline")
}

## five rows of one bandwidth block: the label sits on the first line
band_block <- function(lab, chs, vals) {
  rows <- character(length(chs))
  for (i in seq_along(chs)) {
    end  <- if (i == length(chs)) " \\\\[4pt]" else " \\\\"
    lead <- if (i == 1L) paste0(lab, "\n") else ""
    rows[i] <- paste0(lead, "&", sprintf("%.1f", chs[i]), "&",
                      num_row(vals[[i]]), end)
  }
  paste(rows, collapse = "\n")
}

## ---------------------------------------------------------------------
## 3.  family 1:  H1 linear, a = 0 / a = 1 pairs, T_min, T_Fisher, RP_n, GRP_n
## ---------------------------------------------------------------------
build_h1_rpgrp <- function(path) {
  d  <- read_sim(path)
  tb <- make_index(d, c("label", "example", "a", "pho", "error_type", "estimate_type"))
  n  <- unique(d$n)
  chunks <- character(0)
  for (ex in sort(unique(d$example))) {
    for (er in sort(unique(d$error_type))) {
      blocks <- character(0)
      for (ph in sort(unique(d$pho))) {
        sg   <- pho_tex(ph)
        rows <- character(0)
        for (est in 1:3) {
          nm <- EST_NAME[est]
          lab_min <- if (nm == "Lasso")
            paste0("$\\hat{T}_{\\rm{min, post-Lasso}}^C, \\ \\Sigma = ", sg, "$")
          else
            paste0("$\\hat{T}_{\\rm{min}, post-", nm, "}^C, \\ \\Sigma = ", sg, "$")
          lab_fis <- if (nm == "Lasso")
            paste0("$\\hat{T}^C_{\\rm{Fisher, post-Lasso}}, \\ \\Sigma = ", sg, "$")
          else
            paste0("$\\hat{T}^C_{\\rm{Fisher}, post-", nm, "}, \\ \\Sigma = ", sg, "$")
          for (pp in list(list(lab_min, "PLS_min_power"),
                          list(lab_fis, "PLS_fisher_power"))) {
            v0 <- lookup(tb, pp[[2]], ex, 0, ph, er, est)
            v1 <- lookup(tb, pp[[2]], ex, 1, ph, er, est)
            rows <- c(rows, paste0(pp[[1]], "\n&0.0&", num_row(v0),
                                   " \\\\\n&1.0&", num_row(v1), " \\\\[4pt]"))
          }
        }
        for (key in c("rp_power", "grp_power")) {
          lab <- paste0("$", if (key == "rp_power") "RP_n" else "GRP_n",
                        ", \\ \\Sigma = ", sg, "$")
          v0 <- lookup_first(tb, key, ex, 0, ph, er)
          v1 <- lookup_first(tb, key, ex, 1, ph, er)
          rows <- c(rows, paste0(lab, "\n&0.0&", num_row(v0),
                                 " \\\\\n&1.0&", num_row(v1), " \\\\[4pt]"))
        }
        blocks <- c(blocks, paste(rows, collapse = "\n"))
      }
      tail <- paste0(
        "\n\\hline\n\\end{tabular}}\n\\captionsetup{width=\\wd\\TAB}\n\\centering\n",
        "\\caption{Empirical sizes and powers of $\\hat{T}_{\\rm{min}}^C$, ",
        "$\\hat{T}_{\\rm{Fisher}}^C$,\n",
        "$RP_n$, and $GRP_n$ for the bandwidth $h = n^{-2/9}$ under $H_{11}$ with\n",
        "$\\varepsilon \\sim ", ERR_TEX[er], "$ in Study 1 (Example ", ex,
        ").}\\label{table-H1-", ex, ERR_TAG[er], "}\n",
        "\\usebox{\\TAB}\n\\end{table}\n\\clearpage\n")
      chunks <- c(chunks, paste0(head_h1_rpgrp(n), "\n",
                                 paste(blocks, collapse = "\n\\hline\n"), tail))
    }
  }
  doc_assemble(chunks)
}

## ---------------------------------------------------------------------
## 4.  family 2:  H2 logistic, a = 0 / a = 1 pairs, T_min, T_Fisher, GRP_n
## ---------------------------------------------------------------------
build_h2_grp <- function(path) {
  d  <- read_sim(path)
  tb <- make_index(d, c("label", "example", "a", "pho", "estimate_type"))
  n  <- unique(d$n)
  chunks <- character(0)
  for (ex in sort(unique(d$example))) {
    blocks <- character(0)
    for (ph in sort(unique(d$pho))) {
      sg   <- pho_tex(ph)
      rows <- character(0)
      for (est in 1:3) {
        nm <- EST_NAME[est]
        for (pp in list(
          list(paste0("$\\hat{T}^C_{\\rm{min, post-", nm, "}}, \\ \\Sigma = ", sg, "$"),
               "PLS_min_power"),
          list(paste0("$\\hat{T}^C_{\\rm{Fisher, post-", nm, "}}, \\ \\Sigma = ", sg, "$"),
               "PLS_fisher_power"))) {
          v0 <- lookup(tb, pp[[2]], ex, 0, ph, est)
          v1 <- lookup(tb, pp[[2]], ex, 1, ph, est)
          rows <- c(rows, paste0(pp[[1]], "\n&0.0&", num_row(v0),
                                 " \\\\\n&1.0&", num_row(v1), " \\\\[4pt]"))
        }
      }
      lab <- paste0("$GRP_n, \\ \\Sigma = ", sg, "$")
      v0  <- lookup_first(tb, "grp_power", ex, 0, ph)
      v1  <- lookup_first(tb, "grp_power", ex, 1, ph)
      rows <- c(rows, paste0(lab, "\n&0.0&", num_row(v0),
                             " \\\\\n&1.0&", num_row(v1), " \\\\[4pt]"))
      blocks <- c(blocks, paste(rows, collapse = "\n"))
    }
    tail <- paste0(
      "\n\\hline\n\\end{tabular}}\n\\captionsetup{width=\\wd\\TAB}\n\\centering\n",
      "\\caption{Empirical sizes and powers of $\\hat{T}^C_{\\rm{min}}$, ",
      "$\\hat{T}^C_{\\rm{Fisher}}$, and\n",
      "$GRP_n$ for the bandwidth $h = n^{-2/9}$ under $H_{21}$ in Study 2 ",
      "(Example ", ex, ").}\\label{table-H21-ex", ex, "}\n",
      "\\usebox{\\TAB}\n\\end{table}\n\\clearpage\n")
    chunks <- c(chunks, paste0(head_h2_grp(n), "\n",
                               paste(blocks, collapse = "\n\\hline\n"), tail))
  }
  doc_assemble(chunks)
}

## ---------------------------------------------------------------------
## 5.  family 3:  H1 linear, different bandwidths c_h
## ---------------------------------------------------------------------
build_h1_band <- function(path) {
  d   <- read_sim(path)
  tb  <- make_index(d, c("label", "example", "a", "pho", "c_h",
                         "error_type", "estimate_type"))
  n   <- unique(d$n)
  chs <- sort(unique(d$c_h))
  chunks <- character(0)
  for (ex in sort(unique(d$example))) {
    for (er in sort(unique(d$error_type))) {
      for (est in 1:3) {
        nm <- EST_TEX[est]
        groups <- character(0)
        for (ph in sort(unique(d$pho))) {
          blocks <- character(0)
          for (a in c(0, 1)) {
            for (pp in list(list("min", "PLS_min_power"),
                            list("Fisher", "PLS_fisher_power"))) {
              lab <- paste0("$\\hat{T}^C_{\\rm{", pp[[1]], ",", nm, "}}$, ",
                            "$\\Sigma = ", pho_tex(ph), "$, $a = ", a, "$")
              vals <- lapply(chs, function(ch)
                lookup(tb, pp[[2]], ex, a, ph, ch, er, est))
              blocks <- c(blocks, band_block(lab, chs, vals))
            }
          }
          groups <- c(groups, paste(blocks, collapse = "\n"))
        }
        tail <- paste0(
          "\n\\hline\n\\end{tabular}}\n\\captionsetup{width=\\wd\\TAB}\n\\centering\n",
          "\\caption{Empirical sizes and powers of $\\hat{T}^C_{\\rm{min,", nm, "}}$ and\n",
          "$\\hat{T}^C_{\\rm{Fisher,", nm, "}}$ for different bandwidths and ", nm, "\n",
          "under $H_{1", ex, "}$ with $\\varepsilon \\sim ", ERR_TEX[er],
          "$ in Study 1.}\\label{table-S-ex", ex, "-e", er, "-m", est, "}\n",
          "\\usebox{\\TAB}\n\\end{table}\n\\clearpage\n")
        chunks <- c(chunks, paste0(head_h1_band(n), "\n",
                                   paste(groups, collapse = "\n\\hline\n"), tail))
      }
    }
  }
  doc_assemble(chunks, thetable = "S\\arabic{table}")
}

## ---------------------------------------------------------------------
## 6.  family 4:  H2 logistic, different bandwidths c_h
## ---------------------------------------------------------------------
build_h2_band <- function(path) {
  d   <- read_sim(path)
  tb  <- make_index(d, c("label", "example", "a", "pho", "c_h", "estimate_type"))
  n   <- unique(d$n)
  chs <- sort(unique(d$c_h))
  chunks <- character(0)
  for (ex in sort(unique(d$example))) {
    for (est in 1:3) {
      nm <- EST_TEX[est]
      groups <- character(0)
      for (ph in sort(unique(d$pho))) {
        blocks <- character(0)
        for (a in c(0, 1)) {
          for (pp in list(list("min", "PLS_min_power"),
                          list("Fisher", "PLS_fisher_power"))) {
            ## row labels exactly as in the manuscript's Table S32
            lab <- if (pp[[1]] == "min")
              paste0("$\\hat{T}_{\\rm{min, ", nm, "}}^C$, $\\Sigma = ",
                     pho_tex(ph), "$, $a = ", a, "$")
            else
              paste0("$\\hat{T}^C_{\\rm{Fisher, ", nm, "}}$, $\\Sigma = ",
                     pho_tex(ph), "$, $a = ", a, "$")
            vals <- lapply(chs, function(ch)
              lookup(tb, pp[[2]], ex, a, ph, ch, est))
            blocks <- c(blocks, band_block(lab, chs, vals))
          }
        }
        groups <- c(groups, paste(blocks, collapse = "\n"))
      }
      tail <- paste0(
        "\n\\hline\n\\end{tabular}}\n\\captionsetup{width=\\wd\\TAB}\n\\centering\n",
        "\\caption{Empirical sizes and powers of $\\hat{T}_{\\rm{min, ", nm,
        "}}^C$ and $\\hat{T}_{\\rm{Fisher, ", nm,
        "}}^C$ for different bandwidths and ", nm,
        " under $H_{2", ex, "}$ in Study 2.}\\label{table-bandwidth-H2", ex,
        "m", est, "}\n",
        "\\usebox{\\TAB}\n\\end{table}\n\\clearpage\n")
      chunks <- c(chunks, paste0(head_h2_band(n), "\n",
                                 paste(groups, collapse = "\n\\hline\n"), tail))
    }
  }
  doc_assemble(chunks, thetable = "S\\arabic{table}")
}

## ---------------------------------------------------------------------
## 7.  write the files (and compile the PDFs if LaTeX is available)
## ---------------------------------------------------------------------
## pdflatex from PATH first, then the usual TeX Live / MiKTeX locations
find_latex <- function(cmd = "pdflatex") {
  if (file.exists(cmd)) return(normalizePath(cmd, winslash = "/"))
  p <- Sys.which(cmd)
  if (nzchar(p)) return(unname(p))
  cand <- c(
    Sys.glob("C:/texlive/*/bin/windows/pdflatex.exe"),
    Sys.glob("D:/texlive/*/bin/windows/pdflatex.exe"),
    Sys.glob("D:/LaTeX/texlive/*/bin/windows/pdflatex.exe"),
    Sys.glob("C:/Program Files/MiKTeX/miktex/bin/x64/pdflatex.exe"),
    Sys.glob("C:/Users/*/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdflatex.exe"))
  cand <- cand[file.exists(cand)]
  if (length(cand)) cand[1] else ""
}

latex_engine <- find_latex(latex_cmd)
cat("working directory : ", getwd(), "\n", sep = "")
cat("LaTeX engine      : ",
    if (nzchar(latex_engine)) latex_engine else "not found (only .tex will be written)", "\n",
    sep = "")
if (nzchar(latex_engine) && grepl("[^ -~]", getwd()))
  cat("note: the working directory contains non-ASCII characters;",
      "pdflatex may fail there, run the script from an ASCII path instead.\n")

write_set <- function(name, tex) {
  dir <- out_dir[[name]]
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir, winslash = "/", mustWork = TRUE)
  target <- file.path(dir, "tables.tex")
  con <- file(target, open = "w")
  writeLines(tex, con, sep = "")
  close(con)
  ntab <- lengths(regmatches(tex, gregexpr("\\\\begin\\{table\\}", tex)))
  cat(sprintf("  %-20s %s  (%d tables)\n", name, target, ntab))
  if (build_pdf && nzchar(latex_engine)) {
    old <- setwd(dir)
    ## pdflatex output goes to a *file*: handing it the console of a GUI
    ## (RStudio) or a discarded handle makes pdfTeX die on Windows with
    ## "fflush() failed (Bad file descriptor)".
    outfile <- file.path(dir, "tables.out")
    errfile <- file.path(dir, "tables.err")
    st <- suppressWarnings(system2(latex_engine,
                                   c("-interaction=nonstopmode", "-halt-on-error",
                                     "tables.tex"),
                                   stdout = outfile, stderr = errfile))
    setwd(old)
    if (is.null(st) || is.na(st)) st <- 0L
    log <- file.path(dir, "tables.log")
    if (st == 0) {
      unlink(c(outfile, errfile, file.path(dir, c("tables.aux", "tables.log"))))
      cat("                       pdf ok\n")
    } else {
      cat("                       pdflatex failed (status ", st, "). Last lines:\n", sep = "")
      shown <- character(0)
      for (f in c(outfile, errfile))
        if (file.exists(f) && file.size(f) > 0)
          shown <- c(shown, tail(readLines(f, warn = FALSE), 6))
      if (!length(shown) && file.exists(log))
        shown <- tail(readLines(log, warn = FALSE), 6)
      for (l in shown) cat("                         ", l, "\n", sep = "")
    }
  }
  invisible(target)
}

cat("generate_tables.R\n")
write_set("h1_rpgrp", build_h1_rpgrp(csv_file$h1_rpgrp))
write_set("h2_grp",   build_h2_grp(csv_file$h2_grp))
write_set("h1_bw",    build_h1_band(csv_file$h1_bw))
write_set("h2_bw",    build_h2_band(csv_file$h2_bw))
cat("done.\n")

