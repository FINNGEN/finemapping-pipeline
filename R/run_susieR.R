#!/usr/bin/env Rscript

options(stringsAsFactors = F)

library(argparse)
library(data.table)
library(dplyr)
library(susieR)

load_R <- function(path, snp) {
  if (tools::file_ext(path) == "npz") {
    require(reticulate)
    np <- import("numpy")
    R <- np$load(path)["ld"]
  } else {
    R <- fread(cmd = paste("zcat", path), header = FALSE)
  }

  R <- as.matrix(R)
  if (length(snp) != dim(R)[1]) {
    stop("LD matrix size doesn't match the number of SNPs")
  }

  return(R)
}

susie_bhat_wrapper <- function(beta, se, R, n, L, var_y = 1.0) {
  fitted_bhat <- susie_bhat(
    bhat = beta,
    shat = se,
    R = R,
    n = n,
    var_y = var_y,
    L = L,
    scaled_prior_variance = 0.1,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    standardize = TRUE,
    check_input = FALSE
  )

  variables <-
    summary(fitted_bhat)$vars %>%
      rename(prob = variable_prob) %>%
      arrange(variable) %>%
      mutate(
        mean = susie_get_posterior_mean(fitted_bhat),
        sd = susie_get_posterior_sd(fitted_bhat)
      )
  cs <- summary(fitted_bhat)$cs

  return(list(variables = variables, cs = cs))
}

main <- function(args) {
  df <- fread(args$z)
  n <- args$n_samples
  L <- args$L

  if (!is.null(args$pheno)) {
    pheno <- fread(args$pheno)
    var_y <- var(pheno[, 0])
  } else if (!is.null(args$var_y)) {
    var_y <- args$var_y
  } else {
    var_y <- NULL
  }

  if (!is.null(args$ld)) {
    R <- load_R(args$ld, df$rsid)
    res <- susie_bhat_wrapper(df$beta, df$se, R, n, L, var_y)
    variables <- cbind(df, res$variables[c("mean", "sd", "prob", "cs")])
    cs <- res$cs %>% select(-variable)
    write.table(variables, args$snp, sep = "\t", row.names = F, quote = F)
    write.table(cs, args$cred, sep = "\t", row.names = F, quote = F)
  } else {
    require(Matrix)
    R <- Matirx::.sparseDiagonal(nrow(df))
  }
}

parser <- ArgumentParser()

parser$add_argument("--z", "-z", type = "character", required = TRUE)
parser$add_argument("--ld", type = "character")
parser$add_argument("--out", type = "character")
parser$add_argument("--snp", type = "character")
parser$add_argument("--cred", type = "character")
parser$add_argument("--log", type = "character")
parser$add_argument("--pheno", type = "character")
parser$add_argument("--n-samples", "-n", type = "integer", required = TRUE)
parser$add_argument("--var-y", type = "double", default = 1.0)
parser$add_argument("--L", type = "integer", default = 10)

args <- parser$parse_args()
print(args)

if (is.null(args$out) & any(sapply(list(args$snp, args$cred, args$log), is.null))) {
  stop("Either --out or all of --snp, --cred, and --log should be specified.")
}

if (is.null(args$out)) {
  args$out <- "tmp"
}
if (is.null(args$snp)) {
  args$snp <- paste0(args$out, ".susie.snp")
}
if (is.null(args$cred)) {
  args$cred <- paste0(args$out, ".susie.cred")
}
if (is.null(args$log)) {
  args$log <- paste0(args$out, ".susie.log")
}

sink(file(args$log, open = "w"), type = "message")
print("Analysis started")
main(args)

print("Finished!")
sink()
