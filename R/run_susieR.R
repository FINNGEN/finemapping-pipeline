#!/usr/bin/env Rscript

options(stringsAsFactors = F)

library(argparse)
library(data.table)
library(dplyr)
library(stringr)
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

susie_bhat_wrapper <- function(beta, se, R, n, L, var_y = 1.0, prior_weights = NULL) {
  fitted_bhat <- susie_bhat(
    bhat = beta,
    shat = se,
    R = R,
    n = n,
    var_y = var_y,
    L = L,
    prior_weights = prior_weights,
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

  return(list(
    susie_obj = fitted_bhat,
    variables = variables,
    cs = cs
  ))
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

  if (!is.null(args$prior_weights)) {
    prior_weights <- fread(args$prior_weights)
    df$v = str_c(df$chromosome, df$position, df$allele1, df$allele2, sep = ":")
    prior_weights <- prior_weights$mean_binpred[match(df$v, prior_weights$SNP)]
    if (is.null(prior_weights)) {
      stop("'min_binpred' is missing (--prior-weights).")
    }
    if (any(is.na(prior_weights))) {
      stop("All SNPs should have prior weights (--prior-weights).")
    }
  } else {
    prior_weights <- NULL
  }

  if (!is.null(args$ld)) {
    R <- load_R(args$ld, df$rsid)
  } else {
    warning("No LD file is provided.")
    require(Matrix)
    R <- Matirx::.sparseDiagonal(nrow(df))
  }

  res <- susie_bhat_wrapper(df$beta, df$se, R, n, L, var_y, prior_weights)
  variables <- cbind(df, res$variables[c("mean", "sd", "prob", "cs")])
  cs <- res$cs
  if (!is.null(cs)) {
    cs <- cs %>% mutate(cs_size = length(str_split(variable, ",")[[1]])) %>% select(-variable)
  }
  if (args$write_alpha) {
    alpha <- t(res$susie_obj$alpha)
    colnames(alpha) <- paste0("alpha", seq(ncol(alpha)))
    variables <- cbind(variables, alpha)
  }
  if (args$save_susie_obj) {
    saveRDS(res$susie_obj, file = args$susie_obj)
  }

  write.table(variables, args$snp, sep = "\t", row.names = F, quote = F)
  write.table(cs, args$cred, sep = "\t", row.names = F, quote = F)
}

parser <- ArgumentParser()

parser$add_argument("--z", "-z", type = "character", required = TRUE)
parser$add_argument("--ld", type = "character")
parser$add_argument("--out", type = "character")
parser$add_argument("--snp", type = "character")
parser$add_argument("--cred", type = "character")
parser$add_argument("--log", type = "character")
parser$add_argument("--susie-obj", type = "character")
parser$add_argument("--pheno", type = "character")
parser$add_argument("--n-samples", "-n", type = "integer", required = TRUE)
parser$add_argument("--var-y", type = "double", default = 1.0)
parser$add_argument("--L", type = "integer", default = 10)
parser$add_argument("--prior-weights", type = "character")
parser$add_argument("--write-alpha", action = "store_true")
parser$add_argument("--save-susie-obj", action = "store_true")

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
if (is.null(args$susie_obj)) {
  args$susie_obj <- paste0(args$out, ".susie.rds")
}

sink(file(args$log, open = "w"), type = "message")
print("Analysis started")
main(args)

print("Finished!")
sink()
