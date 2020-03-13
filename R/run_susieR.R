#!/usr/bin/env Rscript

options(stringsAsFactors = F)

library(argparse)
library(data.table)
library(dplyr)
library(stringr)
library(susieR)

load_R <- function(path, snp, dominant=False) {
  if (tools::file_ext(path) == "ld") {
    R <- fread(path, header = FALSE)
  } else if(tools::file_ext(path) %in% c("gz", "bgz")) {
    R <- fread(cmd = paste("zcat", path), header = FALSE)
  } else{
    stop("LD matrix file extension is not supported")
  }

  R <- as.matrix(R)
  if (length(snp) != dim(R)[1]) {
    stop("LD matrix size doesn't match the number of SNPs")
  }
  if (dominant) {
    warning("--dominant is specified. LD matrix is squared.")
    R <- R ** 2
  }

  return(R)
}

compute_yty <- function(beta, se, p, R, n, k) {
  # beta and se should be standarized
  beta_s <- beta * sqrt(2 * p * (1 - p))
  se_s <- se * sqrt(2 * p * (1 - p))

  # Y'Y =  Bj^2 (Xj'Xj) + Var(Bj)(Xj'Xj)(N - k)
  XjtXj <- (n - 1) * diag(R)
  yty <- beta_s ** 2 * XjtXj + se_s ** 2 * XjtXj * (n - k)

  return(median(yty))
}

summarize.susie.cs = function (object, orig_vars,...) {
  if (is.null(object$sets))
    stop("Cannot summarize SuSiE object because credible set information is not available")
  variables = data.frame(cbind(1:length(object$pip), object$pip, -1))
  colnames(variables) = c('variable', 'variable_prob', 'cs')
  rownames(variables) = NULL
  added_vars <- c()
  if (object$null_index > 0) variables = variables[-object$null_index,]
  if (!is.null(object$sets$cs)) {
    cs = data.frame(matrix(NA, length(object$sets$cs), 5))
    colnames(cs) = c('cs', 'cs_log10bf', 'cs_avg_r2', 'cs_min_r2', 'variable')
    for (i in 1:length(object$sets$cs)) {
      if (any(object$sets$cs[[i]] %in% added_vars)) {
        print( paste0("Skipping cs ", i , " as there is an overlap between variants in this cs and previous credible sets"))
        print("Removed cs variants:")
        print(orig_vars[object$sets$cs[[i]],], max=length(object$sets$cs[[i]]))
        next
      } else {
        added_vars <- append(added_vars, object$sets$cs[[i]])
      }
      variables$cs[variables$variable %in% object$sets$cs[[i]]] = object$sets$cs_index[[i]]
      variables[variables$variable %in% object$sets$cs[[i]],"cs_specific_prob"] = object$alpha[object$sets$cs_index[[i]], object$sets$cs[[i]]]

      cs$cs[i] = object$sets$cs_index[[i]]
      cs$cs_log10bf[i] = object$lbf[cs$cs[i]]
      cs$cs_avg_r2[i] = object$sets$purity$mean.abs.corr[i]^2
      cs$cs_min_r2[i] = object$sets$purity$min.abs.corr[i]^2
      cs$variable[i] = paste(object$sets$cs[[i]], collapse=',')
    }
    variables = variables[order(variables$variable_prob, decreasing = T),]
  } else {
    cs = NULL
  }
  return(list(vars=variables, cs=na.omit(cs)))
}

susie_ss_wrapper <- function(df, R, n, L, var_y = 1.0, prior_weights = NULL, min_abs_corr=0.0) {
  beta <- df$beta
  se <- df$se
  fitted_bhat <- susie_suff_stat(
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
    check_input = FALSE,
    min_abs_corr=min_abs_corr
  )
  cs_summary <- summarize.susie.cs(fitted_bhat, df)
  variables <-
    cs_summary$vars %>%
      rename(prob = variable_prob) %>%
      arrange(variable) %>%
      mutate(
        mean = susie_get_posterior_mean(fitted_bhat),
        sd = susie_get_posterior_sd(fitted_bhat)
      )

  cs <- cs_summary$cs
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
    R <- load_R(args$ld, df$rsid, dominant=args$dominant)
  } else {
    warning("No LD file is provided.")
    require(Matrix)
    R <- Matirx::.sparseDiagonal(nrow(df))
  }

  if (args$compute_yty) {
    yty <- compute_yty(df$beta, df$se, df$maf, R, n, args$n_covariates)

    # we scale var_y accordingly, instead of tweaking susie's code
    # https://github.com/stephenslab/susieR/blob/master/R/susie_ss.R#L334
    # susie_ss(XtX = XtX, Xty = Xty, yty = var_y * (n-1), n = n, ...)
    var_y <- yty / (n - 1)
  } else if (!is.null(args$yty)) {
    var_y <- args$yty / (n - 1)
  } else if (!is.null(args$pheno)) {
    pheno <- fread(args$pheno)
    var_y <- var(pheno[, 0])
  } else if (!is.null(args$var_y)) {
    var_y <- args$var_y
  } else {
    var_y <- NULL
  }

  res <- susie_ss_wrapper(df, R, n, L, var_y, prior_weights, min_abs_corr=args$min_cs_corr)
  variables <- cbind(df, res$variables[c("mean", "sd", "prob", "cs")])
  cs <- res$cs
  if (!is.null(cs)) {
    cs <- cs %>% mutate(cs_size = unlist(lapply(str_split(variable, ","), length))) %>% select(-variable)
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
parser$add_argument("--yty", type = "double")
parser$add_argument("--min_cs_corr", default=0.5, type = "double")
parser$add_argument("--compute-yty", action="store_true")
parser$add_argument("--n-covariates", "-k", type = "integer")
parser$add_argument("--prior-weights", type = "character")
parser$add_argument("--write-alpha", action = "store_true")
parser$add_argument("--save-susie-obj", action = "store_true")
parser$add_argument("--dominant", action = "store_true")

args <- parser$parse_args()

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

logfile <- file(args$log, open = "w")
sink(logfile, type = "output", split = TRUE)
sink(logfile, type = "message")
print(args)

if (is.null(args$out) & any(sapply(list(args$snp, args$cred, args$log), is.null))) {
  stop("Either --out or all of --snp, --cred, and --log should be specified.")
}
if ((args$var_y != 1.0 | !is.null(args$yty)) & args$compute_yty) {
  warning("--compute-yty will override the specified --var-y/--yty.")
}
if (args$compute_yty & is.null(args$n_covariates)) {
  stop("--n-covariates is required for --compute-yty.")
}

print("Analysis started")
tryCatch({
  main(args)
}, error = function(e) {
  sink()
  message(as.character(e))
  sink(type="message")
  stop(e)
})
