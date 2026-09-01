#!/usr/bin/env Rscript
# ============================================================================
# Generate a synthetic Bacterial Cellulose (BC) production optimization dataset
# ============================================================================
#
# Design
# ------
# - Full Factorial Design, 5 factors x 5 levels each  ->  5^5 = 3125 conditions.
# - Factors (5): pH, Temperature, RPM, Inoculum volume, Sugar concentration.
# - Response: Bacterial Cellulose yield (g/L).
# - Every experimental condition is run in triplicate.
#
# The yield is simulated from an unknown "true" model (main effects + one
# interaction + quadratic curvature) plus normally distributed measurement noise.
# The true model and its optimum are stored in a separate metadata file so the
# teacher knows the ground truth that students must rediscover.

RNG_SEED <- 42
N_REPLICATES <- 3

FACTORS <- c("pH", "Temp_C", "RPM", "Inoculum_pct", "Sugar_gL")

LEVELS <- list(
  pH            = c(4.0, 4.5, 5.0, 5.5, 6.0),
  Temp_C        = c(25.0, 27.5, 30.0, 32.5, 35.0),
  RPM           = c(50.0, 87.5, 125.0, 162.5, 200.0),
  Inoculum_pct  = c(5.0, 7.5, 10.0, 12.5, 15.0),
  Sugar_gL      = c(20.0, 30.0, 40.0, 50.0, 60.0)
)

LOW  <- sapply(LEVELS, function(x) x[1])
HIGH <- sapply(LEVELS, function(x) x[length(x)])

TRUE_MODEL <- list(
  b0 = 7.5,
  b = c(pH = 0.5, Temp_C = -0.6, RPM = -0.3, Inoculum_pct = 0.5, Sugar_gL = 0.7),
  interactions = c("pH*Temp_C" = -0.45, "Sugar_gL*RPM" = 0.3),
  quadratic = c(pH = -0.7, Temp_C = -0.8, RPM = -0.5, Inoculum_pct = -0.6, Sugar_gL = -0.7),
  sigma = 0.28
)

DATA_DIR <- file.path(dirname(dirname(normalizePath(sys.frame(0)$ofile))), "data")

# ---- Helper functions --------------------------------------------------------

coded_to_real <- function(code, factor) {
  LOW[factor] + (code + 1.0) * (HIGH[factor] - LOW[factor]) / 2.0
}

real_to_coded <- function(real, factor) {
  2.0 * (real - LOW[factor]) / (HIGH[factor] - LOW[factor]) - 1.0
}

true_response <- function(coded_df) {
  m <- TRUE_MODEL
  y <- rep(m$b0, nrow(coded_df))
  for (f in FACTORS) {
    y <- y + m$b[f] * coded_df[[f]]
  }
  y <- y + m$interactions[["pH*Temp_C"]] * coded_df[["pH"]] * coded_df[["Temp_C"]]
  y <- y + m$interactions[["Sugar_gL*RPM"]] * coded_df[["Sugar_gL"]] * coded_df[["RPM"]]
  for (f in FACTORS) {
    y <- y + m$quadratic[f] * coded_df[[f]]^2
  }
  return(y)
}

true_optimum <- function() {
  m <- TRUE_MODEL
  idx <- setNames(seq_along(FACTORS), FACTORS)
  A <- matrix(0, 5, 5)
  c_vec <- numeric(5)
  for (f in FACTORS) {
    A[idx[f], idx[f]] <- 2.0 * m$quadratic[f]
    c_vec[idx[f]] <- -m$b[f]
  }
  for (pair_name in names(m$interactions)) {
    parts <- strsplit(pair_name, "\\*")[[1]]
    f1 <- parts[1]; f2 <- parts[2]
    A[idx[f1], idx[f2]] <- m$interactions[pair_name]
    A[idx[f2], idx[f1]] <- m$interactions[pair_name]
  }
  x <- solve(A, c_vec)
  coded_opt <- setNames(x, FACTORS)
  coded_df <- as.data.frame(as.list(coded_opt))
  y_opt <- true_response(coded_df)
  return(list(coded = coded_opt, y = y_opt))
}

nearest_condition <- function(coded_opt) {
  best <- list()
  for (f in FACTORS) {
    real_val <- coded_to_real(coded_opt[f], f)
    best[[f]] <- LEVELS[[f]][which.min(abs(LEVELS[[f]] - real_val))]
  }
  return(best)
}

# ---- Main --------------------------------------------------------------------

main <- function() {
  set.seed(RNG_SEED)

  # Build full factorial design
  design <- expand.grid(LEVELS)

  # Convert to coded units
  coded <- as.data.frame(design)
  for (f in FACTORS) {
    coded[[f]] <- real_to_coded(design[[f]], f)
  }

  # True response at each condition
  base_y <- true_response(coded)

  # Generate rows with noise
  rows <- vector("list", nrow(design) * N_REPLICATES)
  idx <- 0
  noise <- rnorm(nrow(design) * N_REPLICATES, 0, TRUE_MODEL$sigma)
  noise_idx <- 0
  for (cond_id in seq_len(nrow(design))) {
    for (rep in 1:N_REPLICATES) {
      noise_idx <- noise_idx + 1
      idx <- idx + 1
      rows[[idx]] <- data.frame(
        run_number    = idx,
        condition_id  = cond_id,
        replicate     = rep,
        pH            = design[cond_id, "pH"],
        Temp_C        = design[cond_id, "Temp_C"],
        RPM           = design[cond_id, "RPM"],
        Inoculum_pct  = design[cond_id, "Inoculum_pct"],
        Sugar_gL      = design[cond_id, "Sugar_gL"],
        BC_yield_gL   = round(base_y[cond_id] + noise[noise_idx], 3)
      )
    }
  }

  df <- do.call(rbind, rows)

  # Write CSV
  dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
  csv_path <- file.path(DATA_DIR, "bc_full_factorial_dataset.csv")
  write.csv(df, csv_path, row.names = FALSE)

  # Write hidden truth JSON
  opt <- true_optimum()
  meta <- list(
    design = "5^5 Full Factorial (5 factors x 5 levels, evenly spaced)",
    factors = lapply(FACTORS, function(f) list(levels = LEVELS[[f]])),
    conditions = nrow(design),
    replicates = N_REPLICATES,
    rows = nrow(df),
    true_model = TRUE_MODEL,
    true_optimum_coded = as.list(round(opt$coded, 4)),
    true_optimum_real  = as.list(setNames(
      round(sapply(FACTORS, function(f) coded_to_real(opt$coded[f], f)), 3),
      FACTORS
    )),
    true_optimum_yield_gL = round(opt$y, 3),
    nearest_design_condition = nearest_condition(opt$coded)
  )

  json_path <- file.path(DATA_DIR, "hidden_truth.json")
  writeLines(toJSON(meta, auto_unbox = TRUE, pretty = TRUE), json_path)

  cat(sprintf("Wrote %s  (%d rows)\n", csv_path, nrow(df)))
  cat(sprintf("Conditions: %d (5^5 full factorial), x%d replicates\n", nrow(design), N_REPLICATES))
  cat(sprintf("Yield range: %.2f - %.2f g/L\n", min(df$BC_yield_gL), max(df$BC_yield_gL)))
  cond_means <- tapply(df$BC_yield_gL, df$condition_id, mean)
  cat(sprintf("Per-condition means: %.2f - %.2f g/L\n", min(cond_means), max(cond_means)))
  cat(sprintf("True optimum (coded): %s\n", paste(names(opt$coded), round(opt$coded, 4), sep="=", collapse=", ")))
  cat(sprintf("True optimum (real):  %s\n",
    paste(FACTORS, sprintf("%.3f", sapply(FACTORS, function(f) coded_to_real(opt$coded[f], f))), sep="=", collapse=", ")))
  cat(sprintf("  -> %.3f g/L\n", opt$y))
}

# Run
main()
