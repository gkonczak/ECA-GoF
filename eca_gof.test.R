# Empty Cell Area (ECA) Test ----

.eca_stat <- function(x_sample, q_fun, p_fun, 
                      params = NULL, delta=1/(2*length(x_sample)), x_min = NULL, x_max = NULL) {
  call_q <- function(p) {
    if (is.null(params)) {
      q_fun(p)
    } else {
      do.call(q_fun, c(list(p = p), params))
    }
  }
  call_p <- function(q) {
    if (is.null(params)) {
      p_fun(q)
    } else {
      do.call(p_fun, c(list(q = q), params))
    }
  }
  L <- call_q(delta )
  U <- call_q(1 - delta )
  total_length <- U - L
  if (total_length <= 0) {
    return(list(
      statistic = 0,
      merged_starts = numeric(0),
      merged_ends = numeric(0),
      L = L,
      U = U
    ))
  }
  
  if (length(x_sample) == 0) {
    return(list(
      statistic = 1,
      merged_starts = numeric(0),
      merged_ends = numeric(0),
      L = L,
      U = U
    ))
  }
  p_vals <- call_p(x_sample)

  p_starts <- pmax(0, p_vals - delta)
  p_ends   <- pmin(1, p_vals + delta)
  
  starts <- call_q(p_starts)
  ends   <- call_q(p_ends)
  
  valid <- is.finite(starts) & is.finite(ends) & starts < ends
  starts <- starts[valid]
  ends <- ends[valid]
  
  if (length(starts) == 0) {
    return(list(
      statistic = 1,
      merged_starts = numeric(0),
      merged_ends = numeric(0),
      L = L,
      U = U
    ))
  }

  merged <- .merge_intervals(starts, ends)
  clip_starts <- pmax(merged$starts, L)
  clip_ends   <- pmin(merged$ends, U)
  
  valid_intervals <- clip_starts < clip_ends
  covered_length <- sum(clip_ends[valid_intervals] - clip_starts[valid_intervals])
  
  eca <- (total_length - covered_length) / total_length
  eca <- max(0, min(1, eca))
  
  list(
    statistic = eca,
    merged_starts = clip_starts[valid_intervals],
    merged_ends = clip_ends[valid_intervals],
    L = L,
    U = U
  )
}

.merge_intervals <- function(starts, ends) {
  if (length(starts) == 0) {
    return(list(starts = numeric(0), ends = numeric(0)))
  }
  
  ord <- order(starts)
  starts <- starts[ord]
  ends <- ends[ord]
  
  merged_starts <- starts[1]
  merged_ends <- ends[1]
  
  for (i in seq_along(starts)[-1]) {
    last_end <- merged_ends[length(merged_ends)]
    
    if (starts[i] <= last_end) {
      merged_ends[length(merged_ends)] <- max(last_end, ends[i])
    } else {
      merged_starts <- c(merged_starts, starts[i])
      merged_ends <- c(merged_ends, ends[i])
    }
  }
  
  list(starts = merged_starts, ends = merged_ends)
}


eca_gof.test <- function(x, dist = "norm", params = list(mean = 0, sd = 1),
                     delta = 1 / (2*length(x)), alpha = 0.05, n.sim = 1000) {
  
  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 5) stop("Sample size must be at least 5")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (delta <= 0 || delta >= 1) stop("delta must be between 0 and 1")
  if (n.sim < 100) stop("n.sim must be at least 100")
  
  q_fun <- match.fun(paste0("q", dist))
  p_fun <- match.fun(paste0("p", dist))
  d_fun <- match.fun(paste0("d", dist))
  r_fun <- match.fun(paste0("r", dist))
  
  result_obs <- .eca_stat(x, q_fun, p_fun, params = params, delta = delta)
  eca_observed <- result_obs$statistic
  
  eca_simulated <- numeric(n.sim)
  for (i in 1:n.sim) {
    x_sim <- do.call(r_fun, c(list(n = n), params))
    result_sim <- .eca_stat(x_sim, q_fun, p_fun, params = params, delta = delta)
    eca_simulated[i] <- result_sim$statistic
  }
  
  p_value <- mean(eca_simulated >= eca_observed)
  decision <- if (p_value < alpha) "REJECT H0" else "DO NOT REJECT H0"
  params_str <- paste(names(params), "=", unlist(params), collapse = ", ")
  statistic <- eca_observed
  names(statistic) <- "ECA"
  
  result <- list(
    statistic = statistic,
    p.value   = p_value,
    alpha     = alpha,
    decision  = decision,
    method    = paste0("Empty Cell Area (ECA) Test\n  Distribution: ", 
                       dist, "(", params_str, ")",
                       "\n  alpha = ", alpha,
                       " | Decision: ", decision),
    data.name = DNAME
  )
  class(result) <- "htest"
  return(result)
}


# Examples ----


## Expected: DO NOT REJECT H0

# 1. H0: X ~ N(5, 2)
x1 <- rnorm(10, mean = 5, sd = 2)

eca_gof.test(
  x     = x1,
  dist  = "norm",
  params = list(mean = 5, sd = 2),
  alpha = 0.05,
  n.sim = 2000
)

# 2. H0: X ~ Exp(rate = 3)
x2 <- rexp(10, rate = 3)

eca_gof.test(
  x      = x2,
  dist   = "exp",
  params = list(rate = 3),
  alpha  = 0.05,
  n.sim  = 2000
)

# 3. H0: X ~ Unif(0, 10)
x3 <- runif(10, min = 0, max = 10)

eca_gof.test(
  x      = x3,
  dist   = "unif",
  params = list(min = 0, max = 10),
  alpha  = 0.05,
  n.sim  = 2000
)

# 4. H0: X ~ Gamma(shape = 2, rate = 0.5)

x4 <- rgamma(10, shape = 2, rate = 0.5)

eca_gof.test(
  x      = x4,
  dist   = "gamma",
  params = list(shape = 2, rate = 0.5),
  alpha  = 0.05,
  n.sim  = 2000
)

# 5. H0: X ~ LogNormal(meanlog = 1, sdlog = 0.5)
x5 <- rlnorm(10, meanlog = 1, sdlog = 0.5)

eca_gof.test(
  x      = x5,
  dist   = "lnorm",
  params = list(meanlog = 1, sdlog = 0.5),
  alpha  = 0.05,
  n.sim  = 2000
)

## EXPECTED: REJECT H0

# 6. H0: X ~ N(0, 1)
x6 <- rt(10, df = 1)
eca_gof.test(x6,
             dist = "norm", 
             params = list(mean = 0, sd = 1),
             alpha = 0.05, 
             n.sim = 2000
)

# 7. H0: X ~ N(1, 1)
x7 <- rexp(10, rate = 1)
eca_gof.test(x7,
             dist   = "norm",
             params = list(mean = 1, sd = 1),
             alpha  = 0.05,
             n.sim  = 2000
)

# 8. H0: X ~ N(0, 3)
# Bimodal mixture vs. unimodal normal
n <- 10
component <- rbinom(n, 1, 0.5)
x8 <- ifelse(component == 0,
             rnorm(n, mean = -3, sd = 0.5),
             rnorm(n, mean =  3, sd = 0.5))

eca_gof.test(x8,
             dist   = "norm",
             params = list(mean = 0, sd = 3),
             alpha  = 0.05,
             n.sim  = 2000
)



