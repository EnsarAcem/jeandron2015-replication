################################################################################
# Robustness checks for the reproduction of
#   Jeandron A, et al. (2015) Water Supply Interruptions and Suspected Cholera
#   Incidence: A Time-Series Regression in the Democratic Republic of the Congo.
#   PLoS Med 12(10): e1001893.
#
# HOW TO USE
#   Run the (corrected) main script pmed.1001893.s004.r first, in the same R
#   session, so that these objects exist:
#     data, volcen, spl, model.basis, dres0, model0_ac, attrdl()
#   Then source this file. It prints one line per check:
#     RR (95% CI) for 0 m3 versus the reference volume, cumulated over all lags,
#     and (linear models only) the attributable fraction with its empirical CI.
#
# NOTE: run in R 4.5 (November 2026); its output is reported in Tables 2-4 of the
#   replication report. The dlnm warnings "centering through 'cen' now applied at the
#   prediction stage" (29 per group of neighbourhoods) do not affect the results:
#   crosspred() centres at the value passed in 'cen'.
################################################################################

library(splines); library(tsModel); library(dlnm); library(MASS)
options(na.action = "na.exclude")
set.seed(2015)                      # the empirical CIs of attrdl() are simulated

################################################################################
# PART 1 - SENSITIVITY MODELS A AND E OF THE ORIGINAL TABLE 2
################################################################################
show_rr <- function(label, pred, at = "0")
  cat(sprintf("%-58s RR %.2f (%.2f-%.2f)\n", label,
              pred$allRRfit[at], pred$allRRlow[at], pred$allRRhigh[at]))

# --- Model A, internally consistent (own basis, own spline, own residuals)
splA <- bs(data$date2, degree = 3, df = 17)
modelA.basis <- crossbasis(data$vol1000, lag = 12, argvar = list(df = 1, cen = volcen),
                           arglag = list(fun = "poly", degree = 2, int = F))
modelA    <- glm(cholera ~ modelA.basis + splA + dow + rain, data, family = quasipoisson)
dresA     <- residuals(modelA, "deviance")
modelA_ac <- update(modelA, . ~ . + Lag(dresA, 1:2))
show_rr("Model A (17 df), own residuals [expected 5.98 (4.09-8.76)]",
        crosspred(modelA.basis, modelA_ac, at = 0:6))

# --- Model A with the autoregressive terms taken from the MAIN model (dres0):
#     this is the combination that returns the published 6.66 (4.31-10.29)
modelA_ac0 <- update(modelA, . ~ . + Lag(dres0, 1:2))
show_rr("Model A (17 df), main-model residuals [expected 6.66 (4.31-10.29)]",
        crosspred(modelA.basis, modelA_ac0, at = 0:6))

# --- Model E with Model D residuals (as distributed) and with its own residuals
modelD.basis <- crossbasis(data$vol1000, lag = 12, argvar = list(df = 1, cen = volcen),
                           arglag = list(fun = "poly", degree = 3, int = F))
modelE.basis <- crossbasis(data$vol1000, lag = 12, argvar = list(df = 1, cen = volcen),
                           arglag = list(fun = "poly", degree = 3, int = T))
modelD <- glm(cholera ~ modelD.basis + spl + dow + rain, data, family = quasipoisson)
modelE <- glm(cholera ~ modelE.basis + spl + dow + rain, data, family = quasipoisson)
dresD  <- residuals(modelD, "deviance"); dresE <- residuals(modelE, "deviance")
show_rr("Model E, Model D residuals [expected 2.99 (1.65-5.45)]",
        crosspred(modelE.basis, update(modelE, . ~ . + Lag(dresD, 1:2)), at = 0:6))
show_rr("Model E, own residuals [expected 3.06 (1.68-5.57)]",
        crosspred(modelE.basis, update(modelE, . ~ . + Lag(dresE, 1:2)), at = 0:6))

################################################################################
# PART 2 - ROBUSTNESS CHECKS
################################################################################
# One function fits a specification with the authors' two-step procedure
# (model without AR terms -> deviance residuals -> model with their lags 1-2).
#
#   y        outcome series
#   other    admissions from the other area (stratified models) or NULL
#   tm       matrix (spline of date) or factor (month-by-year indicators)
#   lag      maximum lag
#   argvar   exposure-response basis; default = authors' linear term
#   x        exposure series (volume, or a 0/1 interruption indicator)
#   family   "quasipoisson" or "negbin"
#   ar       include lag 1-2 deviance residuals?
#   rain     "lin" (same-day, authors), "none", or "dl" (distributed lag 0-12)
#   at, cen  contrast reported: exposure 'at' versus 'cen'
#   af       also return attributable fraction (linear volume models only)
fit_spec <- function(y = data$cholera, other = NULL, tm = spl, lag = 12,
                     argvar = list(df = 1, cen = volcen), x = data$vol1000,
                     family = "quasipoisson", ar = TRUE, rain = "lin",
                     at = 0, cen = volcen, af = TRUE, leads = FALSE) {
  # NB: crosspred()/attrdl() find coefficients by the NAME of the basis object,
  # so no other term in the model may have a name starting with "volb"
  volb <- crossbasis(x, lag = lag, argvar = argvar,
                     arglag = list(fun = "poly", degree = 2, int = F))
  dow <- data$dow
  rhs <- "volb + tm + dow"
  if (rain == "lin") { rainv <- data$rain; rhs <- paste(rhs, "+ rainv") }
  if (rain == "dl")  {
    rainb <- crossbasis(data$rain, lag = 12, argvar = list(fun = "lin"),
                        arglag = list(fun = "poly", degree = 2))
    rhs <- paste(rhs, "+ rainb")
  }
  if (!is.null(other)) rhs <- paste(rhs, "+ other")
  if (leads) {                                   # negative-control exposure
    Q   <- Lag(data$vol1000, -(0:12))            # columns = volume at t, t+1, ..., t+12
    leadb <- crossbasis(Q, lag = 12, argvar = list(df = 1, cen = volcen),
                        arglag = list(fun = "poly", degree = 2, int = F))
    rhs <- paste(rhs, "+ leadb")
  }
  f0 <- as.formula(paste("y ~", rhs))
  fitter <- function(f) if (family == "negbin") glm.nb(f) else glm(f, family = quasipoisson)
  m <- fitter(f0)
  if (ar) {
    ac <- Lag(residuals(m, "deviance"), 1:2)
    m  <- fitter(as.formula(paste("y ~", rhs, "+ ac")))
  }
  p   <- crosspred(volb, m, at = at, cen = cen)
  out <- c(RR = unname(p$allRRfit[1]), low = unname(p$allRRlow[1]),
           high = unname(p$allRRhigh[1]), n = sum(!is.na(fitted(m))))
  if (af && is.null(argvar$fun) && all(x == data$vol1000, na.rm = TRUE)) {
    a    <- attrdl(x, volb, y, m, type = "af", cen = cen, range = c(0, cen))
    asim <- attrdl(x, volb, y, m, type = "af", cen = cen, range = c(0, cen), sim = TRUE)
    out  <- c(out, AF = 100 * a, AFlow = 100 * unname(quantile(asim, .025)),
              AFhigh = 100 * unname(quantile(asim, .975)))
  }
  if (leads) {
    pl  <- crosspred(leadb, m, at = 0, cen = volcen)
    out <- c(out, leadRR = unname(pl$allRRfit[1]), leadlow = unname(pl$allRRlow[1]),
             leadhigh = unname(pl$allRRhigh[1]))
  }
  out
}

show <- function(label, r) {
  s <- sprintf("%-52s RR %.2f (%.2f-%.2f) n=%d", label, r["RR"], r["low"], r["high"], r["n"])
  if ("AF" %in% names(r))
    s <- paste0(s, sprintf("  AF %.1f%% (%.1f-%.1f)", r["AF"], r["AFlow"], r["AFhigh"]))
  if ("leadRR" %in% names(r))
    s <- paste0(s, sprintf("  | leads RR %.2f (%.2f-%.2f)", r["leadRR"], r["leadlow"], r["leadhigh"]))
  cat(s, "\n")
}

run_all <- function(y, other = NULL, title = "") {
  cat("\n=====", title, "=====\n")
  show("Original specification", fit_spec(y, other))

  # A. flexibility of the control for season and trend ------------------------
  for (k in c(4, 6, 8, 10, 14, 16, 18, 24)) {
    df <- round(k * 64 / 12)                     # 64 months of data
    show(sprintf("A. Spline of date: %d df/year (%d df)", k, df),
         fit_spec(y, other, tm = bs(data$date2, degree = 3, df = df)))
  }
  show("A. Month-by-year indicators instead of spline",
       fit_spec(y, other, tm = factor(format(data$date2, "%Y-%m"))))

  # B. maximum lag -------------------------------------------------------------
  for (L in c(7, 10, 14, 21))
    show(sprintf("B. Maximum lag %d d", L), fit_spec(y, other, lag = L))

  # C. shape of the exposure-response -----------------------------------------
  show("C. Exposure-response: natural spline, 2 df",
       fit_spec(y, other, argvar = list(fun = "ns", df = 2), af = FALSE))
  show("C. Exposure-response: natural spline, 3 df",
       fit_spec(y, other, argvar = list(fun = "ns", df = 3), af = FALSE))
  show("C. Binary exposure: day with 0 m3 vs other days",
       fit_spec(y, other, x = as.numeric(data$vol1000 == 0), argvar = list(fun = "lin"),
                at = 1, cen = 0, af = FALSE))
  show("C. Binary exposure: day with <1,000 m3 vs other days",
       fit_spec(y, other, x = as.numeric(data$vol1000 < 1), argvar = list(fun = "lin"),
                at = 1, cen = 0, af = FALSE))

  # D. distribution and autocorrelation ---------------------------------------
  show("D. Negative binomial instead of quasi-Poisson", fit_spec(y, other, family = "negbin"))
  show("D. No autoregressive terms", fit_spec(y, other, ar = FALSE))

  # E. precipitation -----------------------------------------------------------
  show("E. Precipitation omitted", fit_spec(y, other, rain = "none"))
  show("E. Precipitation as distributed lag (0-12 d)", fit_spec(y, other, rain = "dl"))

  # F. leave one calendar year out --------------------------------------------
  yr <- format(data$date2, "%Y")
  for (v in 2009:2014) {
    yy <- y; yy[yr == v] <- NA
    show(sprintf("F. Excluding %d", v), fit_spec(yy, other))
  }

  # G. negative-control exposure: volume supplied on the 12 FOLLOWING days ----
  show("G. Leads 1-12 d added (negative control)", fit_spec(y, other, leads = TRUE, af = FALSE))

  # H. definition of 'optimal production' -------------------------------------
  for (pr in c(.50, .75, .90)) {
    q <- unname(quantile(data$vol1000, pr, na.rm = TRUE))
    show(sprintf("H. Reference = %dth percentile (%.0f m3)", round(100 * pr), 1000 * q),
         fit_spec(y, other, argvar = list(df = 1, cen = q), cen = q))
  }
}

run_all(data$cholera,                          title = "All town")
run_all(data$chol_highcov, data$chol_lowcov,   title = "Higher tap water consumption")
run_all(data$chol_lowcov,  data$chol_highcov,  title = "Lower tap water consumption")

################################################################################
# PART 3 - FIGURE 6: RR BY FLEXIBILITY OF THE SPLINE OF DATE (check A)
#   Saves figure6.tif in the working directory. Uses fit_spec() from Part 2.
################################################################################
ks <- c(4, 6, 8, 10, 12, 14, 16, 18, 24)          # df per year; 12 = original
get_rr <- function(y, other = NULL)
  t(sapply(ks, function(k)
    fit_spec(y, other, tm = bs(data$date2, degree = 3, df = round(k * 64 / 12)),
             af = FALSE)[c("RR", "low", "high")]))
res_all  <- get_rr(data$cholera)
res_high <- get_rr(data$chol_highcov, data$chol_lowcov)

tiff("figure6.tif", width = 6.5, height = 3, units = "in", pointsize = 8, res = 300)
par(mfrow = c(1, 2), mar = c(4, 4.5, 2, 1))
panels <- list(res_all, res_high)
titles <- c("All town", "Higher tap water consumption")
for (i in 1:2) {
  r <- panels[[i]]
  plot(ks, r[, "RR"], log = "y", ylim = c(0.5, 16), type = "n", xaxt = "n",
       xlab = "df per year in the spline of date",
       ylab = "RR, 0 vs 4,794 m3 (log scale)", main = titles[i])
  axis(1, at = ks)
  abline(h = 1, lty = 2)
  arrows(ks, r[, "low"], ks, r[, "high"], angle = 90, code = 3, length = 0.03)
  points(ks, r[, "RR"], pch = ifelse(ks == 12, 19, 21), bg = "white")
}
dev.off()
print(round(cbind(df_per_year = ks, all = res_all, high = res_high), 2))
