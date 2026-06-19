## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
inla_available <- bdlnm::check_inla()

## ----setup, message=FALSE-----------------------------------------------------
library(bdlnm)
library(dlnm)
library(splines)
library(utils)

## -----------------------------------------------------------------------------
head(london)

## -----------------------------------------------------------------------------
# Exposure-response and lag-response spline parameters
dlnm_var <- list(
  var_prc = c(10, 75, 90),
  var_fun = "ns",
  lag_fun = "ns",
  max_lag = 21,
  lagnk = 3
)


# Cross-basis parameters
argvar <- list(
  fun = dlnm_var$var_fun,
  knots = quantile(london$tmean, dlnm_var$var_prc / 100, na.rm = TRUE),
  Bound = range(london$tmean, na.rm = TRUE)
)

arglag <- list(
  fun = dlnm_var$lag_fun,
  knots = logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk)
)

# Create crossbasis
cb <- crossbasis(london$tmean, lag = dlnm_var$max_lag, argvar, arglag)

## -----------------------------------------------------------------------------
# Seasonality of mortality time series
seas <- ns(london$date, df = round(8 * length(london$date) / 365.25))

## -----------------------------------------------------------------------------
# Prediction values (equidistant points)
temp <- round(seq(min(london$tmean), max(london$tmean), by = 0.1), 1)
# Ensure it falls inside the range of temperatures after rounding:
temp <- temp[temp >= min(london$tmean) & temp <= max(london$tmean)]

## ----eval = inla_available----------------------------------------------------
mod <- bdlnm(
  mort_75plus ~ cb + factor(dow) + seas,
  data = london,
  family = "poisson",
  sample.arg = list(n = 1000, seed = 5243)
)

## ----eval = inla_available----------------------------------------------------
str(mod, max.level = 1)

## ----eval = inla_available----------------------------------------------------
mmt <- optimal_exposure(mod, exp_at = temp)

str(mmt)

## ----fig.width = 7, eval = inla_available-------------------------------------
plot(
  mmt,
  xlab = "Temperature (ºC)",
  main = paste0("MMT (Median = ", round(mmt$summary[["0.5quant"]], 1), "ºC)")
)

## ----eval = inla_available----------------------------------------------------
cen <- mmt$summary[["0.5quant"]]
cen

## ----eval = inla_available----------------------------------------------------
cpred <- bcrosspred(mod, exp_at = temp, cen = cen)

## ----eval = inla_available----------------------------------------------------
str(cpred)

## ----eval = inla_available----------------------------------------------------
cpred$coefficients |>
  head(c(5, 5))

## ----eval = inla_available----------------------------------------------------
cpred$matRRfit[,, "sample1"] |>
  head()

## ----eval = inla_available----------------------------------------------------
cpred$allRRfit |>
  head(c(5, 5))

## ----eval = inla_available----------------------------------------------------
cpred$matRRfit.summary |>
  head(5)

cpred$allRRfit.summary |>
  head(5)

## ----fig.width = 7, fig.height = 4, eval = inla_available---------------------
plot(
  cpred,
  "overall",
  xlab = "Temperature (ºC)",
  ylab = "Relative Risk",
  col = 4,
  main = "Overall",
  log = "y"
)

## ----fig.width = 7, fig.height = 4, eval = inla_available---------------------
plot(
  cpred,
  "overall",
  xlab = "Temperature (ºC)",
  ylab = "Relative Risk",
  col = 4,
  main = "Overall",
  log = "y",
  ci = "sampling"
)

## ----fig.width = 7, fig.height = 6, eval = inla_available---------------------
plot(
  cpred,
  "3d",
  zlab = "Relative risk",
  col = 4,
  lphi = 60,
  cex.axis = 0.6,
  xlab = "Temperature (ºC)",
  main = "3D graph of temperature effect"
)

## ----fig.width = 7, fig.height = 6, eval = inla_available---------------------
plot(
  cpred,
  "contour",
  xlab = "Temperature (ºC)",
  ylab = "Lag",
  main = "Contour plot"
)

## ----fig.width = 7, fig.height = 6, eval = inla_available---------------------
htemp <- round(quantile(london$tmean, 0.99), 1)
plot(
  cpred,
  "slices",
  ci = "bars",
  type = "p",
  pch = 19,
  exp_at = htemp,
  ylab = "RR",
  main = paste0("Association for a high temperature (", htemp, "ºC)")
)

## ----fig.width = 7, fig.height = 6, eval = inla_available---------------------
plot(
  cpred,
  "slices",
  lag_at = 0,
  col = 4,
  ylab = "RR",
  main = paste0("Association at lag 0")
)

