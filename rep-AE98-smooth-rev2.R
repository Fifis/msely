rm(list = ls())

# Tested with:
# pnd 0.1.2
# smoothemplik 0.0.17


tryCatch(setwd("~/Dropbox/HSE/10/missing/application/angrist-evans-1998/"), error = function(e) return(NULL))
# cores <- if (Sys.info()["nodename"] == "UL0012159") 24 else 14

library(lmtest)
library(sandwich)

library(parallel)
library(smoothemplik)
library(momentfit)
library(pnd)

library(tikzDevice)

# CHANGE THE NUMBER OF CORES DEPENDING ON THE MACHINE
ncores <- if (.Platform$OS.type == "unix") 4 else 1
# Please use all available cores! Each of the final simulations uses all
# non-parallel functions (SEL in chunks, gradients etc.), so the more cores,
# the better; a decent super-computer cluster should have 100+ cores
RcppParallel::setThreadOptions(numThreads = 1)
data.table::setDTthreads(1)


endog.vars <- "MOREKIDS"
# Cases:
# Year: 1980 or 1990
# Sample: 1. All women; 2. Married women
# Sub-sample: 1. Original model, 2. Whites + nonwhites, 3. Whites + blacks, 4. Whites only
# Formula: 1. All variables, 2. Drop BOY1ST
# Dependent variable (DV): one out of 5

pcomb <- expand.grid(year = c(1980, 1990), sample = 1:2, subsample = 1:4, formula = 1:2,
                     depvar = c("WORKEDM", "WEEKSM1", "HOURSWM", "INCOMEM", "FAMINCL"), stringsAsFactors = FALSE)
pcomb <- pcomb[pcomb$depvar == "INCOMEM", ]
p <- which(pcomb$sample == 2 & pcomb$subsample == 4 & pcomb$year == 1980 & pcomb$formula == 1)
print(pcomb[p, ])

# sink("GMM-overid.txt")

dep.var <- pcomb$depvar[p]
print(pcomb[p, ])

# Selecting by year
d <- if (pcomb$year[p] == 1980) readRDS("pums80-twoa.rds") else readRDS("pums90-onea.rds")

# Selecting by sample
if (pcomb$sample[p] == 2) d <- d[d$MSAMPLE == 1, ]

# Selecting by sub-sample
# We are selecting only white women
if (pcomb$year[p] == 1990) d <- d[d$WHITEM == 1, ] else d <- d[d$BLACKM + d$HISPM + d$OTHRACEM == 0, ]
year <- pcomb$year[p]

# Race in 1990: 001 = White, 002 = Black,
# Addressing different ethnicity encoding in different years
# if (max(d$HISPM) > 1) d$HISPM <- as.integer(as.logical(d$HISPM))
# # Selecting only one mother's nationality for a perfect partitioning
# # White & Hispanic --> Hispanic <- 1, White <- 0
# c(sum(d$WHITEM == 1), sum(d$HISPM == 1), sum(d$WHITEM == 1 & d$HISPM == 1))
# if (any(d$WHITEM == 1 & d$HISPM == 1)) d$WHITEM[d$WHITEM == 1 & d$HISPM == 1] <- 0
# # Black & Hispanic --> Black <- 1, Hispanic <- 0
# c(sum(d$BLACKM == 1), sum(d$HISPM == 1), sum(d$BLACKM == 1 & d$HISPM == 1))
# if (any(d$BLACKM == 1 & d$HISPM == 1)) d$HISPM[d$HISPM == 1 & d$BLACKM == 1] <- 0
# if (is.null(d$OTHRACEM)) d$OTHRACEM <- as.integer(d$WHITEM + d$BLACKM + d$HISPM == 0)
# # if (is.null(d$WHITEM))   d$WHITEM <- d$BLACKM + d$HISPM + d$OTHRACEM == 0
# # if (is.null(d$NWHITEM))  d$NWHITEM <- as.integer(!as.logical(d$WHITEM))
# if (pcomb$subsample[p] == 3) d <- d[d$WHITEM == 1 | d$BLACKM == 1, ]
# if (pcomb$subsample[p] == 4) d <- d[d$WHITEM == 1, ]

# A categorical ethnicity variable for proper formulae
# if (pcomb$subsample[p] == 1) {
#   d$ETHN <- factor(d$WHITEM + 2*d$BLACKM + 3*d$HISPM + 4*d$OTHRACEM,
#                    labels = c("WHITE", "BLACK", "HISP", "OTHER"))
# } else if (pcomb$subsample[p] == 2) {
#   d$ETHN <- factor(d$WHITEM + 2*d$NWHITEM, labels = c("WHITE", "NWHITE"))
# } else if (pcomb$subsample[p] == 3) {
#   d$ETHN <- factor(d$WHITEM + 2*d$BLACKM, labels = c("WHITE", "BLACK"))
# }

# Changing the units of measurement for proper scaling
d$INCOMEM <- d$INCOMEM / 1000

# Properly encoding child gender for model matrices
d$CHLDORD <- 1 * d$BOYS2 + 2 * d$GIRLS2 + 3*(1-d$BOYS2)*(1-d$GIRLS2)*d$BOY1ST + 4*(1-d$BOYS2)*(1-d$GIRLS2)*(1-d$BOY1ST)
d$CHLDORD <- factor(d$CHLDORD, labels = c("BB", "GG", "BG", "GB"))

colnames(d)[colnames(d) == "AGEM1"] <- "AGEM"

excl.inst <- c("BOYS2", "GIRLS2")
incl.inst <- c("AGEM", "AGEFSTM", if (pcomb$formula[p] == 1) "BOY1ST" else NULL)


f1l <- formula(paste0(dep.var, " ~ ", paste0(c(endog.vars, incl.inst), collapse = " + ")))
f1r <- formula(paste0(" ~ ", paste0(c(excl.inst, incl.inst), collapse = " + ")))
f2r <- formula(paste0(" ~  CHLDORD * AGEM + AGEFSTM"))
f3r <- formula(paste0(" ~ CHLDORD * (AGEM + AGEFSTM)"))
f4r <- formula(paste0(" ~ (CHLDORD + AGEM + AGEFSTM)^2"))
f5r <- formula(paste0(" ~ (CHLDORD + AGEM + AGEFSTM)^2 + I(AGEM^2) + I(AGEFSTM^2)"))


mmX <- model.matrix(f5r, data = d)[, -1]
colnames(mmX) <- gsub("[^A-Za-z0-9]+", "_", colnames(mmX))
colnames(mmX) <- gsub("^I_", "", colnames(mmX))
colnames(mmX) <- gsub("_$", "", colnames(mmX))
gc()

# Sets of instruments
i1 <- 1:5  # AE98 baseline
i2 <- c(1:5, 8:10)  # Added AGEM*CHLDORDR
i3 <- c(1:5, 8:13)  # Added AGEFSTM*CHLDORDR
i4 <- c(1:5, 8:14)  # Added AGEM*AGEFSTM
i5 <- 1:14          # Added AGEM^2, AGEFSTM^2

# Basic summary
mm <- function(x) c(Mean = mean(x), quantile(x, 1:3/4))
mm(d$INCOMEM)
mm(d$INCOMEM[d$MOREKIDS == 0])
mm(d$INCOMEM[d$MOREKIDS == 1])
t.test(INCOMEM ~ MOREKIDS, data = d)

#################################
# GMM
#################################

# Checking for perfect collinearity in the variables
caret::findLinearCombos(mmX)
# Nothing, good

# Estimator 2: SEL, no missingness

cat("\n\n=======================================\n")

mf.mod1 <- momentModel(f1l, f1r, data = d, vcov = "MDS")
mf.mod2 <- momentModel(f1l, f2r, data = d, vcov = "MDS")
mf.mod3 <- momentModel(f1l, f3r, data = d, vcov = "MDS")
mf.mod4 <- momentModel(f1l, f4r, data = d, vcov = "MDS")
mf.mod5 <- momentModel(f1l, f5r, data = d, vcov = "MDS")

mf.GMM1 <- gmmFit(mf.mod1, type = "iter", method = "BFGS", control = list(REPORT = 1, trace = 6))
mf.GMM2 <- gmmFit(mf.mod2, type = "iter", method = "BFGS", control = list(REPORT = 1, trace = 6))
mf.GMM3 <- gmmFit(mf.mod3, type = "iter", method = "BFGS", control = list(REPORT = 1, trace = 6))
mf.GMM4 <- gmmFit(mf.mod4, type = "iter", method = "BFGS", control = list(REPORT = 1, trace = 6))
mf.GMM5 <- gmmFit(mf.mod5, type = "iter", method = "BFGS", control = list(REPORT = 1, trace = 6))

#####################
### Table 1, GMM part
print(mf.GMM.sum1 <- summary(mf.GMM1)) 
#####################
print(mf.GMM.sum2 <- summary(mf.GMM2))
print(mf.GMM.sum3 <- summary(mf.GMM3))
print(mf.GMM.sum4 <- summary(mf.GMM4))
print(mf.GMM.sum5 <- summary(mf.GMM5))

# round(cbind(TSLS = coef(mf.IV), TSLS.SE = mf.IV.sum@coef[,2], IGMM = coef(mf.GMM), IGMM.SE = mf.GMM.sum@coef[,2]), 6)
# round(cbind(TSLS = coef(mf.IV), TSLS.SE = mf.IV.sum@coef[,2], IGMM = coef(mf.GMM), IGMM.SE = mf.GMM.sum@coef[,2]), 6)
slist <- list(mf.GMM.sum1, mf.GMM.sum2, mf.GMM.sum3, mf.GMM.sum4, mf.GMM.sum5)

round(do.call(cbind, lapply(slist, function(x) x@coef[, 1])), 3)

rm(pcomb, slist, mmX, i1, i2, i3, i4, i5, f1l, f1r, f2r, f3r, f4r, f5r,
   mf.mod2, mf.GMM2, mf.GMM.sum2, mf.mod3, mf.GMM3, mf.GMM.sum3,
   mf.mod4, mf.GMM4, mf.GMM.sum4, mf.mod5, mf.GMM5, mf.GMM.sum5)

# Conditioning on top-coded age
table(d$AGEFSTM)   # We need to merge that one unit
d$AGEFSTM2 <- d$AGEFSTM
d$AGEFSTM2[d$AGEFSTM2 > 32] <- 32
mf.mod98 <- momentModel(INCOMEM ~ MOREKIDS + AGEM + AGEFSTM + BOY1ST,
                        ~ AGEFSTM + factor(AGEFSTM2)*(BOYS2 + GIRLS2 + AGEM + BOY1ST),
                        data = d, vcov = "MDS")
mf.GMM98 <- gmmFit(mf.mod98, type = "iter", method = "BFGS", control = list(REPORT = 1, trace = 6))
print(mf.GMM.sum98 <- summary(mf.GMM98))
ddd <- model.matrix(~ AGEFSTM + factor(AGEFSTM2)*(BOYS2 + GIRLS2 + AGEM + BOY1ST), data = d)
colnames(ddd)

# For the output
round(as.vector(t(mf.GMM.sum98@coef[c(1, 3:5, 2), 1:2])), 3)


##################
# SEL part
##################

# Why continuous variables (like AGE...) matter!
# The endogenous variable seems to be binary
cont.vars <- c("AGEM", "AGEFSTM")
excl.inst <- c("BOYS2", "GIRLS2")
incl.inst <- c("AGEM", "AGEFSTM", "BOY1ST")

interaction0 <- function(x, type = c("integer", "factor", "character")) {
  type <- type[1]
  chr <- if (NCOL(x) > 1) apply(x, 1, paste0, collapse = "!") else as.character(x)
  if (type == "character") return(chr)
  fctr <- factor(chr)
  if (type == "factor") return(fctr) else return(as.integer(fctr))
}

# Because interaction() is MUCH slower with many columns and distinct values
# Checking if included instruments create the same partitioning
d$catnum.RHS   <- interaction0(d[, c(setdiff(incl.inst, cont.vars), endog.vars)]) # (Z, Xin partition)
d$catnum.dscrt <- interaction0(d[, c(setdiff(incl.inst, cont.vars), excl.inst)]) # (Xex1, Xin partition)
d$catnum.incl  <- interaction0(d[, setdiff(incl.inst, cont.vars)]) # (Xin partition)
d$catnum.YZXW   <- interaction0(d[, c(dep.var, endog.vars, cont.vars, incl.inst, excl.inst)])
d$catnum.cont   <- interaction0(d[, cont.vars])
d$catnum.XdZW   <- interaction0(d[, c(endog.vars, setdiff(incl.inst, cont.vars), excl.inst)])

#
catnum.exog <- interaction0(d[, c(cont.vars, incl.inst, excl.inst)])
length(unique(catnum.exog))

d <- d[, unique(c(cont.vars, dep.var, endog.vars, incl.inst, excl.inst, "catnum.dscrt", "catnum.YZXW", "catnum.cont"))]
# The SEL weights should be based on Xin, Xex, though, but the matrix will
# still be blockwise and sparse because there are many dummy variables

n <- nrow(d)
n == sum(complete.cases(d[, c(dep.var, excl.inst, incl.inst, endog.vars)]))


# Moment function without missingness (full sample)
g.complete <- function(theta, data, dep.var, RHSvars, ...) {
  if (length(theta) == length(RHSvars) + 1) {
    RHSvars <- c("(Intercept)", RHSvars)
    if (!all(sort(RHSvars) == sort(names(theta)))) stop("The initial value vector must have the same names as RHSvars!")
  } else stop("The RHS variables do not appear in the names of theta!")
  dataRHS <- data[, setdiff(names(theta), "(Intercept)")]
  Xtheta <- as.numeric(as.matrix(cbind(1, dataRHS)) %*% theta)
  Dg <- data[, dep.var] - Xtheta
  Dg[is.na(data[, dep.var])] <- 0
  return(Dg)
}


YZXW <- aggregate(d$catnum.YZXW, by = list(d$catnum.YZXW), FUN = length)
d$count.YZXW <- YZXW[d$catnum.YZXW, ncol(YZXW)]
uniq.inds <- !duplicated(d$catnum.YZXW)
dw <- d[uniq.inds, ]

theta0 <- coef(mf.GMM98)
RHSvars <- setdiff(names(theta0), "(Intercept)")


# d$res.this.mod <- g.complete(theta = theta0, data = d, dep.var = dep.vars[var.index], RHSvars = RHSvars)
# aggregate(res.this.mod ~ catnum.dscrt, data = d, FUN = function(x) round(c(mean = mean(x), SD = sd(x), n = length(x), Tmean = mean(x) / sd(x) * sqrt(length(x))), 3))

# We are going to work *exclusively* with the lighter data set
res.this.mod <- g.complete(theta = theta0, data = dw, dep.var = dep.var, RHSvars = RHSvars)
boxplot(res.this.mod ~ dw$catnum.dscrt, frame = FALSE, notch = TRUE, main = "Residuals by discrete category")
boxplot(abs(res.this.mod) ~ dw$catnum.dscrt, frame = FALSE, notch = TRUE, main = "Abs. residuals by discrete category")

tab.cont <- table(dw$AGEM, dw$AGEFSTM)
image(sort(unique(dw$AGEM)), sort(unique(dw$AGEFSTM)), log1p(tab.cont), asp = 1, xlab = "Age", ylab = "Age at 1st birth")

tab.cat.discr <- table(dw$catnum.dscrt)
cat("Gains through blocking: up to ", round(100*(1 - sum((tab.cat.discr / sum(tab.cat.discr))^2)), 2), "% faster\n", sep = "")

# Doing the SEL in groups with smoothing
# Creating only the unique observations based on exogenous variables and their counts
cond.vars <- c(incl.inst, excl.inst)
this.count <- aggregate(d[, cond.vars], by = as.list(d[, cond.vars]), FUN = length)
d.uniq.XW <- this.count[, seq_along(cond.vars)]
d.uniq.XW$count <- this.count[, ncol(this.count)]
count.tab <- table(d.uniq.XW$count)
xt <- as.numeric(names(count.tab))
plot(xt, as.numeric(count.tab), main = "Observation per the tiniest discrete cell", ylim = c(0, max(count.tab)))
round(quantile(d.uniq.XW$count, 1:19/20))

tr <- function(x) sqrt(x); trinv <- function(x) x^2
trc <- tr(d.uniq.XW$count)
dnc <- density(trc, adjust = 0.5, from = 0)
plot(dnc$x, dnc$y, type = "l", yaxt = "n", bty = "n", ylab = "", xlab = "", xaxt = "n"); rug(trc)
axis(1, at = pretty(dnc$x, n = 9), labels = trinv(pretty(dnc$x, n = 9)))


##########################
#### SEL

bwSEL <- 1.2

optMixSELfull <- function(rho, theta0, data, pi.var = NULL, D.var = NULL, helper = NULL, helper.by.var = NULL,
                          by, cont.vars, speedup.YZXW = TRUE, bw.mu = c(2.5, 2.5), degree.mu = 0,
                          bw = 1.2, kernel = "triangular", weight.tolerance = 1e-10,
                          parallel = TRUE, cores = 4, cl = NULL, hessian = TRUE,
                          ...) {
  tic0 <- Sys.time()
  n <- n0 <- nrow(data)
  YZXW <- unique(c(dep.var, endog.vars, cont.vars, incl.inst, excl.inst))
  
  if (speedup.YZXW) {
    data$fac.YZXW <- interaction0(data[, YZXW], type = "factor")
    data$int.YZXW <- as.integer(data$fac.YZXW)
    
    counts <- aggregate(int.YZXW ~ fac.YZXW, data = data, FUN = NROW)
    names(counts)[2] <- "Acount"
    data <- merge(data, counts, by = "fac.YZXW", all.x = TRUE)
    uniq.inds <- !duplicated(data$int.YZXW)
    data <- data[uniq.inds, ]
    cat("Average duplicated full observations:", sprintf("%1.2f", n0 / nrow(data)), "\n")
    n <- nrow(data)
  } else cat("Not de-duplicating full observations\n")
  data$cont.id <- as.integer(interaction(data[, cont.vars], drop = TRUE))
  # Sorting by exog. discr, cont., and all
  if (speedup.YZXW) data <- data[order(data[, by], data$cont.id, data$int.YZXW), ] else
    data <- data[order(data[, by], data$cont.id), ]
  
  # Creating global smoothing weights
  data.cont.uniq <- data[!duplicated(data$cont.id), c(cont.vars, "cont.id")]
  w <- kernelWeights(as.matrix(data.cont.uniq[, cont.vars]), bw = bw, kernel = kernel)
  # Nearest neighbours
  if (FALSE) {
    xmat <- as.matrix(data.cont.uniq[, cont.vars])
    xdist <- as.matrix(dist(xmat))
    xdist <- apply(xdist, 1, rank, ties.method = "first") - 1
    w <- kernelFun(xdist / (0.1*nrow(xmat)), kernel = "triangular")
    w <- w / rowSums(w)
    # print(diag(w))
  }
  cat("Unique continuous variable combinations:", nrow(data.cont.uniq), "\n")
  
  data.split <- split(data, data[, by])
  k <- length(data.split)
  bs <- sapply(data.split, nrow)
  cat("Unique discrete variable combinations (blocks):", k, "\n")
  cat("Block sizes:", bs, "\n")
  
  # Counting unique conditioning variable combinations --
  # continuous within the split by discrete
  cont.list <- lapply(data.split, function(x) {
    x.tally <- if (speedup.YZXW) aggregate(x$Acount, by = list(x$cont.id), FUN = sum) else
      aggregate(x$cont.id, by = list(x$cont.id), FUN = length)
    colnames(x.tally) <- c("cont.id", "tally")
    out <- unique(x[, c(cont.vars, "cont.id")])
    out <- merge(out, x.tally, by = "cont.id") # Adding tally to all original continuous variables
    out
  })
  clen <- sapply(cont.list, nrow)
  cat("Approx. avg. continuous duplicates by category:", round(sum(bs)/sum(clen), 1), "\n")
  
  bw.list <- replicate(k, bw.mu, simplify = FALSE)
  SEL <- function(theta, parallel = FALSE, cores = 1) {
    SELi <- function(i) {
      x <- data.split[[i]]
      this.bw.mu <- bw.list[[i]]
      rho.series <- rho(theta = theta, data = x, dep.var = dep.var, RHSvars = RHSvars, pi.var = pi.var, D.var = D.var, helper = helper, helper.by.var = helper.by.var, bw.mu = this.bw.mu, degree.mu = degree.mu)
      for (attempt in 1:10) {
        if (all(is.finite(rho.series))) break
        this.bw.mu  <- this.bw.mu + 1
        rho.series <- rho(theta = theta, data = x, dep.var = dep.var, RHSvars = RHSvars, pi.var = pi.var, D.var = D.var, helper = helper, helper.by.var = helper.by.var, bw.mu = this.bw.mu, degree.mu = degree.mu)
      }
      
      i.xcu <- cont.list[[i]]$cont.id
      n.xcu <- cont.list[[i]]$tally  # Identical cont. X = identical SEL values
      empliklist <- lapply(i.xcu, function(j) {
        wrep <- w[j, x$cont.id]
        if (speedup.YZXW) wrep <- wrep * x$Acount
        ret <- EL0(z = rho.series, ct = wrep, renormalise = TRUE, weight.tolerance = weight.tolerance)
        ret$wsum <- sum(wrep)
        ret$wnonz <- sum(wrep > weight.tolerance)
        return(ret)
      })
      logelr.vec  <- unlist(lapply(empliklist, "[[", "logelr")) # ELR is identical where w_ij are identical
      lambda.vec  <- unlist(lapply(empliklist, "[[", "lam"))
      dlambda.vec <- unlist(lapply(empliklist, "[[", "f.root"))
      exit.vec    <- unlist(lapply(empliklist, "[[", "exitcode"))
      wsum.vec    <- unlist(lapply(empliklist, "[[", "wsum"))
      wnonz.vec    <- unlist(lapply(empliklist, "[[", "wnonz"))
      logsemplik <- sum(logelr.vec * n.xcu)
      attr(logsemplik, "logelr") <- logelr.vec
      attr(logsemplik, "lambda") <- lambda.vec
      attr(logsemplik, "dlambda") <- dlambda.vec
      attr(logsemplik, "exitcode") <- exit.vec
      attr(logsemplik, "cont.id") <- i.xcu
      attr(logsemplik, "tally") <- n.xcu
      attr(logsemplik, "wsum") <- wsum.vec
      attr(logsemplik, "wnonz") <- wnonz.vec
      attr(logsemplik, "attempt") <- attempt
      logsemplik
    }
    
    SEL.list <- if (inherits(cl, "cluster")) parallel::parLapply(cl, 1:k, SELi)
    else if (parallel & cores > 1) parallel::mclapply(1:k, SELi, mc.cores = cores, mc.preschedule = FALSE)
    else lapply(1:k, SELi)
    SELs <- unlist(SEL.list)
    SELsum <- sum(SELs)
    atrs <- c("logelr", "lambda", "dlambda", "exitcode", "cont.id", "tally", "wsum", "wnonz", "attempt")
    for (a in atrs) attr(SELsum, a) <- lapply(SEL.list, attr, a)
    attr(SELsum, "SEL.by") <- SELs
    return(SELsum)
  }
  
  
  # print(microbenchmark::microbenchmark(SEL(theta0, parallel = FALSE), SEL(theta0, parallel = TRUE, cores = 2), SEL(theta0, parallel = TRUE, cores = 4), times = 10))
  
  BFGS.ctrl <- list(trace = 6, REPORT = 1, maxit = 100, reltol = 1e-12, ndeps = pmax(rep(6e-6, length(theta0)), abs(theta0)*6e-6))
  cores2 <- min(cores, k, max(round(k/2), 8)) # For function evaluation
  f <- function(x, verbose = TRUE) {
    out <- -SEL(x, parallel = parallel, cores = cores2) # Not too many due to overhead
    if (verbose) cat("point: c(", paste0(formatC(x, digits=6), collapse = ", "),
                     "); SEL = ", sprintf("%1.8f", -out),
                     if (parallel & cores2 > 1) paste0(" -- ", cores2, " cores") else " -- 1 core", "\n", sep = "")
    return(out)
  }
  gf <- function(x, verbose = TRUE) {
    out <- -pnd::Grad(FUN = function(y) SEL(y, parallel = FALSE), x = x, cores = ncores, cl = cl, elementwise = FALSE, vectorised = FALSE, multivalued = FALSE)
    if (verbose) cat("grad:  c(", paste0(formatC(out, digits=6), collapse = ", "),
                     "); ||g|| = ", sqrt(sum(out^2)),
                     if (parallel & cores > 1) paste0(" -- ", min(length(x) * 4, cores), " cores") else " -- 1 core", "\n", sep = "")
    return(out)
  }
  
  
  tic1 <- Sys.time()
  td1 <- as.numeric(difftime(tic1, tic0, units = "secs"))
  a <- f(theta0)
  if (max(ntry <- unlist(attr(a, "attempt"))) > 1) { # If the bandwidth is small for some units, change is straight away
    for (j in which(ntry > 1)) bw.list[[j]] <- bw.list[[j]] + ntry[j] - 1
    tic1 <- Sys.time()
    cat("! Had to increase bw.mu by up to ", max(ntry) - 1, " in ", sum(ntry > 1), "/", k, " subsets because of gaps.\n", sep = "")
    a <- f(theta0)
  }
  wnonz <- attr(a, "wnonz")
  num.ops <- sum(sapply(wnonz, sum))
  nonz.tot <- sapply(1:k, function(i) mean(wnonz[[i]] / bs[i]))
  nonz.ave <- weighted.mean(nonz.tot, w = bs)
  cat("Share of weights above tol. in blocks: ", round(nonz.ave, 5), ".\n", sep = "")
  cat("Time savings due to all optimisations, ratio: ", round(n0^2 / num.ops), ".\n", sep = "")
  
  tic2 <- Sys.time()
  td2 <- as.numeric(difftime(tic2, tic1, units = "secs"))
  
  ag <- gf(theta0)
  tic3 <- Sys.time()
  td3 <- as.numeric(difftime(tic3, tic2, units = "secs"))
  cat("Data preparation took ", round(td1, 1), " s; one SEL eval takes ", round(td2, 1), " s, grad takes ", round(td3, 1)," s. Optimising...\n", sep = "")
  
  o.bfgs <- optim(theta0, fn = f, gr = gf, method = "BFGS", control = BFGS.ctrl)
  og.bfgs <- -gf(o.bfgs$par)
  tic4 <- Sys.time()
  td4 <- as.numeric(difftime(tic4, tic3, units = "secs"))
  cat("Optimisation took", round(td4, 1), "s.\n")
  
  if (hessian) {
    cat("Computing the Hessian via optimHess...\n")
    h1 <- optimHess(par = o.bfgs$par, fn = f, gr = gf, control = BFGS.ctrl)
    v1 <- solve(h1)
    o.bfgs$hessian <- -h1
  } else v1 <- NULL
  tic5 <- Sys.time()
  td5 <- as.numeric(difftime(tic5, tic4, units = "secs"))
  cat(" took ", round(td5, 1), " s.\n", sep = "")
  
  return(list(o.bfgs = o.bfgs, og.bfgs = og.bfgs, vcov1 = v1, bw.mu = bw.list, degree.mu = degree.mu,
              times = c(data.prep = td1, SEL.eval = td2, grad.eval = td3, optim = td4, hessian = td5)))
}

res0 <- g.complete(theta = theta0, data = dw, dep.var = dep.var, RHSvars = RHSvars)
boxplot(res0 ~ dw$catnum.dscrt)


est.rho <- optMixSELfull(rho = g.complete, theta0 = theta0, data = d,
                         by = "catnum.dscrt", cont.vars = cont.vars, speedup.YZXW = TRUE,
                         hessian = TRUE, parallel = TRUE, cl = NULL, cores = ncores)

round(est.rho$o.bfgs$par, 3)
1 - pchisq(est.rho$o.bfgs$value, df = length(est.rho$o.bfgs$par))
# 2*(315.43931198 - est.rho$o.bfgs$value)


v1 <- solve(-est.rho$o.bfgs$hessian)


SELb <- est.rho$o.bfgs$par



#####################
### Table 1; columns 1-4 = GMM part, columns 5:8 = SELrho part
round(cbind(mf.GMM.sum1@coef[, c(1, 2, 3, 4)],
      Coef = SELb, SE = sqrt(diag(v1)), t = SELb/sqrt(diag(v1)), p = 2*pnorm(-abs(SELb/sqrt(diag(v1)))))[c(1, 3:5, 2), ], 3)
#####################


os <- sapply(ls(), function(x) object.size(get(x)))
sort(os)
suppressWarnings(rm(dw, dnc, res.this.mod, uniq.inds, res0, subset.df.uniq, vhat, ddml1, ddml2,
   SMD, SMDn, SMD0, JSMD, B, J, res2, JMK, Omega, Omegahat, res.smd, Zhat, ii, AVar))
gc()
os <- sapply(ls(), function(x) object.size(get(x)))
sort(os)

################# MISSINGNESS

d$catnum.XdZW <- as.integer(interaction(d[, c(endog.vars, setdiff(incl.inst, cont.vars), excl.inst)], drop = TRUE, lex.order = TRUE))

# Propensity score based on observables
clamp <- function(x) (x - min(x)) / (max(x) - min(x))
# d$prop.score <- with(d, 0.95 - 0.3 * MOREKIDS - 0.1 * clamp(AGEM) - 0.15 * clamp(AGEFSTM))

# SEL_g and SEL_rho

g.Dg_pi <- function(theta, data, dep.var, RHSvars, pi.var, ...) {
  if (length(theta) == length(RHSvars) + 1) {
    RHSvars <- c("(Intercept)", RHSvars)
    if (!all(sort(RHSvars) == sort(names(theta)))) stop("The initial value vector must have the same names as RHSvars!")
  } else stop("The RHS variables do not appear in the names of theta!")
  dataRHS <- data[, setdiff(names(theta), "(Intercept)")]
  Xtheta <- as.numeric(as.matrix(cbind(1, dataRHS)) %*% theta)
  Dg <- data[, dep.var] - Xtheta
  Dg <- Dg / data[, pi.var]
  Dg[is.na(data[, dep.var])] <- 0
  return(Dg)
}

rho.eff <- function(theta, data, dep.var, RHSvars, pi.var, D.var, bw.mu = c(2.5, 2.5), degree.mu = 0,
                    helper = c("gstar", "Ystar", "Dg", "DY"), helper.by.var
) {
  if (is.null(helper)) helper <- "gstar"
  helper <- helper[1]
  if (length(theta) == length(RHSvars) + 1) {
    RHSvars <- c("(Intercept)", RHSvars)
    if (!all(sort(RHSvars) == sort(names(theta)))) stop("The initial value vector must have the same names as RHSvars!")
  } else stop("The RHS variables do not appear in the names of theta!")
  
  dataRHS <- data[, setdiff(names(theta), "(Intercept)")]
  Xtheta <- as.numeric(as.matrix(cbind(1, dataRHS)) %*% theta)
  Y <- data[, dep.var]
  Dg <- Y - Xtheta
  good.obs <- is.finite(Y)
  Dg[!good.obs] <- 0
  
  if (helper == "Ystar") {
    Ystar.hat <- kernelMixedSmooth(x = data[good.obs, cont.vars], xout = data[, cont.vars], by = data[good.obs, helper.by.var], byout = data[, helper.by.var], y = Y[good.obs], bw = bw.mu, degree = degree.mu, kernel = "triangular")
    mu.hat    <- Ystar.hat - Xtheta
  } else if (helper == "gstar") {
    mu.hat    <- kernelMixedSmooth(x = data[good.obs, cont.vars], xout = data[, cont.vars], by = data[good.obs, helper.by.var], byout = data[, helper.by.var], y = (Y - Xtheta)[good.obs], bw = bw.mu, degree = degree.mu, kernel = "triangular")
  } else if (helper == "DY") {
    DY <- Y
    DY[is.na(Y)] <- 0
    Ystar.hat <- kernelMixedSmooth(x = data[good.obs, cont.vars], xout = data[, cont.vars], by = data[good.obs, helper.by.var], byout = data[, helper.by.var], y = DY[good.obs], bw = bw.mu, kernel = "triangular") / data[, pi.var]
    mu.hat <- Ystar.hat - Xtheta
  } else if (helper == "Dg") {
    mu.hat    <- kernelMixedSmooth(x = data[good.obs, cont.vars], xout = data[, cont.vars], by = data[good.obs, helper.by.var], byout = data[, helper.by.var], y = Dg[good.obs], bw = bw.mu, kernel = "triangular") / data[, pi.var]
  }
  
  rho <- (Dg - data[, D.var]*mu.hat) / data[, pi.var] + mu.hat
  attr(rho, "mu") <- mu.hat
  attr(rho, "bw") <- bw.mu
  return(rho)
}


#################################################
# Main analysis
###############################################


# b <- 33
# s <- seeds[1]
# Run one simulation on a fixed seed for all missingness levels in a sequence
# This speeds up the execution by a factor of 3 since the previous optimum
# in a nice starting value for a slightly perturbed data set
doOneSeed <- function(s) {
  wdir <- "ae98rep/"
  if (!dir.exists(wdir)) dir.create(wdir)
  fn.seed <- paste0(wdir, "res-", year, "-", s, ".RData")
  if (file.exists(fn.seed)) {
    load(fn.seed) # Contains result.list
  } else {
    result.list <- vector("list", length(svec))
  }
  print("--------------------------------")
  print(Sys.time())
  cat("Year: ", year, ", seed: ", s, "\n", sep = "")
  ####################
  
  for (k in seq_along(svec)) {
    out <- result.list[[k]]
    if (is.null(out$ED)) { # Nothing done, not even the first element
      out <- list(ED = NULL, # Preparing the output structure
                  bw.pi = NULL, bw.pi.CV = NULL, bw.mu = NULL, bw.mu.CV = NULL,
                  coef = NULL, se = NULL, SELg = NULL, SELr = NULL)
    } else if (!is.null(out$SELr)) {
      cat("Year: ", year, ", seed: ", s, ", k: ", k, "/", length(svec), " already done, nothing to compute.\n", sep = "")
      next # Everything is ready, no processing needed
    }
    # Otherwise, do the steps below that have not been finished
    
    prop <- as.numeric(with(d, 0.99 - apply(cbind(MOREKIDS, AGEM, AGEFSTM), 2, clamp) %*% c(0.15, 0.1, -0.05) * svec[k]))
    prop[prop < 0.05] <- 0.05
    prop[prop > 0.99] <- 0.99
    d$prop.score <- prop
    cat("k=", k, ", propensity min, mean, max: ", paste0(round(c(min(prop), mean(prop), max(prop)), 2), collapse = " "), "\n", sep = "")
    set.seed(s)
    d$unif <- runif(n)
    d$D <- as.numeric(d$unif < d$prop.score)
    
    out$ED <- mean(d$D)
    
    if (is.null(out$bw.pi)) {
      cat("Cross-validating the pi bandwidth...", sep = "")
      bw.pi.grid <- c(1.1, seq(1.5, 4, 0.5))
      # Avoiding LSCV because it has not been optimised for mixed cases yet
      d.dedup <- smoothemplik:::prepareKernel(x = d[, c(cont.vars, "catnum.XdZW")], y = d$D, bw = 1)
      CV.dedup <- function(b) {
        pihat <- kernelMixedSmooth(x = d.dedup$x[, cont.vars], y = d.dedup$y, by = d.dedup$x[, "catnum.XdZW"], weights = d.dedup$weights, LOO = TRUE, kernel = "triangular", bw = b, no.dedup = TRUE)
        weighted.mean((pihat - d.dedup$y)^2, w = d.dedup$weights)
      }
      CV.pi.vals <- sapply(bw.pi.grid, CV.dedup)
      # plot(bw.pi.grid, CV.pi.vals)
      bw.pi0 <- bw.pi.grid[which.min(CV.pi.vals)]
      umat <- rbind(diag(2), -diag(2))
      cmat <- c(1.0999, 1.0999, -apply(d.dedup$x[, cont.vars], 2, function(x) diff(range(x))))
      gCV <- function(x) pnd::Grad(CV.dedup, x = x, acc.order = 2, h = 2e-6, elementwise = FALSE, vectorised = FALSE, multivalued = FALSE)
      bw.pi.CV <- constrOptim(rep(bw.pi0, 2), f = CV.dedup, grad = gCV, ui = umat, ci = cmat, outer.eps = 1e-3, outer.iterations = 5, method = "BFGS", control = list(maxit = 10, reltol = 1e-3))$par
      cat(" optimum at bw.pi = (", paste0(sprintf("%1.1f", bw.pi.CV), collapse = ", "), ")\n", sep = "")
      # summary(d$pihat)
      rm(d.dedup); gc()
      
      out$bw.pi    <- bw.pi.CV
      out$bw.pi.CV <- list(x = bw.pi.grid, y = CV.pi.vals)
      result.list[[k]] <- out
      save(result.list, file = fn.seed, compress = "xz")
    }
    d$pihat  <- kernelMixedSmooth(x = d[, cont.vars], y = d$D, by = d$catnum.XdZW, kernel = "triangular", bw = out$bw.pi, parallel = FALSE, cores = ncores)
    d$pihat[d$pihat < 0.01] <- 0.01
    
    # png("propens-est.png", 400, 400, type = "cairo", pointsize = 18)
    # par(mar = c(4, 4, 0.1, 0.1))
    # plot(d$prop.score, d$pihat, asp = 1, bty = "n", xlab = "True propensity", ylab = "Estimated propensity", xlim = c(0, 1), ylim = c(0, 1), cex = 0.2)
    # abline(0, 1, lty = 2)
    # abline(h = c(0, 1), v = c(0, 1), lty = 3, col = "#00000088")
    # legend("topleft", c("Lo", "Med", "Hi"), title = "Missingness", bty = "n", pch = 16, col = substr(mycols, 1, 7))
    # dev.off()
    
    # Adding more bandwidths for pi
    pi.bws <- mean(out$bw.pi) / c(1.5, 2, 3)  # Choice of bandwidth: see footnote in the Appendix
    names(pi.bws) <- c(15, 20, 30)
    for (j in seq_along(pi.bws)) {
      br <- names(pi.bws)[j]
      pihat <- kernelMixedSmooth(x = d[, cont.vars], y = d$D, by = d$catnum.XdZW, kernel = "triangular", bw = pi.bws[j])
      pihat[pihat < 0.01] <- 0.01
      d[, paste0("pihat", br)] <- pihat
    }
    rm(pihat, br); gc()
    cat("Added estimators of pi to the data frame.\n")
    
    d1 <- d
    d1[d1$D == 0, dep.var] <- NA
    
    g <- function(theta, data, pi = "pihat") {
      W <- as.matrix(cbind(`(Intercept)` = 1, data[, c(incl.inst, excl.inst)]))
      X <- as.matrix(cbind(`(Intercept)` = 1, data[, c(endog.vars, incl.inst)]))
      res <- data[, dep.var] - as.vector(X %*% theta)
      g <- W * (res / data[, pi])
    }
    gp <- function(theta, data, pi = "pihat") {
      W <- as.matrix(cbind(`(Intercept)` = 1, data[, c(incl.inst, excl.inst)]))
      X <- as.matrix(cbind(`(Intercept)` = 1, data[, c(endog.vars, incl.inst)]))
      X <- X / data[, pi]
      -crossprod(W, X) / nrow(X)
    }
    
    # First, 2SLS, then, GMM
    f1 <- formula(paste0(dep.var, " ~ ", paste0(c(endog.vars, incl.inst), collapse = " + ")))
    f2 <- formula(paste0("~ ", paste0(c(excl.inst, incl.inst), collapse = " + ")))
    specIV.miss <- momentModel(f1, f2, data = d1[d1$D == 1, ], vcov = "MDS")
    modIV.miss  <- tsls(specIV.miss)
    coef.IV <- coef(modIV.miss)
    rm(specIV.miss, modIV.miss, prop); gc()
    cat("Estimated 2SLS.\n")
    
    
    # Preparing to store the results in matrices for easy recall
    # The closer the initial value with minimal differences, the better
    
    pi.bws <- c(pi.bws, NA) # CV last
    brs <- names(pi.bws)
    # IGMM for many bandwidths
    # Overwriting the specs to save memory -- doing other bandwidths
    
    if (!is.null(out$coef)) {
      cat("Seed " , s, ", k=", k, ", IGMM estimates exist, skipping.\n", sep = "")
    } else {
      doGMMi <- function(i) {
        br <- names(pi.bws)[i]
        cat("Doing seed", s, "k", k, "GMM", br, "\n")
        specGMM.miss <- momentModel(g = function(theta, data) g(theta, data, pi = paste0("pihat", br)), x = d1[d1$D == 1, ],
                                    theta0 = coef.IV,
                                    grad = function(theta, data) gp(theta, data, pi = paste0("pihat", br)), vcov = "MDS")
        cat("Created moment model, fitting...\n")
        modGMM.miss  <- gmmFit(specGMM.miss, type = "iter", itertol = 1e-6, method = "BFGS", control = list(REPORT = 10, trace = 1))
        cat("Fitted moment model.\n")
        list(coef = coef(modGMM.miss), se = summary(modGMM.miss)@coef[, 2])
      }
      
      # IGMM.list <- mclapply(seq_along(pi.bws), doGMMi, mc.cores = min(ncores, length(pi.bws)))
      IGMM.list <- lapply(seq_along(pi.bws), doGMMi)
      coefs.GMM <- do.call(cbind, lapply(IGMM.list, "[[", "coef"))
      SEs.GMM <- do.call(cbind, lapply(IGMM.list, "[[", "se"))
      colnames(coefs.GMM) <- colnames(SEs.GMM) <- paste0("IGMM", brs)
      cat("Seed " , s, ", k=", k, ", IGMM done at ", as.character(Sys.time()), "\n", sep = "")
      
      out$coef <- coefs.GMM
      out$se <- SEs.GMM
      result.list[[k]] <- out
      save(result.list, file = fn.seed, compress = "xz")
      rm(doGMMi, IGMM.list, coefs.GMM, SEs.GMM); gc()
    }
    
    
    
    diag.msg <- paste0("Seed s=", s, ", loop iteration k=", k)
    
    extrCoef <- function(x) if (!is.character(x)) x$o.bfgs$par else {a <- rep(NA, length(coef.IV)); attr(a, "code") <- x; a}
    extrSE <- function(x) if (!is.character(x)) sqrt(diag(x$vcov1)) else {a <- rep(NA, length(coef.IV)); attr(a, "code") <- x; a}
    
    # SELg for many bandwidths
    if (!is.null(out$SELg) & !is.character(out$SELg[[1]])) {
      cat("Seed " , s, ", k=", k, ", SELg estimates exist, skipping.\n", sep = "")
    } else {
      doSELgi <- function(i) {
        br <- names(pi.bws)[i]
        cat("Doing seed", s, ", missingness iteration k", k, ", SELg, bandwidth divider (without dot)", br, "\n")
        init.val.g <- if (k == 1) out$coef[, paste0("IGMM", br)] else result.list[[k-1]]$coef[, paste0("SELg", br)]
        
        trySEL <- function(initval, bw = c(1.2, 1.2)) tryCatch(optMixSELfull(theta0 = initval, data = d1,
                                                                             pi.var = paste0("pihat", br),
                                                                             rho = g.Dg_pi, by = "catnum.dscrt", cont.vars = cont.vars, kernel = "triangular", bw = bw,
                                                                             hessian = TRUE, parallel = TRUE, cores = ncores, cl = cl), speedup.YZXW = TRUE,
                                                               error = function(e) return(paste0(diag.msg, " -- error in SELg", br, ": ", as.character(unlist(e)))),
                                                               warning = function(e) return(paste0(diag.msg, " -- warning in SELg", br, ": ", as.character(unlist(e)))))
        est.miss.g <- trySEL(init.val.g)
        spanfail <- NULL
        if (is.character(est.miss.g)) {
          spanfail <- TRUE
          if (grepl("Spanning", est.miss.g[1])) {
            # Gradually increasing the SEL bandwidth
            est.miss.g <- trySEL(init.val.g, bw = c(1.2, 2.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(init.val.g, bw = c(2.2, 1.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(init.val.g, bw = c(2.2, 2.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(init.val.g, bw = c(2.2, 3.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(init.val.g, bw = c(3.2, 2.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(init.val.g, bw = c(3.2, 3.2))
          } else {
            est.miss.g <- trySEL(out$coef[, paste0("IGMM", br)])
          }
        }
        if (is.character(est.miss.g)) {
          if (grepl("Spanning", est.miss.g[1])) {
            est.miss.g <- trySEL(out$coef[, paste0("IGMM", br)], bw = c(1.2, 2.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(out$coef[, paste0("IGMM", br)], bw = c(2.2, 1.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(out$coef[, paste0("IGMM", br)], bw = c(2.2, 2.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(out$coef[, paste0("IGMM", br)], bw = c(2.2, 3.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(out$coef[, paste0("IGMM", br)], bw = c(3.2, 2.2))
            if (is.character(est.miss.g)) est.miss.g <- trySEL(out$coef[, paste0("IGMM", br)], bw = c(3.2, 3.2))
          } else {
            est.miss.g <- trySEL(coef.IV)
          }
        }
        return(list(spec = est.miss.g, coef = extrCoef(est.miss.g), se = extrSE(est.miss.g), spanfail = spanfail))
      }
      
      SELg.list  <- lapply(seq_along(pi.bws), doSELgi)
      coefs.SELg <- do.call(cbind, lapply(SELg.list, "[[", "coef"))
      SEs.SELg   <- do.call(cbind, lapply(SELg.list, "[[", "se"))
      colnames(coefs.SELg) <- colnames(SEs.SELg) <- paste0("SELg", brs)
      cat("Seed " , s, ", k=", k, ", SELg done at ", as.character(Sys.time()), "\n", sep = "")
      
      if (all(colnames(coefs.SELg) %in% colnames(out$coef))) { # Something had been done before but an error occurred
        out$coef[, colnames(coefs.SELg)] <- coefs.SELg
        out$se[, colnames(SEs.SELg)] <- SEs.SELg
      } else {
        out$coef <- cbind(out$coef, coefs.SELg)
        out$se <- cbind(out$se, SEs.SELg)
      }
      out$SELg <- lapply(SELg.list, "[[", "spec")
      result.list[[k]] <- out
      save(result.list, file = fn.seed, compress = "xz")
      rm(SELg.list, coefs.SELg, SEs.SELg); gc()
    }
    
    if (is.null(out$bw.mu)) {
      cat("Cross-validating the mu bandwidth...\n", sep = "")
      g1.VS <- g.complete(theta = out$coef[, ncol(out$coef)], data = d1[d1$D == 1, ], dep.var = dep.var, RHSvars = RHSvars)
      # This trick enables huge speed-ups with as little information loss as possible
      # cor(g1.VS, g1r.VS) # ~ 0.96...
      ngr <- 100
      g1q <- unique(unname(quantile(g1.VS, 1:(ngr-1)/ngr)))
      g1g <- cut(g1.VS, breaks = c(-Inf, g1q, Inf), labels = 1:ngr)
      g1r.VS <- ave(g1.VS, g1g)
      bw.mu.grid <- c(1.1, 2:4)
      d.dedup <- smoothemplik:::prepareKernel(x = d1[d1$D == 1, c(cont.vars, "catnum.XdZW")], y = g1r.VS, kernel = "triangular", bw = 1.5)
      CV.dedup <- function(b) {
        muhat <- kernelMixedSmooth(x = d.dedup$x[, cont.vars], y = d.dedup$y,
                                   by = d.dedup$x[, "catnum.XdZW"], weights = d.dedup$weights,
                                   LOO = TRUE, kernel = "triangular", bw = b, no.dedup = TRUE)
        weighted.mean((muhat - d.dedup$y)^2, w = d.dedup$weights)
      }
      # CV.dedup(4)
      CV.mu.vals <- sapply(bw.mu.grid, CV.dedup)
      if (all(!is.finite(CV.mu.vals))) { # Grid too small, division by 0
        d.temp <- split(as.data.frame(d.dedup$x[, cont.vars]), f = factor(d.dedup$x[, "catnum.XdZW"]))
        maxgaps <- lapply(d.temp, function(d) apply(d, 2, function(x) max(diff(sort(unique(x))))))
        for (b in 4:10 + 0.1) {
          bw.mu.grid <- c(bw.mu.grid, b)
          this.CV <- CV.dedup(b)
          CV.mu.vals <- c(CV.mu.vals, this.CV)
          if (is.finite(this.CV)) break
        }
      }
      # plot(bw.mu.grid, CV.mu.vals)
      bw.mu0 <- bw.mu.grid[which.min(CV.mu.vals)]
      gCV <- function(x) pnd::Grad(CV.dedup, x = x, h = 1e-5, elementwise = FALSE, vectorised = FALSE, multivalued = FALSE)
      umat <- rbind(diag(2), -diag(2))
      cmat <- c(1.0999, 1.0999, -apply(d.dedup$x[, cont.vars], 2, function(x) diff(range(x))))
      bw.mu.CV <- constrOptim(rep(bw.mu0, 2), CV.dedup, grad = gCV, ui = umat, ci = cmat, outer.eps = 1e-2, outer.iterations = 3, method = "BFGS", control = list(maxit = 5, reltol = 1e-2))$par
      cat(" optimum at bw.mu = (", paste0(sprintf("%1.1f", bw.mu.CV), collapse = ", "), ")\n", sep = "")
      rm(d.dedup, g1.VS, g1q, g1g, g1r.VS, CV.dedup, gCV, umat, cmat); gc()
      
      out$bw.mu    <- bw.mu.CV / 2 # Similar under-smoothing
      out$bw.mu.CV <- list(x = bw.mu.grid, y = CV.mu.vals)
      result.list[[k]] <- out
      save(result.list, file = fn.seed, compress = "xz")
      rm(bw.mu0, bw.mu.CV)
    }
    
    
    if (!is.null(out$SELr) & !is.character(out$SELr[[1]])) {
      cat("Seed " , s, ", k=", k, ", SELr estimates exist, skipping.\n", sep = "")
    } else {
      doSELri <- function(i) {
        br <- names(pi.bws)[i]
        cat("Doing seed", s, "k", k, "SELr", br, "\n")
        init.val.r <- if (k == 1) out$coef[, paste0("SELg", br)] else result.list[[k-1]]$coef[, paste0("SELr", br)]
        
        trySEL <- function(initval, bw = c(1.2, 1.2)) tryCatch(optMixSELfull(theta0 = initval, data = d1, rho = rho.eff,
                                                                             cont.vars = cont.vars, speedup.YZXW = TRUE, bw.mu = out$bw.mu, D.var = "D",
                                                                             pi.var = paste0("pihat", br),
                                                                             helper = "gstar", helper.by.var = "catnum.XdZW", by = "catnum.dscrt",
                                                                             kernel = "triangular", bw = bw, hessian = TRUE,
                                                                             parallel = TRUE, cores = ncores, cl = cl),
                                                               error = function(e) return(paste0(diag.msg, " -- error in SELr -- ", paste0(as.character(e), collapse = " - "))))
        est.miss.r <- trySEL(init.val.r)
        spanfail <- NULL
        if (is.character(est.miss.r)) {
          spanfail <- TRUE
          if (grepl("Spanning", est.miss.r[1])) {
            est.miss.r <- trySEL(init.val.r, bw = c(1.2, 2.2))
            if (is.character(est.miss.r)) est.miss.r <- trySEL(init.val.r, bw = c(2.2, 1.2))
            if (is.character(est.miss.r)) est.miss.r <- trySEL(init.val.r, bw = c(2.2, 2.2))
          }
        }
        if (is.character(est.miss.r)) est.miss.r <- trySEL(out$coef[, paste0("IGMM", br)])
        if (is.character(est.miss.r)) est.miss.r <- trySEL(out$coef[, paste0("SELg", br)])
        if (is.character(est.miss.r)) est.miss.r <- trySEL(coef.IV)
        list(spec = est.miss.r, coef = extrCoef(est.miss.r), se = extrSE(est.miss.r), spanfail = spanfail)
      }
      
      # SELr for many bandwidths
      SELr.list  <- lapply(seq_along(pi.bws), doSELri)
      coefs.SELr <- do.call(cbind, lapply(SELr.list, "[[", "coef"))
      SEs.SELr   <- do.call(cbind, lapply(SELr.list, "[[", "se"))
      colnames(coefs.SELr) <- colnames(SEs.SELr) <- paste0("SELr", brs)
      cat("Seed " , s, ", k=", k, ", SELr done at ", as.character(Sys.time()), "\n", sep = "")
      
      if (all(colnames(coefs.SELr) %in% colnames(out$coef))) { # Something had been done before but an error occurred
        out$coef[, colnames(coefs.SELr)] <- coefs.SELr
        out$se[, colnames(SEs.SELr)] <- SEs.SELr
      } else {
        out$coef <- cbind(out$coef, coefs.SELr)
        out$se <- cbind(out$se, SEs.SELr)
      }
      out$SELr <- lapply(SELr.list, "[[", "spec")
      result.list[[k]] <- out
      save(result.list, file = fn.seed, compress = "xz")
      rm(SELr.list, coefs.SELr, SEs.SELr); gc()
    }
    
  }
  # Nothing to return after the main loop; everything is in the files
  return(NULL)
}

doOneSeedSafe <- function(s) tryCatch(doOneSeed(s), error = function(e) return(e))

svec <- seq(0, 12, 0.2) # Missingnesss strength multiplier
seeds <- 1:1000
ncores <- 1
RcppParallel::setThreadOptions(numThreads = 1)
data.table::setDTthreads(1)
cl <- NULL  # Do not parallelise gradients -- parallelise over seeds instead

# RUN THIS TO BENCHMARK THE CODE
doOneSeedSafe(1)
# Suggestion: leave this overnight

# RUN THIS IN THE CLUSTER
# res.list[seeds] <- mclapply(seeds, doOneSeedSafe, mc.cores = ncores, mc.preschedule = FALSE)

# Main loop; very slow -- this is the workhorse
# for (s in seeds) {
#   a <- doOneSeedSafe(s)
#   gc()
#   cat("Seed", s, "done completely, congratulations!\n")
# }


clus <- makeCluster(ncores, outfile = "")
clusterExport(clus, setdiff(ls(), "clus"))
clusterEvalQ(clus, {library(smoothemplik); library(momentfit); library(pnd); library(data.table)})
clusterEvalQ(clus, {RcppParallel::setThreadOptions(numThreads = 1); data.table::setDTthreads(1)})
res.list <- parLapply(clus, seeds, doOneSeedSafe)
stopCluster(clus)

stop("Enough is enough.")

#################################################
# End of estimation; beginning of output formatting
################################################

MC <- 1000
fl <- paste0("ae98rep/res-", year, "-", 1:MC, ".RData")
# fl <- paste0("many-sims-2025-04-23/res-", year, "-", 1:MC, ".RData")  # On the authors' machine
res.list <- vector("list", MC)

for (fil in fl) {
  if (file.exists(fil)) {
    sid <- as.integer(gsub(".+-(\\d+)\\..+", "\\1", fil))
    load(fil, verbose = (sid %% 50 == 0))
    res.list[[sid]] <- result.list
    if (sid %% 50 == 0) cat(sid, "\n")
  }
}

# The following must be null after successful completion
bads <- which(sapply(res.list, is.null))
cat(length(bads), "completely unfinished seeds:", bads, "\n")
if (length(bads) > 0) res.list <- res.list[-bads]
# Keep only the first 1000 seeds for exact reproducibility
if (length(res.list) > 1000) res.list <- res.list[1:1000]

# Success diagnostics
ncoef1 <- t(sapply(res.list, function(x) {
  # if (inherits(x, "error")) return(rep(0, length(svec)))
  return(sapply(x, function(y) if (is.matrix(y$coef)) sum(grepl("IGMM", colnames(y$coef))) else 0))
}))
ncoef2 <- t(sapply(res.list, function(x) {
  # if (inherits(x, "error")) return(rep(0, length(svec)))
  return(sapply(x, function(y) if (is.matrix(y$coef)) sum(grepl("SELg", colnames(y$coef))) else 0))
}))
ncoef3 <- t(sapply(res.list, function(x) {
  # if (inherits(x, "error")) return(rep(0, length(svec)))
  return(sapply(x, function(y) if (is.matrix(y$coef)) sum(grepl("SELr", colnames(y$coef))) else 0))
}))
# Must have length one
table(ncoef1)
table(ncoef2)
table(ncoef3)
image(ncoef1)
image(ncoef2)
image(ncoef3)

length(res.list)  # Must be 1000

# Did the estimation stop somewhere between the stages? Not in the final version, buit better check
IGMMnoSELg <- which(apply(ncoef1 & !ncoef2, 1, any)) # Any simulations with IGMM but not SELg?
SELgnoSELr <- which(apply(ncoef2 & !ncoef3, 1, any)) # Any simulations with SELg but not SELr?
if (length(IGMMnoSELg) > 0) print(res.list[[IGMMnoSELg[1]]])
if (length(SELgnoSELr) > 0) print(res.list[[SELgnoSELr[1]]])

# cat(paste0("Rscript --no-save --no-restore common-1980-", SELgnoSELr, ".R"), sep = "\n")

# Generating a script with all these tasks


# Variance of the sample average of D
EDmat <- t(sapply(res.list, \(x) sapply(x, \(y) if (!is.null(y[["ED"]])) y[["ED"]] else NA)))
plot(density(EDmat[, ncol(EDmat)], bw = "SJ", na.rm = TRUE), bty = "n", yaxt = "n", ylab = "", main = "Sample avg. of D for s=smax in 1000 MC reps")
rug(EDmat[, ncol(EDmat)], lwd = 2, col = "#00000022")
round(apply(EDmat, 2, range, na.rm = TRUE), 3)

ps <- colMeans(EDmat, na.rm = TRUE)
fps <- is.finite(ps) # Were some probabilities a complete failure?
if (any(!fps)) {
  xps <- seq_along(ps)
  ps <- unname(predict(lm(ps ~ xps), newdata = list(xps = xps)))
}
rm(fps)
mr <- 1 - ps # missingness rate

plot(svec, mr, main = "Missingness rate as a function of s", bty = "n", ylim = c(0, 1))
abline(v = seq(0, 12, 0.5), h = seq(0, 1, 0.1), lty = 3, col = "#00000066")
round(rbind(svec, mr), 3)
print(rbind(svec, mr)[, 1:24], 2)

ncoefperc <- cbind(colMeans(ncoef1) / max(ncoef1) * 100,
                   colMeans(ncoef2) / max(ncoef2) * 100,
                   colMeans(ncoef2) / max(ncoef3) * 100)
matplot(mr, ncoefperc, main = "Percentage of successful simulations", ylim = c(0, 100), type = "l", bty = "n", lwd = 2, lty = 1)
abline(h = c(0, 100), lty = 3)
text(mr, 90 - (100-ncoefperc)*3, round(ncoefperc), srt = 90, cex = 0.8)

# Timings for the estimator
tms <- do.call(rbind, lapply(res.list, function(x) do.call(rbind, lapply(x, function(y) if (is.list(y$SELr[[1]])) c(y$SELg[[1]]$times, y$SELr[[1]]$times) else rep(NA, 10))) ))
tms <- cbind(seed = rep(1:length(res.list), each = length(res.list[[1]])),
             k = rep(1:length(res.list[[1]]), length(res.list)),
             tms)
tms <- data.frame(est = rep(c("SELg", "SELr"), each = nrow(tms)), rbind(tms[, c(1, 2, 2+1:5)], tms[, c(1, 2, 7+1:5)]))

tms.qs <- aggregate(I(optim+hessian) ~ est*k, data = tms, FUN = \(x) quantile(x, c(0.25, 0.5, 0.75)))
round(rbind(SELg = apply(tms.qs[[3]][tms.qs$est == "SELg", ], 2, median),
            SELr = apply(tms.qs[[3]][tms.qs$est == "SELr", ], 2, median)))

tms.qs <- aggregate(optim ~ est*k, data = tms, FUN = \(x) quantile(x, c(0.05, 0.5, 0.95)))

matplot(mr, tms.qs$optim[tms.qs$est == "SELg", ], main = "Optimisation time", lty = c(2, 1, 2), type = "l",
        las = 1, xlab = "Missingness", ylab = "Seconds per simulation", ylim = range(0, tms.qs$optim, na.rm = TRUE), bty = "n", col = 1)
matplot(mr, tms.qs$optim[tms.qs$est == "SELr", ], lty = c(2, 1, 2), type = "l", col = 2, add = TRUE)
# Ratio
tms <- do.call(rbind, lapply(res.list, function(x) do.call(rbind, lapply(x, function(y) if (is.list(y$SELr[[1]])) c(y$SELg[[1]]$times, y$SELr[[1]]$times) else rep(NA, 10))) ))
tms <- cbind(seed = rep(1:length(res.list), each = length(res.list[[1]])),
             k = rep(1:length(res.list[[1]]), length(res.list)),
             tms)
tms <- data.frame(cbind(tms[, 1:2], tms[, 7+1:5] / tms[, 2+1:5]))
boxplot(optim ~ k, data = tms, main = "Optimisation time ratio (SELr / SELg)", las = 2, xlab = "Missingness", ylim = range(0, tms$optim, na.rm = TRUE), frame = FALSE, xaxt = "n"); abline(h = 1, lty = 2)
axis(1, at = seq_along(ps), labels = round(mr, 2))
boxplot(SEL.eval ~ k, data = tms, main = "SEL evaluation time ratio (SELr / SELg)", las = 2, xlab = "Missingness", ylim = range(0, tms$SEL.eval, na.rm = TRUE), frame = FALSE, xaxt = "n"); abline(h = 1, lty = 2)
axis(1, at = seq_along(ps), labels = round(mr, 2))
boxplot(grad.eval ~ k, data = tms, main = "Gradient of SEL evaluation time ratio (SELr / SELg)", las = 1, xlab = "", ylim = range(0, tms$grad.eval, na.rm = TRUE), frame = FALSE, xaxt = "n"); abline(h = 1, lty = 2)
axis(1, at = seq_along(ps), labels = round(mr, 2))
# Can SELg be slower?
min(tms$optim, na.rm = TRUE)

vs <- colnames(res.list[[1]][[1]]$coef)
b.all <- s.all <- list()
# Merging all the estimators
for (v in vs) {
  # Coefs
  a <- lapply(res.list, function(x) do.call(rbind, lapply(x, function(y) if (!is.null(y)) tryCatch(if (is.null(y$coef)) rep(NA, 5) else y$coef[, v], error = \(e) rep(NA, 5)) else rep(NA, 5))))
  aa <- simplify2array(a)
  dimnames(aa) <- list(paste0("s", sprintf("%1.2f", ps)), colnames(a[[1]]), NULL)
  b.all[[v]] <- aa
  # SE
  a <- lapply(res.list, function(x) do.call(rbind, lapply(x, function(y) if (!is.null(y)) tryCatch(if (is.null(y$se)) rep(NA, 5) else y$se[, v], error = \(e) rep(NA, 5)) else rep(NA, 5))))
  aa <- simplify2array(a)
  dimnames(aa) <- list(paste0("s", sprintf("%1.2f", ps)), colnames(a[[1]]), NULL)
  s.all[[v]] <- aa
  print(v)
}

b.all <- simplify2array(b.all)
s.all <- simplify2array(s.all)

dim(b.all)
dimnames(b.all)

# Including simulations with missingness rate up to 46%, as in the paper
pscond <- ps > 1 - 0.46 # Including all simulations
# Including all simulations, and less dense if two decimal digits coincide
# pscond <- ps > 0 & (!duplicated(round(ps, 2)))

round(ps, 2)
b.all <- b.all[pscond,,,]
s.all <- s.all[pscond,,,]
ps <- ps[pscond]
mr <- 1 - ps # missingness rate

# Example comparison -- Median abs. difference between point estimates: VS and efficient
round(apply(abs((b.all[,,,"SELr"] - b.all[,,,"SELg"]) / b.all[,,,"SELg"]), c(1, 2), median, na.rm = TRUE), 2)
# Effect of the pi bandwidth on the point estimate
round(apply(abs((b.all[,,,"SELr"] - b.all[,,,"SELr30"]) / b.all[,,,"SELr30"]), c(1, 2), median, na.rm = TRUE), 2)

med2 <- rbind(apply(abs((b.all[19,,,"SELr"] - b.all[19,,,"SELr15"]) / b.all[19,,,"SELr"]), 1, median, na.rm = TRUE),
              apply(abs((b.all[19,,,"SELr"] - b.all[19,,,"SELr20"]) / b.all[19,,,"SELr"]), 1, median, na.rm = TRUE),
              apply(abs((b.all[19,,,"SELr"] - b.all[19,,,"SELr30"]) / b.all[19,,,"SELr"]), 1, median, na.rm = TRUE))
print(med2, 2)
median(med2)

# Efficiency gains
round(apply(s.all[,,,"SELg30"] / s.all[,,,"SELr30"], c(1, 2), median, na.rm = TRUE), 3) # Small bandwidth = good

# Looking at non-convergence for the highest missingness rate
dimnames(b.all) # Dimensions: (1) missingness, (2) coef. names, (3) seeds, (4) estimators
colMeans(is.finite(b.all[dim(b.all)[1], "MOREKIDS", , ]))
which(!is.finite(b.all[dim(b.all)[1], "MOREKIDS", , "SELg20"]))
which(!is.finite(b.all[dim(b.all)[1], "MOREKIDS", , "SELr20"]))

cord <- c(paste0("SELg", c("30", "20", "15", "")), paste0("SELr", c("30", "20", "15", "")))
cat("MOREKIDS: median ^theta_SEL", "\n")
print(apply(b.all[, "MOREKIDS", , cord], c(1, 3), median, na.rm = TRUE), 3)
cat("MOREKIDS: median SE of (^theta_SEL)", "\n")
print(apply(s.all[, "MOREKIDS", , cord], c(1, 3), median, na.rm = TRUE), 3)

typ.se <- apply(s.all[, "MOREKIDS", , c("IGMM30", "SELg30", "SELr30")], c(1, 3), median, na.rm = TRUE)
round(typ.se[, "IGMM30"] / typ.se[, "SELg30"], 2)
round(typ.se[, "SELg30"] / typ.se[, "SELr30"], 2) - 1

print(apply(b.all[, "MOREKIDS", , c("IGMM30", "SELg30", "SELr30")], c(1, 3), median, na.rm = TRUE), 3)

# What are the chosen bandwidths?
bw.pi.opt <- t(apply(simplify2array(lapply(res.list, function(x) sapply(x, \(y) if (!is.null(a <- y[["bw.pi"]])) mean(a) else NA))), c(1, 2), median))
bw.mu.opt <- t(apply(simplify2array(lapply(res.list, function(x) sapply(x, \(y) if (!is.null(a <- y[["bw.mu"]])) mean(a) else NA))), c(1, 2), median))
matplot(bw.pi.opt, ylim = range(bw.pi.opt, na.rm = T) + c(-0.01, 0.01), type = "l", lwd = 2, lty = 1, bty = "n", main = "Cross-validated bw_pi")
matplot(bw.mu.opt, ylim = range(bw.mu.opt, na.rm = T) + c(-0.01, 0.01), type = "l", lwd = 2, lty = 1, bty = "n", main = "Cross-validated bw_mu")

median(bw.pi.opt, na.rm = TRUE)
median(bw.mu.opt*2, na.rm = TRUE) # Division was made inside the algorithm

vs <- c("IGMM", "IGMM15", "IGMM20", "IGMM30", "SELg", "SELg15", "SELg20", "SELg30", "SELr", "SELr15", "SELr20", "SELr30")

alpha <- 0.05
cnames <- dimnames(b.all)[[2]]
# The importance of double robustness in the following plot: SELg drifts, SELrho does not
for (i in which(cnames %in% c("MOREKIDS"))) { # Variables
  mycols <- rainbow(length(ps), end = 0.45, v = 0.8, rev = TRUE)
  pdf(paste0("SEL-est-", cnames[i], ".pdf"), 12, 6, pointsize = 14)
  par(mar = c(2, 4, 2, 0) + 0.1, mfrow = c(3, 4))
  olap <- NULL
  for (v in vs) { # Estimators
    b <- b.all[,i,,v]
    s <- s.all[,i,,v]
    xl <- quantile(c(b.all[,i,,] + 2*s.all[,i,,], b.all[,i,,] - 2*s.all[,i,,]), c(alpha/2, 1-alpha/2), na.rm = TRUE)
    plot(NULL, NULL, xlim = xl, ylim = c(0, max(mr)),
         main = if (grepl("[0-9]", v)) paste0(substr(v, 1, 4), " CV / ", as.numeric(gsub("[^0-9]", "", v)) / 10) else paste0(v, " CV"),
         bty = "n",
         xlab = "", ylab = "Non-missingness rate", las = 1)
    abline(v = pretty(xl, n = 10), lty = 3, col = "#00000044")
    abline(v = 0)
    # abline(h = c(0.5, 1), lty = 2)
    xleft  <- b + qnorm(alpha/2)*s
    xright <- b + qnorm(1-alpha/2)*s
    xdat <- cbind(LQ = t(apply(xleft,  1, quantile, probs = c(alpha/2, 0.5, 1-alpha/2), na.rm = TRUE)),
                  ME = t(apply(b,  1, quantile, probs = c(alpha/2, 0.5, 1-alpha/2), na.rm = TRUE)),
                  RQ = t(apply(xright, 1, quantile, probs = c(alpha/2, 0.5, 1-alpha/2), na.rm = TRUE)),
                  SUCC = apply(xleft, 1, function(x) mean(is.finite(x))))
    points(xdat[, 2], mr-0.005, pch = 1, lty = 2, col = "#AA0000CC")
    arrows(xdat[, 1], mr-0.005, xdat[, 3], mr-0.005, code = 3, angle = 90, length = 0.05, lty = 2, col = "#AA0000CC")
    points(xdat[, 8], mr+0.005, pch = 1, lty = 2, col = "#008800CC")
    arrows(xdat[, 7], mr+0.005, xdat[, 9], mr+0.005, code = 3, angle = 90, length = 0.05, lty = 2, col = "#008800CC")
    points(xdat[, 5], mr, pch = 16, lwd = 2)
    arrows(xdat[, 4], mr, xdat[, 6], mr, code = 3, angle = 90, length = 0.05)
    # Percentage of simulations where the coefficient would appear insignificant
    olap <- cbind(olap, rowMeans(xleft < 0 & xright > 0))
    # Skipping two frames
  }
  colnames(olap) <- vs
  print(dimnames(b.all)[[2]][i])
  print(round(olap, 2))
  print("")
  # legend("topright", v, lty = 1:4, title = "Estimator",  bg = "#FFFFFFEE", box.col = "#FFFFFFEE", ncol = 2)
  # legend("topleft", paste0(c("Low", "Med", "High"), " (", round(100-colMeans(d[, paste0("D", 1:3)])*100), "%)"), col = mycols, pch = 15, title = "Missingness", bg = "#FFFFFFEE", box.col = "#FFFFFFEE")
  dev.off()
}

est.names <- c("IGMM30", "SELg30", "SELr30") # These ones to plot

# Paper version
tikz("application-compare-CI.tex", width = 6, height = 2.5, pointsize = 12, verbose = TRUE)
par(mar = c(4, 4, 2, 0.5) + 0.1, mfrow = c(1, 3))
olap <- NULL
xl <- -2 + c(-3.5, 4)
leftend <- rightend <- leftci <- rightci <- CIlen <- NULL
CIends <- NULL
for (v in est.names) { # Estimators
  b <- b.all[,"MOREKIDS",,v]
  s <- s.all[,"MOREKIDS",,v]
  
  # abline(h = c(0.5, 1), lty = 2)
  xleft  <- b + qnorm(0.025)*s
  xright <- b + qnorm(0.975)*s
  xdat <- cbind(LQ = t(apply(xleft,  1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)),
                ME = t(apply(b,  1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)),
                RQ = t(apply(xright, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)))
  leftend <- cbind(leftend, xdat[, 4])
  rightend <- cbind(rightend, xdat[, 6])
  leftci <- cbind(leftci, xdat[, 2])
  rightci <- cbind(rightci, xdat[, 8])
  
  plot(NULL, NULL, xlim = xl, ylim = c(0, 0.5), main = switch(v, IGMM30 = "$\\hat\\gamma_{\\mathrm{GMM,IPW}}$", SELg30 = "$\\hat\\gamma_{\\mathrm{SEL, IPW}}$", SELr30 = "$\\hat\\gamma$"), bty = "n", xlab = "", ylab = "Missingness", yaxt = "n")
  axis(2, at = seq(0, 0.7, 0.1), labels = paste0(seq(0, 70, 10), "\\%"), las = 1)
  # abline(v = -4:1, lty = 3, col = "#00000044")
  abline(v = 0)
  # points(xdat[, 2], mr, pch = 1, lty = 2, col = "#AA0000CC")
  # arrows(xdat[, 1], mr, xdat[, 3], mr, code = 3, angle = 90, length = 0.05, lty = 2, col = "#AA0000CC")
  # for (i in c(1, 3, 7, 9)) lines(xdat[, i], mr, lty = 3, lwd = 2)
  for (i in c(2, 8)) lines(xdat[, i], mr, lty = 2, lwd = 2)
  points(xdat[, 5], mr, pch = 16, lwd = 2)
  arrows(xdat[, 4], mr, xdat[, 6], mr, code = 3, angle = 90, length = 0.05, lwd = 2)
  
  # CI widths
  CIlen <- cbind(CIlen, xdat[, 8] - xdat[, 2])
  CIends[[v]] <- xdat[, c(2,8)]
}
# legend("topright", v, lty = 1:4, title = "Estimator",  bg = "#FFFFFFEE", box.col = "#FFFFFFEE", ncol = 2)
dev.off()
print(CIlen)

gamma.GMM <- coef(mf.GMM1)["MOREKIDS"]
gamma.SEL <- est.rho$o.bfgs$par["MOREKIDS"]


# Appendix, figure B.1; CHANGE pscond to all TRUE TO OBTAIN THE FULL PICTURE
tikz("application-CI-length.tex", width = 2.7, height = 1.7, pointsize = 12, verbose = TRUE)
par(mar = c(2, 2, 0.5, 0.5) + .1)
matplot(1-ps, CIlen[, -2], xlab = "", ylab = "", type = "l", lty = 2:1, bty = "n", lwd = 2, ylim = range(0, CIlen[, -2]),
        col = "#000000", las = 1, xaxt = "n")
axis(1, at = seq(0, 0.7, 0.1), labels = paste0(seq(0, 70, 10), "\\%"), las = 1)
dev.off()


# Ratio of confidence-interval lengths
tikz("application-CI-ratio.tex", width = 2.7, height = 1.7, pointsize = 12, verbose = TRUE)
par(mar = c(2, 2.5, 0.5, 0.5) + .1)
plot(1-ps, CIlen[,1] / CIlen[,3], ylim = c(0, 1.4), xlab = "", ylab = "", type = "l", lty = 1, bty = "n", lwd = 2,
        col = "#000000", las = 1, xaxt = "n")
axis(1, at = seq(0, 0.7, 0.1), labels = paste0(seq(0, 70, 10), "\\%"), las = 1)
abline(h = 1, lty = 2)
dev.off()


colnames(rightend) <- colnames(leftend) <- colnames(rightci) <- colnames(leftci)  <- est.names
round(rightend - leftend, 2)
round(rightci - leftci, 2)

# Table 2, visual form
olap <- NULL
for (a in c(0.10, 0.05, 0.01)) {
  for (v in est.names) { # Confidence intervals
    b <- b.all[,"MOREKIDS",,v]
    s <- s.all[,"MOREKIDS",,v]
    xleft  <- b + qnorm(a/2)*s
    xright <- b + qnorm(1-a/2)*s
    olap <- cbind(olap, rowMeans(xleft < 0 & xright > 0, na.rm = TRUE))
    if (any(!is.finite(b.all)) | any(!is.finite(s.all))) warning("Failure in some simulations")
  }
}
colnames(olap) <- apply(expand.grid(est.names, c(0.10, 0.05, 0.01)), 1, paste0, collapse = "_")
rownames(olap) <- paste0("s", sprintf("%1.2f", mr))
print(olap, 2)
dev.off()
par(mar = c(7, 4, 0, 0)+.1)
image(1:ncol(olap),  mr, round(t(olap)[, nrow(olap):1], 3), ylab = "Missingness", zlim = c(0, 1), xaxt = "n", xlab = "", yaxt = "n")
axis(1, seq_len(ncol(olap)), colnames(olap), las = 2)
axis(2, mr, rev(round(mr, 2)), las = 1)
xp <- rep(1:12, each = length(ps))
yp <- rep(rev(mr), ncol(olap))
for (i in seq_along(xp)) {
  text(xp[i], yp[i], round(as.vector(olap)[i]*100), font = 2)
}

############## TABLE 2
irows <- unique(c(seq(1, length(mr), 3), length(mr)))
for (i in irows) cat(sprintf("%1.0f", 100*(mr[i])), "\\% & ", paste0(gsub("^0\\.", ".", sprintf("%1.2f", olap[i, ])), collapse = " & "), " \\\\\n", sep = "")
# With thousands
for (i in irows) cat(sprintf("%1.0f", 100*(mr[i])), "\\% & ", paste0(gsub("^0\\.", ".", sprintf("%1.3f", olap[i, ])), collapse = " & "), " \\\\\n", sep = "")

colnames(olap) <- c(paste0("10%", est.names), paste0("5%", est.names), paste0("1%", est.names))
rownames(olap) <- paste0("miss", sprintf("%1.2f", mr))
round(olap, 2)

mycols <- c("#d85148",            "#c39e32",            "#64b84e",            "#706bde", "#000000")
doLeg <- function(i, where = 2, pos = "bottomright") if (i == where) legend(pos, est.names, ncol = 2, lty = 1, title = "Estimator",  bg = "#FFFFFFAA", box.col = "#FFFFFFAA", col = mycols, lwd = 2) else return(NULL)

# Median coef
pdf(paste0("SEL-rho-beta.pdf"), 10, 6, pointsize = 14)
par(mar = c(2, 2, 2, 0) + 0.1, mfrow = c(2, 2))
for (i in 2:5) { # Variables
  betas <- apply(b.all[,i,,est.names], c(1, 3), median, na.rm = TRUE) # Median coef for each estimator
  betasl <- apply(b.all[,i,,est.names], c(1, 3), quantile, prob = 0.025, na.rm = TRUE)
  betasu <- apply(b.all[,i,,est.names], c(1, 3), quantile, prob = 0.975, na.rm = TRUE)
  matplot(mr, betas, xlim = c(0, 0.45), ylim = range(0, betasl, betasu), main = paste0(dimnames(b.all)[[2]][i], " point est."), bty = "n", xlab = "Missingness rate", ylab = "", type = "l", col = mycols, lwd = 2.5, lty = 1)
  matplot(mr, betasl, type = "l", col = mycols, lwd = 1, lty = 2, add = TRUE)
  matplot(mr, betasu, type = "l", col = mycols, lwd = 1, lty = 2, add = TRUE)
  abline(h = 0, lty = 3)
  doLeg(i, where = 3)
}
dev.off()

# Median SE
pdf(paste0("SEL-rho-SE.pdf"), 10, 6, pointsize = 14)
par(mar = c(2, 2, 2, 0) + 0.1, mfrow = c(2, 2))
for (i in 2:5) { # Variables
  sms <- apply(s.all[,i,,est.names], c(1, 3), median, na.rm = TRUE) # Median SE for each estimator
  smsl <- apply(s.all[,i,,est.names], c(1, 3), quantile, prob = 0.05, na.rm = TRUE)
  smsu <- apply(s.all[,i,,est.names], c(1, 3), quantile, prob = 0.95, na.rm = TRUE)
  matplot(mr, sms, xlim = c(0, 0.45), ylim = c(0, quantile(c(sms, smsl, smsu), 0.97)), main = paste0(dimnames(s.all)[[2]][i], " SE"), bty = "n", xlab = "Non-missingness rate", ylab = "", type = "l", col = mycols, lwd = 3, lty = 1)
  matplot(mr, smsl, type = "l", col = mycols, lwd = 1, lty = 2, add = TRUE)
  matplot(mr, smsu, type = "l", col = mycols, lwd = 1, lty = 2, add = TRUE)
  abline(h = 0, lty = 3)
  doLeg(i, pos = "topleft")
}
dev.off()

################# FIGURE 1
library(tikzDevice)
tikz("application-compare-SE.tex", width = 5.5, height = 4, pointsize = 12, verbose = TRUE)
par(mar = c(4, 4, 0.5, 0.5) + 0.1, mfrow = c(2, 2))
for (i in c(3:5, 2)) { # Variables
  sms <- apply(s.all[,i,, est.names], c(1, 3), median, na.rm = TRUE) # Median SE for each estimator
  matplot(mr, sms, xlim = c(0, 0.5), ylim = c(0, max(sms)), main = "",
          bty = "n", xlab = "", ylab = "", type = "l", col = "black",
          lwd = 2, lty = 3:1, las = 1, xaxt = "n")
  # abline(h = 0, lty = 3)
  axis(1, seq(0, 0.5, 0.1), paste0(seq(0, 50, 10), "\\%"))
  mtext(paste0("{\\sl ", tolower(dimnames(s.all)[[2]][i]), "}"), 3, line = -1.5, cex = 5/6)
  mtext("Missingness", 1, line = 2, cex = 5/6)
  if (i == 2) print(round(sms[, 2] / sms[, 3], 2))
}
# legend("bottomright", c(expression(SE(hat(theta)[2*SLS])), expression(SE(hat(theta)[VS])), expression(SE(hat(theta)))), bty = "n", lwd = 2, lty = c(1, 2, 4))
# legend("bottomright", c("$\\mathrm{SE}(\\hat\\theta_{\\mathrm{2SLS}})$", "$\\mathrm{SE}(\\hat\\theta_{\\mathrm{VS}})$", "$\\mathrm{SE}(\\hat\\theta)$"), bty = "n", lwd = 2, lty = c(1, 2, 4))
dev.off()

# Presentation version
cairo_pdf(paste0("SEL-rho-SE-pres.pdf"), 6, 4, pointsize = 12)
par(family = "Fira Sans")
cols <- c("#E66100", "#5D3A9B")
par(mar = c(3, 2.5, 2, 0.5) + 0.1, mfrow = c(2, 2))
for (i in 2:5) { # Variables
  sms <- apply(s.all[,i,,est.names], c(1, 3), median, na.rm = TRUE) # Median SE for each estimator
  this.sms <- sms[, c("IGMM30", "SELr30")]
  matplot(mr, this.sms, xlim = c(0, 0.5), ylim = c(0, max(sms)), main = switch(i, "", "I(more than 2 kids)", "Mother's current age", "Mother's age at first birth", "I(first-born is a boy)"), bty = "n", xlab = "", ylab = "", type = "l", col = cols, lwd = 3, lty = 2:1, las = 1, xaxt = "n", yaxt = "n")
  axis(1, seq(0, 0.4, 0.1), labels = paste0(seq(0, 0.4, 0.1)*100, "%"))
  axis(2, at = pretty(c(0, max(sms)), n =  if (i == 2) 4 else 3), las = 1)
  mtext("Missingness", 1, line = 1.75, cex = 0.7)
  abline(h = switch(i, NULL, seq(0, max(this.sms), 0.2), seq(0, max(this.sms), 0.005), seq(0, max(this.sms), 0.01), seq(0, max(this.sms), 0.02)), lty = 3, col = "#00000044")
  if (i == 5) legend("bottomright", c("GMM", "Our estimator"), col = cols, bg = "#FFFFFFBB", box.col = "#FFFFFFBB", lwd = 2, lty = 2:1)
}
dev.off()

pdf(paste0("SEL-rho-tstat.pdf"), 10, 6, pointsize = 14)
par(mar = c(2, 2, 2, 0) + 0.1, mfrow = c(2, 2))
t.all <- b.all / s.all
for (i in 2:5) { # Variables
  tsts <- apply(t.all[,i,,-5], c(1, 3), median, na.rm = TRUE) # Median t for each estimator
  tstsl <- apply(t.all[,i,,-5], c(1, 3), quantile, prob = 0.05, na.rm = TRUE)
  tstsu <- apply(t.all[,i,,-5], c(1, 3), quantile, prob = 0.95, na.rm = TRUE)
  matplot(mr, tsts, xlim = c(0, 0.45), ylim = (yl <- range(0, tsts, tstsl, tstsu)), main = paste0(dimnames(t.all)[[2]][i], " t-stat"), bty = "n", xlab = "Non-missingness rate", ylab = "", type = "l", col = mycols, lwd = 3, lty = 1)
  matplot(mr, tstsl, type = "l", col = mycols, lwd = 1, lty = 2, add = TRUE)
  matplot(mr, tstsu, type = "l", col = mycols, lwd = 1, lty = 2, add = TRUE)
  abline(h = 0, lty = 3)
  doLeg(i, where = 3)
}
dev.off()


# Finding the seed where it looks the closest
mycols2 <- rainbow(5, end = 0.4, v = 0.7, alpha = 0.8, rev = TRUE)

dev.mat <- t(b.all[, "MOREKIDS", ,"IGMM30"] - b.all[, "MOREKIDS", , "SELr30"])
pdf("TSLS-SELrho.pdf", 6, 4.5)
plot(NULL, NULL, xlim = c(-2, 2), ylim = c(0.5, 1), main = expression("MOREKIDS" ~ (2*SLS - SEL[rho]) ~ "in simulations"), bty = "n", xlab = "Difference between estimators", ylab = "Non-missingness")
for (i in 1:nrow(dev.mat)) lines(dev.mat[i, ], ps, col = "#00000022")
dev1 <- rowSums(abs(dev.mat)^2)

best.inds <- order(dev1) # Best seeds for all ED
for (i in 1:5) lines(dev.mat[best.inds[i], ], ps, lwd = 3, col = mycols2[i])
abline(v = 0, lty = 2, lwd = 2)
legend("topright", paste0("Seed = ", best.inds[1:5]), col = mycols2, lwd = 3, bty = "n")
dev.off()


### TABLE in the Appendix, B.1 -- change pscond to all TRUE and re-run!
dimnames(b.all)
gammas.SEL <- t(b.all[c("s0.75", "s0.61", "s0.45", "s0.30"), "MOREKIDS", , "SELr30"])
gammas.GMM <- t(b.all[c("s0.75", "s0.61", "s0.45", "s0.30"), "MOREKIDS", , "IGMM30"])

# Relative bias
tab.bias <- rbind(GMM = colMeans(gammas.GMM - gamma.GMM),
      SEL = colMeans(gammas.SEL - gamma.SEL))
# MSE
tab.mse <- rbind(GMM = colMeans((gammas.GMM - gamma.GMM)^2),
      SEL = colMeans((gammas.SEL - gamma.SEL)^2))
# Gains - MSE ratio minus 1
tab.gains <- tab.mse[1,] / tab.mse[2,] - 1
round(rbind(tab.bias, tab.mse, tab.gains), 3)

# Hausman test: MOREKIDS
dimnames(b.all)
num   <- t(b.all[, "MOREKIDS", , "IGMM"])   - t(b.all[, "MOREKIDS", , "SELr30"])  # b_Ineff - b_Eff
denom <- sqrt(t(s.all[, "MOREKIDS", , "IGMM"])^2 - t(s.all[, "MOREKIDS", , "SELr30"])^2)  # sqrt(VarIneff - VarEff)
sum(!is.finite(denom))  # One bad case
which(!is.finite(denom), arr.ind = TRUE)
denom[452, 32] <- (denom[452, 31] + denom[452, 33]) / 2
tstat <- num / denom
pval <- 2 * pnorm(-abs(tstat))

brks <- seq(0, 1, 0.01)
h1 <- hist(pval[, "s0.90"], breaks = brks, main = "10% missing")
h2 <- hist(pval[, "s0.75"], breaks = brks, main = "25% missing")
h3 <- hist(pval[, "s0.61"], breaks = brks, main = "40% missing")
h4 <- hist(pval[, "s0.45"], breaks = brks, main = "55% missing")
h5 <- hist(pval[, "s0.30"], breaks = brks, main = "70% missing")

plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 10), main = "Distribution of Hausman test p-values",
     bty = "n", ylab = "p-value density")
lines(h1$mids, h1$counts/10, lwd = 2, col = 1)
lines(h2$mids, h2$counts/10, lwd = 2, col = 2)
lines(h3$mids, h3$counts/10, lwd = 2, col = 3)
lines(h4$mids, h4$counts/10, lwd = 2, col = 4)
lines(h5$mids, h5$counts/10, lwd = 2, col = 5)
legend("topright", c("10%", "25%", "40%", "55%", "70%"), bty = "n", lwd = 2, col = 1:5, title = "Missingness")

ks.test(pval[, "s0.45"], punif)

hist(pval[, "s0.50"], breaks = seq(0, 1, 0.02), freq = FALSE, main = "P-value of Hausman test at 50% missingness",
     xlab = "p-value")
abline(h = 1, lty = 2, lwd = 2)

png("qq-plot-pvalue-hausman-miss10.png", 512, 512, type = "cairo", pointsize = 16)
qqplot(ppoints(nrow(pval)), pval[, "s0.90"], col = "#00000044",  conf.level = 0.95,
       asp = 1, main = "Q-Q plot + 95% confidence band", bty = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Theoretical quantiles of Uniform[0, 1]", ylab = "Quantiles of p-value empirical distribution")
abline(a = 0, b = 1, lwd = 3, col = "#FFFFFFFF")
abline(a = 0, b = 1, lty = 2, lwd = 1.5)
abline(h = seq(0, 1, 0.2), v = seq(0, 1, 0.2), lty = 3, col = "#00000033")
legend("topleft", "Missingness: 10%")
dev.off()

png("qq-plot-pvalue-hausman-miss25.png", 512, 512, type = "cairo", pointsize = 16)
qqplot(ppoints(nrow(pval)), pval[, "s0.75"], col = "#00000044",  conf.level = 0.95,
       asp = 1, main = "Q-Q plot + 95% confidence band", bty = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Theoretical quantiles of Uniform[0, 1]", ylab = "Quantiles of p-value empirical distribution")
abline(a = 0, b = 1, lwd = 3, col = "#FFFFFFFF")
abline(a = 0, b = 1, lty = 2, lwd = 1.5)
abline(h = seq(0, 1, 0.2), v = seq(0, 1, 0.2), lty = 3, col = "#00000033")
legend("topleft", "Missingness: 25%")
dev.off()

png("qq-plot-pvalue-hausman-miss40.png", 512, 512, type = "cairo", pointsize = 16)
qqplot(ppoints(nrow(pval)), pval[, "s0.61"], col = "#00000044",  conf.level = 0.95,
       asp = 1, main = "Q-Q plot + 95% confidence band", bty = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Theoretical quantiles of Uniform[0, 1]", ylab = "Quantiles of p-value empirical distribution")
abline(a = 0, b = 1, lwd = 3, col = "#FFFFFFFF")
abline(a = 0, b = 1, lty = 2, lwd = 1.5)
abline(h = seq(0, 1, 0.2), v = seq(0, 1, 0.2), lty = 3, col = "#00000033")
legend("topleft", "Missingness: 40%")
dev.off()

png("qq-plot-pvalue-hausman-miss50.png", 512, 512, type = "cairo", pointsize = 16)
qqplot(ppoints(nrow(pval)), pval[, "s0.50"], col = "#00000044",  conf.level = 0.95,
       asp = 1, main = "Q-Q plot + 95% confidence band", bty = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Theoretical quantiles of Uniform[0, 1]", ylab = "Quantiles of p-value empirical distribution")
abline(a = 0, b = 1, lwd = 3, col = "#FFFFFFFF")
abline(a = 0, b = 1, lty = 2, lwd = 1.5)
abline(h = seq(0, 1, 0.2), v = seq(0, 1, 0.2), lty = 3, col = "#00000033")
legend("topleft", "Missingness: 50%")
dev.off()

png("qq-plot-pvalue-hausman-miss55.png", 512, 512, type = "cairo", pointsize = 16)
qqplot(ppoints(nrow(pval)), pval[, "s0.56"], col = "#00000044",  conf.level = 0.95,
       asp = 1, main = "Q-Q plot + 95% confidence band", bty = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Theoretical quantiles of Uniform[0, 1]", ylab = "Quantiles of p-value empirical distribution")
abline(a = 0, b = 1, lwd = 3, col = "#FFFFFFFF")
abline(a = 0, b = 1, lty = 2, lwd = 1.5)
abline(h = seq(0, 1, 0.2), v = seq(0, 1, 0.2), lty = 3, col = "#00000033")
legend("topleft", "Missingness: 55%")
dev.off()

png("qq-plot-pvalue-hausman-miss70.png", 512, 512, type = "cairo", pointsize = 16)
qqplot(ppoints(nrow(pval)), pval[, "s0.30"], col = "#00000044",  conf.level = 0.95,
       asp = 1, main = "Q-Q plot + 95% confidence band", bty = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Theoretical quantiles of Uniform[0, 1]", ylab = "Quantiles of p-value empirical distribution")
abline(a = 0, b = 1, lwd = 3, col = "#FFFFFFFF")
abline(a = 0, b = 1, lty = 2, lwd = 1.5)
abline(h = seq(0, 1, 0.2), v = seq(0, 1, 0.2), lty = 3, col = "#00000033")
legend("topleft", "Missingness: 70%")
dev.off()


