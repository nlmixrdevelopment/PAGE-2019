## ------------------------------------------------------------------------
## Load libraries
library(ggplot2)
library(nlmixr)

str(theo_sd)

ggplot(theo_sd, aes(TIME, DV)) + geom_line(aes(group=ID), col="red") + scale_x_continuous("Time (h)") + scale_y_continuous("Concentration") + labs(title="Theophylline single-dose", subtitle="Concentration vs. time by individual")



## ------------------------------------------------------------------------
one.cmt <- function() {
    ini({
        tka <- .45 # log Ka
        tcl <- 1 # log Cl
        tv <- 3.45 # log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.err <- 0.7
    })
    model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.err)
    })
}


## ------------------------------------------------------------------------
fit <- nlmixr(one.cmt, theo_sd, est="nlme")
print(fit)


## ------------------------------------------------------------------------
one.compartment <- function() {
    ini({
        tka <- .45 # log Ka
        tcl <- 1 # log Cl
        tv <- 3.45 # log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.err <- 0.7
    })
    model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.err)
    })
}


## ------------------------------------------------------------------------
fit <- nlmixr(one.compartment, theo_sd, est="saem", control=saemControl(print=500))

fit


## ------------------------------------------------------------------------
## Example 1 -- Fix tka to 0.5 and re-estimate.
one.ka.0.5 <- fit %>%
    update(tka=fix(0.5)) %>%
    nlmixr(est="saem", control=saemControl(print=500))

one.ka.0.5


## ------------------------------------------------------------------------
## Example 2 -- Fix tka to model estimated value and re-estimate.
one.ka.0.5 <- fit %>%
    update(tka=fix) %>%
    nlmixr(est="saem", control=saemControl(print=500))

one.ka.0.5


## ------------------------------------------------------------------------
## Example 3 -- Change tka to 0.7 in orginal model function and then estimate
one.ka.0.7 <- fit %>%
    update(tka=0.7) %>%
    nlmixr(est="saem", control=saemControl(print=500))

one.ka.0.7


## ------------------------------------------------------------------------
## Example 4 -- Remove eta.ka on ka
noEta <- fit %>%
    update(ka <- exp(tka)) %>%
    nlmixr(est="saem", control=saemControl(print=500))

noEta


## ------------------------------------------------------------------------
## Example 5 -- Add eta.ka on ka
addBackKa <- noEta %>%
    update({ka <- exp(tka + bsv.ka)}) %>%
    update(bsv.ka=0.1) %>%
    nlmixr(est="saem", control=saemControl(print=500))

addBackKa


## ------------------------------------------------------------------------
addBackKa$omega


## ------------------------------------------------------------------------
## Example 6 -- Add covariate
## Note currently cov is needed as a prefix so nlmixr knows this is an
## estimated parameter not a parameter
wt70 <- fit %>% ## FIXME could be based on if it finds the covarite in the last used nlmixr data.
    update({cl <- exp(tcl + eta.cl)*(WT/70)^covWtPow}) %>%
    update(covWtPow=fix(0.75)) %>%
    nlmixr(est="saem", control=saemControl(print=500))

wt70


## ------------------------------------------------------------------------
## Example 7 -- Change residual error
## Since there are 0 observations in the data, these are changed to 0.0150 to show proportional error change.
d <- theo_sd
d$DV[d$EVID == 0 & d$DV == 0] <- 0.0150

addPropModel <- fit %>%
    update({cp ~ add(add.err)+prop(prop.err)}) %>%
    update(prop.err=0.1) %>%
    nlmixr(d,est="saem", control=saemControl(print=500))

addPropModel


## ------------------------------------------------------------------------
## Add Diagnostic Plots
library(xpose.nlmixr)
addCwres(fit)
xpdb <- xpose_data_nlmixr(fit)
dv_vs_pred(xpdb)
dv_vs_ipred(xpdb)
res_vs_pred(xpdb)
res_vs_idv(xpdb)
prm_vs_iteration(xpdb)
absval_res_vs_idv(xpdb, res = 'IWRES')
absval_res_vs_pred(xpdb, res = 'IWRES')
ind_plots(xpdb, nrow=3, ncol=4)
res_distrib(xpdb)
vpc.ui(fit, show=list(obs_dv=T))

