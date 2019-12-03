library(nlmixr)
library(ggplot2)
library(xpose.nlmixr)
#-----------------------------------------------------------
#Exploratory plots
df = theo_sd

ggplot(df, aes(x=TIME, y=DV, color=factor(ID))) +
  geom_point() +
    geom_line() +
    theme(legend.position = "none")

ggplot(df, aes(x=TIME, y=DV, color = factor(ID))) +
  geom_point() + geom_line() +
  stat_smooth(aes(color = NULL), size = 1.3) +
      scale_y_log10() +
    theme(legend.position = "none")

#----------------------------------------------------------
# Theophylline model using ODEs
m2 <- function() {
  ini({
    tka <- .45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    ## After the differential equations are defined
    kel <- cl / v;
    d/dt(depot) = -ka*depot;
    d/dt(centr) = ka*depot-kel*centr;
    ## And the concentration is then calculated
    cp = centr / v;
    ## Last, nlmixr is told that the plasma concentration follows
    ## a additive error (estimated by the parameter add.err)
    cp ~ add(add.sd)
  })
}

fit2 <- nlmixr(m2, theo_sd, est="saem", table=tableControl(cwres=TRUE, npde=TRUE))

print(fit2)

#fit2 <- fit2 %>% addCwres() # In case this was not specified under model fit prior to estimation one could add this here to the results

#GoF by xpose
xpdb <- xpose_data_nlmixr(fit2)

dv_vs_pred(xpdb) +
  ylab("Observed Theophylline Concentrations (ng/mL)") +
  xlab("Population Predicted Theophylline Concentrations (ng/mL)")
dv_vs_ipred(xpdb) +
  ylab("Observed Theophylline Concentrations (ug/mL)") +
  xlab("Individual Predicted Theophylline Concentrations (ng/mL)")
res_vs_pred(xpdb) +
  ylab("Conditional Weighted Residuals") +
  xlab("Population Predicted Theophylline Concentrations (ng/mL)")
res_vs_idv(xpdb) +
  ylab("Conditional Weighted Residuals") +
  xlab("Time (h)")
prm_vs_iteration(xpdb)
absval_res_vs_idv(xpdb, res = 'IWRES') +
  ylab("Individual Weighted Residuals") +
  xlab("Time (h)")
absval_res_vs_pred(xpdb, res = 'IWRES')  +
  ylab("Individual Weighted Residuals") +
  xlab("Population Predicted Theophylline Concentrations (ng/mL)")
ind_plots(xpdb, nrow=3, ncol=4) +
  ylab("Predicted and Observed Theophylline concentrations (ng/mL)") +
  xlab("Time (h)")
res_distrib(xpdb) +
  ylab("Density") +
  xlab("Conditional Weighted Residuals")

p1 <- nlmixr::vpc(fit1,nsim=500, show=list(obs_dv=T),
            ylab = "Theophylline Concentrations (ng/mL)", xlab = "Time (h)")

p2 <- nlmixr::vpc(fit1,nsim=500, show=list(obs_dv=T), pred_corr=TRUE,
                  ylab = "Prediction-Corrected Theophylline Concentrations (ng/mL)", xlab = "Time (h)")

library(gridExtra)

grid.arrange(p1,p2)
