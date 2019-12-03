library(nlmixr)
library(ggplot2)
library(xpose.nlmixr)

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
# Theophylline model using linCmt and WT Allometric scaling as Covariate on CL
theo_sd$lnWT70 = log(theo_sd$WT/70)
m4 <- function() {
  ini({
    tka <- .5
    tcl <- -3.2
    tv <- -1
    beta.wt <- 0.75
    eta.ka ~ 1
    eta.cl ~ 2
    eta.v ~ 1
    add.sd <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + beta.wt * lnWT70 + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

fit4 <- nlmixr(m4, theo_sd, est="saem", table=tableControl(cwres=TRUE, npde=TRUE))

print(fit4)

#fit4 <- fit4 %>% addCwres() # In case this was not specified under model fit prior to estimation one could add this here to the results

#GoF by xpose
xpdb <- xpose_data_nlmixr(fit4)

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

library(gridExtra)
p1 <- nlmixr::vpc(fit1,nsim=500, show=list(obs_dv=T),
                  ylab = "Theophylline Concentrations (ng/mL)", xlab = "Time (h)")

p2 <- nlmixr::vpc(fit1,nsim=500, show=list(obs_dv=T),pred_corr=TRUE,
                  ylab = "Prediction-Corrected Theophylline Concentrations (ng/mL)", xlab = "Time (h)")

grid.arrange(p1, p2)
