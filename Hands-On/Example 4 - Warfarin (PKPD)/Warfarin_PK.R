## nlmixr course PAGE 2019
## Example courtesy of Tomoo Funaki and Nick Holford
## Some changes from Matthew Fidler and Wenping Wang
## Load libraries

library(nlmixr)
library(dplyr)

d <- warfarin %>%
  filter(dvid == "cp")

pk.warfarin <- function() {
  ini({
    tktr <- log(0.0001)
    tka <- log(1)
    tcl <- log(0.1)
    tv <- log(1)
    
    eta.ktr ~ 0.1
    eta.ka ~ 0.1
    eta.cl ~ 0.1
    eta.v ~ 0.1
    prop.err <- 1
    pkadd.err <- 0.00002
  })
  model({
    ktr <- exp(tktr + eta.ktr)
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    
    d/dt(depot) = -ktr * depot
    d/dt(gut) =  ktr * depot -ka * gut
    d/dt(center) =  ka * gut - cl / v * center
    
    cp = center / v
    cp ~ prop(prop.err) + add(pkadd.err)
  })
}

fit <- nlmixr(pk.warfarin, d, "saem")
