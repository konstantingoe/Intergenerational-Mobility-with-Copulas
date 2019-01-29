# Preparing R -------------------------------------------------------------

# Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  broom,
  car,
  effects,
  foreign,
  ggeffects,
  haven,
  influence.ME,
  interplot,
  lme4,
  lmtest,
  lubridate,
  lm.beta,
  pastecs,
  xtable,
  data.table,
  plyr,
  plm,
  pglm,
  readstata13,
  readxl,
  rlang,
  sandwich,
  sjlabelled,
  sjmisc,
  sjPlot,
  sjstats,
  stargazer,
  strucchange,
  rio,
  gridExtra,
  cowplot,
  rgl,
  plotly,
  copula,
  VineCopula,
  CDVine,
  distrEx,
  statip,
  grid,
  scatterplot3d,
  dplyr,
  psych,
  Hmisc,
  MASS,
  fitdistrplus,
  reshape2,
  # those beneath always load last
  ggplot2,
  tidyverse)


#options(scipen = 100)
options(digits = 2)
