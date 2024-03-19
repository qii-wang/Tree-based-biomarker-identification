# load the R functions
source("RealData_cox.R")

# read the data
dat = read.csv("surv_example.csv")

# x is the explanatory variables
# categorical variables should be converted to factor variables
x = select(dat, -time, -status)
x$var17 = factor(x$var17)
x$var18 = factor(x$var18)
x$var21 = factor(x$var21)

# event time and status
eventtime = dat$time
status = dat$status

# 10 different functions for biomarker selection
rf_impurity_cox(x, eventtime, status)
rf_corrected_impurity_cox(x, eventtime, status)

# this method may return NA values for variable importance or p-value
# this is due to small sample size and large dimensionality
# I am still thinking about the solution to this issue
rf_permutation_importance_cox(x, eventtime, status)

cif_cox(x, eventtime, status)
dart20_cox(x, eventtime, status)
dart50_cox(x, eventtime, status)
dart200_cox(x, eventtime, status)
rf_min_depth_cox(x, eventtime, status)
XGBoost_cox(x, eventtime, status)
tree_cox(x, eventtime, status)