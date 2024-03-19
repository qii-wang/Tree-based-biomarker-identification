# load the R functions
# one of the packages (bartMachine) relies on the installation of Java
# so you may need to install java first http://www.java.com
source("RealData_lm.R")

# read the data
dat = read.csv("linreg_example.csv")

# x is the explanatory variables
# categorical variables should be converted to factor variables
x = select(dat, -outcome)
x$var17 = factor(x$var17)
x$var18 = factor(x$var18)
x$var21 = factor(x$var21)

# continuous outcome
y = dat$outcome

# 11 different functions for biomarker selection
rf_impurity_lm(x, y)
rf_corrected_impurity_lm(x, y)
rf_permutation_importance_lm(x, y)
cif_lm(x, y)
dart20_lm(x, y)
dart50_lm(x, y)
dart200_lm(x, y)
bart_lm(x, y)
rf_min_depth_lm(x, y)
XGBoost_lm(x, y)
tree_lm(x, y)