# load the R functions
# one of the packages (bartMachine) relies on the installation of Java
# so you may need to install java first http://www.java.com
source("RealData_lm.R")
set.seed(512)

# read the data
dat = read.csv("shiny.ZY.cytokine.anc.csv")

# x is the explanatory variables
# categorical variables should be converted to factor variables
x = select(dat, -Y)
# x$var17 = factor(x$var17)
# x$var18 = factor(x$var18)
# x$var21 = factor(x$var21)

# continuous outcome
y = dat$Y

# 11 different functions for biomarker selection
tmp3 = rf_permutation_importance_lm(x, y, alpha = 0.1)
tmp4 = cif_lm(x, y, alpha = 0.1)
tmp6 = dart50_lm(x, y)
tmp8 = bart_lm(x, y, alpha = 0.1)
tmp9 = rf_min_depth_lm(x, y)
tmp10 = XGBoost_lm(x, y, alpha = 0.1)

rf_res = tmp3$selection_result
crf_res = tmp4$selection_result
dart_res = tmp6$selection_result
bart_res = tmp8$selection_result[, c("xnames", "vip", "pval_local", "selection_local")]
colnames(bart_res)[2:4] = c("importance", "pval", "selected")
rf_min_depth_res = tmp9$selection_result
xgboost_res = tmp10$selection_result_gain

cbind(rf_res[order(rf_res$pval),],
      rf_min_depth_res[order(rf_min_depth_res$variable_depth),],
      crf_res[order(crf_res$pval),],
      xgboost_res[order(xgboost_res$pval),],
      bart_res[order(bart_res$pval),],
      dart_res[order(-dart_res$importance),]) %>% write_csv("ENAR_real_data_results_seed4096.csv")

