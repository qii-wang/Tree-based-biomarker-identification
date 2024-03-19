if(!require(ranger)){
  install.packages("ranger")
}

if(!require(party)){
  install.packages("party")
}

if(!require(BART)){
  install.packages("BART")
}

if(!require(randomForestSRC)){
  install.packages("randomForestSRC")
}

if(!require(tidyverse)){
  install.packages("tidyverse")
}

if(!require(xgboost)){
  install.packages("xgboost")
}

if(!require(rpart)){
  install.packages("rpart")
}

if(!require(bartMachine)){
  install.packages("bartMachine")
}

library(ranger)
library(party)
library(BART)
library(randomForestSRC)
library(tidyverse)
library(xgboost)
library(rpart)
library(bartMachine)

rf_impurity_lm = function(x, y, num_permutations = 100, alpha = 0.05){
  df = cbind(y, x)
  p = ncol(x)
  t_start = proc.time()
  rf_model = ranger(y~.,
                    data = df,
                    mtry = ceiling(sqrt(p)),
                    importance = "impurity",
                    write.forest = FALSE, oob.error = FALSE)
  pval_rf_permu = importance_pvalues(rf_model,
                                     method = "altmann",
                                     formula = y~.,
                                     data = df, num.permutations = num_permutations)
  runtime = (proc.time() - t_start)[[3]]
  
  xnames = colnames(x)
  imp = pval_rf_permu[, "importance"]
  pval = pval_rf_permu[, "pvalue"]
  selection_result = data.frame(xnames = xnames, importance = imp,
                                pval = pval, selected = as.integer(pval < alpha))
  
  return(list(selection_result = selection_result, runtime = runtime,
              method = "RF Impurity", num_permutations = num_permutations, alpha = alpha))
  
}

rf_corrected_impurity_lm = function(x, y, num_permutations = 100, alpha = 0.05){
  df = cbind(y, x)
  p = ncol(x)
  t_start = proc.time()
  rf_model = ranger(y~.,
                    data = df,
                    mtry = ceiling(sqrt(p)),
                    importance = "impurity_corrected",
                    write.forest = FALSE, oob.error = FALSE)
  pval_rf_permu = importance_pvalues(rf_model,
                                     method = "altmann",
                                     formula = y~.,
                                     data = df, num.permutations = num_permutations)
  runtime = (proc.time() - t_start)[[3]]
  
  xnames = colnames(x)
  imp = pval_rf_permu[, "importance"]
  pval = pval_rf_permu[, "pvalue"]
  selection_result = data.frame(xnames = xnames, importance = imp,
                                pval = pval, selected = as.integer(pval < alpha))
  
  return(list(selection_result = selection_result, runtime = runtime,
              method = "RF Corrected Impurity", num_permutations = num_permutations, alpha = alpha))
  
}

rf_permutation_importance_lm = function(x, y, num_permutations = 100, alpha = 0.05){
  df = cbind(y, x)
  p = ncol(x)
  t_start = proc.time()
  rf_model = ranger(y~.,
                    data = df,
                    mtry = ceiling(sqrt(p)),
                    importance = "permutation",
                    write.forest = FALSE)
  
  pval_rf_permu = importance_pvalues(rf_model,
                                     method = "altmann",
                                     formula = y~.,
                                     data = df, num.permutations = num_permutations)
  runtime = (proc.time() - t_start)[[3]]
  
  xnames = colnames(x)
  imp = pval_rf_permu[, "importance"]
  pval = pval_rf_permu[, "pvalue"]
  selection_result = data.frame(xnames = xnames, importance = imp,
                                pval = pval, selected = as.integer(pval < alpha))
  
  return(list(selection_result = selection_result, runtime = runtime,
              method = "RF Permutation Importance", num_permutations = num_permutations, alpha = alpha))
  
}


cif_lm = function(x, y, num_permutations = 100, alpha = 0.05,
                   ntree = 100){
  df = cbind(y, x)
  p = ncol(x)
  n = nrow(x)
  mtry = ceiling(sqrt(p))
  formu = y~.
  
  t_start = proc.time()
  forest = cforest(formu, data = df,
                   controls = cforest_unbiased(mtry = mtry, ntree = ntree))
  invisible(capture.output(obs_varimp <- varimp(forest, pre1.0_0 = TRUE)))
  selection = names(obs_varimp)
  perm_mat = matrix(NA, ncol = length(selection), nrow = num_permutations, 
                    dimnames = list(1:num_permutations, selection))
  perm_dat = df
  for (pp in 1:num_permutations) {
    perm_id = sample(1:n, n, replace = FALSE)
    perm_dat[, "y"] = perm_dat[perm_id, "y"]
    perm_forest = cforest(formu, data = perm_dat, 
                          controls = cforest_unbiased(mtry = mtry, ntree = ntree))
    invisible(capture.output(perm_mat[pp, ] <- varimp(perm_forest, pre1.0_0 = T)))}
  pval = sapply(selection, function(x) sum(perm_mat[, x] >= obs_varimp[x]) / num_permutations)
  
  runtime = (proc.time() - t_start)[[3]]
  
  xnames = colnames(x)
  selection_result = data.frame(xnames = xnames, importance = obs_varimp,
                                pval = pval, selected = as.integer(pval < alpha))
  
  return(list(selection_result = selection_result, runtime = runtime,
              method = "Conditional Inference Forest", num_permutations = num_permutations,
              alpha = alpha, ntree = ntree))
  
}


dart20_lm = function(x, y, nskip = 2000,
                      ndpost = 5000, threshold = 0.5){
  
  df = cbind(y, x)
  p = ncol(x)
  ntree = 20
  t_start = proc.time()
  invisible(capture.output(dart <- wbart(x.train = df[, c(-1)], y.train = df$y,
                                         sparse = TRUE, ntree = ntree, 
                                         printevery = 1e6, ndpost = ndpost,
                                         nskip = nskip)))
  runtime = (proc.time() - t_start)[[3]]
  dart_varcount = dart$varcount
  
  varcount = matrix(0, nrow = nrow(dart_varcount), ncol = p)
  cnt = 0
  for(i in 1:p){
    if(is.factor(x[, i])){
      for(j in 1:length(levels(x[, i]))){
        cnt = cnt + 1
        varcount[, i] = varcount[, i] + dart_varcount[, cnt]
      }
    }
    else{
      cnt = cnt + 1
      varcount[, i] = dart_varcount[, cnt]
    }
  }
  
  mpvip_dart = apply(varcount, MARGIN = 2, function(x){
    mean(x>0)
  })
  xnames = colnames(x)
  selection_result = data.frame(xnames = xnames, importance = mpvip_dart,
                                selected = as.integer(mpvip_dart > threshold))
  
  return(list(selection_result = selection_result, runtime = runtime,
              method = "DART 20",
              nskip = nskip, ndpost = ndpost,
              threshold = threshold))
  
}


dart50_lm = function(x, y, nskip = 2000,
                     ndpost = 5000, threshold = 0.5){
  
  df = cbind(y, x)
  p = ncol(x)
  ntree = 50
  t_start = proc.time()
  invisible(capture.output(dart <- wbart(x.train = df[, c(-1)], y.train = df$y,
                                         sparse = TRUE, ntree = ntree, 
                                         printevery = 1e6, ndpost = ndpost,
                                         nskip = nskip)))
  runtime = (proc.time() - t_start)[[3]]
  dart_varcount = dart$varcount
  
  varcount = matrix(0, nrow = nrow(dart_varcount), ncol = p)
  cnt = 0
  for(i in 1:p){
    if(is.factor(x[, i])){
      for(j in 1:length(levels(x[, i]))){
        cnt = cnt + 1
        varcount[, i] = varcount[, i] + dart_varcount[, cnt]
      }
    }
    else{
      cnt = cnt + 1
      varcount[, i] = dart_varcount[, cnt]
    }
  }
  
  mpvip_dart = apply(varcount, MARGIN = 2, function(x){
    mean(x>0)
  })
  xnames = colnames(x)
  selection_result = data.frame(xnames = xnames, importance = mpvip_dart,
                                selected = as.integer(mpvip_dart > threshold))
  
  return(list(selection_result = selection_result, runtime = runtime,
              method = "DART 50",
              nskip = nskip, ndpost = ndpost,
              threshold = threshold))
  
}

dart200_lm = function(x, y, nskip = 2000,
                     ndpost = 5000, threshold = 0.5){
  
  df = cbind(y, x)
  p = ncol(x)
  ntree = 200
  t_start = proc.time()
  invisible(capture.output(dart <- wbart(x.train = df[, c(-1)], y.train = df$y,
                                         sparse = TRUE, ntree = ntree, 
                                         printevery = 1e6, ndpost = ndpost,
                                         nskip = nskip)))
  runtime = (proc.time() - t_start)[[3]]
  dart_varcount = dart$varcount
  
  varcount = matrix(0, nrow = nrow(dart_varcount), ncol = p)
  cnt = 0
  for(i in 1:p){
    if(is.factor(x[, i])){
      for(j in 1:length(levels(x[, i]))){
        cnt = cnt + 1
        varcount[, i] = varcount[, i] + dart_varcount[, cnt]
      }
    }
    else{
      cnt = cnt + 1
      varcount[, i] = dart_varcount[, cnt]
    }
  }
  
  mpvip_dart = apply(varcount, MARGIN = 2, function(x){
    mean(x>0)
  })
  xnames = colnames(x)
  selection_result = data.frame(xnames = xnames, importance = mpvip_dart,
                                selected = as.integer(mpvip_dart > threshold))
  
  return(list(selection_result = selection_result, runtime = runtime,
              method = "DART 200",
              nskip = nskip, ndpost = ndpost,
              threshold = threshold))
  
}


bart_lm = function(x, y, num_permutations = 100, alpha = 0.05,
                   nburn = 2000, npost = 5000, ntree = 20){
  invisible(capture.output(bart <- bartMachine(X = x, y = y, num_burn_in = nburn, 
                                               num_iterations_after_burn_in = npost, 
                                               num_trees = ntree, run_in_sample = FALSE, verbose = FALSE)))
  
  t_start = proc.time()
  invisible(capture.output(bart_var_sel <- var_selection_by_permute(bart, num_reps_for_avg = 10, 
                                                                    num_permute_samples = num_permutations,
                                                                    num_trees_for_permute = ntree, 
                                                                    alpha = alpha, plot = FALSE)))
  runtime = (proc.time() - t_start)[[3]]
  
  bart_vip = bart_var_sel$var_true_props_avg
  permute_mat = bart_var_sel$permute_mat
  xnames = names(bart_vip)
  tmp = data.frame(feature = xnames)
  
  pval_local = data.frame(feature = names(bart_vip), pval = sapply(1:ncol(permute_mat), function(j){
    mean(permute_mat[, j] >= bart_vip[j])
  }))
  pval_local = left_join(tmp, pval_local, by = "feature")
  
  vip_max = apply(permute_mat, MARGIN = 1, max)
  pval_global_max = data.frame(feature = names(bart_vip), pval = sapply(bart_vip, function(x){
    mean(vip_max >= x)
  }))
  pval_global_max = left_join(tmp, pval_global_max, by = "feature")
  
  mk = colMeans(permute_mat)
  sk = apply(permute_mat, MARGIN = 2, sd)
  C_star = (bart_vip - mk)/sk
  pval_global_se = data.frame(feature = names(bart_vip), pval = sapply(C_star, function(C){
    1 - mean(sapply(1:nrow(permute_mat), function(x){
      all(permute_mat[x, ] <= mk + sk * C)
    }))
  }))
  pval_global_se = left_join(tmp, pval_global_se, by = "feature")
  
  selection_result = data.frame(xnames = xnames,
                                vip=bart_vip,
                                pval_local=pval_local$pval,
                                pval_global_max=pval_global_max$pval,
                                pval_global_se=pval_global_se$pval,
                                selection_local=ifelse(pval_local$pval < alpha, 1, 0),
                                selection_global_max=ifelse(pval_global_max$pval < alpha, 1, 0),
                                selection_global_se=ifelse(pval_global_se$pval < alpha, 1, 0))
  return(list(selection_result = selection_result, runtime = runtime,
              method = "BART 20",
              num_permutations = num_permutations, alpha = alpha,
              nburn = nburn, npost = npost,
              ntree = ntree))
}


rf_min_depth_lm = function(x, y){
  df = cbind(y, x)
  t_start = proc.time()
  invisible(capture.output(rf_min_depth_res <- var.select(formula = y~., data = df, 
                                                          method = "md")))
  runtime = (proc.time() - t_start)[[3]]
  
  tmp = data.frame(feature = colnames(x))
  varselect_full = rf_min_depth_res$varselect
  varselect_full$feature = rownames(varselect_full)
  varselect_full$selection = ifelse(varselect_full$feature %in% rf_min_depth_res$topvars, 1, 0)
  varselect_full = left_join(tmp, varselect_full, by = "feature")
  
  xnames = tmp$feature
  selection_result = data.frame(xnames = xnames, variable_depth=varselect_full$depth,
                                selected = varselect_full$selection)
  
  return(list(selection_result = selection_result, runtime = runtime,
              method = "Random Forest Minimal Depth"))
  
}



xgb_permu_p_val = function(x_df, label, bst, obs_imp_df, nrounds, num_permutations, early_stopping = TRUE){
  sample_size = nrow(x_df)
  tmp = data.frame(Feature = colnames(x_df))
  obs_imp_df = left_join(tmp, obs_imp_df, by = "Feature")
  obs_gain_imp = ifelse(is.na(obs_imp_df$Gain), 0, obs_imp_df$Gain)
  obs_cover_imp = ifelse(is.na(obs_imp_df$Cover), 0, obs_imp_df$Cover)
  
  extreme_gain_cnt = 0
  extreme_cover_cnt = 0
  for(i in 1:num_permutations){
    shuffle_label = sample(label, size = sample_size, replace = FALSE)
    shuffle_xgb_df = xgb.DMatrix(data = x_df, label = shuffle_label)
    if(early_stopping){
      shuffle_bst = xgboost(data = shuffle_xgb_df, max.depth = bst$params$max_depth, nrounds = bst$best_iteration, objective = bst$params$objective, verbose = 0)
    }
    else{
      shuffle_bst = xgboost(data = shuffle_xgb_df, max.depth = bst$params$max_depth, nrounds = nrounds, objective = bst$params$objective, verbose = 0)
    }
    shuffle_imp_df = xgb.importance(model = shuffle_bst)
    shuffle_imp_df = left_join(tmp, shuffle_imp_df, by = "Feature")
    shuffle_gain_imp = ifelse(is.na(shuffle_imp_df$Gain), 0, shuffle_imp_df$Gain)
    shuffle_cover_imp = ifelse(is.na(shuffle_imp_df$Cover), 0, shuffle_imp_df$Cover)
    
    extreme_gain_cnt = extreme_gain_cnt + as.integer(shuffle_gain_imp >= obs_gain_imp)
    extreme_cover_cnt = extreme_cover_cnt + as.integer(shuffle_cover_imp >= obs_cover_imp)
  }
  gain_res = data.frame(Feature = colnames(x_df), Gain = obs_gain_imp, Gain_pval = (extreme_gain_cnt + 1) / (num_permutations + 1))
  cover_res = data.frame(Feature = colnames(x_df), Cover = obs_cover_imp, Cover_pval = (extreme_cover_cnt + 1) / (num_permutations + 1))
  return(list(gain = gain_res, cover = cover_res))
}


XGBoost_lm = function(x, y, num_permutations = 100, alpha = 0.05,
                       max_depth = 6, early_stopping_rounds = 10){
  df = cbind(y, x)
  p = ncol(x)
  label = y
  x_df = model.matrix(~., data = x)[, -1]
  xgb_df = xgb.DMatrix(data = x_df, label = label)
  
  t_start = proc.time()
  bst = xgboost(data = xgb_df, max.depth = max_depth, nrounds = 500, objective = "reg:squarederror", 
                early_stopping_rounds = early_stopping_rounds, verbose = 0)
  obs_imp_df = xgb.importance(model = bst)
  
  res = xgb_permu_p_val(x_df, label, bst, obs_imp_df, nrounds = 500,
                        num_permutations = num_permutations, early_stopping = TRUE)
  runtime = (proc.time() - t_start)[[3]]
  
  res$gain$selection = ifelse(res$gain$Gain_pval < alpha, 1, 0)
  res$cover$selection = ifelse(res$cover$Cover_pval < alpha, 1, 0)
  
  xnames = colnames(x_df)
  selection_result_gain = data.frame(xnames = xnames, importance = res$gain$Gain,
                                     pval = res$gain$Gain_pval,
                                     selected = res$gain$selection)
  
  selection_result_cover = data.frame(xnames = xnames, importance = res$cover$Cover,
                                      pval = res$cover$Cover_pval,
                                      selected = res$cover$selection)
  
  return(list(selection_result_gain = selection_result_gain, selection_result_cover = selection_result_cover,
              runtime = runtime,
              method = "XGBoost", num_permutations = num_permutations, alpha = alpha,
              max_depth = max_depth, early_stopping_rounds = early_stopping_rounds))
  
}


cp_select_1se = function(big_tree) {
  min_xerror = min(big_tree$cptable[, "xerror"])
  min_xerror_st = big_tree$cptable[which.min(big_tree$cptable[, "xerror"]), "xstd"]
  for(i in 1:nrow(big_tree$cptable)) {
    if(big_tree$cptable[i, "xerror"] < min_xerror + min_xerror_st){
      return(big_tree$cptable[i, "CP"])}
  }
}

cp_select_min = function(big_tree){
  return(big_tree$cptable[which.min(big_tree$cptable[, "xerror"]), "CP"])
}


tree_lm = function(x, y, cv_num = 5){
  df = cbind(y, x)
  p = ncol(x)
  xnames = colnames(x)
  tstart_fit = proc.time()
  tree_md = rpart(y~., data = df, 
                  control = rpart.control(xval = cv_num, maxsurrogate = 0, 
                                          cp = 0, minsplit = 2, minbucket = 1))
  time_fit = (proc.time() - tstart_fit)[[3]]
  
  tstart_1se = proc.time()
  tree_pruned_1se = prune(tree_md, cp = cp_select_1se(tree_md))
  runtime_1se = ((proc.time() - tstart_1se)[[3]]) + time_fit
  
  tstart_min = proc.time()
  tree_pruned_min = prune(tree_md, cp = cp_select_min(tree_md))
  runtime_min = ((proc.time() - tstart_min)[[3]]) + time_fit
  
  imp_tree = tree_md$variable.importance
  imp_var_names = names(imp_tree)
  imp_all = rep(0, p)
  names(imp_all) = xnames
  for(q in 1:length(imp_all)){
    if(names(imp_all[q]) %in% imp_var_names){
      imp_all[q] = imp_tree[which(imp_var_names == names(imp_all[q]))]} 
    else {
      imp_all[q]=0
    }
  }
  
  imp_tree_1se = tree_pruned_1se$variable.importance
  imp_var_names_tree_1se = names(imp_tree_1se)
  imp_all_tree_1se = rep(0, p)
  names(imp_all_tree_1se) = xnames
  for(q in 1:length(imp_all_tree_1se)){
    if(names(imp_all_tree_1se[q]) %in% imp_var_names_tree_1se){
      imp_all_tree_1se[q] = imp_tree_1se[which(imp_var_names_tree_1se == names(imp_all_tree_1se[q]))]} 
    else {
      imp_all_tree_1se[q]=0
    }
  }
  
  imp_tree_min = tree_pruned_min$variable.importance
  imp_var_names_tree_min = names(imp_tree_min)
  imp_all_tree_min = rep(0, p)
  names(imp_all_tree_min) = xnames
  if(length(imp_tree_min) != 0){
    for(q in 1:length(imp_all_tree_min)){
      if(names(imp_all_tree_min[q]) %in% imp_var_names_tree_min){
        imp_all_tree_min[q] = imp_tree_min[which(imp_var_names_tree_min == names(imp_all_tree_min[q]))]} 
      else {
        imp_all_tree_min[q]=0
      }
    }
  }
  
  selection_result = data.frame(xnames = xnames, importance_no_pruning=imp_all,
                                importance_1se = imp_all_tree_1se,
                                importance_min = imp_all_tree_min)
  return(list(selection_result = selection_result, runtime_1se = runtime_1se,
              runtime_min = runtime_min,
              method = "Decision Tree",
              cv_num = cv_num))
}
