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

if(!require(survival)){
  install.packages("survival")
}

if(!require(rpart)){
  install.packages("rpart")
}

library(ranger)
library(party)
library(BART)
library(randomForestSRC)
library(tidyverse)
library(xgboost)
library(survival)
library(rpart)

rf_impurity_cox = function(x, eventtime, status, num_permutations = 100, alpha = 0.05){
  df = cbind(eventtime, status, x)
  p = ncol(x)
  t_start = proc.time()
  rf_model = ranger(Surv(eventtime, status)~.,
                    data = df,
                    mtry = ceiling(sqrt(p)),
                    importance = "impurity",
                    write.forest = FALSE, oob.error = FALSE)
  pval_rf_permu = importance_pvalues(rf_model,
                                     method = "altmann",
                                     formula = Surv(eventtime, status)~.,
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

rf_corrected_impurity_cox = function(x, eventtime, status, num_permutations = 100, alpha = 0.05){
  df = cbind(eventtime, status, x)
  p = ncol(x)
  t_start = proc.time()
  rf_model = ranger(Surv(eventtime, status)~.,
                    data = df,
                    mtry = ceiling(sqrt(p)),
                    importance = "impurity_corrected",
                    write.forest = FALSE, oob.error = FALSE)
  pval_rf_permu = importance_pvalues(rf_model,
                                     method = "altmann",
                                     formula = Surv(eventtime, status)~.,
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

rf_permutation_importance_cox = function(x, eventtime, status, num_permutations = 100, alpha = 0.05){
  df = cbind(eventtime, status, x)
  p = ncol(x)
  t_start = proc.time()
  rf_model = ranger(Surv(eventtime, status)~.,
                    data = df,
                    mtry = ceiling(sqrt(p)),
                    importance = "permutation",
                    write.forest = FALSE)
  
  pval_rf_permu = importance_pvalues(rf_model,
                                     method = "altmann",
                                     formula = Surv(eventtime, status)~.,
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


cif_cox = function(x, eventtime, status, num_permutations = 100, alpha = 0.05,
                   ntree = 100){
  df = cbind(eventtime, status, x)
  p = ncol(x)
  n = nrow(x)
  mtry = ceiling(sqrt(p))
  formu = Surv(eventtime, status)~.
  
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
    perm_dat[, c("eventtime", "status")] = perm_dat[perm_id, c("eventtime", "status")]
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


dart20_cox = function(x, eventtime, status, nskip = 2000,
                      ndpost = 5000, keepevery = 10,
                      threshold = 0.5){
  
  df = cbind(eventtime, status, x)
  p = ncol(x)
  ntree = 20
  t_start = proc.time()
  invisible(capture.output(dart <- surv.bart(x.train = select(df, -eventtime, -status),
                                             times = df$eventtime, delta = df$status,
                                             sparse = TRUE, ntree = ntree, 
                                             printevery = 1e6, ndpost = ndpost,
                                             keepevery = keepevery, nskip = nskip)))
  runtime = (proc.time() - t_start)[[3]]
  dart_varcount = dart$varcount[, -1]
  
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
              nskip = nskip, ndpost = ndpost, keepevery = keepevery,
              threshold = threshold))
  
}


dart50_cox = function(x, eventtime, status, nskip = 2000,
                      ndpost = 5000, keepevery = 10,
                      threshold = 0.5){
  
  df = cbind(eventtime, status, x)
  p = ncol(x)
  ntree = 50
  t_start = proc.time()
  invisible(capture.output(dart <- surv.bart(x.train = select(df, -eventtime, -status),
                                             times = df$eventtime, delta = df$status,
                                             sparse = TRUE, ntree = ntree, 
                                             printevery = 1e6, ndpost = ndpost,
                                             keepevery = keepevery, nskip = nskip)))
  runtime = (proc.time() - t_start)[[3]]
  dart_varcount = dart$varcount[, -1]
  
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
              nskip = nskip, ndpost = ndpost, keepevery = keepevery,
              threshold = threshold))
  
}

dart200_cox = function(x, eventtime, status, nskip = 2000,
                      ndpost = 5000, keepevery = 10,
                      threshold = 0.5){
  
  df = cbind(eventtime, status, x)
  p = ncol(x)
  ntree = 200
  t_start = proc.time()
  invisible(capture.output(dart <- surv.bart(x.train = select(df, -eventtime, -status),
                                             times = df$eventtime, delta = df$status,
                                             sparse = TRUE, ntree = ntree, 
                                             printevery = 1e6, ndpost = ndpost,
                                             keepevery = keepevery, nskip = nskip)))
  runtime = (proc.time() - t_start)[[3]]
  dart_varcount = dart$varcount[, -1]
  
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
              nskip = nskip, ndpost = ndpost, keepevery = keepevery,
              threshold = threshold))
  
}


rf_min_depth_cox = function(x, eventtime, status){
  df = cbind(eventtime, status, x)
  t_start = proc.time()
  invisible(capture.output(rf_min_depth_res <- var.select(formula = Surv(eventtime, status)~., data = df, 
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


XGBoost_cox = function(x, eventtime, status, num_permutations = 100, alpha = 0.05,
                       max_depth = 6, early_stopping_rounds = 10){
  df = cbind(eventtime, status, x)
  p = ncol(x)
  label = ifelse(df$status == 1, df$eventtime, -df$eventtime)
  x_df = model.matrix(~., data = x)[, -1]
  xgb_df = xgb.DMatrix(data = x_df, label = label)
  
  t_start = proc.time()
  bst = xgboost(data = xgb_df, max.depth = max_depth, nrounds = 500, objective = "survival:cox", 
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


tree_cox = function(x, eventtime, status, cv_num = 5){
  df = cbind(eventtime, status, x)
  p = ncol(x)
  xnames = colnames(x)
  tstart_fit = proc.time()
  tree_md = rpart(Surv(eventtime, status)~., data = df, 
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
