library(tidyverse)
library(Matrix)
library(phyloseq)
library(GUniFrac)
library(quantreg)
library(parallel)
library(lmerTest)


source("./Code/util.R")

numCores = 1 # you can set up multiple cores for parallelization
nsim = 1000 # you can modify the number of simulations
level = 0.05 
job.id = 1  # customize job.id according to the `opt` data frame below

print(job.id)

opt = expand.grid(
  model = c('homo'),
  error = c('rnorm'),
  setting = c('null', 'alt'),
  m = c(100, 200, 500),
  ni = c(5, 10, 20),
  stringsAsFactors = F
) 

m = opt$m[job.id]
ni = opt$ni[job.id]
model = opt$model[job.id]
error = opt$error[job.id]
setting = opt$setting[job.id]




#149 taxa on 531 samples and 4 metadata
load(file = "./Data/genus.Rdata") # # Please contact authors for data

## rarify
if (!"otu_tab_rff.rds" %in% list.files("./Data/")){
  set.seed(2022)
  otu.tab = as.matrix(Dat[, c(2:150)])
  otu.names = colnames(otu.tab)
  colnames(otu.tab) = otu.names
  otu.tab.rff = GUniFrac::Rarefy(otu.tab)$otu.tab.rff
  saveRDS(otu.tab.rff,
          "./Data/otu_tab_rff.rds")
} else{
  otu.tab.rff = readRDS("./Data/otu_tab_rff.rds")
}


## step 1: logit reg
if (!"logit_coef.rds" %in% list.files("./Data/")){
  logit_coef = lapply(1:ncol(otu.tab.rff), function(j){
    df = data.frame(
      y = otu.tab.rff[, j],
      HBP = Dat$HBP30, 
      PA = scale(Dat$I18TOTAL),
      age = scale(Dat$ex9_age),
      dietscore = scale(Dat$dietscore2)) %>%
      mutate(ybin = as.integer(y > 0))
    logit.fit = glm(ybin ~ ., 
                    data =df %>% dplyr::select(-y), 
                    family = binomial(link = "logit"))
    logit.coef = coef(logit.fit)
  })
  names(logit_coef) = colnames(otu.tab.rff)
  saveRDS(logit_coef, "./Data/logit_coef.rds")
} else{
  logit_coef = readRDS("./Dataa/logit_coef.rds")
}


## step 2: quantile reg
if (!"quantile_coef.rds" %in% list.files("./Data")){
  taus = seq(0.01, 0.99, by = 0.01)
  quantile_coef = lapply(1:ncol(otu.tab.rff), function(j){
    df_nz = data.frame(
      y = dither(otu.tab.rff[, j], type = "right", value = 1),
      HBP = Dat$HBP30, 
      PA = scale(Dat$I18TOTAL),
      age = scale(Dat$ex9_age),
      dietscore = scale(Dat$dietscore2)) %>%
      filter(y > 0)
    quantile.fit = rq(y ~ ., data = df_nz, tau = taus)
    quantile.coef = coef(quantile.fit); coef.name = rownames(quantile.coef)
    quantile.coef = lapply(1:nrow(quantile.coef), function(i) approxfun(taus,quantile.coef[i,]))
    names(quantile.coef) = coef.name
    quantile.coef
  })
  names(quantile_coef) = colnames(otu.tab.rff)
  saveRDS(quantile_coef, "./Data/quantile_coef.rds" )
} else{
  quantile_coef = readRDS("./Data/quantile_coef.rds")
}



df = data.frame(
  HBP = Dat$HBP30, 
  PA = Dat$I18TOTAL,
  age = Dat$ex9_age,
  dietscore = Dat$dietscore2)

zerorate = apply(otu.tab.rff, 2, function(x) mean(x == 0))


set.seed(2022)
taus = c(0.1, 0.25, 0.5, 0.75, 0.9)

simu_func = function(i, model = c('homo'), 
                     error = c('rnorm'),
                     setting = c('null', 'alt')
){
  print(i)
  
  model = match.arg(model)
  error = match.arg(error)
  setting = match.arg(setting)
  
  ## step 1: sampling
  resample_idx = sample(nrow(otu.tab.rff), size = m, replace = F)
  
  ## step 2: generate covariates and outcome
  df_longi = lapply(resample_idx, function(id){
    HBP = rep(df$HBP[id],each = ni)
    age = seq(from = df$age[id], to = df$age[id] + 0.1 * (ni-1), by = 0.1)
    PA = rnorm(ni, mean = df$PA[id], sd = 1)
    DS = round(rnorm(ni, mean = df$dietscore[id], sd = 1))
    data.frame(1, HBP, PA, age, dietscore = DS)
  }) %>% do.call(rbind, .)
  
  df_longi$PA = scale(df_longi$PA)
  df_longi$age = scale(df_longi$age)
  df_longi$dietscore = scale(df_longi$dietscore)
  
  
  X = df_longi$HBP
  Z = as.matrix(df_longi %>% dplyr::select(-HBP))
  study = factor(rep(c(1:m),each=ni))
  
  if (setting == 'null'){
    
    ### generate ybin
    hi = rnorm(m, 0, 1)
    h = rep(hi, each = ni)
    ybin = lapply(1:ncol(otu.tab.rff), function(j){
      logit_coef_j = logit_coef[[j]]
      logit_coef_j['HBP'] = 0
      eta = as.matrix(df_longi) %*% logit_coef_j + h
      eta_j = log((1-zerorate[j])/zerorate[j])
      eta = scale(eta) + eta_j
      probs = exp(eta)/(1+ exp(eta))
      probs[is.nan(probs)] = 1 # if exp(eta) = Inf
      ybin = sapply(probs, function(prob) rbinom(1, 1, prob = prob))
    }) %>% do.call(cbind, .)
    colnames(ybin) = colnames(otu.tab.rff)
    
    ### generate y
    u = runif(m*ni, 0, 1)
    hi = rnorm(m, 0, 1)
    h = rep(hi, each = ni)
    
    yij = lapply(1:ncol(otu.tab.rff), function(j){
      func_j = quantile_coef[[j]]
      quantile_coef_j = lapply(func_j, function(f) f(u)) %>%
        do.call(cbind, .)
      quantile_coef_j[,'HBP'] = 0
      yij = round(rowSums(as.matrix(df_longi) * quantile_coef_j) + h)
    }) %>% do.call(cbind, .)
    yij[yij < 0] = 0
    colnames(yij) = colnames(otu.tab.rff)
    yij = yij * ybin
    
    
  }
  
  
  if (setting == 'alt'){
    ## step 1: sampling
    resample_idx = sample(nrow(otu.tab.rff), size = m, replace = F)
    
    
    hi = rnorm(m, 0, 1)
    h = rep(hi, each = ni)
    ybin = lapply(1:ncol(otu.tab.rff), function(j){
      logit_coef_j = logit_coef[[j]]
      eta = as.matrix(df_longi) %*% logit_coef_j + h
      eta_j = log((1-zerorate[j])/zerorate[j])
      eta = scale(eta) + eta_j
      probs = exp(eta)/(1+ exp(eta))
      probs[is.nan(probs)] = 1 # if exp(eta) = Inf
      ybin = sapply(probs, function(prob) rbinom(1, 1, prob = prob))
    }) %>% do.call(cbind, .)
    colnames(ybin) = colnames(otu.tab.rff)
    
    ### generate y
    u = runif(m*ni, 0, 1)
    hi = rnorm(m, 0, 1)
    h = rep(hi, each = ni)
    
    yij = lapply(1:ncol(otu.tab.rff), function(j){
      func_j = quantile_coef[[j]]
      quantile_coef_j = lapply(func_j, function(f) f(u)) %>%
        do.call(cbind, .)
      yij = round(rowSums(as.matrix(df_longi) * quantile_coef_j) + h)
    }) %>% do.call(cbind, .)
    yij[yij < 0] = 0
    colnames(yij) = colnames(otu.tab.rff)
    yij = yij * ybin
  }
  
  # quantile test
  result = lapply(1:ncol(yij), function(t){
    print(paste("taxa", t))
    yijt = yij[, t]
    yijt = dither(yijt, type = "right", value = 1)
    ybint = ybin[, t]
    
    na_idx = which(is.na(yijt))
    
    result = suppressMessages(try(
      zero_inf_ind_rank_score_test(X1m = X[-na_idx],
                                   Z[-na_idx, ],
                                   study[-na_idx],
                                   y = yijt[-na_idx],
                                   ybin = ybint[-na_idx],
                                   taus = taus)
    ))
    ## linear regression and linear mixed model
    lmfit = lm(yijt ~ X + Z - 1)
    pval_lm = summary(lmfit)$coefficients['X', 'Pr(>|t|)']
    
    lmmfit = suppressMessages(lmerTest::lmer(yijt ~ X + Z - 1 + (1 | study)))
    pval_lmm = summary(lmmfit)$coefficients['X', 'Pr(>|t|)']
    result_bench = c(pval_lm, pval_lmm); names(result_bench) = c('lm','lmm')
    result = c(result, result_bench)
    return(result)
  })
  result_old = do.call(rbind, result)
  
  ## remove some rare taxa
  otu.tab = yij
  metadata = df_longi
  metadata$study = study
  rownames(otu.tab) = rownames(metadata) = rownames(ybin) = 1:nrow(otu.tab)
  
  missing_idx = sapply(1:ncol(otu.tab), function(j) length(which(is.na(otu.tab[, j]))))
  remove_taxa = which(missing_idx >= ni * m * 0.9)
  if(length(remove_taxa) > 0){
    otu.tab = otu.tab[, -remove_taxa] %>% na.omit()
  } else{
    otu.tab = otu.tab %>% na.omit()
  }
  
  metadata = metadata[as.numeric(rownames(otu.tab)), ]
  
  ## existing methods
  pval_linda = suppressWarnings(
    suppressMessages(
      #hush(
        LinDA_test(otu.tab, metadata)$pval
      #)
    )
  )
  
  pval_ma = suppressWarnings(
    suppressMessages(
      #hush(
        MaAsLin2_test(otu.tab, metadata, output = paste0('./Scratch/tmp_', job.id))$pval
      #)
    )
  )
  
  # pval_ldm = suppressWarnings(
  #   suppressMessages(
  #     hush(
  #       LDM_test(otu.tab, metadata)$pval
  #     )
  #   )
  # )
  
  pval_nb = suppressWarnings(
    suppressMessages(
      #hush(
        NBZIMM_test(otu.tab, metadata, method = 'nb')$pval
      #)
    )
  )
  pval_zinb = suppressWarnings(
    suppressMessages(
      #hush(
        NBZIMM_test(otu.tab, metadata, method = 'zinb')$pval
      #)
    )
  )
  pval_zig_count = suppressWarnings(
    suppressMessages(
      #hush(
        NBZIMM_test(otu.tab, metadata, method = 'zig_count')$pval
      #)
    )
  )
  pval_zig_prop = suppressWarnings(
    suppressMessages(
      #hush(
        NBZIMM_test(otu.tab, metadata, method = 'zig_prop')$pval
      #)
    )
  )
  
  
  result = cbind(pval_linda, pval_ma, #pval_ldm, 
                 pval_nb, pval_zinb, pval_zig_count, pval_zig_prop)
  colnames(result) = c('LinDA', 'MaAsLin2', 'NB', 'ZINB', 'ZIG_count', 'ZIG_prop')
  result = cbind(result_old, result)
  
  
  ## result summary
  power = apply(result, 2, function(x) mean(x < level, na.rm = T))
  
  
  print(power)
  print("----------------")
  return(power)
  
}
result = mclapply(1:nsim, simu_func, 
                  model = model, 
                  error = error,
                  setting = setting, 
                  mc.cores = numCores)

dir.create("./Results/sim2_complete")
saveRDS(result, paste0("./Results/sim2_complete/power_",job.id,".rds"))


