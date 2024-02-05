library(tidyverse)
library(Matrix)
library(phyloseq)
library(GUniFrac)
library(quantreg)
library(parallel)
library(lmerTest)

source('./Code/util.R')


numCores = 1 # you can set up multiple cores for parallelization
nsim = 1000 # you can modify the number of simulations
level = 0.05 
job.id = 1  # customize job.id according to the `opt` data frame below

print(job.id)

opt = expand.grid(
  model = c('homo'),
  error = c('rnorm','rchisq'),
  setting = c('null', 'alt'),
  sigma = 1, #c(0.5, 1, 2),
  m = c(100, 200, 500),
  ni = c(3, 5, 10, 20),
  stringsAsFactors = F
) 


m = opt$m[job.id]
ni = opt$ni[job.id]
model = opt$model[job.id]
error = opt$error[job.id]
setting = opt$setting[job.id]
sigma = opt$sigma[job.id]



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

set.seed(2022)
taus = c(0.1, 0.25, 0.5, 0.75, 0.9)


random_err = function(n, error = 'rnorm'){
  if (error == 'rnorm'){
    return(rnorm(n))
  } else if (error == 'rchisq'){
    return(rchisq(n, df = 2))
  } else if (error == 'rt'){
    return(rt(n, df = 2))
  } else if (error == 'rcauchy'){
    return(rcauchy(n))
  }
  return(0)
}

nonHBP_idx = which(Dat$HBP30 == 0)
HBP_idx = which(Dat$HBP30 == 1)
idx1 = 22
idx2 = 51
idx3 = 35
idx4 = 42

idx = c(idx1, idx2, idx3, idx4) # four representative taxa

simu_func = function(i, model = c('homo'), 
                     error = c('rnorm','rchisq','rt','rcauchy'),
                     setting = c('null', 'alt'),
                     n_taxa = 50
){
  print(i)
  
  model = match.arg(model)
  error = match.arg(error)
  setting = match.arg(setting)
  
  if (setting == 'null'){
    yij = lapply(idx, function(t){
      
      y = otu.tab.rff[, t]
      
      ## resampling from non-HBP sample only
      
      otu_tab = sapply(1:n_taxa, function(i){
        resample_tau = runif(m, 0, 1)
        y = quantile(y[nonHBP_idx], probs = resample_tau, type = 3) + 1
        names(y) = NULL
        
        ## step 2: generate longitudinal data
        yij = lapply(y, function(yi){
          yij = exp(log(yi) + random_err(ni, error))
          yij = ifelse(yij < 1, 0, round(yij) - 1)
        }) %>% do.call(c, .)
        
        return(yij)
      })
      colnames(otu_tab) = paste0('taxa',t, '_rep', 1:ncol(otu_tab))
      return(otu_tab)
    }) %>% do.call(cbind, .)
    study = factor(rep(c(1:m),each=ni))
    Z = matrix(1, nrow = nrow(yij), ncol = 1); colnames(Z) = 'Z';
    X = rep(rbinom(m, 1, 0.5), each = ni)  
    ybin = (yij > 0) * 1
    
  }
  
  
  if (setting == 'alt'){
    yij = lapply(idx, function(t){
      y = otu.tab.rff[, t]
      ## resampling from both HBP and non-HBP samples
        otu_tab = sapply(1:n_taxa, function(i){
        resample_tau1 = runif(m/2, 0, 1)
        resample_tau2 = runif(m/2, 0, 1)
        y1 = quantile(y[nonHBP_idx], probs = resample_tau1, type = 3) + 1
        y2 = quantile(y[HBP_idx], probs = resample_tau2, type = 3) + 1
        y = c(y1, y2); names(y) = NULL
        
        ## step 2: generate longitudinal data
        yij = lapply(y, function(yi){
          yij = exp(log(yi) + random_err(ni, error))
          yij = ifelse(yij < 1, 0, round(yij) - 1)
        }) %>% do.call(c, .)
        
        #ybin = as.integer(yij > 0)
        return(yij)
      })
      colnames(otu_tab) = paste0('taxa',t, '_rep', 1:ncol(otu_tab))
      return(otu_tab)
    }) %>% do.call(cbind, .)
    
    study = factor(rep(c(1:m),each=ni))
    Z = matrix(1, nrow = nrow(yij), ncol = 1); colnames(Z) = 'Z'; 
    X = c(rep(0, m/2 * ni), rep(1, m/2 * ni))
    ybin = (yij > 0) * 1
    
  }
  
  otu.tab = yij
  metadata = data.frame(X, Z, study)
  rownames(otu.tab) = rownames(metadata) = rownames(ybin) = 
    paste0("subject", rep(1:m, each=ni), "_visit", rep(1:ni, m))
  
  ## quantile test 
  result = lapply(1:ncol(otu.tab), function(t){
    #print(paste('taxa', t))
    yijt = otu.tab[, t]
    yijt = dither(yijt, type = "right", value = 1)
    ybint = ybin[, t]
    
    
    result = suppressMessages(try(
      zero_inf_ind_rank_score_test(X1m = X, 
                                   Z, 
                                   study, 
                                   y = yijt, 
                                   ybin = ybint, 
                                   taus = taus)
    ))
    if ('try-error' %in% class(result)){
      result = rep(NA, 24)
    }
    ## linear regression and linear mixed model
    lmfit = lm(yijt ~ X + Z - 1)
    pval_lm = summary(lmfit)$coefficients['X', 'Pr(>|t|)']
    
    lmmfit = suppressMessages(lmerTest::lmer(yijt ~ X + Z - 1 + (1 | study)))
    pval_lmm = summary(lmmfit)$coefficients['X', 'Pr(>|t|)']
    result_bench = c(pval_lm, pval_lmm); names(result_bench) = c('lm','lmm')
    result = c(result, result_bench)
    return(result)
  })
  
  pval = do.call(rbind, result)
  
  
  ## existing methods
  pval_linda = suppressWarnings(
    suppressMessages(
      #hush(
        LinDA_test_unadj(otu.tab, metadata)$pval
      #)
    )
  )
  
  pval_ma = suppressWarnings(
    suppressMessages(
      #hush(
        MaAsLin2_test_unadj(otu.tab, metadata, output = paste0('./tmp_', job.id))$pval
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
        NBZIMM_test_unadj(otu.tab, metadata, method = 'nb')$pval
      #)
    )
  )
  pval_zinb = suppressWarnings(
    suppressMessages(
      #hush(
        NBZIMM_test_unadj(otu.tab, metadata, method = 'zinb')$pval
      #)
    )
  )
  pval_zig_count = suppressWarnings(
    suppressMessages(
      #hush(
        NBZIMM_test_unadj(otu.tab, metadata, method = 'zig_count')$pval
      #)
    )
  )
  pval_zig_prop = suppressWarnings(
    suppressMessages(
      #hush(
        NBZIMM_test_unadj(otu.tab, metadata, method = 'zig_prop')$pval
      #)
    )
  )
  
  
  pval_combo = cbind(pval_linda, pval_ma, #pval_ldm, 
                     pval_nb, pval_zinb, pval_zig_count, pval_zig_prop)
  colnames(pval_combo) = c('LinDA', 'MaAsLin2', 'NB', 'ZINB', 'ZIG_count', 'ZIG_prop')
  pval = cbind(pval, pval_combo)
  rownames(pval) = colnames(otu.tab)
  print(str(pval))
  return(pval)
  
}
result = mclapply(1:nsim, simu_func, 
                  model = model, 
                  error = error,
                  setting = setting, 
                  mc.cores = numCores)

dir.create("./Results")
saveRDS(result, paste0("./Results/pval_",job.id,".rds"))


