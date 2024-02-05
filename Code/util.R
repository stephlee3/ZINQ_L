library(lme4)
library(tidyverse)
library(quantreg)
library(CompQuadForm)
library(expm)
library(FactoMineR)
library(Matrix)
library(MASS)
library(LDM)
library(MicrobiomeStat)
library(Maaslin2)
library(NBZIMM)


## force suppress output from function: 
## ref: https://stackoverflow.com/questions/2723034/suppress-output-of-a-function
hush=function(code){ 
  sink("/dev/null") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

score_func = function(tau, u){
  tau - ifelse(u < 0, 1, 0)
}

LDM_test = function(otu.tab, metadata, n.rej.stop = 100){
  metadata_longi = metadata %>%
    group_by(study) %>%
    mutate(time = 1:n()) %>%
    ungroup()
  
  min_samples = metadata_longi %>% group_by(study) %>%
    summarise(count = n()) %>% pull(count) %>% min()
  
  metadata_longi = metadata_longi %>%
    mutate(idx = 1:nrow(.)) %>%
    filter(time <= min_samples) %>%
    as.data.frame()
  
  otu.tab.cut = otu.tab[metadata_longi$idx,]

  metadata_longi$otu.tab.cut = otu.tab.cut
  
  res.ldm = ldm(
    formula = otu.tab.cut| (PA + age + dietscore) ~ HBP,
    data = metadata_longi,
    seed = 34794,
    cluster.id = study,
    perm.within.type = "none",
    perm.between.type = "free",
    n.rej.stop = n.rej.stop
  )
  pval = res.ldm$p.otu.omni
  pval.adj = res.ldm$q.otu.omni
  return(list(pval = pval,
              pval.adj = pval.adj))
}

LinDA_test = function(otu.tab, metadata){
  res.linda = linda(
    t(otu.tab), 
    metadata, 
    formula = '~HBP + PA + age + dietscore + (1|study)',
    feature.dat.type = 'count',
    prev.filter = 0, 
    is.winsor = TRUE, 
    outlier.pct = 0.03,
    p.adj.method = "BH", 
    alpha = 0.05)
  summary.tab = res.linda$output$HBP
  return(list(pval = summary.tab$pvalue, 
              pval.adj = summary.tab$padj))
}

LinDA_test_unadj = function(otu.tab, metadata){
  res.linda = linda(
    t(otu.tab), 
    metadata, 
    formula = '~X + (1|study)',
    feature.dat.type = 'count',
    prev.filter = 0, 
    is.winsor = TRUE, 
    outlier.pct = 0.03,
    p.adj.method = "BH", 
    alpha = 0.05)
  summary.tab = res.linda$output$X
  return(list(pval = summary.tab$pvalue, 
              pval.adj = summary.tab$padj))
}

MaAsLin2_test = function(otu.tab, metadata, output = './tmp'){
  if(!dir.exists(output)){
    dir.create(output)
  }
  
  res.ma = Maaslin2(
    input_data = otu.tab, 
    input_metadata = metadata, 
    output = output, 
    #transform = "AST",
    fixed_effects = c('HBP','X1','PA', 'age', 'dietscore'),
    random_effects = c('study'),
    #normalization = 'NONE',
    plot_heatmap = F,
    plot_scatter = F,
    standardize = FALSE)
  summary.tab = res.ma$results %>% filter(name == 'HBP') 
  summary.tab = summary.tab[match(colnames(otu.tab), summary.tab$feature), ]
  
  pval = summary.tab %>% pull(pval)
  qval = summary.tab %>% pull(qval)
  return(list(pval = pval,
              pval.adj = qval))
}

MaAsLin2_test_unadj = function(otu.tab, metadata, output = './tmp'){
  if(!dir.exists(output)){
    dir.create(output)
  }
  
  res.ma = Maaslin2(
    input_data = otu.tab, 
    input_metadata = metadata, 
    output = output, 
    #transform = "AST",
    fixed_effects = c('X'),
    random_effects = c('study'),
    #normalization = 'NONE',
    plot_heatmap = F,
    plot_scatter = F,
    standardize = FALSE)
  summary.tab = res.ma$results %>% filter(name == 'X') 
  summary.tab = summary.tab[match(colnames(otu.tab), summary.tab$feature), ]
  
  pval = summary.tab %>% pull(pval)
  qval = summary.tab %>% pull(qval)
  return(list(pval = pval,
              pval.adj = qval))
}

NBZIMM_test = function(otu.tab, metadata, method = 'nb'){
  N = rowSums(otu.tab)
  metadata$N = N
  if (method == 'nb'){
    fit = mms(y = otu.tab, 
              fixed = ~  HBP+ PA + age + dietscore + stats::offset(log(N)), 
              random = ~ 1 | study, 
              data = metadata,
              min.p = 0, 
              method = 'nb')
  }
  
  if (method == 'zinb'){
    fit = mms(y = otu.tab, 
              fixed = ~  HBP+ PA + age + dietscore + stats::offset(log(N)), 
              random = ~ 1 | study, 
              data = metadata,
              min.p = 0, 
              zi_fixed = ~ 1, #HBP + PA + age + dietscore,
              zi_random = NULL, #~ 1 | study,
              method = 'zinb')
  }
  
  if (method == 'zig_count'){
    otu.tab.log = log2(otu.tab + 1)
    fit = mms(y = otu.tab.log, 
              fixed = ~  HBP+ PA + age + dietscore + stats::offset(log(N)), 
              random = ~ 1 | study, 
              data = metadata,
              min.p = 0, 
              zi_fixed = ~1, #~ HBP + PA + age + dietscore,
              zi_random = NULL,
              method = 'zig')
    
  }
  
  if (method == 'zig_prop'){
    otu.tab.prop = t(apply(otu.tab, 1, function(x) x/sum(x)))
    otu.tab.asin = asin(sqrt(otu.tab.prop))
    fit = mms(y = otu.tab.asin, 
              fixed = ~  HBP+ PA + age + dietscore, 
              random = ~ 1 | study, 
              data = metadata,
              min.p = 0, 
              zi_fixed = ~1, #~ HBP + PA + age + dietscore,
              zi_random = NULL, #~ 1 | study,
              method = 'zig')
    
  }
  
  res = NBZIMM::fixed(fit)$dist
  rn = rownames(res)
  res = res[grepl('HBP', rn), ]
  
  ## handle missing taxa
  taxa = res$responses
  pval = pval.adj = rep(NA, ncol(otu.tab))
  pval[match(taxa, colnames(otu.tab))] = res$pvalue
  pval.adj[match(taxa,colnames(otu.tab))] = res$padj
  return(list(pval = pval,
              pval.adj = pval.adj))
  
}

NBZIMM_test_unadj = function(otu.tab, metadata, method = 'nb'){
  N = rowSums(otu.tab)
  metadata$N = N
  if (method == 'nb'){
    fit = mms(y = otu.tab, 
              fixed = ~  X + stats::offset(log(N)), 
              random = ~ 1 | study, 
              data = metadata,
              min.p = 0, 
              method = 'nb')
  }
  
  if (method == 'zinb'){
    fit = mms(y = otu.tab, 
              fixed = ~  X + stats::offset(log(N)), 
              random = ~ 1 | study, 
              data = metadata,
              min.p = 0, 
              zi_fixed = ~ 1, #HBP + PA + age + dietscore,
              zi_random = NULL, #~ 1 | study,
              method = 'zinb')
  }
  
  if (method == 'zig_count'){
    otu.tab.log = log2(otu.tab + 1)
    fit = mms(y = otu.tab.log, 
              fixed = ~ X + stats::offset(log(N)), 
              random = ~ 1 | study, 
              data = metadata,
              min.p = 0, 
              zi_fixed = ~1, #~ HBP + PA + age + dietscore,
              zi_random = NULL,
              method = 'zig')
    
  }
  
  if (method == 'zig_prop'){
    otu.tab.prop = t(apply(otu.tab, 1, function(x) x/sum(x)))
    otu.tab.asin = asin(sqrt(otu.tab.prop))
    fit = mms(y = otu.tab.asin, 
              fixed = ~  X , 
              random = ~ 1 | study, 
              data = metadata,
              min.p = 0, 
              zi_fixed = ~1, #~ HBP + PA + age + dietscore,
              zi_random = NULL, #~ 1 | study,
              method = 'zig')
    
  }
  
  res = NBZIMM::fixed(fit)$dist
  rn = rownames(res)
  res = res[grepl('X', rn), ]
  
  ## handle missing taxa
  taxa = res$responses
  pval = pval.adj = rep(NA, ncol(otu.tab))
  pval[match(taxa, colnames(otu.tab))] = res$pvalue
  pval.adj[match(taxa,colnames(otu.tab))] = res$padj
  return(list(pval = pval,
              pval.adj = pval.adj))
  
}


NBZIMM_indiv_test = function(y, metadata, method = 'nb'){
  if (method == 'nb'){
    fit = glmm.nb(
              fixed = y ~  X + Z - 1, #+ stats::offset(log(N)), 
              random = ~ 1 | study,
              data = metadata
              )
  }
  
  if (method == 'zinb'){
    fit = glmm.zinb(
              fixed = y ~  X + Z - 1, #+ stats::offset(log(N)),
              random = ~ 1 | study, 
              data = metadata, 
              zi_fixed = ~ 1, #HBP + PA + age + dietscore,
              zi_random = NULL #~ 1 | study,
              )
  }
  
  if (method == 'zig_count'){
    y = log2(y + 1)
    fit = lme.zig(
              fixed = y ~  X + Z - 1, #+ stats::offset(log(N)),
              random = ~ 1 | study, 
              data = metadata, 
              zi_fixed = ~1, #~ HBP + PA + age + dietscore,
              zi_random = NULL
             )
    
  }
  
  if (method == 'zig_prop'){
    y = asin(sqrt(y)) # 
    fit = lme.zig(
              fixed = y ~  X + Z - 1, #+ stats::offset(log(N)),
              random = ~ 1 | study, 
              meatadata = metadata, 
              zi_fixed = ~1, #~ HBP + PA + age + dietscore,
              zi_random = NULL #~ 1 | study,
              )
    
  }
  
  res = NBZIMM::fixed(fit)$dist
  rn = rownames(res)
  pval = res[grepl('X', rn), 'pvalue']
  return(pval)
}




ZIBR_test = function(otu.tab, metadata){
  metadata_longi = metadata %>%
    group_by(study) %>%
    mutate(time = 1:n()) %>%
    ungroup()
  
  min_samples = metadata_longi %>% group_by(study) %>%
    summarise(count = n()) %>% pull(count) %>% min()
  
  metadata_longi = metadata_longi %>%
    mutate(idx = 1:nrow(.)) %>%
    filter(time <= min_samples)
  
  
  otu.tab.prop = t(apply(otu.tab, 1, function(x) x/sum(x)))[metadata_longi$idx,]
  
  beta_cov = logit_cov = metadata_longi %>% dplyr::select(HBP, PA, age, dietscore)
  
  res = sapply(1:ncol(otu.tab.prop), function(j){
    y = otu.tab.prop[, j]
    zibr.fit <- zibr(logistic.cov = logit_cov, 
                     beta.cov = beta_cov, 
                     Y = y, 
                     subject.ind = as.numeric(metadata_longi$study),
                     time.ind = metadata_longi$time)
    zibr.fit$joint.p['HBP']
  })
  return(res)
}


#####

QRank_test = function(X1m,Z, study, y, ybin, taus = c(0.1, 0.25, 0.5, 0.75, 0.9)){
  qfit = QRank(gene = y, snp = X1m, cov = Z, tau = taus)
  pval = qfit$composite.pvalue
  return(pval)
}

logit_mixed_test = function(X1m, Z, study, y, ybin){
  logit_model = try(glmer(ybin ~ -1 + X1m + Z + (1|study),
                          family = binomial(link = "logit")), silent = T)
  
  if ('try-error' %in% class(logit_model)){
    logit_pval = NA
    names(logit_pval) = 'logit_Wald'
    logit_pval_LRT = NA
    names(logit_pval_LRT) = 'logit_LRT'
  } else{
    logit_summ = summary(logit_model)
    logit_pval = logit_summ$coefficients['X1m', 'Pr(>|z|)']
    names(logit_pval) = 'logit_Wald'
    
    logit_model_null = try(glmer(ybin ~ -1 + Z + (1|study),
                                 family = binomial(link = "logit")), silent = T)
    logit_pval_LRT = anova(logit_model, logit_model_null)$`Pr(>Chisq)`[2]
    names(logit_pval_LRT) = 'logit_LRT'
  }
  pval = c(logit_pval, logit_pval_LRT)
  return(pval)
}


minP_combine_test = function(pval,
                             taus = c(0.1, 0.25, 0.5, 0.75, 0.9),
                             Sigma.hat,
                             M = 10000){
  width = length(taus)
  df = 1
  t.obs = min(pval)
  qmin.quantile = rep(qchisq(1-t.obs, df= df), width)
  
  eSig = eigen(Sigma.hat,symmetric = T)
  eSig.value = eSig$values
  if (any(eSig.value < -1e-10)){
    eSig.vector = eSig$vectors[,which(eSig.value>1e-10)]
    Lambda = as.matrix(diag(eSig.value[which(eSig.value>1e-10)]))
    Sigma.hat = eSig.vector %*% Lambda %*% t(eSig.vector)
  }
  beta.sim = MASS::mvrnorm(n=M, mu=rep(0, width*df), Sigma=Sigma.hat)
  prob.quantile = mean(apply(beta.sim, 1, function(z){
    obs.sim = z^2/diag(Sigma.hat)
    all(obs.sim < qmin.quantile)
  }) )
  min_pval = 1-(1-t.obs)*prob.quantile
  return(min_pval)
  
}

cauchy_combine_test = function(zerorate, 
                               taus = c(0.1, 0.25, 0.5, 0.75, 0.9),
                               pval,
                               method = c('unequal',
                                          'unequal_clip',
                                          'equal',
                                          'equal_clip')){
  method = match.arg(method)
  
  if (method == 'unequal'){
    w = taus*(taus <= 0.5) + (1-taus)*(taus > 0.5)
    w = c(zerorate, w / sum(w) * (1-zerorate))
  } 
  
  if (method == 'unequal_clip'){
    taus_clip = c(0.25, 0.5, 0.75)
    w = taus_clip*(taus_clip <= 0.5) + (1-taus_clip)*(taus_clip > 0.5)
    w = c(zerorate, 0, w / sum(w) * (1 - zerorate),0)
  }
  
  logit_pval = pval[1]
  
  if (method == 'equal'){
    if (is.na(logit_pval)){
      w = c(0, rep(1,5)/5)
    } else{
      w = rep(1,6)/6
    }
  }
  
 
  if (method == 'equal_clip'){
    if (is.na(logit_pval)){
      w = c(0, 0, rep(1,3)/3, 0)
    } else{
      w = c(1/4, 0, rep(1/4, 3), 0)
    }
  }
  
  stats_cauchy = sum(w*tan((0.5-pval)*pi))
  pval_cauchy = 1 - pcauchy(stats_cauchy)
  return(pval_cauchy)

}


zero_inf_ind_rank_score_test = function(X1m,Z,study,y, ybin, taus){
  m = length(unique(study))
  N = length(X1m)
  ni = as.vector(table(study))
  L = sum(ni * (ni - 1))
  q = ncol(Z)
  
  nonzero_idx = which(ybin > 0)
  #print(length(nonzero_idx))
  X1m_tilde = X1m * ybin
  Z_tilde = Z * ybin
  
  ## GLMM - logit reg
  logit_model = try(glmer(ybin ~ -1 + X1m + Z + (1|study),
                          family = binomial(link = "logit")), silent = T)
  
  
  if ('try-error' %in% class(logit_model)){
    logit_pval = NA
    names(logit_pval) = 'logit_Wald'
    logit_pval_LRT = NA
    names(logit_pval_LRT) = 'logit_LRT'
  } else{
    logit_summ = summary(logit_model)
    logit_pval = logit_summ$coefficients['X1m', 'Pr(>|z|)']
    names(logit_pval) = 'logit_Wald'
    
    logit_model_null = try(glmer(ybin ~ -1 + Z + (1|study),
                                 family = binomial(link = "logit")), silent = T)
    logit_pval_LRT = anova(logit_model, logit_model_null)$`Pr(>Chisq)`[2]
    names(logit_pval_LRT) = 'logit_LRT'
  }
    
  
  
  
  
  
  X1m_star = try(X1m_tilde - Z_tilde %*% solve(crossprod(Z_tilde), crossprod(Z_tilde, X1m_tilde)))
  if ('try-error' %in% class(X1m_star)){
    ZtZ= crossprod(Z_tilde)
    eZ = eigen(ZtZ,symmetric = T)
    eZ.value = eZ$values
    if (any(abs(eZ.value) < 1e-10)){
      n = length(which(abs(eZ.value) > 1e-10))
      eZ.vector = eZ$vectors[,which(abs(eZ.value)>1e-10), drop = F]
      if (n == 1){
        Lambda_inv = matrix(eZ.value[which(abs(eZ.value)>1e-10)]^(-1), 1, 1)
      } else{
        Lambda_inv = as.matrix(diag(eZ.value[which(abs(eZ.value)>1e-10)]^(-1)))
      }
      ZtZ_inv = eZ.vector %*% Lambda_inv %*% t(eZ.vector)
    }
    X1m_star = X1m_tilde - Z_tilde %*% ZtZ_inv %*% crossprod(Z_tilde, X1m_tilde)
  }
  df = data.frame(y, ybin, X1m = X1m, Z, study, idx = 1:length(y)) %>% filter(ybin > 0)
  nonzero_ni = as.vector(table(df$study))
  nonzero_L = sum(nonzero_ni * (nonzero_ni-1))
  
  qformula = as.formula(paste0("y ~ -1 + ", paste0(colnames(Z), collapse = '+')))
  qmf = model.frame(qformula, data= df)
  qfit = rq(qformula, tau = taus, data = df)
  gamma = coef(qfit)
  #print(gamma)
  u = df$y - predict(qfit)
  S = sapply(1:length(taus), function(i){
    1/sqrt(N) * sum(X1m_star[nonzero_idx] * score_func(taus[i],u[,i]))
  })
  ## indep
  Q_indep = sapply(1:length(taus), function(i){
    tau = taus[i]
    1/N * sum(X1m_star[nonzero_idx]^2) * tau * (1- tau)
  })
  
  ## combine multiple quantiles
  width = length(taus)
  V0 = matrix(0, ncol=width, nrow=width)
  for (kk in 1:(width-1)){
    for (ll in (kk+1):width){
      V0[kk, ll] = min(taus[kk], taus[ll]) - taus[kk]*taus[ll]
    }
  }
  V0 = V0 + t(V0) + diag(taus*(1 - taus))
  V0 = V0 * sum(X1m_star[nonzero_idx]^2) / N
  
  Q_dep_result = lapply(1:length(taus), function(t){
    tau = taus[t]
    diag_ele = tau*(1-tau)* sum(X1m_star[nonzero_idx]^2)
    
    off_diag_ele_result = lapply(unique(df$study), function(i){
      idx = df$idx[which(df$study == i)]
      #eg = expand.grid(j1 = idx, j2 = idx) %>%
      #  as.data.frame() %>%
      #  filter(j1 != j2)
      
      Xlist = lapply(idx, function(id) X1m_star[id,])
      Xsum = Reduce('+', Xlist)
      t1 = Xsum %*% t(Xsum)  - t(X1m_star[idx,]) %*% X1m_star[idx,]
      
      #tmp = lapply(1:nrow(eg), function(r){
      #  j1 = eg$j1[r]
      #  j2 = eg$j2[r]
      #  X_tilde[j1, ] %*% t(X_tilde[j2, ]) * (u[j1, t] < 0) * (u[j2, t] < 0)
      #}) %>% Reduce('+', .)
      
      
      Xlist_neg_res = lapply(idx, function(id) X1m_star[id,] * (u[which(df$idx == id),t] <0))
      X1m_neg_res = do.call(rbind, Xlist_neg_res)
      Xsum_neg_res = Reduce('+', Xlist_neg_res)
      t2 = Xsum_neg_res %*% t(Xsum) 
      t2m = lapply(1:length(Xlist), function(l){
        Xlist_neg_res[[l]] %*% t(Xlist[[l]])
      }) %>% Reduce('+', .)
      t2 = t2 - t2m
      
      t3 = Xsum_neg_res %*% t(Xsum_neg_res)  - t(X1m_neg_res) %*% X1m_neg_res
      
      ans = list(t1 = t1,
                 t2 = t2,
                 t = tau^2 * t1 - 2 * tau* t2 + t3,
                 Xlist_neg_res = Xlist_neg_res
                 )
      return(ans)
    }) 
    
    off_diag_ele = lapply(off_diag_ele_result, function(i) i$t) %>% Reduce('+', .)
    t1 = lapply(off_diag_ele_result, function(i) i$t1) %>% Reduce('+', .)
    t2 = lapply(off_diag_ele_result, function(i) i$t2) %>% Reduce('+', .)
    
    Xlist_neg_res = lapply(off_diag_ele_result, function(i) i$Xlist_neg_res)

    Q_dep = as.numeric((diag_ele + off_diag_ele)/N)
    
    ans = list(t1 = t1,
               t2 = t2,
               Q_dep = Q_dep,
               Xlist_neg_res = Xlist_neg_res)
  })
  
  Q_dep = sapply(Q_dep_result, function(q) q$Q_dep)
  Xlist_neg_res = lapply(Q_dep_result, function(q) q$Xlist_neg_res)
  t1 = sapply(Q_dep_result, function(q) q$t1)
  t2 = sapply(Q_dep_result, function(q) q$t2)
  
  
  ## combine multiple quantiles
  width = length(taus)
  V = matrix(0, width, width)
  for(i in 1:(width-1)){
    for(j in (i+1):width){
      tau_i = taus[i]
      tau_j = taus[j]
      d0 = (min(tau_i, tau_j) - tau_i * tau_j) * sum(X1m_star[nonzero_idx]^2)
      d1 = tau_i * tau_j * t1[1] 
      d2 = -tau_j * t2[i]-tau_i * t2[j]
      li = Xlist_neg_res[[i]]
      lj = Xlist_neg_res[[j]]
      d3 = lapply(1:length(li), function(k){
        tmp_i = Reduce('+', li[[k]])
        tmp_j = Reduce('+', lj[[k]])
        d3 = tmp_i %*% t(tmp_j)  - t(do.call(rbind, li[[k]])) %*% do.call(rbind, lj[[k]])
      }) %>% Reduce('+', .) 
      V[i, j] = (d0 + d1 + d2 + d3)/N
    }
  }
  V = V + t(V) + diag(Q_dep)
  
  
  Ttau_indep = S^2/Q_indep
  pval_indep = sapply(Ttau_indep, function(tt){
    pchisq(tt, df = 1, lower.tail = F)
  })
  names(pval_indep)= paste0("indep:", taus)
  
  Ttau_dep = S^2/Q_dep
  pval_dep = sapply(Ttau_dep, function(tt){
    pchisq(tt, df = 1, lower.tail = F)
  })
  names(pval_dep)= paste0("dep:", taus)
  
  ## minP combine
  pval_indep_minP = minP_combine_test(pval = c(logit_pval_LRT, pval_indep),
                                      taus = taus,
                                      Sigma.hat = V0,
                                      M = 10000)
  names(pval_indep_minP) = 'indep:minP'
  
  pval_indep_minP_clip = minP_combine_test(pval = c(logit_pval_LRT, pval_indep[c(2:4)]),
                                      taus = taus[c(2:4)],
                                      Sigma.hat = V0[2:4,2:4],
                                      M = 10000)
  names(pval_indep_minP_clip) = 'indep:minP_clip'
  
  
  pval_dep_minP = minP_combine_test(pval = c(logit_pval_LRT, pval_dep),
                                    taus = taus,
                                    Sigma.hat = V,
                                    M = 10000)
  names(pval_dep_minP) = 'dep:minP'
  
  pval_dep_minP_clip = minP_combine_test(pval = c(logit_pval_LRT, pval_dep[c(2:4)]),
                                           taus = taus[c(2:4)],
                                           Sigma.hat = V[2:4,2:4],
                                           M = 10000)
  names(pval_dep_minP_clip) = 'dep:minP_clip'
  
  ## combine
  zerorate = ifelse(is.na(logit_pval_LRT), 0, mean(ybin == 0))
  pval_indep_comb = c(logit_pval_LRT, pval_indep)
  pval_dep_comb = c(logit_pval_LRT, pval_dep)
  
  methods = c('unequal', 'unequal_clip', 'equal', 'equal_clip')
  pval_indep_list = sapply(methods, cauchy_combine_test, pval = pval_indep_comb,
                           zerorate = zerorate, taus= taus)
  names(pval_indep_list) = paste0("indep:cauchy_", methods)

  pval_dep_list = sapply(methods, cauchy_combine_test, pval = pval_dep_comb,
                           zerorate = zerorate, taus= taus)
  names(pval_dep_list) = paste0("dep:cauchy_", methods)
  
  
  res = c(pval_indep, pval_dep, 
          logit_pval, logit_pval_LRT,
          pval_indep_list, pval_dep_list,
          pval_indep_minP, pval_indep_minP_clip,
          pval_dep_minP, pval_dep_minP_clip
  )
  
  return(res)
}



