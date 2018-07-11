

# Functions ---------------------------------------------------------------


# Joint Variance test statistics ------------------------------------------
JVAR.stat = function(V00, C1, C2){
  roi.names = c('bankssts','caudalanteriorcingulate','caudalmiddlefrontal','cuneus',
                'entorhinal','fusiform',
                'inferiorparietal','inferiortemporal','isthmuscingulate',
                'lateraloccipital','lateralorbitofrontal','lingual',
                'medialorbitofrontal','middletemporal','parahippocampal',
                'paracentral','parsopercularis','parsorbitalis','parstriangularis',
                'pericalcarine','postcentral','posteriorcingulate','precentral',
                'precuneus','rostralanteriorcingulate','rostralmiddlefrontal',
                'superiorfrontal','superiorparietal','superiortemporal','supramarginal',
                'frontalpole','temporalpole','transversetemporal')
  roi.names = c(paste("L.",roi.names,sep=""), paste("R.",roi.names,sep=""))
  
    
  JVAR = data.frame(row.names = roi.names, roi.index = 1:length(roi.names))
  JVAR$thick = NA;
  JVAR$area = NA;
  JVAR$volume = NA;
  JVAR$stat = NA;
  
  for(j in 1:66){
    ROI = V00[,c(j,j+66,j+66*2)]
    p.vartest = apply(ROI,2,function(a,C1,C2){var.test(a[C1], a[C2],alternative = "less")$p.value}, C1 = C1, C2 = C2)
    JVAR[j,c("thick", "area", "volume")] = p.vartest; 
    JVAR$stat[j] = 1 - prod(1-p.vartest)
  }
  return(JVAR)
}

permute.JVAR.test = function(V00,C1,obs.stat){
  n = nrow(V00); n1 = length(C1);
  perm.IND = sample.int(n,n,replace=F)
  pC1 = perm.IND[1:n1]
  pC2 = perm.IND[(n1+1):n]
  perm.JVAR = JVAR.stat(V00, pC1, pC2)
  (obs.stat >= perm.JVAR$stat)
}


JVAR.significance.perm = function(V00, JVAR, C1, C2, ncore=ncore,nperm=1000, fn.functions){
  nC1 = length(C1); nC2 = length(C2)
  
  require(foreach); require(doParallel)  
  cl = makeCluster(ncore)
  registerDoParallel(cl)
  
  a = Sys.time()
  x = foreach(i=1:ncore,.combine = "cbind") %dopar% {
    source(fn.functions)
    e = c()
    for(j in 1:nperm){
      set.seed(i*ncore*nperm + j)
      e = cbind(e, permute.JVAR.test(V00,C1,JVAR$stat))  
    }
    return(e)
  }
  stopCluster(cl)
  print(Sys.time()-a)
  
  JVAR$p.val = apply(x,1,mean)
  JVAR = JVAR[order(JVAR$p.val),]
  JVAR$p.fdr = p.adjust(JVAR$p.val, method = 'fdr')
  return(JVAR)
}




# Remove covariates -------------------------------------------------------

remove.covariates = function(a,W){
  scale(resid(lm(a~.,W)))
}



# Testing single SNPs -----------------------------------------------------

test.SNP = function(Z11, SNP, GWAS.Covariates){
  require(car)
  L1 = Z11[,c(GWAS.Covariates,"NSLP","Y")]
  L1$SNP = SNP
  L1 = L1[which(L1$SNP != 0),] - -1
  
  # Main effect model
  L1.main = L1[,c("Y", "SNP", GWAS.Covariates)];
  fit = lm(Y~., L1.main)
  p.main = linearHypothesis(fit, "SNP", test = "Chisq")$Pr[2]
  
  # Interaction model
  L1.intx = L1[,c("Y", "SNP", "NSLP", GWAS.Covariates)];
  L1.intx$NSLP.SNP = L1.intx$NSLP * L1.intx$SNP
  fit = lm(Y~., L1.intx)
  p.intx = linearHypothesis(fit, c("SNP", "NSLP.SNP"), test = "Chisq")$Pr[2]
  
  return(c(p.main = p.main, p.intx = p.intx))
}



# Imputation for missing SES ----------------------------------------------
impute.SES  = function(SES){
  require(Amelia)
  priors.HI = matrix(NA, nrow=0, ncol = 5)
  priors.PE = matrix(NA, nrow=0, ncol = 5)
  for(i in seq(1,nrow(SES))){
    x1 = SES$HI[i]
    x2 = SES$PE[i]
    x3 = SES$OC[i]
    if( is.na(x1) && !is.na(x2) && !is.na(x3)){
      priors.HI = rbind(priors.HI, c(i, 1, 0, x2+x3, 0.99999))
    }
    if( !is.na(x1) && is.na(x2) && !is.na(x3)){
      priors.PE = rbind(priors.PE, c(i, 1, 0, x1+x3, 0.99999))
    }
  }
  make.impute = function(x,v,imp.SES){
    require(Amelia)
    a = x$imputations$imp1[rownames(imp.SES),]
    miss = which(is.na(imp.SES[,v]))
    imp.SES[miss,v] = a[miss,v]
    miss2 = which(is.na(imp.SES[,v]))
    imp.SES[miss2,v] = mean(a[,v], na.rm = T)
    return(imp.SES)
  }
  
  set.seed(1)
  x.HI = amelia(SES, m=1, bound = rbind(c(1,0,Inf), c(2, 0, Inf), c(3, 0, Inf)), pr = priors.HI)
  x.PE = amelia(SES, m=1, bound = rbind(c(1,0,Inf), c(2, 0, Inf), c(3, 0, Inf)), pr = priors.PE)
  
  imp.SES = SES
  imp.SES = make.impute(x.HI, "HI", imp.SES)
  imp.SES = make.impute(x.PE, "PE", imp.SES)
  
  return(imp.SES)
  
}

