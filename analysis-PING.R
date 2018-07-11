# Chintan Mehta
# April 2, 2016
# PING analysis - PERFORM ANALYSIS
# chintan.mehta@yale.edu


# Remove covariates from Cognitive Data -----------------------
COG.Covariates = c("Age", "Sex", "HI", "PE", "PC1", "PC2")
COG = D11[ID,1635:1642]
COG = apply(D11[ID,1635:1642],2,remove.covariates,Z11[,COG.Covariates])

# Remove covariates from Cognitive and Imaging Data -----------------------
Whole.brain = T # Include whole brain measure as covariate?
IMG.Covariates = c("Age", "Sex", "Hand", "Manufacturer", "PC1", "PC2")

IMG = c()
IMG$roi.thick = D11[ID,c(8+1:66)]
IMG$roi.area = D11[ID,c(8+1:66) + 105]
IMG$roi.vol = D11[ID,c(8+1:66) + 105*2]
IMG$Whole = data.frame(row.names=rownames(D11),D11[,8 + 105*1:3])
colnames(IMG$Whole) = c("thick","area","volume")

V00 = data.frame(row.names = ID)
for(j in 1:3){
  W.IMG = Z11[ID, IMG.Covariates]
  if(Whole.brain){
    W.IMG$Whole = IMG$Whole[ID,j]
  }
  if(j==2){
    # For area, control for Household income & Parental education
    # This is shown for PING sample in Noble et al. (nature Neuroscience, 2015)
    W.IMG$HI = Z11[ID, "HI"]
    W.IMG$PE = Z11[ID, "PE"]
  }
  I0 = apply(IMG[[j]], 2, function(s,W0){scale(resid(lm(s~.,W0)))}, W0=W.IMG)
  V00[,colnames(IMG[[j]])] = I0
}


# Full Analysis -----------------------------------------------------------


# Prepare Cognitive Response ----------------------------------------------
COG.PRINCOMP = princomp(COG)
sign(COG.PRINCOMP$loadings[,1:2])
# Comp.1 Comp.2
# TBX_ibam_scr            -1     -1
# TBX_ls                  -1     -1
# TBX_pspac_scr           -1      1
# TBX_VOCAB_THETA         -1     -1
# TBX_reading_score       -1     -1
# TBX_flanker_score       -1      1
# TBX_attention_score     -1      1
# TBX_dccs_score          -1      1

# Average scores (controlled for covariates and scaled to zero mean & variance one) for ibam, list sort, vocab, reading
COG.LEARN = COG[,which(sign(COG.PRINCOMP$loadings[,2]) == sign(COG.PRINCOMP$loadings["TBX_reading_score",2]))]
Z11$Learning = apply(COG.LEARN, 1, mean)



# Getting Neuroimaging risk score -----------------------------------------
## Steps
### 1: Perform Joint variance test statistic for each ROI
### 2: Pick # of significant ROIs to keep for further analysis
### 3: Identify and remove neuroimaging outliers in diagnosed subgroup
### 4: Calculate NS-LP scores

# Get JVAR Test Statistic -------------------------------------------------------
# nperm = # of permutations (multiplied by # of cores) 
#         for getting significance of JVAR test. At least 10000 (more the better)
nperm = 1000; ncore = 10
fn.functions = "functions.R"; source(fn.functions)

# Indices for diagnosed & Undiagnosed subgroups
C1 = which(Z11$LP == 1) 
C2 = which(Z11$LP == 0)

JVAR = JVAR.stat(V00,C1,C2)
JVAR = JVAR.significance.perm(V00, JVAR, C1, C2, ncore=ncore,nperm=nperm, fn.functions)
head(JVAR,12)

# Set threshold for keeping ROIs -----------------------------
## Keep ROIs with p-values less than max.JVAR into further analysis.
### Hard-coding to keep only thickness measurements because they contribute most 
### consistently to significance of most important ROIs
max.JVAR = 0.075

KEEP.ROI.IND = JVAR$roi.index[which(JVAR$p.val < max.JVAR)]
KEEP.VARIABLES = colnames(V00)[as.numeric(KEEP.ROI.IND)]
KEEP.VARIABLES
# [1] "MRI_cort_thick.ctx.lh.temporalpole"             "MRI_cort_thick.ctx.lh.lingual"                 
# [3] "MRI_cort_thick.ctx.rh.middletemporal"           "MRI_cort_thick.ctx.lh.caudalanteriorcingulate" 
# [5] "MRI_cort_thick.ctx.rh.rostralanteriorcingulate" "MRI_cort_thick.ctx.lh.transversetemporal"      

V11 = V00[,KEEP.VARIABLES]

# Get pairwise distances  -------------------------------------------
DIST = as.matrix(dist(V11))

# Find neuroimaging outliers in diagnosed group (C1) ----------------------
outlying.threshold = 0.01

outlying.measure = apply(DIST[C1,],1,function(a,C1,C2){wilcox.test(a[C1],a[C2],alternative = 'less')$p.value},C1=C1,C2=C2)
outlying.C1 = C1[which(outlying.measure > outlying.threshold)]

if(length(outlying.C1) > 0){
  C1.new = C1[-which(outlying.measure > outlying.threshold)]
  C2.new = c(C2, outlying.C1)
}else{
  C1.new = C1;
  C2.new = C2;
}

# Get neuroimaging risk score ---------------------------------------------
get.NSLP.score = function(a,C1.new){
  a = rank(a)
  1/sum(a[C1.new])
}

Z11$NSLP = apply(DIST,1, get.NSLP.score, C1.new)
Z11$NSLP = Z11$NSLP/min(Z11$NSLP)
range(Z11$NSLP)


# Check: ------------------------------------------------------------------
plot(Z11$NSLP, Z11$Learning, main = "NSLP Model of Learning score")
h=summary(lm(Learning~NSLP, Z11))
h$coefficients


# Perform GWAS ------------------------------------------------------------
GWAS.Covariates = c("LP", "PC1", "PC2")
# Name of response variable is assumed to be "Y"
Z11$Y = Z11$Learning


## Set minimum minor allele frequency & HWE p.val
min.MAF = 0.05; min.HWE = 1e-4;

GENO = G0$genotypes[ID,]
GENO.QA = col.summary(GENO)
GENO.QA$HWE.pval = 1 - pchisq((GENO.QA$z.HWE)^2, df=1)
KEEP.SNPs = rownames(GENO.QA)[which(GENO.QA$MAF > min.MAF &  GENO.QA$HWE.pval > min.HWE)]
GENO.KEEPSNPS = GENO[,KEEP.SNPs]
GENO11 = matrix(as.numeric(GENO.KEEPSNPS), nrow = nrow(GENO.KEEPSNPS))
GENO11 = data.frame(GENO11, row.names = rownames(GENO.KEEPSNPS))
colnames(GENO11) = colnames(GENO.KEEPSNPS)

## Prepare Data for GWAS in parallel processing
N = ncol(GENO11)
nbatch = ceiling(N/ncore)

GL = c()
for(j in 1:ncore){
  i1 = 1 + nbatch*(j-1)
  i2 = min(nbatch*j,N)
  GL = c(GL, list(GENO11[,i1:i2]))
  print(i2)
}

require(foreach); require(doParallel)  
cl = makeCluster(ncore)
registerDoParallel(cl)

a = Sys.time()
x = foreach(Gj=iter(GL), .combine = "rbind") %dopar% {
  source("functions.R")
  require(car)
  
  e = c()
  for(j in 1:ncol(Gj)){
    result.j = test.SNP(Z11, Gj[,j], GWAS.Covariates)
    e = rbind(e, result.j)
  }
  qsnp = colnames(Gj)
  e = as.data.frame(e,row.names= qsnp);
  return(e)
}
a1 = Sys.time()-a; print(a1)
stopCluster(cl)



SNP.names = rownames(x)
gwas.results = data.frame(row.names = SNP.names,
                          CHR = G0$map[SNP.names,"chromosome"],
                          POS = G0$map[SNP.names,"position"])
gwas.results[,colnames(x)] = x

gwas.results = gwas.results[order(gwas.results$p.intx),]
head(gwas.results)

write.csv(JVAR, file = "jvar-results.csv")
write.csv(gwas.results, file = "gwas-results.csv")






















































