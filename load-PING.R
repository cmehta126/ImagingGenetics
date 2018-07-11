# Chintan Mehta
# April 1, 2016
# PING analysis - LOAD DATA
# chintan.mehta@yale.edu

require(snpStats)
require(Amelia)
require(car)
require(foreach)
require(doParallel)
source("functions.R")


# Define File Names for PCA, IBD, genotype, and phenotype files -----------
fn.PCA.eigenval = "PING_PCA.eigenval"
fn.PCA.eigenvec = "PING_PCA.eigenvec"
fn.IBD = "PING_IBD.genome"
fn.phenotype = "PING-FULL.csv"
fn.genotype = 'PING_660_final_filtered'
bed = paste(fn.genotype,'.bed',sep='')
bim = paste(fn.genotype,'.bim',sep='')
fam = paste(fn.genotype,'.fam',sep='')

# Load files --------------------------------------------------------------
G0 = read.plink(bed,bim,fam)
D0 = read.csv(file = fn.phenotype, header = T)
IBD = read.table(file = fn.IBD,header = T)
PCA.eigenvec = read.table(file = fn.PCA.eigenvec)

# Get Subject IDs having genotype data. -----------------------------------
mean(rownames(G0$fam) == rownames(G0$genotypes))
rownames(G0$genotypes) <- rownames(G0$fam) <- G0$fam$member
rownames(D0) = D0$SubjID
KEEP.SUBJ = intersect(rownames(D0), rownames(G0$genotypes))
D1 = D0[KEEP.SUBJ,]

# Make Working Dataframe --------------------------------------------------------
Z1 = data.frame(row.names = rownames(D1))

# Set Age, Sex, Hand, Manufacturer, Device -----------------------------------------------------------
Z1$Age = apply(D1[, c("Age_At_NPExam", "Age_At_IMGExam")],1,mean,na.rm = T)
Z1$Age[which(is.na(Z1$Age))] = -1 # Set missing age as -1; subjects will be removed later on
Z1$Sex = (D1$Gender == "M")
Z1$Hand = (D1$FDH_23_Handedness_Prtcpnt == "Right")
Z1$Manufacturer = D1$Manufacturer
Z1$DeviceSerialNumber = D1$DeviceSerialNumber
Z1$LP = (D1$FDH_31_Ever_Diag_Learn_Prob == "Yes")*1 # Ever diagnosed with a learning problem?

# Genotype PC
rownames(PCA.eigenvec) = PCA.eigenvec$V2
Z1[,c("PC1", "PC2")] = PCA.eigenvec[rownames(Z1),c("V3","V4")] 

# SES variables & Impute missing SES --------------------------------------
SES.VAR = c("FDH_3_Household_Income", "FDH_Highest_Education","FDH_Highest_Occupation")
SES = D1[,SES.VAR];
colnames(SES) = c("HI", "PE", "OC")
SES = impute.SES(SES)
Z1[,c("HI", "PE")] = SES[,c("HI", "PE")]



#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
# Steps for Sample INCLUSION ----------------------------------------------
## 1. Identify European Genetic Ancestry (EGA)
## 2. Identify siblings (through IBD) to remove
## 3. Age restriction
## 4. Availability of all cognitive data
## 5. Availability of all MRI data


# Define Sample with European genetic ancestry ----------------------------
# Determine extreme quantiles of PC1 & PC2 of subjects with self-reported white ancestry.
# Threshold on those extreme quantiles

PC1.White = Z1$PC1[which((D1$FDH_2_White=="Yes")==T)]
PC2.White = Z1$PC2[which((D1$FDH_2_White=="Yes")==T)]

quantile(PC1.White,probs=c(0,0.01,0.99,1))
# >          0%          1%         99%        100% 
# > -0.02227200 -0.02163948  0.03600463  0.07933670 

quantile(PC2.White,probs=c(0,0.01,0.99,1))
# >          0%          1%         99%        100% 
# > -0.08547510 -0.05101131  0.01787369  0.01983160 

# Without loss of generality, designate subject as having EGA if:
# PC1 < 0.04
# PC2 > -0.05

ID.EGA.PC1 = rownames(Z1)[which(Z1$PC1 < 0.04)]
ID.EGA.PC2 = rownames(Z1)[which(Z1$PC2 > -0.05)]
ID.EGA = intersect(ID.EGA.PC1, ID.EGA.PC2)

# Verify through plots:
xlim = range(Z1$PC1); ylim = range(Z1$PC2)
par(mfrow=c(1,2), xpd = T, mar = c(5,4,5,4))
plot(xlim, ylim, type = 'n', xlab = "PC1", ylab = "PC2")
points(Z1$PC1[which((D1$FDH_2_White=="Yes")==T)], Z1$PC2[which((D1$FDH_2_White=="Yes")==T)], col='blue', pch = 16)
points(Z1$PC1[which((D1$FDH_2_White=="Yes")==F)], Z1$PC2[which((D1$FDH_2_White=="Yes")==F)], col='red', pch = 17)
mtext(side = 3,"First two PCs by Self-reported White\nBlue = Self-reported white\nRed=Everyone else", line = 1)

plot(xlim, ylim, type = 'n', xlab = "PC1", ylab = "PC2")
points(Z1[ID.EGA,c("PC1","PC2")], col='blue', pch = 16)
points(Z1[-which(rownames(Z1) %in% ID.EGA),c("PC1","PC2")], col='red', pch = 17)
mtext(side = 3,"First two PCs by EGA\nBlue = EGA\nRed=Everyone else", line = 1)
par(mfrow=c(1,1),xpd=F)



# Find Siblings by IBD for removing ------------------------------------------------
max.IBD = 0.20
IBD.RESTRICT = IBD[which(IBD$PI_HAT > max.IBD),]
IBD.RESTRICT$IID1 = as.character(IBD.RESTRICT$IID1)
IBD.RESTRICT$IID2 = as.character(IBD.RESTRICT$IID2)


SIBLINGS = data.frame(row.names = rownames(Z1), 
                      Age = Z1$Age,
                      KEEP = rep(T,nrow(Z1)))

for(j in 1:nrow(IBD.RESTRICT)){
  IIDS = unlist(IBD.RESTRICT[j, c("IID1", "IID2")])
  if(all(IIDS %in% rownames(SIBLINGS))){
    Ages = SIBLINGS[IIDS,"Age"]
    if(Ages[1] > Ages[2]){
      SIBLINGS[IIDS[2], "KEEP"] = F
    }else{
      SIBLINGS[IIDS[1], "KEEP"] = F
    }  
  }
}

ID.SIBLING.KEEP = rownames(SIBLINGS)[which(SIBLINGS$KEEP==T)]

# Include by EGA, nonmissing data (COG + IMG), siblings, and age ------------------------------------
min.age = 8

ID.AGE = rownames(Z1)[which(Z1$Age >= min.age)]
ID.IMG = rownames(Z1)[which(!is.na(D1[,9]))]
ID.COG = rownames(na.omit(D1[,1635:1642]))


ID = intersect(ID.EGA, ID.SIBLING.KEEP)
ID = intersect(ID, ID.AGE)
ID = intersect(ID, ID.IMG)
ID = intersect(ID, ID.COG)

# Working Data Frames -----------------------------------------------------------
D11 = D1[ID,]
Z11 = Z1[ID,]

