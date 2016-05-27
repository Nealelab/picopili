#####################
# Logistic Mixed Model, using GMMAT
# - logistic link function
# - GRM + FID clusters (block diagonal) as variance components
# - arbitrary covariates
# - see GMMAT package for references
#
# by Raymond Walters
# May 2016
#
# Note: currently no check for invariant covariates 
#
#####################

testing <- FALSE

require("GMMAT",lib="~/Rlibs")

# need:
# plink file stem
# grm file name (assumed square, with same FID/IIDs as bfile (in order)
# covar file name (file must have header, with FID and IID)
# outname

# build:
# fam_clust
# list of SNPs
# order of IDs
# output file name

# parse arguments
args = commandArgs(TRUE)
bfile = as.character(args[1])
grmfile = as.character(args[2])
covfile = as.character(args[3])
outstem = as.character(args[4])

if(testing){
	print(date())
	print(args)
}

# load files
GRM <- as.matrix(read.table(grmfile))
fam <- read.table(paste(bfile,".fam",sep=""),stringsAsFactors=FALSE,header=FALSE)
bim <- read.table(paste(bfile,".bim",sep=""),stringsAsFactors=FALSE,header=FALSE)
snplist <- bim[,2]
incov <- as.data.frame(read.table(covfile,header=TRUE))

if(testing){
	print(date())
}

# verify GRM is square, same dim as fam
if(ncol(GRM) != nrow(GRM)){
	stop("GRM matrix must be square. See e.g. 'plink --make-rel square gz', rather than '--make-grm-gz'.")
}
if(nrow(GRM) != nrow(fam)){
	stop("GRM file should match plink fam file (e.g. number, order of individuals)")
}

# match ID order
fam_id <- paste(fam[,1], ":", fam[,2], sep="")
cov_id <- paste(incov$FID, ":", incov$IID, sep="")
idx <- match(fam_id, cov_id)
idna <- is.na(idx)

# build covariate/phenotype matrix
dat <- as.data.frame(matrix(NA,nrow(fam),ncol(incov)-1))
colnames(dat) <- c("Y",colnames(incov[,-c(1,2)]))
dat$Y <- fam[,6]
dat$Y[dat$Y == -9] <- NA
dat$Y <- dat$Y - 1 # switch from 1/2 to 0/1 coding 
dat[!idna,-1] <- incov[idx[!idna],-c(1,2)]

if(testing){
	print(date())
	print(dim(dat))
	print(head(dat))
}

# build fam_clust
# block diagonal matrix, indicating IDs in same FID
fam_clust <- matrix(NA,nrow(fam),nrow(fam))
for(i in 1:nrow(fam)){
	for(j in 1:i){
		if(fam[i,1]==fam[j,1]){
			fam_clust[i,j] <- 1
		}else{
			fam_clust[i,j] <- 0
		}
	}
}
fam_clust[upper.tri(fam_clust)] <- t(fam_clust)[upper.tri(fam_clust)]
diag(fam_clust) <- 1

if(testing){
	print(date())
}

# fit model
f0 <- glmm.wald(fixed = Y~.,
		data=dat,
		kins=list(GRM, fam_clust),
		family=binomial(link="logit"),
		infile=bfile,
		snps=snplist,
		missing.method="omit")

if(testing){
	print(date())
	print(head(f0))
	print(summary(f0$converged))
}

# return NA for non-converged results
f0[!f0$converged | is.na(f0$converged),c("BETA","SE","PVAL")] <- NA

# format output
# flip beta since gmmat returns results for A2 instead of A1
out <- data.frame(CHR=f0$CHR, 
		  SNP=f0$SNP, 
		  POS=f0$POS, 
		  A1=f0$A1, 
		  BETA=-1*f0$BETA, 
		  SE=f0$SE, 
		  CHISQ=(f0$BETA/f0$SE)^2, 
		  P=f0$PVAL, 
		  N=f0$N, 
		  M_COV=ncol(dat)-1)

# return
outfile = paste("gmmatfam.",outstem,".R.txt",sep="")
write.table(format(out,digits=8), file=outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep = ' ')

if(testing){
	print(date())
}

# eof
