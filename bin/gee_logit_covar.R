#####################
# GEE for PLINK family data
# - logistic link function
# - exchangeable errors (= CE model)
# - robust (sandwich) SEs
# - arbitrary covariates
#
# Adapted from scripts by Camelia Minica (cameliaminica.nl; doi:10.1038/ejhg.2014.94)
# by Raymond Walters
# December 2015
#
#
# Note: currently no check for invariant covariates 
#
#####################

require(geepack)

Rplink <- function(PHENO, GENO, CLUSTER, COVAR){
	
	if(max(PHENO)==2){
		y <- PHENO-1
	}else{
		y <- PHENO
	}
	m <- ncol(COVAR)
	
	# index to exclude IDs missing pheno or covariate
	excl <- is.na(PHENO) | apply(COVAR,1,function(a) any(is.na(a)) )
	
	# verify cluster is sorted
	if( length(levels(as.factor(CLUSTER))) !=  length(rle(CLUSTER)$values) ){
		ord <- order(CLUSTER)
	}else{
		ord <- 1:length(CLUSTER)
	}
	
	f1 <- function(snp){

		# init
		beta <- NA
		se <- NA
		chi <- NA
		p <- NA

		# index exclude missing SNP, pheno, or covariate 
		idx <- !is.na(snp) & !excl
		n <- sum(idx)
		
		# combine sorting plus exclusions
		xord <- ord[idx[ord]]

		# return NAs if SNP or pheno is invariant, or if empty cells
		facx <- as.factor(snp[xord])
		facy <- as.factor(y[xord])
		if(any(table(facx, facy) == 0) || length(levels(facx))==1 || length(levels(facy))==1){
			beta <- NA
			se <- NA
			chi <- NA
			p <- NA
		}else{
		# fit GEE		
			fit <- geeglm( y[xord] ~ snp[xord] + COVAR[xord,], 
						id=CLUSTER[xord], 
						family=binomial(link="logit"), 
						corstr="exchangeable",
						control=geese.control(maxit=100))
						# maxiter=100,
						# na.action=na.omit)
			summ <- summary(fit)
		
			beta <- coef(summ)$"Estimate"[2]
			se <- coef(summ)$"Std.err"[2]
			chi <- coef(summ)$"Wald"[2]
			# p <- 2*pnorm(abs(z),lower=FALSE)
			p <- coef(summ)$"Pr(>|W|)"[2]
			# n <- summ$nobs (from gee rather than geepack)
		}
		
		# TODO: potentially useful later:
		# summ$working.correlation in gee, or
		# summ$corr in geepack
		
		out <- c(beta, se, chi, p, n, m)
		 
		return(c( 6, out ))
	}
	
	apply(GENO, 2, f1)
}

