# return positive integers in matrix x
# 1st column, values
# 2nd & 3rd columns, array index
find.pos.integer <- function(x) {
    logM = (x == round(x))*(x>0)
	tmp = which(logM==TRUE,arr.ind = TRUE)
	return(cbind(x[tmp],tmp))
}


# vector to matrix
v2m <- function(vec,jk,corr= T){
	p = max(jk)
	M = matrix(0,p,p)
	M[jk] = vec
	M = M + t(M)
	if(corr==TRUE){
		diag(M) = 1
	}else{
		diag(M) = diag(M)/2
	}
	return(M)
}

# get population correlation / covariance matrices for each study
#-----------------------------------------------------------------
Get.Pi.Med <- function(n.Pi,Free,Fixed,Theta,Vtheta,max.iter){

	par.names = names(Theta)
	SEM.par.names = c('B','Ve')
	B = Fixed[['B']]
	Ve = Fixed[['Ve']]
	p = nrow(B)

	# any random parameters?
	id = which(diag(Vtheta)>0)
	if(length(id)>0){
		Theta.r = Theta[id]
		Vtheta.r = Vtheta[id,id]
		if(sum(eigen(Vtheta.r)$values>0)<nrow(Vtheta.r)){
			print('The between-studies covariance matrix of SEM parameters is not positive definite')
			break
		}
		Theta.ri = mvrnorm(n.Pi*max.iter,Theta.r,Vtheta.r)
	}

	
	Pi.all = array(NA,dim = c(p,p,n.Pi))
	count = 0
	for(si in 1:n.Pi){
	for(j in 1:max.iter){
		count = count + 1
		Thetai = Theta
		if(length(id)>0){
			Thetai[id] = Theta.ri[count,]
		}
		
		idf = find.pos.integer(Free[['B']])
		B[idf[,-1]] = Thetai[idf[,1]]
		I_B = diag(1,p) - B
		inv.I_B = try(solve(I_B))

		idf = find.pos.integer(Free[['Ve']])
		Ve[idf[,-1]] = Thetai[idf[,1]]
		
		# only correct for the current mediation model
		Ve[2,2] = 1-B[2,1]^2
		Ve[3,3] = 1-B[3,1]^2
		Ve[4,4] = 1-inv.I_B[4,1]^2-B[4,2]^2*Ve[2,2]-B[4,3]^2*Ve[3,3]-2*B[4,2]*B[4,3]*Ve[3,2]
		
		Pi = inv.I_B%*%Ve%*%t(inv.I_B)
		
		crit1 = sum(eigen(Ve)$values>0)==nrow(Ve)
		crit2 = inherits(inv.I_B,'try-error')==FALSE
		crit3 = sum(diag(Pi)==1)==p
		crit4 = sum(eigen(Pi)$values>0)==nrow(Pi)

		if(crit1*crit2*crit3*crit4){break}
	}
		Pi.all[,,si] = Pi
	}
	return(Pi.all)
}



# get population correlation / covariance matrices for each study
#-----------------------------------------------------------------
Get.Pi.CFA <- function(n.Pi,Free,Fixed,Theta,Vtheta,max.iter){

	par.names = names(Theta)
	SEM.par.names = c('Lambda','Phi')
	L = Fixed[['Lambda']]
	Phi = Fixed[['Phi']]
	p = nrow(L)

	# any random parameters?
	id = which(diag(Vtheta)>0)
	if(length(id)>0){
		Theta.r = Theta[id]
		Vtheta.r = Vtheta[id,id]
		if(sum(eigen(Vtheta.r)$values>0)<nrow(Vtheta.r)){
			print('The between-studies covariance matrix of SEM parameters is not positive definite')
			break
		}
		Theta.ri = mvrnorm(n.Pi*max.iter,Theta.r,Vtheta.r)
	}
	
	Pi.all = array(NA,dim = c(p,p,n.Pi))
	count = 0
	for(si in 1:n.Pi){
	for(j in 1:max.iter){
		count = count + 1
		Thetai = Theta
		if(length(id)>0){
			Thetai[id] = Theta.ri[count,]
		}
		
		idf = find.pos.integer(Free[['Lambda']])
		L[idf[,-1]] = Thetai[idf[,1]]

		idf = find.pos.integer(Free[['Phi']])
		Phi[idf[,-1]] = Thetai[idf[,1]]
		
		Pi = L%*%Phi%*%t(L)
		diag(Pi) = 1
		
		crit1 = sum(eigen(Phi)$values>0)==nrow(Phi)
		crit2 = sum(eigen(Pi)$values>0)==nrow(Pi)
		crit3 = sum(abs(L)>1)==0

		if(crit1*crit2*crit3){break}
	}
		Pi.all[,,si] = Pi
	}
	return(Pi.all)
}

# conceal missing values
Make.Missing <- function(R,method,missing,miss.rate.s,miss.rate.v,p,j,k,vil){
	Nstudy = nrow(R)
	pp	= ncol(R)
	if(is.null(missing)){
		return(R)
	}else if(missing =='MCAR'){
		studyid	<- 1:Nstudy
		miss.id	<- studyid[rbinom(Nstudy,1,miss.rate.s)==1] 
		M	<- matrix(rbinom(p*Nstudy,1,miss.rate.v),Nstudy,p)
		for(studyi in miss.id){
			miss.pos = vil[M[studyi,]==1,]
			if(method == 'covariance'){
				R[studyi,miss.pos] = NA
			}else{
				R[studyi,miss.pos[-which(miss.pos== pp+1)]] = NA
			}
		}
		return(R)
	}
}



# impute missing values in covariance / correlation matrices of each study
# to obtain a rough estimate of the covariance matrix of covariance / correlation matrix
# weighted average correlation
Mimpute <- function(R,N,missing){
	if(is.null(missing)){
		return(R)
	}else{
		na.pos = which(is.na(R),arr.ind = TRUE)
		mu.N = mean(N)
		Rbar = apply(R,2,mean,na.rm = TRUE)# Becker's mean r
		
		for(coli in unique(na.pos[,2])){
			id = na.pos[(na.pos[,2] == coli),1]
			R[id,coli] = Rbar[coli]
		}
		return(R)
	}
}

# change the coordinating system of a vectorized matrix to the coordinating system of 
# the original matrix
# e.g., from vS to S, the former uses one coordinate (vil), whereas the latter uses two (j,k).
Get.vi2jk <- function(p,diag.incl=FALSE,byrow=FALSE){
	A = matrix(1,p,p)
	if(diag.incl ==FALSE){
		pp = p*(p-1)/2
		vi2jk <- matrix(NA,pp,3)
		vi2jk[,3] <- 1:pp
		if(byrow == FALSE){
			vi2jk[,1:2] <- which(lower.tri(A)==1,arr.ind = TRUE)
		}else{
			vi2jk[,1:2] <- which(upper.tri(A)==1,arr.ind = TRUE)
		}
		colnames(vi2jk) = c('j','k','vi')
	}else{
		pp = p*(p+1)/2
		vi2jk <- matrix(NA,pp,3)
		vi2jk[,3] <- 1:pp
		if(byrow == FALSE){
			vi2jk[,1:2] <- which(lower.tri(A,diag = TRUE)==1,arr.ind = TRUE)
		}else{
			vi2jk[,1:2] <- which(upper.tri(A,diag = TRUE)==1,arr.ind = TRUE)
		}
		colnames(vi2jk) = c('j','k','vi')		
	}
	return(vi2jk)
}

# change the coordinating system of a matrix to the coordinating system of 
# the corresponding vectorized matrix
# e.g., from S to vS, the former uses two coordinates (j,k), whereas the latter uses only one (vil).
Get.jk2vi <- function(vi2jk,p,diag.incl=FALSE){
	jk2vi = matrix(0,p,p)
	jk2vi[vi2jk[,1:2]] = vi2jk[,3]
	if(diag.incl){
		jk2vi = jk2vi + t(jk2vi)
		diag(jk2vi) = diag(jk2vi)/2
	}else{
		pp = p*(p-1)/2
		jk2vi = jk2vi + t(jk2vi) + diag(rep(pp+1,p))
	}
	return(jk2vi)
}

# compute the covariance matrix of correlation matrix
# based on Steiger (1980)
Corr.Cov <- function(vRN,j,k,vil){
	nvR	<- length(vRN)-1
	vR	<- c(vRN[1:nvR],1)
	N	<- vRN[nvR+1]
	NvR.cov <- matrix(NA,nvR,nvR)

	for(vi in 1:nvR){
		NvR.cov[vi,vi] <- (1-(vR[vi])^2)^2
	}
	for(vi in 1:(nvR-1)){
	for(vj in (vi+1):nvR){
		NvR.cov[vi,vj] <- ((vR[vil[j[vi],j[vj]]]-vR[vi]*vR[vil[k[vi],j[vj]]])*(vR[vil[k[vi],k[vj]]]-vR[vil[k[vi],j[vj]]]*vR[vj])
		 +(vR[vil[j[vi],k[vj]]]-vR[vil[j[vi],j[vj]]]*vR[vj])*(vR[vil[k[vi],j[vj]]]-vR[vi]*vR[vil[j[vi],j[vj]]])
		 +(vR[vil[j[vi],j[vj]]]-vR[vil[j[vi],k[vj]]]*vR[vj])*(vR[vil[k[vi],k[vj]]]-vR[vi]*vR[vil[j[vi],k[vj]]])
		 +(vR[vil[j[vi],k[vj]]]-vR[vi]*vR[vil[k[vi],k[vj]]])*(vR[vil[j[vj],k[vi]]]-vR[vil[k[vi],k[vj]]]*vR[vj]))/2
		NvR.cov[vj,vi] <- NvR.cov[vi,vj]
	}
	}

	vR.cov = NvR.cov/(N)
	vR.cov = as.matrix(nearPD(vR.cov,posd.tol = 1e-5)$mat)
	return(vR.cov)
}

Vj <- function(vR,vR.impute,vR.bar,N,method,pp,Nstudy,j,k,vil){

	mu.N = mean(N)
	# Vj are computed based on mean correlations
	if(method == 'corvariance'){
		S.vR.bar = Cov.Cov(c(vR.bar,mu.N),j,k,vil)
	}else if(method =='correlation'){
		S.vR.bar = Corr.Cov(c(vR.bar,mu.N),j,k,vil)
	}
	inv.S.vR.bar = solve(S.vR.bar)
	tau.vR = array(NA,dim = c(Nstudy,pp,pp))
	S.vR = array(NA,dim = c(Nstudy,pp,pp))
	for(i in 1:Nstudy){
		S.vR[i,,]<- S.vR.bar/N[i]*mu.N
		tau.vR[i,,] <- inv.S.vR.bar/mu.N*N[i]
	}	
	return(list(S.vR = S.vR,tau.vR = tau.vR))
}




# main function to run simulation
# InputL: input information for data analysis
#		p: number of variables per study
#		f:	number of factors
#		prm: character vector specifying parameters to be recorded
#		prm.dims: number of parameters per each parameter component
#		model.DGen: file specifying settings for data generation
#		model.fname: model file for OpenBUGS to run data analysis
#		out: file names for output
# method: covariance or correlation?
# missing: "MCAR"
# miss.rate.s, miss.rate.v: the missing rate for studies, the missing rate for variables within studies
Sim.R2OpenBugs <- function(InputL,method,missing = NULL,
	miss.rate.s = 0.3,miss.rate.v = 0.05,
	begin.sim = 1,stop.sim = 100,seed,
	print = TRUE){

	# unpack input information
	attach(InputL$sim)
	attach(InputL$model)
	attach(InputL$out)
	attach(readRDS(model.DGen))
	out.names = paste('M',1:3,sep='')
	
	pp 	= p*(p-1)/2
	vi2jk	= Get.vi2jk(p)
	j	= vi2jk[,1]
	k	= vi2jk[,2]
	vil	= Get.jk2vi(vi2jk,p)
	data<-list("Nstudy","N","Ninv","mu.N",'p',"pp","j","k",
			'ind1','ind2',
			"vR","tau.vR","mu.vR.psi",'Z',
			'alpha.prior.vE','beta.prior.vE')

	# define variables in the model file for OpenBUGS
	mu.vR.psi	= rep(0,pp)
	ind1 = (j>(p+1)/2)*(k<(p+1)/2)# Do these two indicators belong to the same factor ?
	ind2 = apply(vi2jk[,1:2],1,function(x) as.numeric(sum(x%in%c(1,2,4))==1))
	Z <- matrix(0,pp,p+1)
	for(vi in 1:pp){
		Z[vi,c(j[vi],k[vi])] = NA
	}	
	Z[,p+1] = NA
	V.par <- matrix(0,p+1,p+1)
	diag(V.par) <- NA 
	
	set.seed(seed)
	for(simj in 1:stop.sim){
		print(paste('simj = ',simj,sep=''))

		N <- rzinb(n = Nstudy, k = 0.8, lambda = NperStudy*0.2, omega = 0)
		N <- N + NperStudy*0.8		# sample size for each study
		Ninv <- 1/N
		mu.N = mean(N)
		vR = array(NA,dim = c(Nstudy,pp))
		Pi.list = Get.Pi(Nstudy,Free,Fixed,Values,max.iter,method)		# population covariance matrix for each study

		for(studyi in 1:Nstudy){
			Ri = cor(mvrnorm(N[studyi],rep(0,p),Pi.list[,,studyi]))
			vR[studyi,] = Ri[vi2jk[,1:2]]	
		} 

		vR = Make.Missing(vR,method,missing,miss.rate.s,miss.rate.v,p,j,k,vil)	# generate missing values
			
		df.prelim = 100*pp/mu.N+pp
		alpha.prior.vE = (df.prelim-pp+1)/2
		beta.prior.vE = alpha.prior.vE*(0.3/mu.N)
	
		vR.bar = apply(vR,2,mean,na.rm = TRUE)
		vR.impute = Mimpute(vR,N,missing)
		Stau.vR <- Vj(vR,vR.impute,vR.bar,N,method,'pooled',pp,Nstudy,j,k,vil)
		tau.vR <- Stau.vR$tau.vR
		S.vR <- Stau.vR$S.vR
		
		vR.inits = vR.impute
		vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA
		vR.psi.inits = matrix(0,Nstudy,pp)
		
		if(simj > (begin.sim -1)){
		for(mi in c(1,3)){
			
			np = length(T.Values.l[[3]])
			est <- T.Values.l[[3]]
			
			if((mi == 3)*(est[2]==0)){
				est[c(2,9:14)] = 0.001
			}
			inits.tmp = c(CFA.initials(est,inits.prm.l[[3]],prm.l[[3]],prm.dims.l[[3]]),
				list(vR = vR.inits,vR.rep = vR.impute,vR.psi = vR.psi.inits))
			inits <- list(inits.tmp)
			#print(data)
			#print(work.d)
				
			fit2 = bugs(data,inits,prm.l[[3]],model.fname[[mi]],
			n.chains=1,n.iter=n.iter[2],DIC = TRUE,OpenBUGS.pgm = openbugs.d,
			n.burnin=n.burnin[2],n.thin = 1,#debug = TRUE,
			saveExec = TRUE,working.directory = work.d)

			for(tryi in 1:10){
				fit2.coda <- as.mcmc.list(fit2)
				del.id = na.omit(match(c('ppp'),varnames(fit2.coda)))
				tmp.conv = geweke.diag(fit2.coda[,-del.id])[[1]]$z
				print(tmp.conv)
				print(fit2,3)	
				if(sum((abs(tmp.conv)>1.96),na.rm = TRUE)==0){
					break
				}else{	
					fit2 = bugs(data,inits,prm.l[[3]],model.fname[[mi]],
						n.chains=1,n.iter=n.iter[2]-n.burnin[2]+1000,DIC = TRUE,
						n.burnin=1000,n.thin = 1,OpenBUGS.pgm = openbugs.d,
						restart = TRUE,saveExec = TRUE,working.directory = work.d) 

				}
			}

			organize.results(simj,tryi,tmp.conv,prm.l[[3]],
				prm.dims.l[[3]],T.Values.l[[3]],fit2,out.names[mi],TRUE)
			
						

		}
		}

	}
	ob.names = ls()
	ob.names = ob.names[-which(ob.names=='fit2')]
	rm(list = ob.names )
}



