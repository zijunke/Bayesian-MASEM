library(matrixcalc)
library(MASS)
library(coda)
library(Matrix)
library(R2OpenBUGS)

##-------------------------------------------------------
# 2-factor CFA: Correct Model
#--------------------------------------------------------
missing = 'MCAR'
method = 'correlation'
wd = 'D:/Research/180307A/MASEM/Example/'
source(paste(wd,'RCode/RealData.R',sep=''))
work.d = paste(wd,'Results/CFA/',sep='')

# prepare data for OpenBUGS
#-------------------------------------------------
data<-list("Nstudy","N","Ninv","mu.N",'p',"pp",
	"j","k",'vil','ind1',
	"vR","tau.vR","mu.vR.psi",'V.par','Z','WI',
	'alpha.prior.vE','beta.prior.vE')

fn = paste(wd,'Data/CFA.dat',sep='')
da = read.table(fn,header = T)

Nstudy = nrow(da)
N = da[,2]
Ninv <- 1/N
mu.N = mean(N)
p = 6
pp 	= p*(p-1)/2
vi2jk	= Get.vi2jk(p)
j	= vi2jk[,1]
k	= vi2jk[,2]
vil	= Get.jk2vi(vi2jk,p)
vR = as.matrix(da[,-c(1:2)])
vR.bar = apply(vR,2,mean,na.rm = TRUE)
vR.impute = Mimpute(vR,N,missing)
Stau.vR <- Vj(vR,vR.impute,vR.bar,N,method,pp,Nstudy,j,k,vil)

tau.vR <- Stau.vR$tau.vR
S.vR <- Stau.vR$S.vR
mu.vR.psi	= rep(0,pp)

ind1 = (j>(p+1)/2)*(k<(p+1)/2)
V.par = matrix(0,7,7)
diag(V.par) = NA
V.par[3,6] = V.par[6,3] = NA
Z <- matrix(0,pp,p+1)
for(vi in 1:pp){
	Z[vi,c(j[vi],k[vi])] = NA
}	
Z[,p+1] = NA

df.prelim = 100*pp/mu.N+pp
alpha.prior.vE = (df.prelim-pp+1)/2
beta.prior.vE = alpha.prior.vE*(0.3/mu.N)
WI = diag(1,2)

# prepare initial values for OpenBUGS
#-------------------------------------------------
vR.inits = vR.impute
vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA

initsl <- list(list(mu.rho=0,L0.raw=rep(.6,2),mu.L=c(.6,.6,NA,.6,.6,NA),
	sd.rho = 0.1,tau.Lraw = diag(100,2),sd.L=c(.1,.1,NA,.1,.1,NA),tau.R = 100,
	vR.psi = matrix(0,Nstudy,pp),vR = vR.inits,
	vR.rep = vR.impute,rho = rep(0,Nstudy),
	Lraw = matrix(0.6,Nstudy,2),xi = rep(1,2)))

# parameters to save; model file name
#-------------------------------------------------
prm = c('mu.rho','mu.L','sd.rho','sd.L','cor.L','ppp')
model.fn = paste(wd,'Model/CFA.txt',sep='')

fit = bugs(data,initsl,prm,model.fn,
	n.chains=1,n.iter=60000,DIC = TRUE,
	n.burnin=30000,n.thin = 1,debug = TRUE,
	saveExec = TRUE,working.directory = work.d)

setwd(work.d)
for(tryi in 1:20){
	#fit.coda <- as.mcmc.list(fit)
	fit.coda = read.openbugs(stem="",thin = 1)
	del.id = na.omit(match(c('ppp'),varnames(fit.coda)))
	tmp.conv = geweke.diag(fit.coda[,-del.id])[[1]]$z
	print(tmp.conv)
	print(summary(fit.coda),3)	
	if(sum((abs(tmp.conv)>1.96),na.rm = TRUE)==0){
		break
	}else{	
		fit = bugs(data,initsl,prm,model.fn,
			n.chains=1,n.iter=30001,DIC = TRUE,
			n.burnin=1,n.thin = 1,
			restart = TRUE,saveExec = TRUE,working.directory = work.d)
	}
}



##-------------------------------------------------------
# 2-factor CFA: Factor loading misspecified
#--------------------------------------------------------
missing = 'MCAR'
method = 'correlation'
wd = 'D:/Research/180307A/MASEM/Example/'
source(paste(wd,'RCode/RealData.R',sep=''))
work.d = paste(wd,'Results/CFA/',sep='')

# prepare data for OpenBUGS
#-------------------------------------------------
data<-list("Nstudy","N","Ninv","mu.N",'p',"pp","j","k",'ind2',
	"vR","tau.vR","mu.vR.psi",'Z','WI','V.par',
	'alpha.prior.vE','beta.prior.vE')

fn = paste(wd,'Data/CFA.dat',sep='')
da = read.table(fn,header = T)

Nstudy = nrow(da)
N = da[,2]
Ninv <- 1/N
mu.N = mean(N)
p = 6
pp 	= p*(p-1)/2
vi2jk	= Get.vi2jk(p)
j	= vi2jk[,1]
k	= vi2jk[,2]
vil	= Get.jk2vi(vi2jk,p)
vR = as.matrix(da[,-c(1:2)])
vR.bar = apply(vR,2,mean,na.rm = TRUE)
vR.impute = Mimpute(vR,N,missing)
Stau.vR <- Vj(vR,vR.impute,vR.bar,N,method,pp,Nstudy,j,k,vil)

tau.vR <- Stau.vR$tau.vR
S.vR <- Stau.vR$S.vR
mu.vR.psi	= rep(0,pp)

ind2 = apply(vi2jk[,1:2],1,function(x) as.numeric(sum(x%in%c(1,2,4))==1))
Z <- matrix(0,pp,p+1)
for(vi in 1:pp){
	Z[vi,c(j[vi],k[vi])] = NA
}	
Z[,p+1] = NA
df.prelim = 100*pp/mu.N+pp
alpha.prior.vE = (df.prelim-pp+1)/2
beta.prior.vE = alpha.prior.vE*(0.3/mu.N)

V.par = matrix(0,7,7)
diag(V.par) = NA
V.par[3,6] = V.par[6,3] = NA
WI = diag(1,2)

# prepare initial values for OpenBUGS
#-------------------------------------------------
vR.inits = vR.impute
vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA

initsl <- list(list(mu.rho=0,L0.raw=rep(.6,2),mu.L=c(.6,.6,NA,.6,.6,NA),
	sd.rho = 0.1,tau.Lraw = diag(100,2),sd.L=c(.1,.1,NA,.1,.1,NA),tau.R = 100,
	vR.psi = matrix(0,Nstudy,pp),vR = vR.inits,
	vR.rep = vR.impute,rho = rep(0,Nstudy),
	Lraw = matrix(0.6,Nstudy,2),xi = rep(1,2)))

# parameters to save; model file name
#-------------------------------------------------
prm = c('mu.rho','mu.L','sd.rho','sd.L','ppp')
model.fn = paste(wd,'Model/CFAW.txt',sep='')

fit = bugs(data,initsl,prm,model.fn,
	n.chains=1,n.iter=60000,DIC = TRUE,
	n.burnin=30000,n.thin = 1,debug = TRUE,
	saveExec = TRUE,working.directory = work.d)

setwd(work.d)
for(tryi in 1:20){
	fit.coda <- as.mcmc.list(fit)
	#fit.coda = read.openbugs(stem="",thin = 10)
	del.id = na.omit(match(c('ppp'),varnames(fit.coda)))
	tmp.conv = geweke.diag(fit.coda[,-del.id])[[1]]$z
	print(tmp.conv)
	print(fit,3)	
	if(sum((abs(tmp.conv)>1.96),na.rm = TRUE)==0){
		break
	}else{	
		fit = bugs(data,initsl,prm,model.fn,
			n.chains=1,n.iter=30001,DIC = TRUE,
			n.burnin=1,n.thin = 1,
			restart = TRUE,saveExec = TRUE,working.directory = work.d)
	}
}


