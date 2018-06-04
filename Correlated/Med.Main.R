library(matrixcalc)
library(MASS)
library(coda)
library(Matrix)
library(R2OpenBUGS)

##-------------------------------------------------------
# Mediation: Correct Model
#--------------------------------------------------------
missing = 'MCAR'
method = 'correlation'
wd = 'D:/Research/180307A/MASEM/Example/'
source(paste(wd,'RCode/RealData.R',sep=''))
work.d = paste(wd,'Results/Med/',sep='')

# prepare data for OpenBUGS
#-------------------------------------------------
data<-list("Nstudy","N","Ninv","mu.N",'p',"pp","j","k",
	"vR","tau.vR","mu.vR.psi",'Z','ZV','WI',
	'alpha.prior.vE','beta.prior.vE')

fn = paste(wd,'Data/Med.dat',sep='')
da = read.table(fn,header = T)

Nstudy = nrow(da)
N = da[,2]
Ninv <- 1/N
mu.N = mean(N)
p = 4
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
Z <- matrix(0,pp,pp)
Z[1,1] = 1
Z[2,3] = 1
Z[3,1:4] = NA;Z[3,5] = 1
Z[4,c(1,3)] = NA;Z[4,6] = 1	
Z[5,c(1,3:6)] = NA;Z[5,2] = 1
Z[6,c(1:3,5:6)] = NA;Z[6,4] = 1
ZV <- cbind(matrix(NA,pp,2),matrix(0,pp,4))
df.prelim = 100*pp/mu.N+pp
alpha.prior.vE = (df.prelim-pp+1)/2
beta.prior.vE = alpha.prior.vE*(0.3/mu.N)
WI = diag(1,2)

# prepare initial values for OpenBUGS
#-------------------------------------------------
vR.inits = vR.impute
vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA

initsl <- list(list(theta0.raw = rep(0,2),a2 = 0,b2 = 0,cp = 0,g = 0,
	iQ = diag(100,2),xi = rep(0.1,2),tau.R = 100,
	vR.psi = matrix(0,Nstudy,pp),vR = vR.inits,
	vR.rep = matrix(0,Nstudy,pp),thetai.raw = matrix(0,Nstudy,2)))

# parameters to save; model file name
#-------------------------------------------------
prm = c('a10','b10','a2','b2','cp','g',
		'sd.a1','sd.b1','rho.a1b1','ppp')
model.fn = paste(wd,'Model/Med.txt',sep='')

fit = bugs(data,initsl,prm,model.fn,
	n.chains=1,n.iter=200000,DIC = TRUE,
	n.burnin=100000,n.thin = 1,debug = TRUE,
	saveExec = TRUE,working.directory = work.d)

setwd(work.d)
for(tryi in 1:20){
	#fit.coda <- as.mcmc.list(fit)
	fit.coda = read.openbugs(stem="",thin = 10)
	del.id = na.omit(match(c('ppp'),varnames(fit.coda)))
	tmp.conv = geweke.diag(fit.coda[,-del.id])[[1]]$z
	print(tmp.conv)
	print(summary(fit.coda),3)	
	if(sum((abs(tmp.conv)>1.96),na.rm = TRUE)==0){
		break
	}else{	
		fit = bugs(data,initsl,prm,model.fn,
			n.chains=1,n.iter=100001,DIC = TRUE,
			n.burnin=1,n.thin = 1,
			restart = TRUE,saveExec = TRUE,working.directory = work.d)
	}
}



##-------------------------------------------------------
# Mediation: Misspecified Model
# Sequential mediaiton
# The true model is parallel mediation
#--------------------------------------------------------
missing = 'MCAR'
method = 'correlation'
wd = 'D:/Research/180307A/MASEM/Example/'
source(paste(wd,'RCode/RealData.R',sep=''))
work.d = paste(wd,'Results/Med/',sep='')

# prepare data for OpenBUGS
#-------------------------------------------------
data<-list("Nstudy","N","Ninv","mu.N",'p',"pp","j","k",
	"vR","tau.vR","mu.vR.psi",'Z','WI',
	'alpha.prior.vE','beta.prior.vE')

fn = paste(wd,'Data/Med.dat',sep='')
da = read.table(fn,header = T)

Nstudy = nrow(da)
N = da[,2]
Ninv <- 1/N
mu.N = mean(N)
p = 4
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

Z <- matrix(0,pp,4)
Z[1,1] = 1
Z[2,1:2] = NA
Z[3,1:3] = NA;Z[3,4] = 1
Z[4,2] = 1	
Z[5,1:4] = NA
Z[6,c(1:2,4)] = NA;Z[6,3] = 1

df.prelim = 100*pp/mu.N+pp
alpha.prior.vE = (df.prelim-pp+1)/2
beta.prior.vE = alpha.prior.vE*(0.3/mu.N)
WI = diag(1,4)

# prepare initial values for OpenBUGS
#-------------------------------------------------
vR.inits = vR.impute
vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA

initsl <- list(list(theta0.raw = rep(0,4),
	iQ = diag(100,4),xi = rep(0.1,4),tau.R = 100,
	vR.psi = matrix(0,Nstudy,pp),vR = vR.inits,
	vR.rep = matrix(0,Nstudy,pp),thetai.raw = matrix(0,Nstudy,4)))

# parameters to save; model file name
#-------------------------------------------------
prm = c('a0','b10','b20','cp0','sd.par',
		'cor.par','ppp')
model.fn = paste(wd,'Model/MedW.txt',sep='')

fit = bugs(data,initsl,prm,model.fn,
	n.chains=1,n.iter=200000,DIC = TRUE,
	n.burnin=100000,n.thin = 1,debug = TRUE,
	saveExec = TRUE,working.directory = work.d)

for(tryi in 1:20){
	#fit.coda <- as.mcmc.list(fit)
	fit.coda = read.openbugs(stem="",thin = 10)
	del.id = na.omit(match(c('ppp'),varnames(fit.coda)))
	tmp.conv = geweke.diag(fit.coda[,-del.id])[[1]]$z
	print(tmp.conv)
	print(fit,3)	
	if(sum((abs(tmp.conv)>1.96),na.rm = TRUE)==0){
		break
	}else{	
		fit = bugs(data,initsl,prm,model.fn,
			n.chains=1,n.iter=100001,DIC = TRUE,
			n.burnin=1,n.thin = 1,
			restart = TRUE,saveExec = TRUE,working.directory = work.d)
	}
}




