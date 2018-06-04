library(matrixcalc)
library(MASS)
library(coda)
library(Matrix)
library(R2OpenBUGS)
library(ZIM)

#-------------------------------------------------------
# Mediation 
#-------------------------------------------------------
wd = 'D:/Research/180307A/MASEM/Example/'
source(paste(wd,'RCode/RealData.R',sep=''))

Nstudy = 50
mu.N = 200
p = 4
method = 'correlation'
missing = 'MCAR'
miss.rate.s = .8
miss.rate.v = .2

set.seed(176353920)

Theta = c(0.1,0.1,0.3,-0.3,0.1,-0.1)
Vtheta = diag(0,6)
Vtheta[1:2,1] = c(0.3^2,0)
Vtheta[1:2,2] = c(0,0.3^2) 

Bfree = matrix(0,p,p)
Bfree[2:4,1] = c(1,3,5)
Bfree[4,2:3] = c(2,4) 
Bfixed = matrix(0,p,p)
Bfixed[lower.tri(Bfixed)] = NA
Bfixed[3,2] = 0

Vefree = diag(0,p)
Vefree[3,2] = Vefree[2,3] = 6
Vefixed = diag(1,p)
Vefixed[3,2] = Vefixed[2,3] = NA
diag(Vefixed[2:4,2:4]) = NA
Free = list(B = Bfree,Ve = Vefree)
Fixed = list(B = Bfixed,Ve = Vefixed)

Pi.list = Get.Pi.Med(Nstudy,Free,Fixed,Theta,Vtheta,100)

N <- rzinb(n = Nstudy, k = 0.8, lambda = mu.N*0.2, omega = 0)
N <- N + mu.N*0.8		# sample size for each study

pp 	= p*(p-1)/2
vi2jk	= Get.vi2jk(p)
j	= vi2jk[,1]
k	= vi2jk[,2]
vil	= Get.jk2vi(vi2jk,p)
vR = array(NA,dim = c(Nstudy,pp))
for(studyi in 1:Nstudy){
	Ri = cor(mvrnorm(N[studyi],rep(0,p),Pi.list[,,studyi]))
	vR[studyi,] = Ri[vi2jk[,1:2]]	
} 
vR = Make.Missing(vR,method,missing,miss.rate.s,miss.rate.v,p,j,k,vil)	# generate missing values

da = cbind(1:Nstudy,N,vR)
colnames(da) = c('id','N',paste('vR',1:pp,sep=''))

wd = 'D:/Research/180307A/MASEM/Example/Data/'
fn = paste(wd,'Med.dat',sep='')
write.table(format(da,width = 8),fn,quote = F,col.names = T,row.names = F)


#-------------------------------------------------------
# Bivariate CFA
#-------------------------------------------------------
wd = 'D:/Research/180307A/MASEM/Example/'
source(paste(wd,'RCode/RealData.R',sep=''))

Nstudy = 50
mu.N = 200
p = 6
f = 2
method = 'correlation'
missing = 'MCAR'
miss.rate.s = .5
miss.rate.v = .3

set.seed(176353920)

Theta = c(.7,.6,.5,.7,.6,.5,.3)
Vtheta = diag(c(.01,.02,.04,.01,.02,.04,.05))

Lfree = matrix(0,p,f)
Lfree[1:3,1] = 1:3
Lfree[4:6,2] = 4:6 
Lfixed = matrix(0,p,f)
Lfixed[1:3,1] = NA
Lfixed[4:6,2] = NA 

Phifree = diag(0,f)
Phifree[2,1] = Phifree[1,2] = 7
Phifixed = diag(1,f)
Phifixed[2,1] = Phifixed[1,2] = NA
Free = list(Lambda = Lfree,Phi = Phifree)
Fixed = list(Lambda = Lfixed,Phi = Phifixed)

Pi.list = Get.Pi.CFA(Nstudy,Free,Fixed,Theta,Vtheta,100)

N <- rzinb(n = Nstudy, k = 0.8, lambda = mu.N*0.2, omega = 0)
N <- N + mu.N*0.8		# sample size for each study

pp 	= p*(p-1)/2
vi2jk	= Get.vi2jk(p)
j	= vi2jk[,1]
k	= vi2jk[,2]
vil	= Get.jk2vi(vi2jk,p)
vR = array(NA,dim = c(Nstudy,pp))
for(studyi in 1:Nstudy){
	Ri = cor(mvrnorm(N[studyi],rep(0,p),Pi.list[,,studyi]))
	vR[studyi,] = Ri[vi2jk[,1:2]]	
} 
vR = Make.Missing(vR,method,missing,miss.rate.s,miss.rate.v,p,j,k,vil)	# generate missing values

da = cbind(1:Nstudy,N,vR)
colnames(da) = c('id','N',paste('vR',1:pp,sep=''))

wd = 'D:/Research/180307A/MASEM/Example/Data/'
fn = paste(wd,'CFA.dat',sep='')
write.table(format(da,width = 8),fn,quote = F,col.names = T,row.names = F)
