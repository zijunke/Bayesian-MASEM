#install.packages('Ryacas')
library(Ryacas)
# Specify the B matrix
# B: the matrix of regression coefficients
B = matrix(0,4,4)
B[2,1] = 'a1'
B[3,1] = 'a2'
B[4,1] = 'cp'
B[4,2] = 'b1'
B[4,3] = 'b2'
B = ysym(B)
# Specify the Ve matrix
# Ve: the residual covariance matrix
Ve = matrix(0,4,4)
Ve[2,3] = Ve[3,2] = 'g'
Ve[1,1] = 1
Ve[2,2] = 's1m1'
Ve[3,3] = 's2m2'
Ve[4,4] = 's2m3'
Ve = ysym(Ve)

# Compute the model-implied correlation matrix
# P: model-implied correlation matrix
I_B = diag(1,4)-B
inv.I_B = solve(I_B)
P = inv.I_B %*% Ve %*% t(inv.I_B)
P