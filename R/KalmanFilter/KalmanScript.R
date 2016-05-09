# Kalman filters ----------------------------------------------------------

library(data.table)
library(ggplot2)
library(foreach)
library(iterators)

library(mnormt)
library(dplyr)
library(data.table)
library(ggplot2)

xhat <- c(0.5, -0.5)
P_0 <- matrix(c(1, 0.5, 0.5, 0.5), ncol=2)
y <- c(2.5, -2)

A <- matrix(c(1.5, 0,
              0, -0.3), ncol=2)

G = diag(2)

Q <- 0.5 * P_0 
R <- 0.8 * P_0

PMIKalmanStep <- function(x,P, A, G, Q, R, y)
{
  K1 <- P %*% t(G) %*% solve(G %*% P %*% t(G) + R)
  
  # Step 1
  x_f <- x + K1 %*% (y - G %*% x)
  P_f <- P - K1 %*% G %*% P
  
  # Step 2
  #K2 <- A %*% P %*% t(G) %*% solve(G%*% P %*% t(G) + R)
  K2 <- A %*% K1
  x_new <- A %*% x + K2 %*% (y - G %*% x)
  P_new <- A %*% P %*% t(A) - K2 %*% G %*% P %*% t(A) + Q
  return(list(x_f = x_f, P_f = P_f, K = K, x_new = x_new, P_new = P_new)) 
}

#KalmanFilter <- function(x0, P0, A, G, Q, R, ydata)
#{

# Step 1
xhatf <- xhat + P_0 %*% t(G) %*% solve(G %*% P_0 %*% t(G) + R) %*% (y - G %*% xhat)
Pf <- P_0 - P_0 %*% t(G) %*% solve(G %*% P_0 %*% t(G) + R) %*% G %*% P_0

# Step 2
K <- A %*% P_0 %*% t(G) %*% solve(G%*% P_0 %*% t(G) + R)
xhatnew <- A %*% xhat + K %*% (y - G %*% xhat)
Sigmanew <- A %*% P_0 %*% t(A) - K %*% G %*% P_0 %*% t(A) + Q


#model <- PMIKalmanStep(xhat,P_0,A,G,Q,R,y)
#model$x_new

x1 <- seq(-5,5,0.03)
x2 <- seq(-5,5,0.03)
grid <- data.table(expand.grid(x1,x2))
setnames(grid, c("x","y"))
z1 <- dmnorm(grid,xhat,varcov = P_0)
z2 <- dmnorm(grid,c(2.3, -1.9),varcov = R)
z3 <- dmnorm(grid,c(xhatf),varcov = Pf)
z4 <- dmnorm(grid,c(xhatnew),varcov = Sigmanew)

grid <- data.table(grid,z1,z2,z3,z4)

grid[,p1:=z1]
grid[,p2:=max(z1,z2), by=1:nrow(grid)]
grid[,p3:=max(z1,z2,z3), by=1:nrow(grid)]
grid[,p4:=max(z1,z2,z3,z4), by=1:nrow(grid)]

grid %>% ggplot(aes(x=x,y=y)) +geom_tile(aes(fill=p1))+scale_fill_distiller(palette="Spectral") +theme_light()

g1 <- grid %>% ggplot(aes(x=x,y=y)) +geom_tile(aes(fill=p1))+scale_fill_distiller(palette="Spectral") +theme_light()
g2 <- grid %>% ggplot(aes(x=x,y=y)) +geom_tile(aes(fill=p2))+scale_fill_distiller(palette="Spectral") +theme_light()
g3 <- grid %>% ggplot(aes(x=x,y=y)) +geom_tile(aes(fill=p3))+scale_fill_distiller(palette="Spectral") +theme_light()
g4 <- grid %>% ggplot(aes(x=x,y=y)) +geom_tile(aes(fill=p4))+scale_fill_distiller(palette="Spectral") +theme_light()


# model <- PMIKalmanStep(xhat,P_0,A,G,Q,R,y)
# model$x_new
# 
# ydata <- as.numeric(Nile)
# Size <- nrow(ydata)
# xdata <- rep(0,Size)
# 
# xt <- matrix
# 
# i <- 1
# 
# A <- diag(1)
# P <- diag(1)
# G <- diag(1)
# Q <- diag(1)
# R <- diag(1)
# 
# xt[,,]
# Pt[,,]
# filterted <- foreach(i=iter(1:Size),.combine = rbind) %do%
# {
#   y <- as.matrix(ydata[i])
#   x <- ifelse((i-1)==0,0, as.matrix(xdata[i-1]))
#   P <- ifelse((i-1)==0,P,Pt[,,i-1])
#   
#   result <- PMIKalmanStep(x,P, A, G, Q, R, y)
#   
#   xdata[i] <- result$x_new
#   P <- result$P_new
#   
#   xt[,,i] <- result$x_new
#   Pt[,,i] <- result$P_new
#    
#   print(i)
#   i<- i+1
#   return(x)
# }
# 
# 
# 
# # N <- 1000
# # sigma_e <- 2
# # sigma_eta <- 1
# # 
# # e_t <- rnorm(N, 0, sigma_e)
# # eta_t <- rnorm(N, 0, sigma_eta)
# # 
# # y_t <- cumsum(eta_t)
# # y_t <- y_t + e_t
# # 
# # plot(y_t)
# # 
# # library(rucm)
# # 
# # modelNile <- ucm(Nile~0, data = Nile, level = T)
# # plot(Nile)
# # 
# # modelNile
# # lines(modelNile$s.level)
# # 
# # a1 <- 0
# # P1 <- 10^7 
# # sigmasq_e <- 15099
# # sigmasq_eta <- 1469.1
# # 
# # P <- function(t)
# # {
# #   if(t==1) return(P1)
# #   else
# #   {
# #     P_tmp <- P(t-1)
# #     return(P_tmp-P_tmp^2/(P_tmp+sigmasq_e)+sigmasq_eta)
# #   }
# # }
# # 
# # a <- function(t,y)
# # {
# #   if(t==1) return(a1)
# #   else
# #   {
# #     a_tmp <- a(t-1,y)
# #     P_tmp <- P(t-1)
# #     return(a_tmp+P_tmp/(P_tmp+sigmasq_e)*(y[t-1]-a_tmp))
# #   }
# # }
# # 
# # L <- function(t) return(sigmasq_e/(P(t)+sigmasq_e))
# # 
# # P2 <- function(t) 
# # {
# #   if(t==1) return(L(1)*P(1)+sigmasq_eta)
# #   else
# #     return(L(t-1)*P(t-1)+sigmasq_eta)
# # }
# # 
# # VarV<- function(t)
# # {
# #   return(P(t)+sigmasq_e)
# # }
# # 
# # time <- 1:length(Nile)
# # 
# # a_wrap <- function(t) return(a(t,Nile))
# # 
# # NilePred_1 <- unlist(lapply(time, a_wrap))
# # VarAlfa <- unlist(lapply(time, P))
# # VarAlfa2 <- unlist(lapply(time, P2))
# # 
# # df <- data.table(time, Nile,NilePred_1,VarAlfa,VarAlfa2)
# # 
# # ggplot(df[2:100,]) + geom_line(aes(x=time, y = Nile, color = "Original")) + geom_line(aes(x=time, y = NilePred_1, color = "Prediction"))
# # ggplot(df[2:100,])+geom_line(aes(x=time, y=VarAlfa))
# 
