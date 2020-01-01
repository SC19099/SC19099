## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
	echo = TRUE,
	message = FALSE,
	warning = FALSE
)

## ----eval=FALSE----------------------------------------------------------
#  function(m){
#    n <- nrow(m)
#    G <- matrix(nrow = n,ncol = n)
#    for (i in 1:n) {
#      for (j in 1:n) {
#        G[i,j] <- m[i,] %*% m[j,]
#      }
#    }
#    return(G)
#  }

## ----eval=FALSE----------------------------------------------------------
#  library(SC19099)
#  data(data)
#  g <- Gram(state.x77)
#  print(g)

## ----eval=FALSE----------------------------------------------------------
#  function(p,pr){
#    mat <- matrix(0,nrow = p,ncol = p)
#    for (i in 1:(p-1)) {
#      for (j in (i+1):p) {
#        prob <- pr[i,j]
#        lineornot <- sample(c(0,1),1,prob = c(1-prob,prob))
#        if(lineornot == 1){
#          mat[i,j] <- 1
#        }
#      }
#    }
#    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
#  }

## ----eval=FALSE----------------------------------------------------------
#  library(SC19099)
#  p <- 100
#  x <- runif(p,0,1)
#  y <- runif(p,0,1)
#  node <- cbind(x,y)
#  distance <- as.matrix(dist(node))
#  pr <- matrix(0,nrow = p,ncol = p)
#  for (i in 1:(p-1)) {
#  for (j in (i+1):p) {
#   d <- distance[i,j]
#   pr[i,j] <- dnorm(d * sqrt(p))
#   }}
#  e <- edgemat(p,pr)
#  print(e)

## ------------------------------------------------------------------------
knitr::kable(head(airquality))

## ------------------------------------------------------------------------
height <- c(165,168,170,172,175,178,180,182,183,185)
weight <- c(60,55,64,65,70,70,71,75,80,78)
lm.hw <- lm(weight~height)
plot(lm.hw)

## ----warning=FALSE-------------------------------------------------------
rayleigh <- function(n,sigma){
 u <- runif(n)
 x <- (-2 * sigma ^ 2 * log(1 - u)) ^ 0.5 #F(x)=1-exp(-x^2/2/sigma^2),x>=0
 return(x)
}

## ------------------------------------------------------------------------
r_sam1 <- rayleigh(1000,1)
r_sam2 <- rayleigh(1000,5)
r_sam3 <- rayleigh(1000,8)

## ------------------------------------------------------------------------
hist(rayleigh(1000,1),prob = TRUE,,main = expression(f(x)==x*exp(-x^2/(2*sigma^2))/sigma^2~~~sigma==1))
y <- seq(0, 100, 0.1)
sigma <- 1
lines(y,y*exp(-y^2/(2*sigma^2))/sigma^2)
hist(rayleigh(1000,5),prob = TRUE,,main = expression(f(x)==x*exp(-x^2/(2*sigma^2))/sigma^2~~~sigma==5))
y <- seq(0, 100, 0.1)
sigma <- 5
lines(y,y*exp(-y^2/(2*sigma^2))/sigma^2)
hist(rayleigh(1000,8),prob = TRUE,,main = expression(f(x)==x*exp(-x^2/(2*sigma^2))/sigma^2~~~sigma==8))
y <- seq(0, 100, 0.1)
sigma <- 8
lines(y,y*exp(-y^2/(2*sigma^2))/sigma^2)

## ----message=FALSE, warning=FALSE----------------------------------------
Mix <- function(n,p){
X1 <- rnorm(n,0,1)
X2 <- rnorm(n,3,1)
r <- sample(c(0,1),n,replace=TRUE,prob=c(1-p,p))
M <- r*X1+(1-r)*X2
return(M)
}

## ------------------------------------------------------------------------
M1 <- Mix(1000,0.75)
hist(M1,main = expression(p==0.75))
M2 <- Mix(1000,0.25)
hist(M2,main = expression(p==0.25))
M3 <- Mix(1000,0.5)
hist(M3,main = expression(p==0.5))

## ------------------------------------------------------------------------
W_sam <- function(sigma,n){
  p <- nrow(sigma)
  L <- chol(sigma) #对sigma进行Cholesky分解
  A <- matrix(nrow = p,ncol = p) #生成p*p的空矩阵
  A[upper.tri(A)] <- 0 #将上三角用0填满
  c <- numeric(p)
  for (i in 1:p) {
    c[i] <- rchisq(1,n - i + 1) 
  }
  diag(A) <- c #用卡方分布随机数填满对角线
  A[lower.tri(A)] <- rnorm((p ^ 2 - p) / 2) #用标准正态随机数填满下三角
  S <- L%*%A%*%t(A)%*%t(L) #Bartlett分解
  return(S)
}

## ------------------------------------------------------------------------
sigma <- 6 * diag(6)
n <- 10
W_sam(sigma,n)

## ------------------------------------------------------------------------
m <- 1e4 #the number of random number
x <- runif(m, min=0, max=pi / 3) 
theta.hat <- mean(pi / 3 * sin(x)) #the estimate
print(c(theta.hat,1 / 2)) #compare the estimate with the theoretical value

## ------------------------------------------------------------------------
m <- 1e4 
set.seed(1)
x <- runif(m) 
thetas <- exp(-x) / (1 + x^2)
theta.hat <- mean(thetas) #the estimate 
var1 <- var(thetas) / m #the variance without variance reduction

MC <- numeric(10000)
for (i in 1:10000){
y <- runif(m / 2)
thetas.new <- (exp(-y) / (1 + y^2) + exp(-(1-y)) / (1 + (1-y)^2)) * 0.5
MC[i] <- mean(thetas.new) 
}
theta.hat.new <- mean(MC) #the estimate using antithetic variables
print(theta.hat.new) #the new estimate
var2 <- var(MC) #the variance with variance reduction
print(c(var1,var2,(var1-var2) / var1)) #show the approximate reduction 

## ------------------------------------------------------------------------
q <- numeric(6)
a <- c(0,.2,.4,.6,.8,1)
for (i in 1:6){
  q[i] <- -log(1 - (1-exp(-1))*a[i]) #generate the numbers dividing the integral equally
}
print(q) 

m <- 2000
theta <- numeric(5)
var <- numeric(5)
for (j in 1:5) {
  set.seed(123)
  U <- runif(m)
  X <- -log(exp(-q[j]) - (1-exp(-1))*U/5) #generate the random number of 5*f_3(x),(q[j],q[j+1])
  thetas <- (1-exp(-1)) / (5*(1+X^2))
  var[j] <- var(thetas) / m #the variance of the integral estimate of subinterval
  theta[j] <- mean(thetas) #the integral of subinterval 
}
theta.hat <- sum(theta) #the esitimate of the intergral
std <- sqrt(sum(var)) #the estimated standard error 
print(c(theta.hat,std))
std.old <- 0.0970314 #from the textbook
print((std.old-std) / std.old) #the percentage reduction of the standard error

## ------------------------------------------------------------------------
n <- 20 #the number of sample data 
alpha <- .05 
set.seed(1)
CL <- replicate(1000, expr = {
  x <- rchisq(n,2) #generate the sample data from χ2(2)
  LCL <- mean(x) + qt(alpha/2,n-1)*sd(x)/sqrt(n) #the lower confidence limit
  UCL <- mean(x) - qt(alpha/2,n-1)*sd(x)/sqrt(n) #the upper confidence limit
  c(LCL,UCL)
} )
fre <- 0 
for (i in 1:1000) {
  if(CL[1,i] < 2 && 2 < CL[2,i]) #the theoretical mean of χ2(2) is 2
    fre <- fre + 1
}
print(CP <- fre / 1000) #the coverage probability of the t-interval

#compute the coverage probability for variance if random samples come from χ2(2)
UCL_var <- replicate(1000, expr = {
y <- rchisq(n,2)
(n-1) * var(y) / qchisq(alpha, df = n-1)
} )
print(mean(UCL_var > 4))#the theoretical variance of χ2(2) is 4

## ------------------------------------------------------------------------
library(moments) #use the package to compute the sample skewness
set.seed(2)
n <- 1000
sk <- replicate(1000,expr = {
  x <- rnorm(n)
  skewness(x)
})
q <- c(0.025,0.05,0.95,0.975)
quantiles <- c(quantile(sk,q[1]),quantile(sk,q[2]),quantile(sk,q[3]),quantile(sk,q[4])) 
# estimate the 0.025, 0.05, 0.95, and 0.975 quantiles of the skewness
print(quantiles)

mean_sk <- mean(sk)
sd_sk <- sd(sk)
std_estimate <- numeric(4)
for (i in 1:4){
  Xq <- quantiles[i]
  std_estimate[i] <- sqrt(q[i] * (1-q[i]) / (n*dnorm(Xq,mean = 0,sd = 1)^2)) 
}
print(std_estimate)#compute the standard error of the estimates from (2.14)

quantile_large <- c(qnorm(0.025,mean = 0,sd = sqrt(6/n)),qnorm(0.05,mean = 0,sd = sqrt(6/n)),qnorm(0.95,mean = 0,sd = sqrt(6/n)),qnorm(0.975,mean = 0,sd = sqrt(6/n)))
# the quantiles of the large sample approximation √b1 ≈ N(0, 6/n)
print(quantile_large)

## ------------------------------------------------------------------------
library(energy)
set.seed(123)
alpha <- .05
n <- 30
m <- 2500 #try smaller m for a trial run
test1 <- test2 <- numeric(m)
para1 <- para2 <- seq(1,10,1)
power1 <- power2 <- numeric(10)

#Estimate the power of the skewness test of normality against symmetric Beta(α, α) distributions
for (j in 1:10) {
  for (i in 1:m) {
    x <- rbeta(n,para1[j],para1[j]) 
    test1[i] <- as.integer(
    shapiro.test(x)$p.value <= alpha)
  }
  power1[j] <- mean(test1)
}
par(mfcol=c(1,2))
plot(para1,power1,main = "The power curve against Beta(α, α)",xlab = "α",ylab = "power",type = "l",col = "red")

#Estimate the power of the skewness test of normality against symmetric t(ν) distributions
for (j in 1:10) {
  for (i in 1:m) {
    y <- rt(n,para2[j]) 
    test2[i] <- as.integer(
    shapiro.test(y)$p.value <= alpha)
  }
  power2[j] <- mean(test2)
}
plot(para2,power2,main = "The power curve against t(ν)",xlab = "ν",ylab = "power",type = "l",col = "red")
print(power1)
print(power2)

## ------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 1
sigma <- 100
m <- 10000 #number of replicates
p1 <- p2 <- p3 <- numeric(m) #storage for p-values

set.seed(12345)
for (j in 1:m) {
  x1 <- rchisq(n,1)
  ttest1 <- t.test(x1, alternative = "greater", mu = mu0)
  p1[j] <- ttest1$p.value
}
p.hat1 <- mean(p1 < alpha)

for (j in 1:m) {
  x2 <- runif(n,0,2)
  ttest2 <- t.test(x2, alternative = "greater", mu = mu0)
  p2[j] <- ttest2$p.value
}
p.hat2 <- mean(p2 < alpha)

for (j in 1:m) {
  x3 <- rexp(n,1)
  ttest3 <- t.test(x3, alternative = "greater", mu = mu0)
  p3[j] <- ttest3$p.value
}
p.hat3 <- mean(p3 < alpha)
print(c(alpha,p.hat1,p.hat2,p.hat3))

## ------------------------------------------------------------------------
library(bootstrap)
library(corrplot)
library(kableExtra)
set.seed(0)
pairs(~mec+vec+alg+ana+sta,data = scor,main = "Scatterplot Matrix")

cor <- cor(scor)
print(cor)
corrplot(cor)

r12 <- function(x, i) {
#want correlation of columns 1 and 2
cor(x[i,1], x[i,2])
}
r34 <- function(x, i) {
#want correlation of columns 1 and 2
cor(x[i,3], x[i,4])
}
r35 <- function(x, i) {
#want correlation of columns 1 and 2
cor(x[i,3], x[i,5])
}
r45 <- function(x, i) {
#want correlation of columns 1 and 2
cor(x[i,4], x[i,5])
}
library(boot) #for boot function
t12 <- boot(data = scor, statistic = r12, R = 2000)$t
std12 <- sd(t12)
t34 <- boot(data = scor, statistic = r34, R = 2000)$t
std34 <- sd(t34)
t35 <- boot(data = scor, statistic = r35, R = 2000)$t
std35 <- sd(t35)
t45 <- boot(data = scor, statistic = r45, R = 2000)$t
std45 <- sd(t45)
print(c(std12,std34,std35,std45))

df <- data.frame(std12,std34,std35,std45)
colnames(df) <- c("ρ^12", "ρ^34", "ρ^35","ρ^45")
kable(df, escape = FALSE) %>% kable_styling(position = "center")

## ------------------------------------------------------------------------
library(boot)
library(kableExtra)
set.seed(123)
sk <- function(x,i) {
m3 <- mean((x[i] - mean(x[i]))^3)
m2 <- mean((x[i] - mean(x[i]))^2)
return( m3 / m2^1.5 )
}

cou.n <- cou.b <- cou.p <- 0
for (i in 1:100) {
  x <- rnorm(50)
  true.sk <- 0
  boot.obj <- boot(x,statistic = sk,R = 2000)
  ci <- boot.ci(boot.obj,type = c("norm","basic","perc"))
  ci.normal <- ci$normal
  ci.basic <- ci$basic
  ci.perc <- ci$percent
  if(ci.normal[2] < true.sk && true.sk < ci.normal[3])
  {cou.n <- cou.n + 1}
  if(ci.basic[4] < true.sk && true.sk < ci.basic[5])
  {cou.b <- cou.b + 1}
  if(ci.perc[4] < true.sk && true.sk < ci.perc[5])
  {cou.p <- cou.p + 1}
}
normal.fre.n <- cou.n / 100
normal.fre.b <- cou.b / 100
normal.fre.p <- cou.p / 100

cou.n <- cou.b <- cou.p <- 0
for (i in 1:100) {
  y <- rchisq(50,5)
  true.sk <- sqrt(8 / 5)
  boot.obj <- boot(y,statistic = sk,R = 2000)
  ci <- boot.ci(boot.obj,type = c("norm","basic","perc"))
  ci.normal <- ci$normal
  ci.basic <- ci$basic
  ci.perc <- ci$percent
  if(ci.normal[2] < true.sk && true.sk < ci.normal[3])
  {cou.n <- cou.n + 1}
  if(ci.basic[4] < true.sk && true.sk < ci.basic[5])
  {cou.b <- cou.b + 1}
  if(ci.perc[4] < true.sk && true.sk < ci.perc[5])
  {cou.p <- cou.p + 1}
}
chisq.fre.n <- cou.n / 100
chisq.fre.b <- cou.b / 100
chisq.fre.p <- cou.p / 100

df <- data.frame(x = c(normal.fre.n, chisq.fre.n),
                 y = c(normal.fre.b, chisq.fre.b),
                 z = c(normal.fre.p, chisq.fre.p))
colnames(df) <- c("norm", "basic", "perc")
rownames(df) <- c("N(0,1)", "χ2(5)")
kable(df, escape = FALSE) %>% row_spec(1:2, bold = T) %>% kable_styling(position = "center")

## ------------------------------------------------------------------------
library(bootstrap)
library(kableExtra)

data(scor,package = "bootstrap")
n <- nrow(scor)
lamda <- eigen(cov(scor))$values 
theta.hat <- lamda[1] / sum(lamda) 
theta.jack <- numeric(n) 
for (i in 1:n)
{lamda.jack <- eigen(cov(scor[-i,]))$values 
theta.jack[i] <- lamda.jack[1] / sum(lamda.jack) 
}
bias.jack <- (n - 1) * (mean(theta.jack) - theta.hat) #the jackknife estimate of bias
se.jack <- sqrt((n - 1) * mean((theta.jack - mean(theta.jack)) ^ 2))  #the jackknife estimate of standard error 

df <- data.frame(theta.hat,bias.jack,se.jack)
colnames(df) <- c("theta.hat", "bias", "standard error")
kable(df, escape = FALSE) %>% kable_styling(position = "center")

## ----message=FALSE-------------------------------------------------------
library(DAAG)
library(kableExtra)
attach(ironslag)

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)

for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]

J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1

J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2

J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3

J4 <- lm(y ~ x + I(x^2) + I(x^3))
yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] +
J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
e4[k] <- magnetic[k] - yhat4
}

df <- data.frame(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))
colnames(df) <- c("mean(e1^2)", "mean(e2^2)", "mean(e3^2)","mean(e4^2)")
kable(df, escape = FALSE) %>% kable_styling(position = "center")

a <- seq(10, 40, .1) #sequence for plotting fits
par(mar=c(1,1,1,1))

L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)

L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)

L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)

L4 <- lm(log(magnetic) ~ log(chemical))
plot(log(chemical), log(magnetic), main="Log-Log", pch=16)
logyhat4 <- L4$coef[1] + L4$coef[2] * log(a)
lines(log(a), logyhat4, lwd=2)

adj1 <- summary(L1)$adj.r.squared
adj2 <- summary(L2)$adj.r.squared
adj3 <- summary(L3)$adj.r.squared
adj4 <- summary(L4)$adj.r.squared

df <- data.frame(adj1,adj2,adj3,adj4)
colnames(df) <- c("adj1", "adj2", "adj3","adj4")
kable(df, escape = FALSE) %>% kable_styling(position = "center")

## ------------------------------------------------------------------------
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 10000

set.seed(1234)
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)

# compute the maximum numebers of extreme points m1 and m2
log(.025) / log(n1 / (n1 + n2))
log(.025) / log(n2 / (n1 + n2))
m1 <- 4
m2 <- 7

# original statistic
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer((outx > m1) | (outy > m2)))
}

R <- 9999
z <- c(x,y)
K <- n1 + n2
reps <- numeric(R)
t0 <- count5test(x,y)
for (i in 1:R) {
  k <- sample(K, size = n1, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k]
  X <- x1 - mean(x1)
  Y <- y1 - mean(y1)
  reps[i] <- count5test(x1, y1)
}

# compute alphahat
alphahat <- mean(c(t0, reps) > t0)
print(alphahat)

## ------------------------------------------------------------------------
library(Ball)
library(mvtnorm)
library(boot)
library(ggplot2)
# distance correlation function
dCov <- function(x, y) {
  x <- as.matrix(x); y <- as.matrix(y)
  n <- nrow(x); m <- nrow(y)
  if (n != m || n < 2) stop("Sample sizes must agree")
  if (! (all(is.finite(c(x, y)))))
  stop("Data contains missing or infinite values")
  Akl <- function(x) {
  d <- as.matrix(dist(x))
  m <- rowMeans(d); M <- mean(d)
  a <- sweep(d, 1, m); b <- sweep(a, 2, m)
  b + M
  }
A<- Akl(x); B <- Akl(y)
sqrt(mean(A * B))
}
ndCov2 <- function(z, ix, dims) {
#dims contains dimensions of x and y
p <- dims[1]
q <- dims[2]
d <- p + q
x <- z[ , 1:p] #leave x as is
y <- z[ix, -(1:p)] #permute rows of y
return(nrow(z) * dCov(x, y)^2)
}
# generate sample
n<-seq(from=10,to=100,by=10)
# loop
k<-100
# significant level
alpha<-0.05
pow_dCor_Model1<-pow_ball_Model1<-pow_dCor_Model2<-pow_ball_Model2<-numeric(length(n))
for (j in 1:length(n)) {

  #storage of temp data
  p_ball1<-numeric(k)
  dcor1<-numeric(k)
  p_ball2<-numeric(k)
  dcor2<-numeric(k)
  dcor1<-dcor2<-p_ball1<-p_ball2<-numeric(k)
  for (i in 1:k) {
    set.seed(i)
    # the function "rmvnorm" is used to 
    # generate the multidimensional normal data
    X<-rmvnorm(n[j],rep(0,2),diag(1,2))
    err<-rmvnorm(n[j],rep(0,2),diag(1,2))
    Y1<-(X/4)+err
    Y2<-(X/4)*err
    Z1<-cbind(X,Y1)
    Z2<-cbind(X,Y2)
    t1<-bcov.test(X,Y2,R=99)
    p_ball1[i]<-t1$p.value
    boot.obj1<-boot(data=Z1,statistic=ndCov2,R=99,sim="permutation",dims=c(2, 2))
    temp1<-c(boot.obj1$t0, boot.obj1$t)
    dcor1[i]<-mean(temp1>=temp1[1])
    
    t2<-bcov.test(X,Y2,R=99)
    p_ball2[i]<-t2$p.value
    boot.obj2<-boot(data=Z2,statistic=ndCov2,R=99,sim="permutation",dims=c(2, 2))
    temp2<-c(boot.obj2$t0, boot.obj2$t)
    dcor2[i]<-mean(temp2>=temp2[1])
    }
  pow_dCor_Model1[j]<-mean(dcor1<alpha)
  pow_ball_Model1[j]<-mean(p_ball1<alpha)
  pow_dCor_Model2[j]<-mean(dcor2<alpha)
  pow_ball_Model2[j]<-mean(p_ball2<alpha)  
}
dat<-data.frame(pow_dCor_Model1,pow_ball_Model1,pow_dCor_Model2,pow_ball_Model2)
# the red one is distance correlation test and the blue one is ballcovariance test
ggplot(dat,aes(n))+geom_point(y=pow_dCor_Model1,fill="white")+geom_line(y=pow_dCor_Model1,colour="yellow")+geom_point(y=pow_ball_Model1,fill="white")+geom_line(y=pow_ball_Model1,colour="green")

ggplot(dat,aes(n))+geom_point(y=pow_dCor_Model2,fill="white")+geom_line(y=pow_dCor_Model2,colour="yellow")+geom_point(y=pow_ball_Model2,fill="white")+geom_line(y=pow_ball_Model2,colour="green")

## ------------------------------------------------------------------------
library(kableExtra)
library(GeneralizedHyperbolic)
rw.Metropolis <- function(sigma,x0,N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= exp(abs(x[i-1]) - abs(y))){
x[i] <- y 
k <- k + 1}
else {
x[i] <- x[i-1]

} }
return(list(x=x, k=k))
}

set.seed(123)
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 0
rw1 <- rw.Metropolis(sigma[1],x0,N)
rw2 <- rw.Metropolis(sigma[2],x0,N)
rw3 <- rw.Metropolis(sigma[3],x0,N)
rw4 <- rw.Metropolis(sigma[4],x0,N)

x1 <- rw1$x
x2 <- rw2$x
x3 <- rw3$x
x4 <- rw4$x
hh <- c(qskewlap(0.025),qskewlap(0.975))
par(mar=c(1,1,1,1))
index <- 1:2000
plot(index, x1, type="l", main="σ=0.05", ylab="x1",ylim=c(-4,4) )
abline(h=hh)
plot(index, x2, type="l", main="σ=0.5", ylab="x2")
abline(h=hh)
plot(index, x3, type="l", main="σ=2", ylab="x3")
abline(h=hh)
plot(index, x4, type="l", main="σ=16", ylab="x4")
abline(h=hh)
df <- data.frame(rw1$k / N, rw2$k / N, rw3$k / N, rw4$k / N) #acceptance rates
colnames(df) <- c("σ=0.05","σ=0.5","σ=2","σ=16")
kable(df, escape = FALSE) %>% kable_styling(bootstrap_options = "striped",position = "center")

## ------------------------------------------------------------------------
a <- log(exp(10))
b <- exp(log(10))
a == b
isTRUE(all.equal(a,b))
print(a-b)

## ------------------------------------------------------------------------
f <- function(a,k) 1-pt(sqrt(a^2*k/(k+1-a^2)),k)
k <- 1000
g <- function(x) f(x,k-1) - f(x,k)

a <- seq(0,sqrt(k),.1)
plot(a,g(a),type = 'l')
abline(h = 0)
uniroot(g,c(1,5))$root

f1 <- function(k) 2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))
ck <- function(a,k) sqrt(a^2*k/(k+1-a^2))
g2 <- function(u,k) (1+u^2/k)^(-(k+1)/2)
k <- 1000

fun <- function(a) f1(k)*integrate(function(u) {g2(u,k)}, 0, ck(a,k))$value - f1(k-1)*integrate(function(u) {g2(u,k-1)}, 0, ck(a,k-1))$value 
uniroot(fun,c(1,5))$root

## ------------------------------------------------------------------------
nA <- 28
nB <- 24
nOO <- 41
nAB <- 70

theta0 <- c(0.1,0.1)
l <- numeric(1000)
for (j in 1:1000) {
  E <- function(theta) {
    p <- theta[1]
    q <- theta[2]
    r <- 1-p-q
    p0 <- theta0[1]
    q0 <- theta0[2]
    r0 <- 1-p0-q0
    return(2*nA*(log(p)+r0/(p0+2*r0)*log(2*r/p))+2*nB*(log(q)+r0/(q0+2*r0)*log(2*r/q))+2*nOO*log(r)+nAB*log(2*p*q))
  }
  Obj <- function(theta) -E(theta)
  optim <- optim(c(0.1,0.1), Obj)
  theta0 <- optim$par
  l[j] <- -optim$value
}
print(theta0)

plot(l[1:20], type = 'l', xlab = 'iterations', ylab = 'log-maximum likelihood values')

## ----warning=FALSE-------------------------------------------------------
formulas <- list( mpg ~ disp, mpg ~ I(1 / disp), mpg ~ disp + wt, mpg ~ I(1 / disp) + wt )
fit1 <- fit2 <- list(4)

#for loops
for (i in 1:4) {
  fit1[[i]] <- lm(formulas[[i]], mtcars)
}
print(fit1)

#lapply()
fit2 <- lapply(formulas, function(x) lm(x,mtcars))
print(fit2)

## ----warning=FALSE-------------------------------------------------------
set.seed(1)
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
n <- length(bootstraps)
fit3 <- fit4 <- fit5 <- list(n)

#a for loop
for (i in 1:n) {
  fit3[[i]] <- lm(mpg ~ disp, bootstraps[[i]])
}
print(fit3)

#lapply()
fit4 <- lapply(1:n, function(x) lm(mpg ~ disp, bootstraps[[x]]))
print(fit4)

#lapply() without an anonymous function
fit5 <- lapply(bootstraps, lm, formula = mpg ~ disp)
print(fit5)

## ------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
r.squared1 <- list(4)
r.squared2 <- list(10)

#Question 1
r.squared1 <- lapply(fit1, rsq)
print(r.squared1)

#Question 2
r.squared2 <- lapply(fit3, rsq)
print(r.squared2)

## ------------------------------------------------------------------------
set.seed(12)
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)

#sapply()
sapply(trials, function(x) x$p.value)

#"[["
sapply(trials, '[[', 'p.value')

## ----message=FALSE, warning=FALSE, eval=FALSE----------------------------
#  library(parallel)
#  boot_df <- function(x) x[sample(nrow(x), rep = T), ]
#  rsquared <- function(mod) summary(mod)$r.squared
#  boot_lm <- function(i) {
#    rsquared(lm(mpg ~ wt + disp, data = boot_df(mtcars)))
#  }
#  system.time(sapply(1:1e5, boot_lm))
#  system.time(unlist(mclapply(1:1e5, boot_lm, mc.cores = 4)))

## ----warning=FALSE-------------------------------------------------------
#the function written before
library(kableExtra)
library(GeneralizedHyperbolic)
rw.Metropolis <- function(sigma,x0,N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= exp(abs(x[i-1]) - abs(y))){
x[i] <- y 
k <- k + 1}
else {
x[i] <- x[i-1]

} }
return(list(x=x, k=k))
}

#the function using Rcpp
library(Rcpp)
cppFunction('List rw_Metropolis(double sigma, double x0, int N) {
NumericVector x(N);
x[0] = x0;
int k = 0;
for (int i = 1;i < N;i++) {
double u = runif(1)[0];
double y = rnorm(1, x[i-1], sigma)[0];
if (u <= exp(abs(x[i-1]) - abs(y))){
x[i] = y;
k = k + 1;}
else 
x[i] = x[i-1];
}
List result = List::create(x,k);
return(result);
}')

#generate random samples
set.seed(123)
N <- 1000
sigma <- 1
x0 <- 0
sample1 <- rw.Metropolis(sigma,x0,N)$x
sample2 <- rw_Metropolis(sigma,x0,N)[[1]]

#qq plot
library(car)
qqplot(sample1, sample2, xlab = "the samples using R",
       ylab = "the samples using Rcpp")
x <- seq(-4,4,.01)
lines(x,x,col = "red")

#Campare the computation time
library(microbenchmark)
ts <- microbenchmark(rw1 = rw.Metropolis(sigma,x0,N),rw2 = rw_Metropolis(sigma,x0,N))
summary(ts)[,c(1,3,5,6)]

