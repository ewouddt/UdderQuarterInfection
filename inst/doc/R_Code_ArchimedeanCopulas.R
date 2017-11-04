###########################
### ARCHIMEDEAN COPULAS ###
###########################



data("udderquarterinfection",envir=environment())
udder <- udderquarterinfection

LLH_temp <- rep(NA,4)
names(LLH_temp) <- paste0("model",0:3)
estimate_temp <- matrix(NA,nrow=12,ncol=4,dimnames = list(
  c(paste0("lambda",1:4),paste0("rho",1:4),"beta",paste0("theta",0:2)),
  paste0("model",0:3)
))

OUT <- list(
  onestageparametric=list(estimate=estimate_temp,stderror=NULL,LLH=LLH_temp),
  twostageparametric=list(estimate=estimate_temp,stderror=NULL,LLH=LLH_temp),
  twostagesemiparametric=list(estimate=estimate_temp,stderror=NULL,LLH=LLH_temp)
)


#=======================================#
# Generate variables for input in models#
#=======================================#

t <- c()
for (i in 1:nrow(udder)){
  if (udder$status[i]==1){t[i] <- ((udder$left[i]+udder$right[i])/2)*4/365.25}
  else {t[i] <- udder$right[i]*4/365.25}}
udder$t <- t

Q1<-udder[udder$quarter=="LF",];Q3<-udder[udder$quarter=="LR",]
Q2<-udder[udder$quarter=="RF",];Q4<-udder[udder$quarter=="RR",]

t1 <- Q1$t;t2 <- Q2$t;t3 <- Q3$t;t4 <- Q4$t
c <- udder$status;c1 <- Q1$status;c2 <- Q2$status;c3 <- Q3$status;c4 <- Q4$status;
C <- cbind(c1,c2,c3,c4)

Heif <- as.numeric(udder$lactation>1)
Heif1 <- as.numeric(Q1$lactation>1);Heif2 <- as.numeric(Q2$lactation>1)
Heif3 <- as.numeric(Q3$lactation>1);Heif4 <- as.numeric(Q4$lactation>1)

K <- nrow(table(udder$cowid)); #number of clusters

udder$cluster <- rep(1:K,each=4);cluster<-udder$cluster
d <- c1+c2+c3+c4; #number of events for each cluster

ProcentCensoring <- 1-mean(c)

#=================================#
# Two-stage parametric estimation #
#=================================#

# =========================================================================
# first stage: estimate baseline and covariate effect (4 baselines, 1 beta)
# =========================================================================

loglik.par.stage1.2A <- function(p,exp=FALSE){
  if (exp==TRUE){
    lambda1 <- exp(p[1])
    rho1 <-    exp(p[2])
    lambda2 <- exp(p[3])
    rho2 <-    exp(p[4])
    lambda3 <- exp(p[5])
    rho3 <-    exp(p[6])
    lambda4 <- exp(p[7])
    rho4 <-    exp(p[8])
    beta <- p[9]}
  else{
    lambda1 <- p[1]
    rho1 <-    p[2]
    lambda2 <- p[3]
    rho2 <-    p[4]
    lambda3 <- p[5]
    rho3 <-    p[6]
    lambda4 <- p[7]
    rho4 <-    p[8]
    beta <-    p[9]}

  s1 <- exp(-lambda1*t1^rho1*exp(beta*Heif1)); #u_i1
  s2 <- exp(-lambda2*t2^rho2*exp(beta*Heif2));
  s3 <- exp(-lambda3*t3^rho3*exp(beta*Heif3));
  s4 <- exp(-lambda4*t4^rho4*exp(beta*Heif4));
  S <- cbind(s1,s2,s3,s4)

  f1 <- lambda1*rho1*t1^(rho1-1)*exp(beta*Heif1)*s1; #-du_i1/dy_i1
  f2 <- lambda2*rho2*t2^(rho2-1)*exp(beta*Heif2)*s2;
  f3 <- lambda3*rho3*t3^(rho3-1)*exp(beta*Heif3)*s3;
  f4 <- lambda4*rho4*t4^(rho4-1)*exp(beta*Heif4)*s4;
  F <- cbind(f1,f2,f3,f4)

  G <- C*log(F)
  H <- (1-C)*log(S)

  loglik<-sum(G+H)
  -sum(loglik)}

res2A <- nlm(loglik.par.stage1.2A,c(log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),0.5),exp=TRUE,hessian=TRUE)

lambda1.2A <- exp(res2A$estimate[1])
rho1.2A    <- exp(res2A$estimate[2])
lambda2.2A <- exp(res2A$estimate[3])
rho2.2A    <- exp(res2A$estimate[4])
lambda3.2A <- exp(res2A$estimate[5])
rho3.2A    <- exp(res2A$estimate[6])
lambda4.2A <- exp(res2A$estimate[7])
rho4.2A    <- exp(res2A$estimate[8])
beta.2A    <- res2A$estimate[9]

# WRITING OUTPUT AWAY:
OUT$onestageparametric$estimate[,"model0"] <- c(lambda1.2A,lambda2.2A,lambda3.2A,lambda4.2A,rho1.2A,rho2.2A,rho3.2A,rho4.2A,beta.2A,NA,NA,NA)

S1 <- exp(-lambda1.2A*t1^rho1.2A*exp(beta.2A*Heif1)); #u_i1
S2 <- exp(-lambda2.2A*t2^rho2.2A*exp(beta.2A*Heif2));
S3 <- exp(-lambda3.2A*t3^rho3.2A*exp(beta.2A*Heif3));
S4 <- exp(-lambda4.2A*t4^rho4.2A*exp(beta.2A*Heif4));

DS1 <- -lambda1.2A*rho1.2A*t1^(rho1.2A-1)*exp(beta.2A*Heif1)*S1; #-du_i1/dy_i1
DS2 <- -lambda2.2A*rho2.2A*t2^(rho2.2A-1)*exp(beta.2A*Heif2)*S2;
DS3 <- -lambda3.2A*rho3.2A*t3^(rho3.2A-1)*exp(beta.2A*Heif3)*S3;
DS4 <- -lambda4.2A*rho4.2A*t4^(rho4.2A-1)*exp(beta.2A*Heif4)*S4;

# MODEL 0: no clustering      !! only 1st stage !!
# ######################

#joint survival function
SS <- S1*S2*S3*S4

#first order partial derivatives
dS1 <- DS1*S2*S3*S4
dS2 <- S1*DS2*S3*S4
dS3 <- S1*S2*DS3*S4
dS4 <- S1*S2*S3*DS4

#second order partial derivatives
d2S12 <- DS1*DS2*S3*S4
d2S13 <- DS1*S2*DS3*S4
d2S14 <- DS1*S2*S3*DS4
d2S23 <- S1*DS2*DS3*S4
d2S24 <- S1*DS2*S3*DS4
d2S34 <- S1*S2*DS3*DS4

#third order partial derivatives
d3S123 <- DS1*DS2*DS3*S4
d3S124 <- DS1*DS2*S3*DS4
d3S134 <- DS1*S2*DS3*DS4
d3S234 <- S1*DS2*DS3*DS4

#fourth order partial derivatives
d4S1234 <- DS1*DS2*DS3*DS4

terms3 <- log(SS^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
                (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
                (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
                (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
                (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
                d2S12^(c1*c2*(1-c3)*(1-c4))*
                d2S13^(c1*(1-c2)*c3*(1-c4))*
                d2S14^(c1*(1-c2)*(1-c3)*c4)*
                d2S23^((1-c1)*c2*c3*(1-c4))*
                d2S24^((1-c1)*c2*(1-c3)*c4)*
                d2S34^((1-c1)*(1-c2)*c3*c4)*
                (-d3S123)^(c1*c2*c3*(1-c4))*
                (-d3S124)^(c1*c2*(1-c3)*c4)*
                (-d3S134)^(c1*(1-c2)*c3*c4)*
                (-d3S234)^((1-c1)*c2*c3*c4)*
                d4S1234^(c1*c2*c3*c4))

LL0 <- sum(terms3)
# LL0 #-4939.196

# WRITING OUTPUT AWAY:
OUT$onestageparametric$LLH["model0"] <- LL0

# Second stage: estimate association
# ==================================

# #######################################################
# MODEL 3: Parent copula with two different child copulas
# #######################################################

likelihood.2st.model3 <- function(p,exp=FALSE){

  if (exp==TRUE){

    theta1 <- exp(p[1])
    theta2 <- exp(p[2])
    theta3 <- exp(p[3])
  }

  else{  #default

    theta1 <- p[1]
    theta2 <- p[2]
    theta3 <- p[3]
  }

  A <- -1+(-1+S1^(-theta2)+S2^(-theta2))^(theta1/theta2)+(-1+S3^(-theta3)+S4^(-theta3))^(theta1/theta3)
  B12 <- -1+S1^(-theta2)+S2^(-theta2)
  B34 <- -1+S3^(-theta3)+S4^(-theta3)
  C1 <- S1^(-theta2-1)*DS1
  C2 <- S2^(-theta2-1)*DS2
  C3 <- S3^(-theta3-1)*DS3
  C4 <- S4^(-theta3-1)*DS4

  #joint survival function
  S <- A^(-1/theta1)

  #first order partial derivatives
  dS1 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C1
  dS2 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C2
  dS3 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C3
  dS4 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C4

  #second order partial derivatives
  d2S12 <- A^(-1/theta1-2)*B12^(theta1/theta2-2)*C1*C2*((1+theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d2S13 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C3
  d2S14 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C4
  d2S23 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C3
  d2S24 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C4
  d2S34 <- A^(-1/theta1-2)*B34^(theta1/theta3-2)*C3*C4*((1+theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #third order partial derivatives
  d3S123 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C3*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S124 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C4*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S134 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C1*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)
  d3S234 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C2*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #fourth order partial derivatives
  d4S1234 <- (1+theta1)*A^(-1/theta1-4)*B12^(theta1/theta2-2)*B34^(theta1/theta3-2)*C1*C2*C3*C4*
    ((1+2*theta1)*(1+3*theta1)*B12^(theta1/theta2)*B34^(theta1/theta3)+(1+2*theta1)*A*((-theta1+theta3)*B12^(theta1/theta2)+(-theta1+theta2)*B34^(theta1/theta3))+(-theta1+theta2)*(-theta1+theta3)*A^2)

  terms<-log(S^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
               (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
               (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
               (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
               (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
               d2S12^(c1*c2*(1-c3)*(1-c4))*
               d2S13^(c1*(1-c2)*c3*(1-c4))*
               d2S14^(c1*(1-c2)*(1-c3)*c4)*
               d2S23^((1-c1)*c2*c3*(1-c4))*
               d2S24^((1-c1)*c2*(1-c3)*c4)*
               d2S34^((1-c1)*(1-c2)*c3*c4)*
               (-d3S123)^(c1*c2*c3*(1-c4))*
               (-d3S124)^(c1*c2*(1-c3)*c4)*
               (-d3S134)^(c1*(1-c2)*c3*c4)*
               (-d3S234)^((1-c1)*c2*c3*c4)*
               d4S1234^(c1*c2*c3*c4))
  -sum(terms)}

res3 <- nlm(likelihood.2st.model3,log(c(0.5,0.5,0.5)),exp=TRUE)
LL3 <- -(res3$minimum)
# LL3 # -3965.733

# WRITING OUTPUT AWAY:
OUT$twostageparametric$estimate[,"model3"] <- OUT$onestageparametric$estimate[,"model0"]
OUT$twostageparametric$estimate[c("theta0","theta1","theta2"),"model3"] <- exp(res3$estimate)
OUT$twostageparametric$LLH["model3"] <- LL3

# MODEL 2: Parent copula with two identical child copulas
# #######################################################

likelihood.2st.model2 <- function(p,exp=FALSE){

  if (exp==TRUE){

    theta1 <- exp(p[1])
    theta2 <- exp(p[2])
    theta3 <- theta2
  }

  else{  #default

    theta1 <- p[1]
    theta2 <- p[2]
    theta3 <- theta2
  }

  A <- -1+(-1+S1^(-theta2)+S2^(-theta2))^(theta1/theta2)+(-1+S3^(-theta3)+S4^(-theta3))^(theta1/theta3)
  B12 <- -1+S1^(-theta2)+S2^(-theta2)
  B34 <- -1+S3^(-theta3)+S4^(-theta3)
  C1 <- S1^(-theta2-1)*DS1
  C2 <- S2^(-theta2-1)*DS2
  C3 <- S3^(-theta3-1)*DS3
  C4 <- S4^(-theta3-1)*DS4

  #joint survival function
  S <- A^(-1/theta1)

  #first order partial derivatives
  dS1 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C1
  dS2 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C2
  dS3 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C3
  dS4 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C4

  #second order partial derivatives
  d2S12 <- A^(-1/theta1-2)*B12^(theta1/theta2-2)*C1*C2*((1+theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d2S13 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C3
  d2S14 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C4
  d2S23 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C3
  d2S24 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C4
  d2S34 <- A^(-1/theta1-2)*B34^(theta1/theta3-2)*C3*C4*((1+theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #third order partial derivatives
  d3S123 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C3*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S124 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C4*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S134 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C1*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)
  d3S234 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C2*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #fourth order partial derivatives
  d4S1234 <- (1+theta1)*A^(-1/theta1-4)*B12^(theta1/theta2-2)*B34^(theta1/theta3-2)*C1*C2*C3*C4*
    ((1+2*theta1)*(1+3*theta1)*B12^(theta1/theta2)*B34^(theta1/theta3)+(1+2*theta1)*A*((-theta1+theta3)*B12^(theta1/theta2)+(-theta1+theta2)*B34^(theta1/theta3))+(-theta1+theta2)*(-theta1+theta3)*A^2)

  terms<-log(S^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
               (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
               (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
               (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
               (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
               d2S12^(c1*c2*(1-c3)*(1-c4))*
               d2S13^(c1*(1-c2)*c3*(1-c4))*
               d2S14^(c1*(1-c2)*(1-c3)*c4)*
               d2S23^((1-c1)*c2*c3*(1-c4))*
               d2S24^((1-c1)*c2*(1-c3)*c4)*
               d2S34^((1-c1)*(1-c2)*c3*c4)*
               (-d3S123)^(c1*c2*c3*(1-c4))*
               (-d3S124)^(c1*c2*(1-c3)*c4)*
               (-d3S134)^(c1*(1-c2)*c3*c4)*
               (-d3S234)^((1-c1)*c2*c3*c4)*
               d4S1234^(c1*c2*c3*c4))
  -sum(terms)}

res2 <- nlm(likelihood.2st.model2,log(c(0.5,0.5)),exp=TRUE)
LL2 <- -(res2$minimum)
# LL2 #-3965.847

# WRITING OUTPUT AWAY:
OUT$twostageparametric$estimate[,"model2"] <- OUT$onestageparametric$estimate[,"model0"]
OUT$twostageparametric$estimate[c("theta0","theta1"),"model2"] <- exp(res2$estimate)
OUT$twostageparametric$LLH["model2"] <- LL2

# MODEL 1: one level of clustering
# ################################

likelihood.2st.model1 <- function(p,exp=FALSE){

  if (exp==TRUE){

    theta1 <- exp(p[1])
    theta2 <- theta1
    theta3 <- theta1
  }

  else{  #default

    theta1 <- p[1]
    theta2 <- theta1
    theta3 <- theta1
  }

  A <- -1+(-1+S1^(-theta2)+S2^(-theta2))^(theta1/theta2)+(-1+S3^(-theta3)+S4^(-theta3))^(theta1/theta3)
  B12 <- -1+S1^(-theta2)+S2^(-theta2)
  B34 <- -1+S3^(-theta3)+S4^(-theta3)
  C1 <- S1^(-theta2-1)*DS1
  C2 <- S2^(-theta2-1)*DS2
  C3 <- S3^(-theta3-1)*DS3
  C4 <- S4^(-theta3-1)*DS4

  #joint survival function
  S <- A^(-1/theta1)

  #first order partial derivatives
  dS1 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C1
  dS2 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C2
  dS3 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C3
  dS4 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C4

  #second order partial derivatives
  d2S12 <- A^(-1/theta1-2)*B12^(theta1/theta2-2)*C1*C2*((1+theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d2S13 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C3
  d2S14 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C4
  d2S23 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C3
  d2S24 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C4
  d2S34 <- A^(-1/theta1-2)*B34^(theta1/theta3-2)*C3*C4*((1+theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #third order partial derivatives
  d3S123 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C3*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S124 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C4*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S134 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C1*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)
  d3S234 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C2*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #fourth order partial derivatives
  d4S1234 <- (1+theta1)*A^(-1/theta1-4)*B12^(theta1/theta2-2)*B34^(theta1/theta3-2)*C1*C2*C3*C4*
    ((1+2*theta1)*(1+3*theta1)*B12^(theta1/theta2)*B34^(theta1/theta3)+(1+2*theta1)*A*((-theta1+theta3)*B12^(theta1/theta2)+(-theta1+theta2)*B34^(theta1/theta3))+(-theta1+theta2)*(-theta1+theta3)*A^2)

  terms<-log(S^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
               (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
               (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
               (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
               (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
               d2S12^(c1*c2*(1-c3)*(1-c4))*
               d2S13^(c1*(1-c2)*c3*(1-c4))*
               d2S14^(c1*(1-c2)*(1-c3)*c4)*
               d2S23^((1-c1)*c2*c3*(1-c4))*
               d2S24^((1-c1)*c2*(1-c3)*c4)*
               d2S34^((1-c1)*(1-c2)*c3*c4)*
               (-d3S123)^(c1*c2*c3*(1-c4))*
               (-d3S124)^(c1*c2*(1-c3)*c4)*
               (-d3S134)^(c1*(1-c2)*c3*c4)*
               (-d3S234)^((1-c1)*c2*c3*c4)*
               d4S1234^(c1*c2*c3*c4))
  -sum(terms)}

res1 <- nlm(likelihood.2st.model1,log(c(0.5,0.5)),exp=TRUE)
LL1 <- -(res1$minimum)
# LL1 #-3973.198

# WRITING OUTPUT AWAY:
OUT$twostageparametric$estimate[,"model1"] <- OUT$onestageparametric$estimate[,"model0"]
OUT$twostageparametric$estimate[c("theta0"),"model2"] <- exp(res2$estimate[1])
OUT$twostageparametric$LLH["model1"] <- LL1




#=====================================#
# Two-stage semiparametric estimation #
#=====================================#

# first stage: estimate baseline and covariate effect (4 baselines, 1 beta)
# =========================================================================

position <- udder$quarter

PH_Z2 <- coxph(Surv(t,c)~Heif+strata(position),data=udder)
beta2 <- as.numeric(PH_Z2$coef[1])
sebeta2 <- sqrt(PH_Z2$var)

# beta2;sebeta2 #0.4086354(0.05095114)

#get estimated survival
fit0 <- survfit(PH_Z2,newdata=data.frame(Heif=0));S0 <- fit0$surv

res <- summary(fit0,censored=TRUE) #hieruit de surv halen per stratum!

S01 <- res$surv[res$strata=="LF"]
S02 <- res$surv[res$strata=="RF"]
S03 <- res$surv[res$strata=="LR"]
S04 <- res$surv[res$strata=="RR"]

t01 <- res$time[res$strata=="LF"]
t02 <- res$time[res$strata=="RF"]
t03 <- res$time[res$strata=="LR"]
t04 <- res$time[res$strata=="RR"]

T1 <- t[udder$quarter=="LF"]
T2 <- t[udder$quarter=="RF"]
T3 <- t[udder$quarter=="LR"]
T4 <- t[udder$quarter=="RR"]
S1 <- c()
S2 <- c()
S3 <- c()
S4 <- c()
ZZ <- vector("list",length=K)
C <- vector("list",length=K)

#cluster <- udder$KNR lukt niet omdat de KNR's niet van 1-1196 gaan
cluster <- rep(1:K,each=4)
Z <- Heif

for(i in 1:K){
  ZZ[[i]] <- Z[cluster==i]
  C[[i]] <- c[cluster==i]
  S1[i] <- (S01[t01==T1[i]])^exp(beta2*ZZ[[i]][1])
  S2[i] <- (S02[t02==T2[i]])^exp(beta2*ZZ[[i]][2])
  S3[i] <- (S03[t03==T3[i]])^exp(beta2*ZZ[[i]][3])
  S4[i] <- (S04[t04==T4[i]])^exp(beta2*ZZ[[i]][4])}

DS1 <- -1;
DS2 <- -1;
DS3 <- -1;
DS4 <- -1;

# Second stage: estimate association
# ==================================

# MODEL 3: Parent copula with two different child copulas
# #######################################################

likelihood.2st.model3 <- function(p,exp=FALSE){

  if (exp==TRUE){

    theta1 <- exp(p[1])
    theta2 <- exp(p[2])
    theta3 <- exp(p[3])
  }

  else{  #default

    theta1 <- p[1]
    theta2 <- p[2]
    theta3 <- p[3]
  }

  A <- -1+(-1+S1^(-theta2)+S2^(-theta2))^(theta1/theta2)+(-1+S3^(-theta3)+S4^(-theta3))^(theta1/theta3)
  B12 <- -1+S1^(-theta2)+S2^(-theta2)
  B34 <- -1+S3^(-theta3)+S4^(-theta3)
  C1 <- S1^(-theta2-1)*DS1
  C2 <- S2^(-theta2-1)*DS2
  C3 <- S3^(-theta3-1)*DS3
  C4 <- S4^(-theta3-1)*DS4

  #joint survival function
  S <- A^(-1/theta1)

  #first order partial derivatives
  dS1 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C1
  dS2 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C2
  dS3 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C3
  dS4 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C4

  #second order partial derivatives
  d2S12 <- A^(-1/theta1-2)*B12^(theta1/theta2-2)*C1*C2*((1+theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d2S13 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C3
  d2S14 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C4
  d2S23 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C3
  d2S24 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C4
  d2S34 <- A^(-1/theta1-2)*B34^(theta1/theta3-2)*C3*C4*((1+theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #third order partial derivatives
  d3S123 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C3*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S124 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C4*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S134 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C1*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)
  d3S234 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C2*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #fourth order partial derivatives
  d4S1234 <- (1+theta1)*A^(-1/theta1-4)*B12^(theta1/theta2-2)*B34^(theta1/theta3-2)*C1*C2*C3*C4*
    ((1+2*theta1)*(1+3*theta1)*B12^(theta1/theta2)*B34^(theta1/theta3)+(1+2*theta1)*A*((-theta1+theta3)*B12^(theta1/theta2)+(-theta1+theta2)*B34^(theta1/theta3))+(-theta1+theta2)*(-theta1+theta3)*A^2)

  terms<-log(S^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
               (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
               (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
               (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
               (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
               d2S12^(c1*c2*(1-c3)*(1-c4))*
               d2S13^(c1*(1-c2)*c3*(1-c4))*
               d2S14^(c1*(1-c2)*(1-c3)*c4)*
               d2S23^((1-c1)*c2*c3*(1-c4))*
               d2S24^((1-c1)*c2*(1-c3)*c4)*
               d2S34^((1-c1)*(1-c2)*c3*c4)*
               (-d3S123)^(c1*c2*c3*(1-c4))*
               (-d3S124)^(c1*c2*(1-c3)*c4)*
               (-d3S134)^(c1*(1-c2)*c3*c4)*
               (-d3S234)^((1-c1)*c2*c3*c4)*
               d4S1234^(c1*c2*c3*c4))
  -sum(terms)}

res3 <- nlm(likelihood.2st.model3,log(c(0.5,0.5,0.5)),exp=TRUE)
LL3 <- -(res3$minimum)
# LL3 # -396.6259

# WRITING OUTPUT AWAY:
OUT$twostagesemiparametric$estimate[,"model3"] <- OUT$onestageparametric$estimate[,"model0"]
OUT$twostagesemiparametric$estimate[c("theta0","theta1","theta2"),"model3"] <- exp(res3$estimate)
OUT$twostagesemiparametric$LLH["model3"] <- LL3

# MODEL 2: Parent copula with two identical child copulas
# #######################################################

likelihood.2st.model2 <- function(p,exp=FALSE){

  if (exp==TRUE){

    theta1 <- exp(p[1])
    theta2 <- exp(p[2])
    theta3 <- theta2
  }

  else{  #default

    theta1 <- p[1]
    theta2 <- p[2]
    theta3 <- theta2
  }

  A <- -1+(-1+S1^(-theta2)+S2^(-theta2))^(theta1/theta2)+(-1+S3^(-theta3)+S4^(-theta3))^(theta1/theta3)
  B12 <- -1+S1^(-theta2)+S2^(-theta2)
  B34 <- -1+S3^(-theta3)+S4^(-theta3)
  C1 <- S1^(-theta2-1)*DS1
  C2 <- S2^(-theta2-1)*DS2
  C3 <- S3^(-theta3-1)*DS3
  C4 <- S4^(-theta3-1)*DS4

  #joint survival function
  S <- A^(-1/theta1)

  #first order partial derivatives
  dS1 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C1
  dS2 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C2
  dS3 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C3
  dS4 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C4

  #second order partial derivatives
  d2S12 <- A^(-1/theta1-2)*B12^(theta1/theta2-2)*C1*C2*((1+theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d2S13 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C3
  d2S14 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C4
  d2S23 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C3
  d2S24 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C4
  d2S34 <- A^(-1/theta1-2)*B34^(theta1/theta3-2)*C3*C4*((1+theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #third order partial derivatives
  d3S123 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C3*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S124 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C4*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S134 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C1*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)
  d3S234 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C2*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #fourth order partial derivatives
  d4S1234 <- (1+theta1)*A^(-1/theta1-4)*B12^(theta1/theta2-2)*B34^(theta1/theta3-2)*C1*C2*C3*C4*
    ((1+2*theta1)*(1+3*theta1)*B12^(theta1/theta2)*B34^(theta1/theta3)+(1+2*theta1)*A*((-theta1+theta3)*B12^(theta1/theta2)+(-theta1+theta2)*B34^(theta1/theta3))+(-theta1+theta2)*(-theta1+theta3)*A^2)

  terms<-log(S^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
               (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
               (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
               (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
               (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
               d2S12^(c1*c2*(1-c3)*(1-c4))*
               d2S13^(c1*(1-c2)*c3*(1-c4))*
               d2S14^(c1*(1-c2)*(1-c3)*c4)*
               d2S23^((1-c1)*c2*c3*(1-c4))*
               d2S24^((1-c1)*c2*(1-c3)*c4)*
               d2S34^((1-c1)*(1-c2)*c3*c4)*
               (-d3S123)^(c1*c2*c3*(1-c4))*
               (-d3S124)^(c1*c2*(1-c3)*c4)*
               (-d3S134)^(c1*(1-c2)*c3*c4)*
               (-d3S234)^((1-c1)*c2*c3*c4)*
               d4S1234^(c1*c2*c3*c4))
  -sum(terms)}

res2 <- nlm(likelihood.2st.model2,log(c(0.5,0.5)),exp=TRUE)
LL2 <- -(res2$minimum)
# LL2 #-396.6815

# WRITING OUTPUT AWAY:
OUT$twostagesemiparametric$estimate[,"model2"] <- OUT$onestageparametric$estimate[,"model0"]
OUT$twostagesemiparametric$estimate[c("theta0","theta1"),"model2"] <- exp(res2$estimate)
OUT$twostagesemiparametric$LLH["model2"] <- LL2

# MODEL 1: one level of clustering
# ################################

likelihood.2st.model1 <- function(p,exp=FALSE){

  if (exp==TRUE){

    theta1 <- exp(p[1])
    theta2 <- theta1
    theta3 <- theta1
  }

  else{  #default

    theta1 <- p[1]
    theta2 <- theta1
    theta3 <- theta1
  }

  A <- -1+(-1+S1^(-theta2)+S2^(-theta2))^(theta1/theta2)+(-1+S3^(-theta3)+S4^(-theta3))^(theta1/theta3)
  B12 <- -1+S1^(-theta2)+S2^(-theta2)
  B34 <- -1+S3^(-theta3)+S4^(-theta3)
  C1 <- S1^(-theta2-1)*DS1
  C2 <- S2^(-theta2-1)*DS2
  C3 <- S3^(-theta3-1)*DS3
  C4 <- S4^(-theta3-1)*DS4

  #joint survival function
  S <- A^(-1/theta1)

  #first order partial derivatives
  dS1 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C1
  dS2 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C2
  dS3 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C3
  dS4 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C4

  #second order partial derivatives
  d2S12 <- A^(-1/theta1-2)*B12^(theta1/theta2-2)*C1*C2*((1+theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d2S13 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C3
  d2S14 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C4
  d2S23 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C3
  d2S24 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C4
  d2S34 <- A^(-1/theta1-2)*B34^(theta1/theta3-2)*C3*C4*((1+theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #third order partial derivatives
  d3S123 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C3*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S124 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C4*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S134 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C1*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)
  d3S234 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C2*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #fourth order partial derivatives
  d4S1234 <- (1+theta1)*A^(-1/theta1-4)*B12^(theta1/theta2-2)*B34^(theta1/theta3-2)*C1*C2*C3*C4*
    ((1+2*theta1)*(1+3*theta1)*B12^(theta1/theta2)*B34^(theta1/theta3)+(1+2*theta1)*A*((-theta1+theta3)*B12^(theta1/theta2)+(-theta1+theta2)*B34^(theta1/theta3))+(-theta1+theta2)*(-theta1+theta3)*A^2)

  terms<-log(S^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
               (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
               (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
               (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
               (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
               d2S12^(c1*c2*(1-c3)*(1-c4))*
               d2S13^(c1*(1-c2)*c3*(1-c4))*
               d2S14^(c1*(1-c2)*(1-c3)*c4)*
               d2S23^((1-c1)*c2*c3*(1-c4))*
               d2S24^((1-c1)*c2*(1-c3)*c4)*
               d2S34^((1-c1)*(1-c2)*c3*c4)*
               (-d3S123)^(c1*c2*c3*(1-c4))*
               (-d3S124)^(c1*c2*(1-c3)*c4)*
               (-d3S134)^(c1*(1-c2)*c3*c4)*
               (-d3S234)^((1-c1)*c2*c3*c4)*
               d4S1234^(c1*c2*c3*c4))
  -sum(terms)}

res1 <- nlm(likelihood.2st.model1,log(c(0.5)),exp=TRUE)
LL1 <- -(res1$minimum)
# LL1 #-404.38

# WRITING OUTPUT AWAY:
OUT$twostagesemiparametric$estimate[,"model1"] <- OUT$onestageparametric$estimate[,"model0"]
OUT$twostagesemiparametric$estimate[c("theta0"),"model1"] <- exp(res1$estimate)
OUT$twostagesemiparametric$LLH["model1"] <- LL1


#=====================================#
# One-stage parametric estimation #
#=====================================#

# MODEL 3: Parent copula with two different child copulas
# #######################################################

likelihood.1st.model3 <- function(p,exp=FALSE){

  if (exp==TRUE){

    lambda1 <- exp(p[1])
    rho1    <- exp(p[2])
    lambda2 <- exp(p[3])
    rho2    <- exp(p[4])
    lambda3 <- exp(p[5])
    rho3    <- exp(p[6])
    lambda4 <- exp(p[7])
    rho4    <- exp(p[8])
    beta    <- p[9] #beta can be negative
    theta1  <- exp(p[10])
    theta2  <- exp(p[11])
    theta3  <- exp(p[12])
  }

  else{  #default

    lambda1 <- p[1]
    rho1    <- p[2]
    lambda2 <- p[3]
    rho2    <- p[4]
    lambda3 <- p[5]
    rho3    <- p[6]
    lambda4 <- p[7]
    rho4    <- p[8]
    beta    <- p[9] #beta can be negative
    theta1  <- p[10]
    theta2  <- p[11]
    theta3  <- p[12]
  }

  S1 <- exp(-lambda1*t1^rho1*exp(beta*Heif1));
  S2 <- exp(-lambda2*t2^rho2*exp(beta*Heif2));
  S3 <- exp(-lambda3*t3^rho3*exp(beta*Heif3));
  S4 <- exp(-lambda4*t4^rho4*exp(beta*Heif4));

  DS1 <- -lambda1*rho1*t1^(rho1-1)*exp(beta*Heif1)*S1;
  DS2 <- -lambda2*rho2*t2^(rho2-1)*exp(beta*Heif2)*S2;
  DS3 <- -lambda3*rho3*t3^(rho3-1)*exp(beta*Heif3)*S3;
  DS4 <- -lambda4*rho4*t4^(rho4-1)*exp(beta*Heif4)*S4;

  A <- -1+(-1+S1^(-theta2)+S2^(-theta2))^(theta1/theta2)+(-1+S3^(-theta3)+S4^(-theta3))^(theta1/theta3)
  B12 <- -1+S1^(-theta2)+S2^(-theta2)
  B34 <- -1+S3^(-theta3)+S4^(-theta3)
  C1 <- S1^(-theta2-1)*DS1
  C2 <- S2^(-theta2-1)*DS2
  C3 <- S3^(-theta3-1)*DS3
  C4 <- S4^(-theta3-1)*DS4

  #joint survival function
  S <- A^(-1/theta1)

  #first order partial derivatives
  dS1 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C1
  dS2 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C2
  dS3 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C3
  dS4 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C4

  #second order partial derivatives
  d2S12 <- A^(-1/theta1-2)*B12^(theta1/theta2-2)*C1*C2*((1+theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d2S13 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C3
  d2S14 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C4
  d2S23 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C3
  d2S24 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C4
  d2S34 <- A^(-1/theta1-2)*B34^(theta1/theta3-2)*C3*C4*((1+theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #third order partial derivatives
  d3S123 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C3*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S124 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C4*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S134 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C1*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)
  d3S234 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C2*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #fourth order partial derivatives
  d4S1234 <- (1+theta1)*A^(-1/theta1-4)*B12^(theta1/theta2-2)*B34^(theta1/theta3-2)*C1*C2*C3*C4*
    ((1+2*theta1)*(1+3*theta1)*B12^(theta1/theta2)*B34^(theta1/theta3)+(1+2*theta1)*A*((-theta1+theta3)*B12^(theta1/theta2)+(-theta1+theta2)*B34^(theta1/theta3))+(-theta1+theta2)*(-theta1+theta3)*A^2)

  terms<-log(S^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
               (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
               (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
               (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
               (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
               d2S12^(c1*c2*(1-c3)*(1-c4))*
               d2S13^(c1*(1-c2)*c3*(1-c4))*
               d2S14^(c1*(1-c2)*(1-c3)*c4)*
               d2S23^((1-c1)*c2*c3*(1-c4))*
               d2S24^((1-c1)*c2*(1-c3)*c4)*
               d2S34^((1-c1)*(1-c2)*c3*c4)*
               (-d3S123)^(c1*c2*c3*(1-c4))*
               (-d3S124)^(c1*c2*(1-c3)*c4)*
               (-d3S134)^(c1*(1-c2)*c3*c4)*
               (-d3S234)^((1-c1)*c2*c3*c4)*
               d4S1234^(c1*c2*c3*c4))
  -sum(terms)}

res3 <- nlm(likelihood.1st.model3,c(log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),0.5,log(0.5),log(0.5),log(0.5)),exp=TRUE)
LL3 <- -(res3$minimum)
# LL3 #-3962.918

# WRITING OUTPUT AWAY:
OUT$onestageparametric$estimate[,"model3"] <- c(exp(res3$estimate[c(1,3,5,7,2,4,6,8)]),res3$estimate[9],exp(res3$estimate[10:12]))
OUT$onestageparametric$LLH["model3"] <- LL3



# MODEL 2: Parent copula with two identical child copulas
# #######################################################

likelihood.1st.model2 <- function(p,exp=FALSE){

  if (exp==TRUE){

    lambda1 <- exp(p[1])
    rho1    <- exp(p[2])
    lambda2 <- exp(p[3])
    rho2    <- exp(p[4])
    lambda3 <- exp(p[5])
    rho3    <- exp(p[6])
    lambda4 <- exp(p[7])
    rho4    <- exp(p[8])
    beta    <- p[9] #beta can be negative
    theta1  <- exp(p[10])
    theta2  <- exp(p[11])
    theta3  <- theta2
  }

  else{  #default

    lambda1 <- p[1]
    rho1    <- p[2]
    lambda2 <- p[3]
    rho2    <- p[4]
    lambda3 <- p[5]
    rho3    <- p[6]
    lambda4 <- p[7]
    rho4    <- p[8]
    beta    <- p[9] #beta can be negative
    theta1  <- p[10]
    theta2  <- p[11]
    theta3  <- theta2
  }

  S1 <- exp(-lambda1*t1^rho1*exp(beta*Heif1));
  S2 <- exp(-lambda2*t2^rho2*exp(beta*Heif2));
  S3 <- exp(-lambda3*t3^rho3*exp(beta*Heif3));
  S4 <- exp(-lambda4*t4^rho4*exp(beta*Heif4));

  DS1 <- -lambda1*rho1*t1^(rho1-1)*exp(beta*Heif1)*S1;
  DS2 <- -lambda2*rho2*t2^(rho2-1)*exp(beta*Heif2)*S2;
  DS3 <- -lambda3*rho3*t3^(rho3-1)*exp(beta*Heif3)*S3;
  DS4 <- -lambda4*rho4*t4^(rho4-1)*exp(beta*Heif4)*S4;

  A <- -1+(-1+S1^(-theta2)+S2^(-theta2))^(theta1/theta2)+(-1+S3^(-theta3)+S4^(-theta3))^(theta1/theta3)
  B12 <- -1+S1^(-theta2)+S2^(-theta2)
  B34 <- -1+S3^(-theta3)+S4^(-theta3)
  C1 <- S1^(-theta2-1)*DS1
  C2 <- S2^(-theta2-1)*DS2
  C3 <- S3^(-theta3-1)*DS3
  C4 <- S4^(-theta3-1)*DS4

  #joint survival function
  S <- A^(-1/theta1)

  #first order partial derivatives
  dS1 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C1
  dS2 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C2
  dS3 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C3
  dS4 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C4

  #second order partial derivatives
  d2S12 <- A^(-1/theta1-2)*B12^(theta1/theta2-2)*C1*C2*((1+theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d2S13 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C3
  d2S14 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C4
  d2S23 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C3
  d2S24 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C4
  d2S34 <- A^(-1/theta1-2)*B34^(theta1/theta3-2)*C3*C4*((1+theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #third order partial derivatives
  d3S123 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C3*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S124 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C4*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S134 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C1*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)
  d3S234 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C2*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #fourth order partial derivatives
  d4S1234 <- (1+theta1)*A^(-1/theta1-4)*B12^(theta1/theta2-2)*B34^(theta1/theta3-2)*C1*C2*C3*C4*
    ((1+2*theta1)*(1+3*theta1)*B12^(theta1/theta2)*B34^(theta1/theta3)+(1+2*theta1)*A*((-theta1+theta3)*B12^(theta1/theta2)+(-theta1+theta2)*B34^(theta1/theta3))+(-theta1+theta2)*(-theta1+theta3)*A^2)

  terms<-log(S^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
               (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
               (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
               (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
               (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
               d2S12^(c1*c2*(1-c3)*(1-c4))*
               d2S13^(c1*(1-c2)*c3*(1-c4))*
               d2S14^(c1*(1-c2)*(1-c3)*c4)*
               d2S23^((1-c1)*c2*c3*(1-c4))*
               d2S24^((1-c1)*c2*(1-c3)*c4)*
               d2S34^((1-c1)*(1-c2)*c3*c4)*
               (-d3S123)^(c1*c2*c3*(1-c4))*
               (-d3S124)^(c1*c2*(1-c3)*c4)*
               (-d3S134)^(c1*(1-c2)*c3*c4)*
               (-d3S234)^((1-c1)*c2*c3*c4)*
               d4S1234^(c1*c2*c3*c4))
  -sum(terms)}

res2 <- nlm(likelihood.1st.model2,c(log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),0.5,log(0.5),log(0.5)),exp=TRUE)
LL2 <- -(res2$minimum)
# LL2 #-3962.946

# WRITING OUTPUT AWAY:
OUT$onestageparametric$estimate[,"model2"] <- c(exp(res2$estimate[c(1,3,5,7,2,4,6,8)]),res2$estimate[9],exp(res2$estimate[10:11]),NA)
OUT$onestageparametric$LLH["model2"] <- LL2

# MODEL 1: one level of clustering
# ################################

likelihood.1st.model1 <- function(p,exp=FALSE){

  if (exp==TRUE){

    lambda1 <- exp(p[1])
    rho1    <- exp(p[2])
    lambda2 <- exp(p[3])
    rho2    <- exp(p[4])
    lambda3 <- exp(p[5])
    rho3    <- exp(p[6])
    lambda4 <- exp(p[7])
    rho4    <- exp(p[8])
    beta    <- p[9] #beta can be negative
    theta1  <- exp(p[10])
    theta2  <- theta1
    theta3  <- theta1
  }

  else{  #default

    lambda1 <- p[1]
    rho1    <- p[2]
    lambda2 <- p[3]
    rho2    <- p[4]
    lambda3 <- p[5]
    rho3    <- p[6]
    lambda4 <- p[7]
    rho4    <- p[8]
    beta    <- p[9] #beta can be negative
    theta1  <- p[10]
    theta2  <- theta1
    theta3  <- theta1
  }

  S1 <- exp(-lambda1*t1^rho1*exp(beta*Heif1));
  S2 <- exp(-lambda2*t2^rho2*exp(beta*Heif2));
  S3 <- exp(-lambda3*t3^rho3*exp(beta*Heif3));
  S4 <- exp(-lambda4*t4^rho4*exp(beta*Heif4));

  DS1 <- -lambda1*rho1*t1^(rho1-1)*exp(beta*Heif1)*S1;
  DS2 <- -lambda2*rho2*t2^(rho2-1)*exp(beta*Heif2)*S2;
  DS3 <- -lambda3*rho3*t3^(rho3-1)*exp(beta*Heif3)*S3;
  DS4 <- -lambda4*rho4*t4^(rho4-1)*exp(beta*Heif4)*S4;

  A <- -1+(-1+S1^(-theta2)+S2^(-theta2))^(theta1/theta2)+(-1+S3^(-theta3)+S4^(-theta3))^(theta1/theta3)
  B12 <- -1+S1^(-theta2)+S2^(-theta2)
  B34 <- -1+S3^(-theta3)+S4^(-theta3)
  C1 <- S1^(-theta2-1)*DS1
  C2 <- S2^(-theta2-1)*DS2
  C3 <- S3^(-theta3-1)*DS3
  C4 <- S4^(-theta3-1)*DS4

  #joint survival function
  S <- A^(-1/theta1)

  #first order partial derivatives
  dS1 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C1
  dS2 <- A^(-1/theta1-1)*B12^(theta1/theta2-1)*C2
  dS3 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C3
  dS4 <- A^(-1/theta1-1)*B34^(theta1/theta3-1)*C4

  #second order partial derivatives
  d2S12 <- A^(-1/theta1-2)*B12^(theta1/theta2-2)*C1*C2*((1+theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d2S13 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C3
  d2S14 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C1*C4
  d2S23 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C3
  d2S24 <- (1+theta1)*A^(-1/theta1-2)*B12^(theta1/theta2-1)*B34^(theta1/theta3-1)*C2*C4
  d2S34 <- A^(-1/theta1-2)*B34^(theta1/theta3-2)*C3*C4*((1+theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #third order partial derivatives
  d3S123 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C3*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S124 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-2)*B34^(theta1/theta3-1)*C1*C2*C4*((1+2*theta1)*B12^(theta1/theta2)+(-theta1+theta2)*A)
  d3S134 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C1*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)
  d3S234 <- (1+theta1)*A^(-1/theta1-3)*B12^(theta1/theta2-1)*B34^(theta1/theta3-2)*C2*C3*C4*((1+2*theta1)*B34^(theta1/theta3)+(-theta1+theta3)*A)

  #fourth order partial derivatives
  d4S1234 <- (1+theta1)*A^(-1/theta1-4)*B12^(theta1/theta2-2)*B34^(theta1/theta3-2)*C1*C2*C3*C4*
    ((1+2*theta1)*(1+3*theta1)*B12^(theta1/theta2)*B34^(theta1/theta3)+(1+2*theta1)*A*((-theta1+theta3)*B12^(theta1/theta2)+(-theta1+theta2)*B34^(theta1/theta3))+(-theta1+theta2)*(-theta1+theta3)*A^2)

  terms<-log(S^((1-c1)*(1-c2)*(1-c3)*(1-c4))*
               (-dS1)^(c1*(1-c2)*(1-c3)*(1-c4))*
               (-dS2)^((1-c1)*c2*(1-c3)*(1-c4))*
               (-dS3)^((1-c1)*(1-c2)*c3*(1-c4))*
               (-dS4)^((1-c1)*(1-c2)*(1-c3)*c4)*
               d2S12^(c1*c2*(1-c3)*(1-c4))*
               d2S13^(c1*(1-c2)*c3*(1-c4))*
               d2S14^(c1*(1-c2)*(1-c3)*c4)*
               d2S23^((1-c1)*c2*c3*(1-c4))*
               d2S24^((1-c1)*c2*(1-c3)*c4)*
               d2S34^((1-c1)*(1-c2)*c3*c4)*
               (-d3S123)^(c1*c2*c3*(1-c4))*
               (-d3S124)^(c1*c2*(1-c3)*c4)*
               (-d3S134)^(c1*(1-c2)*c3*c4)*
               (-d3S234)^((1-c1)*c2*c3*c4)*
               d4S1234^(c1*c2*c3*c4))
  -sum(terms)}

res1 <-nlm(likelihood.1st.model1,c(log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),log(0.5),0.5,log(0.5)),exp=TRUE)
LL1 <- -(res1$minimum)
# LL1 #-3970.024

# WRITING OUTPUT AWAY:
OUT$onestageparametric$estimate[,"model1"] <- c(exp(res1$estimate[c(1,3,5,7,2,4,6,8)]),res1$estimate[9],exp(res1$estimate[10]),NA,NA)
OUT$onestageparametric$LLH["model1"] <- LL1


# Add Model 0 to Semi-parametric and parametric model
OUT$twostageparametric$estimate[,"model0"] <- OUT$onestageparametric$estimate[,"model0"]
OUT$twostageparametric$LLH["model0"] <- OUT$onestageparametric$LLH["model0"]
OUT$twostagesemiparametric$estimate[,"model0"] <- OUT$onestageparametric$estimate[,"model0"]
OUT$twostagesemiparametric$LLH["model0"] <- OUT$onestageparametric$LLH["model0"]

# Temporarily Remove stderror
OUT <- lapply(OUT,FUN=function(x){x[-2]})


print(OUT)
