#' @title Udder Quarter Infection Data
#'
#' @description This package is built up around a data set consisting of clustered infection times in the four udder quarters of cows.
#' Due to the special structure of these data, i.e., quadruples of four clustered event times, thye have been used to develop
#' new techniques in the area of multivariate survival data (see references). This package avails the data on which the analyses
#' presented in the papers were based, and also, where possible, provides the R-code to fit the proposed models.
#'
#'
#' @author Hans Laevens
#' @author Luc Duchateau
#' @author Klara Goethals
#' @author Ewoud de Troyer
#' @author Paul Janssen
#'
#' @references Laevens, H., Deluyker, H., Schukken, Y.H., De Meulemeester, L., Vandermeersch, R., De Muelenaere, E. and De Kruif, A. (1997). Influence of parity and stage of lactation on the somatic cell count in bacteriologically negative dairy cows. Journal of Dairy Science 80, 3219-3226.
#' @references Massonnet, G., Janssen, P. and Duchateau, L. (2009). Modelling udder infection data using copula models for quadruples. Journal of Statistical Planning and Inference 139, 3865-3877.
#' @references Goethals, K., Ampe, B., Berkvens, D., Laevens, H., Janssen, P. and Duchateau, L. (2009). Modeling interval-censored, clustered cow udder quarter infection times through the shared gamma frailty model. Journal of Agricultural Biological and Environmental Statistics 14, 1-14.
#' @references Janssen, P. and Duchateau, L. (2011). Comments on: Inference in multivariate Archimedean copula models. Test 20, 271-275.
#' @references Ampe, B., Goethals, K., Laevens, H. and Duchateau, L. (2012). Investigating clustering in interval-censored udder quarter infection times in dairy cows using a gamma frailty model. Preventive Veterinary Medicine 106, 251-257.
#' @references Kuhn, E., Goethals, K., El-Nouty, C. and Duchateau, L. (2016). Assessing the correlation structure in cow udder quarter infection times through extensions of the correlated frailty model. Journal of Agricultural Biological and Environmental Statistics 21, 601-618.
#' @references Geerdens, C., Claeskens, G., Janssen, P. (2016). Copula based flexible modeling of associations between clustered event times. Lifetime data analysis 22, 363-381.
#' @references Verbelen, R., Antonio, K. and Claeskens, G. (2016). Multivariate mixtures of Erlangs for density estimation under censoring. Lifetime data analysis 22, 429-455.
#' @references Prenen, L., Duchateau, L. and Braekers, R. (2018). Investigating the correlation structure of quadrivariate udder infection times through hierarchical Archimedean copulas. Accepted for publication in Lifetime Data Analysis.
#'
#' @docType package
#' @name UdderQuarterInfectionData
NULL



## UDDERQUARTER DATA ##
#' @title Udder Quarter Infection Data
#'
#' @description
#' The udder quarter infection data contains the times to infection of the individual udder quarters of 100 cows during one lactation period.
#' The 'cowid' variable contains the cow number. Each cow has 4 different quarters, with the variable 'quarter' referring to the specific quarter,
#' i.e., RF: right front,  RR: right rear, LF: left front,  LR: left rear.
#' The lactation number, i.e., the number of calvings that the cow has experienced, is given by the variable 'lactation'.
#' If at least one infection occurs in a quarter during the lactation period, the 'status' variable is set to 1,
#' and the interval [t1,t2] in which an infection takes place is given by the variables 'left' and 'right',
#' where 'right' denotes the time with the first positive result, and 'left' the last time (with a negative result) before the first positive result.
#' Only the first infection is contained in the data set.
#' If no infection occurs, the 'status' variable is set to 1, and the last interval [t1,t2] in which the quarter was observed is put in the variables 'left' and 'right'.,

#' @format A dataframe with 4784 rows and 6 variables.
#' @name udderquarterinfection
#' @examples
#' \dontrun{
#'
#' ## SHORT DATASET ANALYSIS ##
#'
#' library(UdderQuarterInfectionData)
#' data("udderquarterinfection")
#'
#' upper=udderquarterinfection$right/91.31 # Divide by 91.31 to avoid small lambdas
#' lower=udderquarterinfection$left/91.31
#' fail=udderquarterinfection$status
#' cluster=udderquarterinfection$cowid
#' X=udderquarterinfection$lactation>1
#'
#'
#' #Calculate the number of clusters.
#' clusternames <- levels(as.factor(cluster))
#' ncluster <- length(clusternames)
#'
#' # Create a udderquarterinfection with the variables cluster, the lower bound, the upper bound,
#' # the censoring indicator and the covariate.
#' udderquarterinfectionint <- as.matrix(cbind(cluster,lower,upper,fail,X))
#'
#' # create subsets for right-censored and interval-censored observations
#' cendata <- udderquarterinfectionint[udderquarterinfectionint[,4]==0,]
#' intdata <- udderquarterinfectionint[udderquarterinfectionint[,4]==1,]
#'
#' # Create a list of signs that corresponds to the n_ik (here restricted to 4 events)
#' signs <- list(1,c(1,-1))
#' for(i in 3:5) signs[[i]] <- kronecker(signs[[i-1]],c(1,-1))
#'
#' # Function to calculate the loglikelihood per cluster
#' CalcLogLikClust <- function(i,x)
#' {
#'   theta <- x[1]; lambda <- x[2]; gamma <- x[3]; beta <- x[4]
#'   cenX <- cendata[cendata[,1]==clusternames[i],5]
#'   intL <- intdata[intdata[,1]==clusternames[i],2]
#'   nevents <- length(intL)
#'   crossprod <- 1
#'   if(nevents>0){ #comment: if there are no events, crosspod=1
#'     intX <- intdata[intdata[,1]==clusternames[i],5]
#'     # Calculate R*_ij and L*_ij
#'     intRster <- lambda*(intdata[intdata[,1]==clusternames[i],3]^gamma)*exp(intX*beta)
#'     intLster <- lambda*(intL^gamma)*exp(intX*beta)
#'     # Calculate the vector p_i
#'     crossprod <- c(exp(intLster[1]),exp(intRster[1]))
#'     if(nevents>1){
#'       for(ik in 2:nevents) crossprod<-kronecker(crossprod,c(exp(intLster[ik]),exp(intRster[ik])))
#'     }
#'   }
#'   # Loglikelihood for 1 cluster
#'   return(
#'   log(1/(theta^(1/theta))*sum((1/
#'   ((sum(lambda*(as.vector(cendata[cendata[,1]==clusternames[i],3])^gamma)
#'   *exp(cenX*beta))+1/theta+log(crossprod))^(1/theta)))*signs[[nevents+1]]))
#'   )
#' }
#'
#' # Calculate full marginal loglikelihood (formula 5)
#' CalcLogLik <- function(x)
#' {
#'   -sum(sapply(1:ncluster,CalcLogLikClust,x=x))
#' }
#'
#' # Maximising the full marginal loglikelihood to obtain parameter estimates
#' init <- c(1,1,1,1)
#' print(results <- nlm(CalcLogLik,init,print.level=2, hessian=TRUE)) # Can take a while!
#' # $minimum
#' # [1] 5670.491
#' #
#' # $estimate
#' # [1] 3.7967246 0.1201593 1.9672298 0.8590531
#' #
#' # $gradient
#' # [1]  0.0002924871  0.0017653292 -0.0005460029  0.0003265086
#' #
#' # $hessian
#' # [,1]       [,2]      [,3]       [,4]
#' # [1,]   23.22965  -117.7682 -39.93813  -10.10561
#' # [2,] -117.76825 15471.4753 567.24283 1228.87332
#' # [3,]  -39.93813   567.2428 664.76359   24.63047
#' # [4,]  -10.10561  1228.8733  24.63047  147.76479
#' #
#' # $code
#' # [1] 1
#' #
#' # $iterations
#' # [1] 22
#'
#' # Calculate covariance matrix
#' covmatr <- solve(results$hessian)
#' #             [,1]          [,2]          [,3]         [,4]
#' # [1,] 0.049281911  0.0001242730  0.0027853686  0.001872592
#' # [2,] 0.000124273  0.0001982213 -0.0001015391 -0.001623066
#' # [3,] 0.002785369 -0.0001015391  0.0017306214  0.000746460
#' # [4,] 0.001872592 -0.0016230660  0.0007464600  0.020269244
#'
#' }
NULL



