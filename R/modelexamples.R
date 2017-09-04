globalVariables("udderquarterinfection")


#' @importFrom stats nlm
#' @importFrom utils data
NULL


#'@export
#'@title Gamma Frailty Interval Censoring
#'@description Application of the Gamma Frailty Interval Censoring Model on the Udder Quarter Infection Data Set. For more information see Details.
#'@param print.level Parameter of \code{\link[stats]{nlm}} (default=2): this argument determines the level of printing which is done during the minimization process. The default value of 0 means that no printing occurs, a value of 1 means that initial and final details are printed and a value of 2 means that full tracing information is printed.
#'
#'@details
#'ADD MODEL DETAILS HERE
#'
#'@author ?
#'@references ?
#'
#'@return Returns a list with the NLM result in \code{nlm} and the covariance matrix in \code{covmat}.
#'
#'@section R Code for Model :
#'The source R code for this model can found:
#'\itemize{
#'\item in the \code{doc/Models_R_Code.R} file in the package installation folder.
#'\item by accessing the function by calling \code{Gamma_Frailty_Interval_Censoring} (without brackets) or \code{getAnywhere("Gamma_Frailty_Interval_Censoring")}.
#'}
Gamma_Frailty_Interval_Censoring <- function(print.level=2){
  data("udderquarterinfection",envir=environment())

  upper=udderquarterinfection$right/91.31 # Divide by 91.31 to avoid small lambdas
  lower=udderquarterinfection$left/91.31
  fail=udderquarterinfection$status
  cluster=udderquarterinfection$cowid
  X=udderquarterinfection$lactation>1


  #Calculate the number of clusters.
  clusternames <- levels(as.factor(cluster))
  ncluster <- length(clusternames)

  # Create a udderquarterinfection with the variables cluster, the lower bound, the upper bound,
  # the censoring indicator and the covariate.
  udderquarterinfectionint <- as.matrix(cbind(cluster,lower,upper,fail,X))

  # create subsets for right-censored and interval-censored observations
  cendata <- udderquarterinfectionint[udderquarterinfectionint[,4]==0,]
  intdata <- udderquarterinfectionint[udderquarterinfectionint[,4]==1,]

  # Create a list of signs that corresponds to the n_ik (here restricted to 4 events)
  signs <- list(1,c(1,-1))
  for(i in 3:5) signs[[i]] <- kronecker(signs[[i-1]],c(1,-1))

  # Function to calculate the loglikelihood per cluster
  CalcLogLikClust <- function(i,x)
  {
    theta <- x[1]; lambda <- x[2]; gamma <- x[3]; beta <- x[4]
    cenX <- cendata[cendata[,1]==clusternames[i],5]
    intL <- intdata[intdata[,1]==clusternames[i],2]
    nevents <- length(intL)
    crossprod <- 1
    if(nevents>0){ #comment: if there are no events, crosspod=1
      intX <- intdata[intdata[,1]==clusternames[i],5]
      # Calculate R*_ij and L*_ij
      intRster <- lambda*(intdata[intdata[,1]==clusternames[i],3]^gamma)*exp(intX*beta)
      intLster <- lambda*(intL^gamma)*exp(intX*beta)
      # Calculate the vector p_i
      crossprod <- c(exp(intLster[1]),exp(intRster[1]))
      if(nevents>1){
        for(ik in 2:nevents) crossprod<-kronecker(crossprod,c(exp(intLster[ik]),exp(intRster[ik])))
      }
    }
    # Loglikelihood for 1 cluster
    return(
      log(1/(theta^(1/theta))*sum((1/
                                     ((sum(lambda*(as.vector(cendata[cendata[,1]==clusternames[i],3])^gamma)
                                           *exp(cenX*beta))+1/theta+log(crossprod))^(1/theta)))*signs[[nevents+1]]))
    )
  }

  # Calculate full marginal loglikelihood (formula 5)
  CalcLogLik <- function(x)
  {
    -sum(sapply(1:ncluster,CalcLogLikClust,x=x))
  }

  # Maximising the full marginal loglikelihood to obtain parameter estimates
  init <- c(1,1,1,1)
  results <- nlm(CalcLogLik,init,print.level=print.level, hessian=TRUE) # Can take a while!
  # $minimum
  # [1] 5670.491
  #
  # $estimate
  # [1] 3.7967246 0.1201593 1.9672298 0.8590531
  #
  # $gradient
  # [1]  0.0002924871  0.0017653292 -0.0005460029  0.0003265086
  #
  # $hessian
  # [,1]       [,2]      [,3]       [,4]
  # [1,]   23.22965  -117.7682 -39.93813  -10.10561
  # [2,] -117.76825 15471.4753 567.24283 1228.87332
  # [3,]  -39.93813   567.2428 664.76359   24.63047
  # [4,]  -10.10561  1228.8733  24.63047  147.76479
  #
  # $code
  # [1] 1
  #
  # $iterations
  # [1] 22

  # Calculate covariance matrix
  covmatr <- solve(results$hessian)
  #             [,1]          [,2]          [,3]         [,4]
  # [1,] 0.049281911  0.0001242730  0.0027853686  0.001872592
  # [2,] 0.000124273  0.0001982213 -0.0001015391 -0.001623066
  # [3,] 0.002785369 -0.0001015391  0.0017306214  0.000746460
  # [4,] 0.001872592 -0.0016230660  0.0007464600  0.020269244


  return(list(
    nlm=results,
    covmat=covmatr
  ))
}
