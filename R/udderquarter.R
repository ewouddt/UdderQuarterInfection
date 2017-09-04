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
#'
#' The following example analyses of this data set are currently available in the package:
#' \itemize{
#' \item \code{\link{Gamma_Frailty_Interval_Censoring}}
#' }
#'
#' @format A dataframe with 4784 rows and 6 variables.
#' @name udderquarterinfection
#' @examples
#' \dontrun{
#' ## SHORT DATASET ANALYSIS ##
#'
#' library(UdderQuarterInfectionData)
#' data("udderquarterinfection")
#'
#' Gamma_Frailty_Interval_Censoring()
#' # $nlm
#' # $nlm$minimum
#' # [1] 5670.491
#' #
#' # $nlm$estimate
#' # [1] 3.7967246 0.1201593 1.9672298 0.8590531
#' #
#' # $nlm$gradient
#' # [1]  0.0002924871  0.0017653292 -0.0005460029  0.0003265086
#' #
#' # $nlm$hessian
#' # [,1]       [,2]      [,3]       [,4]
#' # [1,]   23.22965  -117.7682 -39.93813  -10.10561
#' # [2,] -117.76825 15471.4753 567.24283 1228.87332
#' # [3,]  -39.93813   567.2428 664.76359   24.63047
#' # [4,]  -10.10561  1228.8733  24.63047  147.76479
#' #
#' # $nlm$code
#' # [1] 1
#' #
#' # $nlm$iterations
#' # [1] 22
#' #
#' #
#' # $covmat
#' # [,1]          [,2]          [,3]         [,4]
#' # [1,] 0.049281911  0.0001242730  0.0027853686  0.001872592
#' # [2,] 0.000124273  0.0001982213 -0.0001015391 -0.001623066
#' # [3,] 0.002785369 -0.0001015391  0.0017306214  0.000746460
#' # [4,] 0.001872592 -0.0016230660  0.0007464600  0.020269244
#' }
NULL



