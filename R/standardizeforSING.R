#' Standardization with double centered and column scaling
#'
#' @param data input matrix with n x px.
#' @param dif.tol the value for the threshold of scaling
#' @param max.iter default value = 10
#'
#' @return standardized matrix with n x px.
#' @export
#'
#' @examples
#' spmwm = 3*matrix(rnorm(100000),nrow=100)+1
#' dim(spmwm)
#' apply(spmwm,1,mean) # we want these to be 0
#' apply(spmwm,2,mean) # we want these to be 0
#' apply(spmwm,2,sd) # we want each of these variances to be 1
#'
#' spmwm_cp=standard(spmwm)
#' max(abs(apply(spmwm_cp,1,mean)))
#' max(abs(apply(spmwm_cp,2,mean)))
#' max(abs(apply(spmwm_cp,2,sd)-1))
standard <- function(data,dif.tol=1e-03,max.iter=10){
  row_mean_max = max(abs(apply(data,1,mean)))
  col_mean_max = max(abs(apply(data,2,mean)))
  col_sd_max= max(abs(apply(data,2,sd)-1))
  n=0
  while(n<=max.iter & max(row_mean_max,col_mean_max,col_mean_max)>= dif.tol) {
    data=scale(data) # centering and scaling for each column
    data=t(scale(t(data),center = T,scale = F))
    n=n+1
  }
  return(data)
}



