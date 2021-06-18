#' Derives the range of a vector
#'
#' @param x vector with unknown range
#'
#' @return the range of the vector
#'
#' @examples
#' custom_vector = c(5,2,7,9,4)
#' ran(custom_vector)
#'
#' @export
ran<-function(x){
  rangex<-range(x)
  diff<-rangex[2]-rangex[1]

  diff
}

## derive ordering of the items in the list by their range
#' Computes the vector of the highest ranges within the vector and returns it as a vector of decreasing size
#'
#' @param beta_list a two dimensional vector of which the ordered range is of intersest
#'
#' @return vector of ordered ranges
#'
#' @examples
#'
#' @export
ran.list<-function(beta_list){
  ran_beta<-numeric(length(beta_list))
  for(j in 1:length(beta_list)){
    ran_beta[j]<-ran(beta_list[[j]])
  }
  indexvec<-order(ran_beta,decreasing=TRUE)

  indexvec
}

## derive ordering of the items in the list by their median locations
#' Computes the increasing order of items based on median.
#'
#' @param beta_list a two dimensional vector of which the ordered median is of interest
#'
#' @return vector of increasing medians of
#'
#' @examples
#'
#' @export
med.list<-function(beta_list){
  med_beta<-numeric(length(beta_list))
  for(j in 1:length(beta_list)){
    med_beta[j]<-median(beta_list[[j]])
  }
  indexvec<-order(abs(med_beta-median(med_beta)),decreasing=FALSE)

  indexvec
}
