#' Function that computes the distance between the two functions estimated by the samples
#'
#' @param sample1 a vector of probabilties for the support points
#' @param sample2 a vector of probabilties for the support points
#' @param samp noclue
#' @param rancor noclue
#' @param gridsize noclue
#'
#' @return dist L2 - distance between the two functions estimated by the samples
#' @examples
#' #l2dist(s)
#' @export
l2dist<-function(sample1, sample2, samp=FALSE,rancor=0.05,gridsize=NA){
  l2distance<-function(sample1, sample2){
    ############################
    ## input:
    ## sample1    a vector of probabilties for the support points
    ## sample2    a vector of probabilties for the support points
    ### !!! both samples must be evaluated at the same grid support points !!!
    #########################################
    ## output:
    ## dist   L2 - distance between the two functions estimated by the samples
    l2<-sqrt(sum((sample1-sample2)^2/sample2)/length(sample1))
    return(l2)
  }

  if(samp){
  ran1<-range(sample1)
  ran2<-range(sample2)
  dens1<-density((sample1),from=(min(ran2[1],ran1[1])-rancor),to=(max(ran2[2],ran1[2])+rancor),kernel="epanechnikov",n=min(length(sample1),2^15),adjust=2.34/0.9)
  dens2<-density((sample2),from=(min(ran2[1],ran1[1])-rancor),to=(max(ran2[2],ran1[2])+rancor),kernel="epanechnikov",n=min(length(sample1),2^15),adjust=2.34/0.9)
  ## identify occurrences of bothsamples in density of sample 1
  val1in2.ind<-findInterval(c(sample1,sample2),signif(dens1$x,digits=5))
    val1in2.y<-dens1$y[val1in2.ind] #posterior21.y
  ## identify occurrences of both samples in density of sample 2
  val2in1.ind<-findInterval(c(sample1,sample2),signif(dens2$x,digits=5))
  val2in1.y<-dens2$y[val2in1.ind] #posterior12.y

  distance<-l2distance(val2in1.y,val1in2.y)
  }else{
    #dens1<-density((sample1),from=min(sample1-rancor),to=max(sample1+rancor),kernel="epanechnikov",n=min(length(sample1),2^15),adjust=2.34/0.9)
    ## sample2 are the exact density values at positions sample1
    #val1in2.ind<-findInterval(sample1,signif(dens1$x,digits=5))
    #val1in2.y<-dens1$y[val1in2.ind] #posterior21.y
    if(is.na(gridsize)) {dens1<-ks:::kdde(sample1,eval.points=sample2[[1]])}else{dens1<-ks:::kdde(sample1,eval.points=sample2[[1]],gridsize = gridsize)}
    distance<-l2distance(dens1$estimate,sample2[[2]])  #val1in2.y
  }

  return(distance)
}
