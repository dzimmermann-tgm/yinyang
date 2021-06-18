#' A method that does stuff
#'
#' @param sample_list parameter
#' @param prior_list parameter
#' @param algo parameter
#' @param method parameter
#' @param prob_old parameter
#' @param burn.in parameter
#' @param bw.adj parameter
#' @param rancor paramete
#' @param bw.adj parameter
#' @param rancor parameter
#' @param neps parameter
#' @param check parameter
#' @param nsam parameter
#'
#' @return
#'
#' @examples
#'
#' @export
yin.yang.sampling<-function(sample_list=NULL,prior_list=NULL,algo="multinom",method='sequential',prob_old="heuristic",burn.in=2500,bw.adj=1,rancor=1,neps=500,check=FALSE,nsam=NULL){
  ## Metropolis Hastings resampling scheme
  ##########################################
  ## sample_list ... list of subset's posterior samples
  ##                  each element must be a numeric vector of samples from the subset's posterior
  ## prior_list  ... list of prior values of the samples' draws
  ## algo        ... variant of the algorithm
  ##                 "multinom"  simple Yin Yang algorithm
  ##                 "MH"  MH Yin Yang algorithm
  ## method      ... method for caluculating the yin yang resampling;
  ##                 possible values: 'tree' if number of sample = 2^n;
  ##                 'sequential' for sequential resampling evaluation
  ## prob_old    ... probability to draw from the 'old' yin sample;
  ##                 if numeric between 0 and 1: (1-prob_old) = probability to draw from the 'new' yang sample
  ##                  either numeric value of probability to draw from the old sample
  ##                   or method to determine probability to draw from the 'old' yin sample;
  ##                 the probability is calculated based on the difference of the variances,
  ##                   "heuristic"   based on samples' variances  (deafult)
  ##                   "exact"   based on simulations to recover exact probabilities
  ##                 (1-prob_old) = probability to draw from the 'new' yang sample
  ## burn.in     ... length of burn-in to initially cut off every MH resampling run and the input samples
  ##                  caveat: add 'artificial' burn-in to inputs, if you correct manually in advance
  ## bw.adj      ... manual adjustment to kernal density estimator's bandwidth
  ## rancor      ... range correction to make range slightly large and density estimation positive value on the outer rims
  ## neps        ... number of steps for approximating the optimal weights
  ## check       ... indicator whether one wants to check the merging behaviour, then summaries for each step of merging are displayed
  ## nsam        ... number of samples to produce
  ######################################################################
  ## Output:
  ## list(merged=sam_list_new,acceptance=acceptance)
  ## merged ... either list of merged output samples ('sequential')
  ##            or final merged sample ('tree')
  ## acceptance ... list of acceptance probabilities (for MH Yin Yang)
  ## ###################################################################
  if(is.null(sample_list)){
    stop(gettextf("please provide list of samples to combine", domain = NA))
  }
  if(is.null(prior_list)){
    stop(gettextf("please provide list of prior values", domain = NA))
  }
  nsample<-length(sample_list)
  nprior<-length(prior_list)
  if(nsample!=nprior){
    stop(gettextf("length of sample list %d is not equal length of prior list %d",nsample,nprior, domain = NA))
  }
  if(is.numeric(prob_old)){
  if(prob_old<=0||prob_old>=1){
      stop(gettextf("probability for yin sample %d must lie between 0 and 1",prob_old, domain = NA))}else{
        prob_old=max(min(prob_old,0.95),0.05)  # to guarantee switching between the samples
      }
}
if(is.null(nsam)){
  nsam<-length(sample_list[[1]])
}
  ## MedMed
medmed<-function(x){
    median(abs(x-median(x)))
  }

  ## Metropolis Hastings yin yang sampler
resampleMH<-function(sample1,sample2,prior1, prior2,prob_old,method,weight=NULL){
    sample1<-sample1[!is.na(sample1)]
    sample2<-sample2[!is.na(sample2)]
    # check lengths and set to equal length
    len_ss<-min(length(sample1),length(sample2))
    sample1<-sample1[1:len_ss]
    sample2<-sample2[1:len_ss]
    ## ranges of the samples
    ran1<-range(sample1)
    ran2<-range(sample2)
    ## estimate densities within the full range of both samples
    ## (min of min to max of maxs)
    dens1<-density((sample1),from=(min(ran2[1],ran1[1])-rancor),to=(max(ran2[2],ran1[2])+rancor),kernel="epanechnikov",n=min(length(sample1),2^15),adjust=bw.adj*2.34/0.9)
    dens2<-density((sample2),from=(min(ran2[1],ran1[1])-rancor),to=(max(ran2[2],ran1[2])+rancor),kernel="epanechnikov",n=min(length(sample1),2^15),adjust=bw.adj*2.34/0.9)
    ## identify occurrences of sample 1 in density of sample 2
    val1in2.ind<-findInterval(sample1,signif(dens2$x,digits=5))
    val1in2.x<-dens2$x[val1in2.ind] #posterior21.x
    val1in2.y<-dens2$y[val1in2.ind] #posterior21.y
    ## identify occurrences of sample 2 in density of sample 1
    val2in1.ind<-findInterval(sample2,signif(dens1$x,digits=5))
    val2in1.x<-dens1$x[val2in1.ind] #posterior12.x
    val2in1.y<-dens1$y[val2in1.ind] #posterior12.y
    ## initialise new sample1 as empty vector
    cons_len<-min(len_ss,length(val2in1.x),length(val2in1.x),length(val2in1.x),length(val2in1.x))
    sam1<-numeric(cons_len+burn.in)
    # write out the new prior for the merged sample
    prior<-numeric(cons_len+burn.in)
    # indicator which sample the last value has been proposed from
    ind.sam.old<-sample(2,1)
    # initialise as sample 1 -> problem with acceptance
    # randomly initialise instead
    # indices which value of the chain to propose from sample 1 or 2
    # go through each chain sequentially
    index.1<-1
    index.2<-1
    ## draw first 'old' value from common part of 1, 2 -> initialisation
    ind.old.1<-1
    ind.old.2<-1
    # initialise new sample and prior
    if(ind.sam.old==1){
      ind.old.1<-2
      sam1[1]<-sample1[ind.old.1] #dens1$x[ind.old.1]
      prior[1]<-prior1[ind.old.1]
      index.1<-ind.old.1
    }else{
      ind.old.2<-2
      sam1[1]<-sample2[ind.old.2] #dens2$x[ind.old.2]
      prior[1]<-prior2[ind.old.2]
      index.2<-ind.old.2
    }
    ##sample dependent acceptance rate
    if(is.numeric(prob_old)){
      prob_old_1<-prob_old
      prob_old_2<-(1-prob_old)
      ###############################
      ##       prob_old_1<-min(prob_old*(var(sample2)/var(sample1)),0.95)
      ##       prob_old_2<-min(prob_old*(var(sample1)/var(sample2)),0.95)
      #       # calculate sample based optimal weights
      #       exp.yin<-mean(sample1)
      #       exp.yang<-mean(sample2)
      #       exp.full<-0.5*exp.yin+0.5*exp.yang
      #       for(i in 1:neps){
      #         g.yin<-sample1*val1in2.y/prior2 # correct for likelihood
      #         g.yang<-sample2*val2in1.y/prior1 # correct for likelihood
      #         # calculate variances
      #         var.yin<-sum((g.yin-exp.full)^2)/cons_len  #cons_len1
      #         var.yang<-sum((g.yang-exp.full)^2)/cons_len  #cons_len2
      #         # calculate weights relative to total variance
      #         prob_old_1<-var.yin/sum(var.yin,var.yang)
      #         prob_old_2<-1-prob_old_1
      #         exp.full<-prob_old_1*exp.yin+prob_old_2*exp.yang
      #       }
    }
    if(prob_old=="heuristic"){
      prob_old_1<-max(0.05,min(0.95,weights_vec[k-1]/(weights_vec[k]+weights_vec[k-1])))
    }
    if(prob_old=="exact"){
      syy.yang=mean(val1in2.y/prior1)
      syy.yin=mean(val2in1.y/prior2)
      w.yang<-(val1in2.y/prior1)
      w.yin<-(val2in1.y/prior2)
      wyy<-data.frame(fir=w.yin[1:(round(len_ss/2))],sec=w.yin[(round(len_ss/2)+1):len_ss])
      Ayinyin<-2/(len_ss*syy.yin)*sum(min(with(wyy,pmin(fir,sec))))
      wyy<-data.frame(fir=w.yang[1:round(len_ss/2)],sec=w.yang[(round(len_ss/2)+1):len_ss])
      Ayangyang<-2/(len_ss*syy.yang)*sum(min(with(wyy,pmin(fir,sec))))
      if((Ayinyin-Ayangyang)< (-0.000005)){prob_old_1=0}else{if((Ayinyin-Ayangyang)> (0.000005)){prob_old_1=1}else{prob_old_1=0.5}}
    }
    prob_old_2<-(1-prob_old_1)
    #cat(prob_old_1,prob_old_2,fill=TRUE)
    #acceptance<-numeric(cons_len+burn.in)
    acceptance<-0
    for(ind in 1:(cons_len+burn.in)){
      # randomise whether to draw from sample 1 or 2
      # increasing the limit above 0.5  (currently 0.7)
      # makes the algorithm more conservative in the sense of
      # staying closer to 'old' sample (sample 1 or the condensed
      # sample of all previous resampling steps)
      if(runif(1)<=prob_old_1){## same sample as before
        ind.sam.new<-1
        index.1<-(index.1+1)
        if(index.1>len_ss){index.1<-index.1-sample(len_ss,size=1)}
        proposal<-sample1[index.1]
        propprior<-prior1[index.1]
      }else{# take other sample for proposal
        ind.sam.new<-2
        index.2<-(index.2+1)
        if(index.2>len_ss){index.2<-index.2-sample(len_ss,size=1)}
        proposal<-sample2[index.2]
        propprior<-prior2[index.2]
      }
      sum.ind<-(ind.sam.old+2*ind.sam.new-2)
      ## switch between 4 cases
      switch(sum.ind,
{# old=new=1
  logpost<-log(val1in2.y[index.1])-log(val1in2.y[ind.old.1])
  logprio<-log(prior1[ind.old.1])-log(propprior)
  accprob<-logpost+logprio
},
{ # old=2,new=1
  logpost<-log(val1in2.y[index.1])-log(val2in1.y[ind.old.2])
  logprio<-log(prior2[ind.old.2])-log(propprior)
  logind<-log(prob_old_1)-log(prob_old_2)
  accprob<-logpost+logprio#+logind
},
{ # old=1,new=2
  logpost<-log(val2in1.y[index.2])-log(val1in2.y[ind.old.1])
  logprio<-log(prior1[ind.old.1])-log(propprior)
  logind<-log(prob_old_2)-log(prob_old_1)
  accprob<-logpost+logprio#+logind
},
{ # old=2,new=2
  logpost<-log(val2in1.y[index.2])-log(val2in1.y[ind.old.2])
  logprio<-log(prior2[ind.old.2])-log(propprior)
  accprob<-logpost+logprio
}
      )
####################################
## Metropolis Hastings acceptance/rejection
if(is.na(accprob)){accprob<-(-Inf)}
accprob<-min(accprob,0)
#acceptance[ind]<-accprob
if(log(runif(1))<accprob){  ## accept new value
  sam1[ind]<-proposal
  prior[ind]<-propprior
  acceptance<-acceptance+1
  # update index
  if(ind.sam.new==1){ind.old.1<-index.1}else{
    ind.old.2<-index.2}
  ind.sam.old<-ind.sam.new
}else{ ## keep old value
  if(ind>1){
      ## all but first sample
      sam1[ind]<-sam1[ind-1]
      prior[ind]<-prior[ind-1]
    } # for first sample it has already been initialised
  }
    }
sample.out<-(sam1[!is.na(sam1)])[(burn.in+1):(cons_len+burn.in)]

#if(method=="tree"){
    return(list(sample.out,acceptance,prior))
 # }
#if(method=="sequential"){
#    return(list(sample.out,acceptance))}
}

resampleyinyang<-function(sample1,sample2,prior1,prior2,prob_old,method,weight=NULL){
    len_sam<-min(length(sample1),length(sample2))
    sample1<-sample1[1:len_sam]
    sample2<-sample2[1:len_sam]
    prior1<-prior1[1:len_sam]
    prior2<-prior2[1:len_sam]
    # calculate sample ranges
    # common support necessary for the method
    ran1<-range(sample1)  # max, min
    ran2<-range(sample2)
    # posterior density of sample 1 on the support of sample 1
    # if the common support is small, then there will be many 0 weights
    #post1<-density((sample1),from=(min(ran2[1],ran1[1])-rancor),to=(max(ran2[2],ran1[2])+rancor),kernel="epanechnikov",n=min(length(sample1),2^15),bw=bw.nrd(sample2),adjust=bw.adj*2.34/0.9)
    #post2<-density((sample2),from=(min(ran2[1],ran1[1])-rancor),to=(max(ran2[2],ran1[2])+rancor),kernel="epanechnikov",n=min(length(sample1),2^15),bw=bw.nrd(sample1),adjust=bw.adj*2.34/0.9)
    post2<-density(sample2,n=min(length(sample1),2^15),from=ran1[1]-rancor,to=ran1[2]+rancor,kernel="epanechnikov",adjust=bw.adj*2.34/0.9,na.rm = TRUE)
    post1<-density(sample1,n=min(length(sample1),2^15),from=ran2[1]-rancor,to=ran2[2]+rancor,kernel="epanechnikov",adjust=bw.adj*2.34/0.9,na.rm = TRUE)
    # fit to values of sample 1
    ind2<-findInterval(sample1,post2$x)
    posterior21.y<-post2$y[ind2]
    posterior21.x<-post2$x[ind2]
    ind1<-findInterval(sample2,post1$x)
    posterior12.y<-post1$y[ind1]
    posterior12.x<-post1$x[ind1]
    ## calculate the weights
    # log ratio of posterior density of sample 2 on support of
    # sample 1 over the prior on sample 2
    lw_star2<-log(posterior21.y)-log(prior2)
    lw_star1<-log(posterior12.y)-log(prior1)
    # cut off values too small or large to get rid of Inf
    # (division by numerical 0)
    lw_star2[is.infinite(lw_star2)]<-sign(lw_star2[is.infinite(lw_star2)])*20
    lw_star1[is.infinite(lw_star1)]<-sign(lw_star1[is.infinite(lw_star1)])*20
    nlw_star2<-pmax(pmin((lw_star2),20),-20)  #-max(lw_star2)
    nlw_star1<-pmax(pmin((lw_star1),20),-20)  #-max(lw_star1)
    # new weights are calculated from standardised log ratios
    # (sum up to 1)
    weights2<-exp(nlw_star2)/sum(exp(nlw_star2) )
    weights1<-exp(nlw_star1)/sum(exp(nlw_star1) )
    #cat(summary(weights),fill=TRUE)
    ## once we have the weights, we draw from the multinomial distribution
    ## over sample 1 with weight vector weights
    sample1.new<-rmultinom(1,size=(len_sam),prob=weights2)
    sample2.new<-rmultinom(1,size=(len_sam),prob=weights1)
    # !!! only values present in the first sample will ever be drawn !!!
    new.sample1<-numeric()
    new.sample2<-numeric()
    new.prior1<-numeric()
    new.prior2<-numeric()
    ### determine probability of yin sample
    if(prob_old=="exact"){
#       syy.yang=mean(posterior21.y/prior2)
#       syy.yin=mean(posterior12.y/prior1)
#       var.yang=1/length(sample2)*sum((posterior21.y/prior2)^2)-syy.yang^2
#       var.yin=1/length(sample1)*sum((posterior12.y/prior1)^2)-syy.yin^2
#       rhoyin=bayesMCClust::calcNumEff(sample1)
#       rhoyang=bayesMCClust::calcNumEff(sample2)
#       prob_old_1<-(rhoyang*var.yang)/(rhoyang*var.yang+rhoyin*var.yin)
#       prob_old_2<-1-prob_old_1
     exp.yin<-mean(sample1)
     exp.yang<-mean(sample2)
     exp.full<-0.5*exp.yin+0.5*exp.yang
      for(i in 1:neps){
        g.yin<-sample1*posterior21.y/prior2 # correct for likelihood
        g.yang<-sample2*posterior12.y/prior1 # correct for likelihood
        # calculate variances
        var.yin<-sum((g.yin-exp.full)^2)/len_sam  #cons_len1
        var.yang<-sum((g.yang-exp.full)^2)/len_sam  #cons_len2
        # calculate weights relative to total variance
        prob_old_1<-var.yin/sum(var.yin,var.yang)
        prob_old_2<-1-prob_old_1
        exp.full<-prob_old_1*exp.yin+prob_old_2*exp.yang
      }
    }
    if(prob_old=="heuristic"){
      prob_old_1<-max(0.05,min(0.95,weights_vec[k-1]/(weights_vec[k]+weights_vec[k-1])))
    }
    if(is.numeric(prob_old)){
      prob_old_1<-max(0.05,min(prob_old*(var(sample2)/var(sample1)),0.95))
    }
    for(i in 1:(len_sam)){
      if(sample1.new[i]!=0){
        # repeat each value of vector sample 1 for as many times
        # as drawn from the multinomial distribution
        new.sample1<-c(new.sample1,rep(sample1[i],sample1.new[i]))
        new.prior1<-c(new.prior1,rep(prior1[i],sample1.new[i]))
      }
      if(sample2.new[i]!=0){
        # repeat each value of vector sample 2 for as many times
        # as drawn from the multinomial distribution
        new.sample2<-c(new.sample2,rep(sample2[i],sample2.new[i]))
        new.prior2<-c(new.prior2,rep(prior2[i],sample2.new[i]))
      }
      #if(is.numeric(prob_old)){prob_old_1<-prob_old}
    }
    ind_sam<-(runif((len_sam),min=0,max=1)<prob_old_1)
    sample.out<-c(new.sample1[ind_sam],new.sample2[!ind_sam])
    prior<-c(new.prior1[ind_sam],new.prior2[!ind_sam])
    #if(method=="tree"){
      return(list(sample.out,prior))
    #}
    #if(method=="sequential"){
    #  return(list(sample.out))}
  }
#######################################
### Metropolis-Hastings Yin Yang sampler
#######################################
########################
## method sequential
if(method=="sequential"){
  ## caluculate sample set weights based on sample variances
  ## heuristic weights of samples
  if(prob_old=="heuristic"|is.null(prob_old)){
    weights_vec<-numeric(nsample)
    for(j in 1:nsample){
      weights_vec[j]<-1/var(sample_list[[j]]) # Scott et al.
    }
    weights_vec<-weights_vec/sum(weights_vec) #weights sum to 1
  }
  ## initialise first sample
  sample1<-sample_list[[1]]
  prior1<-prior_list[[1]]
    #prior1<-prior1#[(burn.in+1):length(sample1)]
    #sample1<-sample1#[(burn.in+1):length(sample1)]
  ## write out samples in a list for further analysis
  sam_list<-list()
  acceptance<-list()
  for(k in 2:nsample){
    # sample2 is the following sample
    assign(paste("sam2",sep=""),sample_list[[k]])
    sample2<-as.numeric(sam2)  #[(burn.in+1):length(sam2)]
    assign(paste("prior2",sep=""),prior_list[[k]])
    prior2<-as.numeric(prior2)  #[(burn.in+1):length(sam2)]
    ## trim all input vectors to same length
    #len_ss<-min(length(sample1),length(sample2))
    #sample1<-sample1[1:len_ss]
    #sample2<-sample2[1:len_ss]
    #prior1<-prior1[1:len_ss]
    #prior2<-prior2[1:len_ss]
    sample1<-rep(sample1,length.out = nsam)
    sample2<-rep(sample2,length.out = nsam)
    prior1<-rep(prior1,length.out = nsam)
    prior2<-rep(prior2,length.out = nsam)
    if(algo=="multinom"){
      resam<-resampleyinyang(sample1,sample2,prior1,prior2,prob_old,method,weight=c(weights_vec[k-1],weights_vec[k]))
      sample1<-resam[[1]]
      prior1<-resam[[2]]
      if(check){cat(summary(resam[[1]]),fill = TRUE)}
      sam_list<-c(sam_list,list(resam[[1]]))
      #prio_list_new<-c(prio_list_new,list(resam[[2]]))
    }
    if(algo=="MH"){
        resam<-resampleMH(sample1,sample2,prior1,prior2,prob_old,method,weight=c(weights_vec[k-1],weights_vec[k]))
          ## output
        sample1<-resam[[1]]
        prior1<-resam[[3]]
        if(check){cat(summary(resam[[1]]),fill = TRUE)}
        acceptance<-c(acceptance,list(resam[[2]]))
        sam_list<-c(sam_list,list(resam[[1]]))
      }
  }
  if(algo=="multinom"){ return(list(merged=sam_list))}
  if(algo=="MH"){ return(list(merged=sam_list,acceptance=acceptance))}
}
  ####################
  ## method tree
if(method=="tree"){
    # exponent -> for multiple of 2 obtain number of hierarchies in the tree
    expon<-round(log(nsample,base=2),digits=0)
    ## caluculate sample set weights based on sample variances
    acceptance<-list()
    sam_list_old<-sample_list
    prio_list_old<-prior_list
    for(k in (expon):1){ #run through tree
      ### number of samples =2^expon
      ## write out samples in a list for further analysis
      ## heuristic weights of samples
      if(prob_old=="heuristic"){
        weights_vec<-numeric(length(sam_list_old))
        for(j in 1:length(sam_list_old)){
          weights_vec[j]<-1/var(sam_list_old[[j]]) # Scott et al.
          #weights_vec[j]<-1/sd(sam_list_old[[j]]) # standard dev.
          #weights_vec[j]<-medmed(sam_list_old[[j]]) # robust weights
          ### little difference for symmetric 'Gaussian' cases
        }
        weights_vec<-weights_vec/sum(weights_vec) #weights sum to 1
      }
      sam_list_new<-list()
      prio_list_new<-list()
      if(check){cat(k,fill = TRUE)}
      for(j in seq(1,2^k-1,by=2)){
        if(check){cat(j,fill = TRUE)}
        # read vectors and priors; remove burn-in
        sam1<-sam_list_old[[j]]
        prior1<-prio_list_old[[j]]
        sample1<-sam1  #[(burn.in+1):length(sam1)]
        prior1<-prior1  # [(burn.in+1):length(sam1)]
        assign(paste("sam2",sep=""),sam_list_old[[j+1]])
        assign(paste("prior2",sep=""),prio_list_old[[j+1]])
        sample2<-sam2  #[(burn.in+1):length(sam2)]
        prior2<-prior2  #[(burn.in+1):length(sam2)]
        sample1<-rep(sample1,length.out = nsam)
        sample2<-rep(sample2,length.out = nsam)
        prior1<-rep(prior1,length.out = nsam)
        prior2<-rep(prior2,length.out = nsam)
        if(algo=="multinom"){
            resam<-resampleyinyang(sample1,sample2,prior1,prior2,prob_old,method,weight=c(weights_vec[k-1],weights_vec[k]))
            if(check){cat(summary(resam[[1]]),fill = TRUE)}
            sam_list_new<-c(sam_list_new,list(resam[[1]]))
            prio_list_new<-c(prio_list_new,list(resam[[2]]))
        }
        if(algo=="MH"){
            resam<-resampleMH(sample1,sample2,prior1,prior2,prob_old,method,weight=c(weights_vec[k-1],weights_vec[k]))
            if(check){cat(summary(resam[[1]]),fill = TRUE)}
            sam_list_new<-c(sam_list_new,list(resam[[1]]))
            acceptance<-c(acceptance,list(resam[[2]]))
            prio_list_new<-c(prio_list_new,list(resam[[3]]))
        }
    }
        sam_list_old<-sam_list_new
    if(check){cat(length(sam_list_old),fill = TRUE)}
        prio_list_old<-prio_list_new
  }
    if(algo=="multinom"){ return(list(merged=sam_list_new))}
    if(algo=="MH"){ return(list(merged=sam_list_new,acceptance=acceptance))}
}
}
