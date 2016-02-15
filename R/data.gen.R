# sim.data() : Function to generate left truncated and right censored data using the Cox's model

#' Generate left-truncated (and right-cencored) data from the Cox model.
#'
#' Various baseline survival functions and truncation distribution are
#' available. Censoring rate can be designated through tuning the parameter
#' \code{Cmax}; \code{Cmas = Inf} means no censoring.
#'
#' @param n An integer for sample size.
#' @param b A numeric vector for true regression coefficients.
#' @param distr.T A character for the baseline survival time (T*) distribution
#'   ("exp" or "weibull").
#' @param shape.T The shape parameter for the Weibull distribution of T*.
#' @param scale.T The scale parameter for the Exponential/Weibull distribution
#'   of T*.
#' @param distr.A A character for the baseline truncation time (A*) distribution
#'   ("weibull", "unif" (Length-Biased Sampling) or "binomial").
#' @param shape.A The shape parameter for the Weibull distribution of A*.
#' @param scale.A The scale parameter for the Weibull distribution of A*.
#' @param p.A The success probability for the binomial distribution of A*.
#' @param Cmax The upper bound of the uniform distribution of the censoring time
#'   (C).
#' @param fix.seed An optional random seed for simulation.
#' @return A data.frame containing the observed survival times (Ys), the
#'   observed truncation times (As), the event indicator (Ds) and the covariates
#'   (Zs).
#' @examples
#' sim.data()
#' sim.data(50, distr.A = "binomial", p.A = 0.8, Cmax = 5)
#' @export
sim.data = function(n=200, b = c(1,1),
                    distr.T = "weibull",
                    shape.T = 2, scale.T = 1,
                    distr.A = "weibull",
                    shape.A = 1, scale.A = 5,
                    p.A = 0.3,
                    Cmax = Inf, fix.seed = NA){
	p.trunc = NULL
	p.censr = NULL
	if(!is.na(fix.seed)) set.seed(fix.seed)
	j = 1
	k = 0

	# T = survival time
    # A = Truncation time
    # ZV = Time-Varying covariate
    # ZF = Fixed covariate
	Ts = NULL
	As = NULL
	Zv = NULL
	Zf = NULL

	while(j <= n){
		k = k + 1

		# Time at which the state changes from 0 to 1.
		Zv.j = rexp(1)
        # A continuous covariate
        Zf.j = runif(1,-1,1)

        # Relative risks: before = rr0.j; after = rr1.j
        RR0.j = exp(b[2]*Zf.j)
        RR1.j = exp(b[1] + b[2]*Zf.j)

        eps.j = rexp(1)

        # Generate Random Samples from the Cox's model.
        # First generate a piece-wise exponential r.v., then use H0^-1
        # to transform it.

	    if(distr.T == "exp"){
            t0 = Zv.j * scale.T
	        Ts.j = (min(eps.j, t0 * RR0.j) / RR0.j +
                    max(eps.j - t0 * RR0.j, 0) / RR1.j) / scale.T
	    }else if(distr.T == "weibull"){
            t0 = Zv.j ^ shape.T
	        Ts.j = (min(eps.j, t0 * RR0.j)/RR0.j +
	                    max(eps.j - t0 * RR0.j, 0) / RR1.j)^(1 / shape.T)
	    }
		if(distr.A == "weibull"){
			As.j = rweibull(1, shape.A, scale.A)
		}else if(distr.A == "unif"){
		    As.j = runif(1,0,100)
		}else if(distr.A == "binomial"){
		    As.j = rbinom(1,5,p.A)
		}else{
		    As.j = sample(0:5, 1)
		}

		# Keep only untruncated (A,T)'s
		if(As.j<Ts.j){
			Ts[j] = Ts.j
			As[j] = As.j
			Zv[j] = Zv.j
            Zf[j] = Zf.j
			j = j+1
		}else{
			next
		}
	}

	# C = Censoring time
	# Censoring is on residual lifetime; censoring occurs after truncation.
	if(Cmax != Inf){
		Cs = runif(n, 0, Cmax)
	}else{
        # ghost runs to keep the same random seed (=0/>0 censoring)
	    Cs = runif(n, 0, 1000)
		Cs = rep(Inf,n)
	}

	# Event indicators
	Cs = As + Cs
	Ds = as.numeric(Ts <= Cs)
	Ys = pmin(Ts,Cs)

	tau = quantile(Ys,seq(0.1,0.9,0.1))

	dat = data.frame(Zf,Zv,As,Ys,Ds)

    # Sort the dataset by observed times.
	dat = dat[order(dat$Ys),]

    dat = cbind(ID = 1:n, dat)

	rownames(dat) = NULL
	p.censr = 1 - mean(dat$Ds)
	p.trunc = 1 - n/k
	return(list(dat=dat,p.censr=p.censr,p.trunc=p.trunc,tau=tau))
}

#' Calculate the censoring rate for a survival dataset correspondes to a
#' specific \code{Cmax}.
#' Increasing the sample size and/or the iteration numbers may help make the
#' calculation more precise.
#'
#' @param I An integer iteration number.
#' @inheritParams sim.data
#' @return A decimal number of the censoring rate.
#' @examples
#' sim.pc()
#' sim.pc(1, n = 800, distr.A = "binomial", p.A = 0.8)
#' @export
sim.pc = function(Cmax=10,n=500,I=500,
                  distr.T="weibull",distr.A = "weibull",
                  shape.T=2, scale.T = 1,
                  shape.A = 1, scale.A=5, p.A = 0.3){
    set.seed(2357)
    vpc = NULL
    for(i in 1:I){
        if(i %% 50 == 0)print(i)
        vpc[i] = sim.data(Cmax=Cmax,n=n,LBS=LBS,distr.T=distr.T,distr.A=distr.A,shape.T=shape.T,scale.T=scale.T,shape.A=shape.A,scale.A=scale.A,p.A=p.A)$p.censr
    }
    pc = mean(vpc)
    return(pc)
}

#' Obtain the required Cmax for a specific censoring proportion (the inverse of
#' \code{sim.pc()}).
#'
#' @param pc numeric. Proportion of censoring.
#' @inheritParams sim.pc
#' @param l numeric. The lower bound of the searching interval.
#' @param u numeric. The upper bound of the searching interval.
#' @param M integer. The maximal iteration for \code{uniroot}.
#' @export
get.Cmax = function(pc=0.5,l=0.1,u=10,n=500,I=500,
                    distr.T="weibull",distr.A = "weibull",
                    shape.T=2, scale.T = 1,
                    shape.A = 1, scale.A=5, p.A = 0.3,
                    M=100){

    s = function(cm)sim.pc(cm,n=n,I=I,LBS=LBS,distr.T=distr.T,shape.T=shape.T,scale.A=scale.A)-pc
    rt=try(uniroot(s,c(l,u))$root,silent=T)
    i=1
    while(i<M){
        if(is.numeric(rt)){
            break
        }else if(s(l)<0){
            # in case that the initial interval does not contain the root
            l = l/2
            rt=try(uniroot(s,c(l,u))$root,silent=T)
        }else{
            u = u+(u-l)/2
            rt=try(uniroot(s,c(l,u))$root,silent=T)
        }
        print(paste0(l,u,"\n"))
        i=i+1
    }
    return(Cmax=rt)
}

