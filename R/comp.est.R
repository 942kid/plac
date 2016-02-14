# library(survival)
# library(Rcpp)
# library(RcppEigen)

#' @useDynLib plac
#' @importFrom Rcpp sourceCpp
#' @importFrom survival tmerge coxph basehaz

# sourceCpp("PLAC_TV_fun.cpp")
# sourceCpp("REPLACE_TV_fun.cpp")

#' @export
comp.est = function(dat){

    dat.td = tmerge(dat,dat,id = ID,death = event(Ys,Ds),trans=tdc(Zv))
    dat.td$id = NULL

    # those who had the transplant before enrolment will be Zv = 1
    # all the time during follow-up (the info before A are useless for the conditional approach??)

    dat.td = subset(dat.td, !(Zv <= As & tstart == 0))

    # those who had the transplant before enrolment should have
    # t.start at enrolment for the remaining interval after the above removal
    dat.td$tstart[dat.td$As > dat.td$tstart] = dat.td$As[dat.td$As > dat.td$tstart]

    cox.LTRC = coxph(Surv(tstart,tstop,death) ~ Zf + trans, dat = dat.td)
    b.cox = unname(coef(cox.LTRC))
    h0.cox = basehaz(cox.LTRC,centered=F)$hazard

    # get h0 from H0 to use as initials.
    if(h0.cox[1]!=0){
        H0.cox = c(0,unique(h0.cox))
        h0.cox = c(h0.cox[1],diff(h0.cox[!duplicated(h0.cox)]))
    }else{
        # if the first obs are censored..
        H0.cox = unique(h0.cox)
        h0.cox = diff(h0.cox[!duplicated(h0.cox)])
    }

    # fixed covariate matrix
	ZF = as.matrix(dat[2:(ncol(dat)-4)])
	ZFt = t(ZF)
	# change point of the time-varying covariate
	ZV = dat[[3]]
	# final status (at Xs) of the time-varying covariate(s)
	ZV_ = subset(dat.td, select = trans, tstop == Ys)[[1]]
	ZFV_ = rbind(ZFt, ZV_)
	# number of covariates (including the time-varying covariate)
	p = ncol(ZF) + 1
	# (As, Xs, Ds)
	X = as.matrix(dat[(ncol(dat)-2):ncol(dat)])
	# sample size
	n = nrow(X)
	# ordered observed distinct event times: w_k
	W = unique(dat$Ys[dat$Ds==1])
	# number of w_k
	m = length(unique(dat$Ys[dat$Ds==1]))
	# time-varying covariate process
	IndZ = TvInd1(ZV,W)
	# time-varying covariates process (including time-invariant covariates)
	Za = array(c(rep(ZF,each=m), IndZ), c(m,n,p))
	Z = matrix(aperm(Za, c(3,1,2)), m*p, n)
	# at-risk processes
	Ind1 = SgInd1(X,W)
	# truncation process
	Ind2 = PwInd1(X,W)
	# number of events at each w_k
	Dn = as.numeric(rle(dat$Ys[dat$Ds==1])$lengths)

    # number of w_k
    m = length(unique(dat$Ys[dat$Ds==1]))
    # sample size
    n = nrow(X)
    # number of covariates (including the time-varying covariate)
    p = ncol(ZF) + 1

    newdat = data.frame(rbind(cbind(ZF,ZV),rep(0,p)))
    newdat = subset(newdat,c(rep(F,n),T))
    names(newdat) = c("Zf", "trans")
    summ.cox = summary(survfit(cox.LTRC,newdata=newdat))
    se.H0.cox = c(0,summ.cox$std.err/summ.cox$surv)

    b = b.cox
    h = h0.cox

    plac.fit = PLAC_Tv3(Z,ZFV_,X,W,Ind1,Ind2,Dn,b,h,100)

    beta = rbind(Cox=b.cox, PLAC=plac.fit$b.hat)

    se.beta = rbind(Cox=sqrt(diag(cox.LTRC$var)), PLAC=plac.fit$se.b.hat)

    colnames(beta) = colnames(se.beta) = c(colnames(ZF), "ZV")

    H0 = cbind(Cox=H0.cox, PLAC=plac.fit$H0.hat)
    se.H0 = cbind(Cox=se.H0.cox, PLAC=plac.fit$se.H0.hat)

    return(list(Event.Time = summ.cox$time,
                b=beta,se.b=se.beta,
                H0=H0, se.H0=se.H0,
                sandwich=plac.fit$swe,
                k=plac.fit$iter))
}

#' @export
comp.estR = function(dat){

    dat.td = tmerge(dat,dat,id = ID,death = event(Ys,Ds),trans=tdc(Zv))
    dat.td$id = NULL

    # those who had the transplant before enrolment will be Zv = 1
    # all the time during follow-up (the info before A are useless for the conditional approach??)

    dat.td = subset(dat.td, !(Zv <= As & tstart == 0))

    # those who had the transplant before enrolment should have
    # t.start at enrolment for the remaining interval after the above removal
    dat.td$tstart[dat.td$As > dat.td$tstart] = dat.td$As[dat.td$As > dat.td$tstart]

    cox.LTRC = coxph(Surv(tstart,tstop,death) ~ Zf + trans, dat = dat.td)
    b.cox = unname(coef(cox.LTRC))
    h0.cox = basehaz(cox.LTRC,centered=F)$hazard

    # get h0 from H0 to use as initials.
    if(h0.cox[1]!=0){
        H0.cox = c(0,unique(h0.cox))
        h0.cox = c(h0.cox[1],diff(h0.cox[!duplicated(h0.cox)]))
    }else{
        # if the first obs are censored..
        H0.cox = unique(h0.cox)
        h0.cox = diff(h0.cox[!duplicated(h0.cox)])
    }

    # fixed covariate matrix
	ZF = as.matrix(dat[2:(ncol(dat)-4)])
	ZFt = t(ZF)
	# change point of the time-varying covariate
	ZV = dat[[3]]
	# final status (at Xs) of the time-varying covariate(s)
	ZV_ = subset(dat.td, select = trans, tstop == Ys)[[1]]
	ZFV_ = rbind(ZFt, ZV_)
	# number of covariates (including the time-varying covariate)
	p = ncol(ZF) + 1
	# (As, Xs, Ds)
	X = as.matrix(dat[(ncol(dat)-2):ncol(dat)])
	# sample size
	n = nrow(X)
	# ordered observed distinct event times: w_k
	W = unique(dat$Ys[dat$Ds==1])
	# number of w_k
	m = length(unique(dat$Ys[dat$Ds==1]))
	# time-varying covariate process
	IndZ = TvInd1(ZV,W)
	# time-varying covariates process (including time-invariant covariates)
	Za = array(c(rep(ZF,each=m), IndZ), c(m,n,p))
	Z = matrix(aperm(Za, c(3,1,2)), m*p, n)
	# at-risk processes
	Ind1 = SgInd1(X,W)
	# truncation process
	Ind2 = PwInd1(X,W)
	# number of events at each w_k
	Dn = as.numeric(rle(dat$Ys[dat$Ds==1])$lengths)

    # number of w_k
    m = length(unique(dat$Ys[dat$Ds==1]))
    # sample size
    n = nrow(X)
    # number of covariates (including the time-varying covariate)
    p = ncol(ZF) + 1

    newdat = data.frame(rbind(cbind(ZF,ZV),rep(0,p)))
    newdat = subset(newdat,c(rep(F,n),T))
    names(newdat) = c("Zf", "trans")
    summ.cox = summary(survfit(cox.LTRC,newdata=newdat))
    se.H0.cox = c(0,summ.cox$std.err/summ.cox$surv)

    b = b.cox
    h = h0.cox

    plac.fit = PLAC_TvR(ZFt,ZFV_,Z,X,W,Ind1,Ind2,Dn,b,h,100)

    beta = rbind(Cox=b.cox, PLAC=plac.fit$b.hat)

    se.beta = rbind(Cox=sqrt(diag(cox.LTRC$var)), PLAC=plac.fit$se.b.hat)

    colnames(beta) = colnames(se.beta) = c(colnames(ZF), "ZV")

    H0 = cbind(Cox=H0.cox, PLAC=plac.fit$H0.hat)
    se.H0 = cbind(Cox=se.H0.cox, PLAC=plac.fit$se.H0.hat)

    return(list(Event.Time = summ.cox$time,
                b=beta,se.b=se.beta,
                H0=H0, se.H0=se.H0,
                sandwich=plac.fit$swe,
                k=plac.fit$iter))
}

#' @export
step.L = function(est, t.eval=c(0.25, 0.75)){

    # evaluation of Lambda_hat
    tim = est$Event.Time

    H0.cox = stepfun(tim, est$H0[,1])
    H0.cox.eval = H0.cox(t.eval)

    H0.plac = stepfun(tim,est$H0[,2])
    H0.plac.eval = H0.plac(t.eval)

    L = rbind(Cox=H0.cox.eval, PLAC=H0.plac.eval)

    se.H0.cox = stepfun(tim,est$se.H0[,1])
    se.H0.cox.eval = se.H0.cox(t.eval)

    se.H0.plac = stepfun(tim,est$se.H0[,2])
    se.H0.plac.eval = se.H0.plac(t.eval)

    se.L = rbind(Cox=se.H0.cox.eval, PLAC=se.H0.plac.eval)

    colnames(L) = t.eval
    colnames(se.L) = t.eval

    return(list(TrueL = t.eval,
                L=L, se.L=se.L,
                L.fun=list(Cox=H0.cox, PLAC=H0.plac),
                se.L.fun=list(Cox=se.H0.cox, PLAC=se.H0.plac)))

}
