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

    cox.LTRC = coxph(Surv(tstart,tstop,death) ~ Zf + trans, data = dat.td)
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

#' Calculate the PLAC estimator when a time-dependent indicator presents
#'
#' Both a conditional approach Cox model and a pairwise likelihood augmented estimator are fitted and the corresponding results are returned in a list.
#' @param formula a formula of of the form \code{Surv(A, X, D) ~ Z} including the time-invariate covariates only.
#' @param a data.frame of the LTRC dataset including the responses, time-invariate covariates and the jump times for the time-depnencent covariate.
#' @param id.var a name of the subject id in \code{data}.
#' @param td.var  aname of the time-dependent covariate in the output.
#' @param tjump a name of the jump time variable in \code{data}.
#' @param init.val a list of the initial values of the coefficients and the baseline hazard function for the PLAC estimator.
#' @importFrom survival Surv tmerge coxph
#' @return a list of model fitting results for both conditional approach and the PLAC estimators.
#' @export
plactv = function(formula, data, id.var, td.var, tjump,
                init.val = NULL, ...){
  if( !inherits(formula, "formula") ) stop("A formula is required!")
  # grep the model.frame for later use
  mf = model.frame(formula=formula, data=data)
  # prepare the response and time-invariant covariate matrix
  X = as.matrix(model.response(mf))
  tmp1 <<- X
  resp = gsub("Surv|\\(|\\)| ","",unlist(strsplit(names(mf)[1],",")))
  colnames(X) = resp
  n = nrow(X)
  W = unique(X[X[,3] == 1, 2])
  m = length(W)
  # at-risk processes
  Ind1 = SgInd1(X, W)
  tmp2 <<- Ind1
  # truncation process
  Ind2 = PwInd1(X, W)
  tmp3 <<- Ind2
  # number of events at each w_k
  Dn = as.numeric(rle(X[X[,3] == 1, 2])$lengths)
  tmp4 <<- Dn
  ZF = as.matrix(model.matrix(attr(mf, "terms"), data=mf)[,-1], nrow = n)
  colnames(ZF) = names(mf)[-1]
  # need counting process expansion of the data
  eval(parse(text = paste0("data.count = tmerge(data,data,id =", id.var,
                           ", death = event(", resp[2], ", ", resp[3], "), ",
                           td.var, " = tdc(", tjump, "))")))
  data.count$id = NULL
  eval(parse(text = paste0("data.count = subset(data.count, !(",
                           tjump, " <= ", resp[1],
                           "& tstart == 0))")))

  eval(parse(text = paste0("data.count$tstart[data.count$",
                           resp[1], "> data.count$tstart] = data.count$",
                           resp[1], "[data.count$",
                           resp[1], "> data.count$tstart]")))
  cox.LTRC = coxph(as.formula(paste0("Surv(tstart, tstop, death) ~ ",
                                     paste(c(names(mf)[-1], td.var), collapse = "+"))), data = data.count)
  b.cox = unname(coef(cox.LTRC))
  h.cox = basehaz(cox.LTRC,centered=F)$hazard
  # get h from H
  if(h.cox[1]!=0){
    H.cox = c(0, unique(h.cox))
    h.cox = diff(H.cox)
  }else{
    # if the first obs is censored..
    H.cox = unique(h.cox)
    h.cox = diff(H.cox)
  }
  newdata = model.matrix(cox.LTRC)[FALSE,]
  newdata = data.frame(rbind(newdata, rep(0, ncol(newdata))))
  summ.cox = summary(survfit(cox.LTRC, newdata = newdata))
  se.H0.cox = c(0, summ.cox$std.err/summ.cox$surv)
  # set the initial values for b (coefficients) and h (baseline hazard)
  if( is.null(init.val) ){
    b_0 = b.cox
    h_0 = h.cox
  }else{
    b_0 = init.val$b_0
    h_0 = init.val$h_0
  }
  p = length(b_0)
  ZFt = t(ZF)
  # the jump times of the time-dependent indicator (zeta)
  ZV = data[[tjump]]
  # covariate values at the observed survival times
  eval(parse(text = paste0("ZV_ = subset(data.count, select = ", td.var, ", tstop == ", resp[2],")[[1]]")))
  ZFV_ = rbind(ZFt, ZV_)
  IndZ = TvInd1(ZV, W)
  Za = array(c(rep(ZF, each = m), IndZ), c(m, n, p))
  Z = matrix(aperm(Za, c(3, 1, 2)), m * p, n)

  plac.fit = PLAC_Tv3(Z, ZFV_, X, W, Ind1, Ind2, Dn, b_0, h_0, K = 100)
  beta = rbind(Cox=b.cox, PLAC=plac.fit$b.hat)
  se.beta = rbind(Cox=sqrt(diag(cox.LTRC$var)), PLAC=plac.fit$se.b.hat)
  colnames(beta) = colnames(se.beta) = c(colnames(ZF), "ZV")
  H0 = cbind(Cox=H.cox, PLAC=plac.fit$H0.hat)
  se.H0 = cbind(Cox=se.H0.cox, PLAC=plac.fit$se.H0.hat)

  return(list(Event.Time = summ.cox$time,
              b=beta,se.b=se.beta,
              H0=H0, se.H0=se.H0,
              sandwich=plac.fit$swe,
              k=plac.fit$iter))
}
