

#' Calculate the PLAC estimator when a time-dependent indicator presents
#'
#' Both a conditional approach Cox model and a pairwise likelihood augmented estimator are fitted and the corresponding results are returned in a list.
#' @param formula a formula of of the form \code{Surv(A, X, D) ~ Z} including the time-invariate covariates only.
#' @param a data.frame of the LTRC dataset including the responses, time-invariate covariates and the jump times for the time-depnencent covariate.
#' @param id.var a name of the subject id in \code{data}.
#' @param td.var  a name of the time-dependent covariate in the output.
#' @param type the type of the time-dependent covariate. Either one of \code{c("independent", "post-trunc", "pre-post-trunc")}.
#' @param tjump a name of the jump time variable in \code{data}.
#' @param init.val a list of the initial values of the coefficients and the baseline hazard function for the PLAC estimator.
#' @useDynLib plac
#' @importFrom Rcpp sourceCpp
#' @importFrom survival Surv tmerge coxph
#' @return a list of model fitting results for both conditional approach and the PLAC estimators.
#' @export
plactv = function(formula, data, id.var,
                  td.var, type = "independent", tjump,
                  init.val = NULL, max.iter = 100, ...){
  if( !inherits(formula, "formula") ) stop("A formula is required!")
  # grep the model.frame for later use
  mf = model.frame(formula=formula, data=data)
  # prepare the response and time-invariant covariate matrix
  X = as.matrix(model.response(mf))
  resp = gsub("Surv|\\(|\\)| ","",unlist(strsplit(names(mf)[1],",")))
  colnames(X) = resp
  n = nrow(X)
  W = unique(X[X[,3] == 1, 2])
  m = length(W)
  # at-risk processes
  Ind1 = SgInd1(X, W)
  # truncation process
  Ind2 = PwInd1(X, W)
  # number of events at each w_k
  Dn = as.numeric(rle(X[X[,3] == 1, 2])$lengths)
  ZF = as.matrix(model.matrix(attr(mf, "terms"), data=mf)[,-1], nrow = n)
  colnames(ZF) = names(mf)[-1]
  if( type == "post-trunc" ){
    assign(td.var, rep(0, n))
    eval(parse(text = paste0("ZF = cbind(ZF,", td.var, ")")))
  }
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

  if( type == "independent" ){
    plac.fit = PLAC_Tv3(Z, ZFV_, X, W, Ind1, Ind2, Dn, b_0, h_0, max.iter)
  }else if( type %in% c("post-trunc", "pre-post-trunc") ){
    plac.fit = PLAC_TvR1(ZFt, ZFV_, Z, X, W, Ind1, Ind2, Dn, b_0, h_0, max.iter)
  }

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

#' Calulate the Values of the cumulative Hazard at Fixed Times
#'
#' @param est an object of the class \code{plac.fit}.
#' @param t.eval time points at which the Lambda(t) is evaluated (for both conditional apporach and the PLAC estimator).
#' @return a list containing the estiamtes and SEs of Lambda(t) for both conditional apporach and the PLAC estimator.
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
