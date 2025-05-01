#' @keywords internal
#' @noRd
estimate.rank.aft = function(y, delta, x, Gamma=NULL, b=NULL, beta=NULL,
                             test=TRUE, regularize=FALSE, tol=10^-12) {

  n.data=length(y)
	if ('data.frame' %in% class(x) ) x=as.matrix(x)
	if( !('matrix' %in% class(x)) ) x=matrix(x, ncol=1)
	if (!is.null(Gamma)) {
	  if( !('matrix' %in% class(Gamma)) ) Gamma=matrix(Gamma,nrow=1)
		if (ncol(Gamma)!=ncol(x)) {
		   print('Error:  Gamma and X must have the same number of columns')
		   return()
		}
	}

	if (is.null(b) & !is.null(Gamma) ) {
		print('Error:  both Gamma and b must be specified for single-taxon
		      constrained analysis')
		return()
	}

  #
  #   center x and y
  #
	y=y-mean(y)
	x=x-rep( colMeans(x), each=n.data )

	setup.res=setup(Gamma=Gamma, Lambda=Lambda, x=x, b=b, tol=tol)
	Gamma=setup.res$Gamma
	Lambda=setup.res$Lambda
	n.gamma=setup.res$n.gamma
	n.lambda=setup.res$n.lambda

	x=setup.res$x
	b=setup.res$b
	n.var=ncol(x)

	if (n.lambda==0) {
	  beta.r=as.vector( t(Gamma) %*% b )
		res=list(beta=NULL, beta.r=beta.r, n.gamma=n.gamma, n.lambda=n.lambda,
		         Gamma=Gamma, Lambda=Lambda, b=b, y=y, x=x, delta=delta,
		         regularize=regularize)
		return(res)
	}

	n.param=n.lambda
	if( is.null(beta) ){
	  lambda=rep(0,n.lambda)
	}else{
	  lambda=as.vector( Lambda %*% beta )
	}
	param0 = lambda


  #
  #   estimate parameters
  #

  if( n.var==1 ) {
    res=optim(param0, fn=fcn.J.n, method='Brent', lower=-10, upper=10, y=y,
              delta=delta, x=x, Gamma=Gamma, Lambda=Lambda, b=b,
              n.gamma=n.gamma, n.lambda=n.lambda, regularize=regularize,
              control=list(reltol=10^-12, abstol=10^-12))
  }else{
	  res=optim(param0, gr=grad.J.n, fn=fcn.J.n, method='BFGS', y=y, delta=delta,
	            x=x, Gamma=Gamma, Lambda=Lambda, b=b, n.gamma=n.gamma,
	            n.lambda=n.lambda, regularize=regularize,
	            control=list(reltol=10^-14, abstol=10^-12))
  }
  par=res$par
  value=res$value

	beta=NULL
	beta.r=NULL


  lambda=par[1:n.lambda]
	if (n.gamma==0) {
		beta=as.vector(t(Lambda)%*% lambda)
	}else{
		beta.r=as.vector( t(Gamma)%*% b + t(Lambda)%*% lambda )
	}


	res=list(beta=beta, beta.r=beta.r, value=value, n.gamma=n.gamma,
	         n.lambda=n.lambda, Gamma=Gamma, Lambda=Lambda, y=y, x=x,
	         delta=delta, regularize=regularize, b=b)

	return(res)
}




#' @keywords internal
#' @noRd
setup = function(Gamma, Lambda, x, b, tol=tol) {

  n.var=ncol(x)

  # assign Lambda=Diag(n.var) if Gamma not specified;
  #   else convert Gamma to matrix if necessary and calculate
  #   orthogonal Gamma and Lambda matrices


	if (is.null(Gamma)) {
	  Lambda=diag(n.var)
		n.gamma=0
		n.lambda=n.var
	}else{
	  if( !('matrix' %in% class(Gamma) ) ) {
	    Gamma=matrix(Gamma,nrow=1)
		}
		proj=t(Gamma)%*%Gamma
    proj.eigen=eigen(proj)
    use.gamma=( abs(proj.eigen$values)>tol )
    use.lambda=!use.gamma
		n.gamma=sum(use.gamma)
		n.lambda=sum(use.lambda)
    R= proj.eigen$vectors[,use.gamma]
		L.D=Gamma %*% R
		Gamma=t(R)
		if (n.lambda>0) Lambda=t( proj.eigen$vectors[,use.lambda] )
        b=solve(L.D,b)
	}
	res=list(Gamma=Gamma, Lambda=Lambda, x=x, n.gamma=n.gamma,
	         n.lambda=n.lambda, b=b)
  return(res)
}








#' @keywords internal
#' @noRd
fast.J.n = function(beta, y, x, delta, regularize) {

	n.data=length(y)
	x.beta=colSums( beta*t(x) )
  e=y-x.beta
	ord=order(e)
	e=e[ord]
	delta=delta[ord]
	cum.delta=cumsum(delta)
	r.min=rank(e, ties.method='min')
	r.max=rank(e, ties.method='max')
	r.i.dot=delta*(n.data - 0.5*(r.min+r.max)) + 0.5*delta
	r.dot.i=0.5*(cum.delta[r.min]+cum.delta[r.max]) - 0.5*delta
#    s=n.data*(r.i.dot - r.dot.i)
  s=r.i.dot - r.dot.i

	Q.n = -sum( e*s )
	if (regularize==TRUE) Q.n=Q.n + 0.5*min(0.01*n.data,1)*sum( beta^2 )

	return(Q.n)
}


#' @keywords internal
#' @noRd
fast.W.n = function(beta,y,x,delta, regularize) {

  n.data=length(y)
	x.beta=colSums( beta*t(x) )
  e=y-x.beta
	ord=order(e)
	rnk=rep(0,n.data)
	rnk[ord]=1:n.data
	e=e[ord]
	delta=delta[ord]
	cum.delta=cumsum(delta)
	r.min=rank(e, ties.method='min')
	r.max=rank(e, ties.method='max')
	r.i.dot=delta*(n.data - 0.5*(r.min+r.max)) + 0.5*delta
	r.dot.i=0.5*(cum.delta[r.min]+cum.delta[r.max]) - 0.5*delta
    s=r.i.dot - r.dot.i
#	s=n.data*s[rnk]
	s=s[rnk]

	W.n = colSums( s*x )
	if (regularize==TRUE) W.n=W.n + min(0.01*n.data,1)*beta

	res=list(W.n=W.n, s=s)
	return(res)
}




#' @keywords internal
#' @noRd
fcn.J.n = function(par, y, delta, x, Gamma, Lambda, n.gamma, n.lambda, b, regularize) {

	n.p=length(par)
	W.n=rep(0,n.p)

  if (is.null(Gamma)) {
		beta=par
	}else if (is.null(Lambda)) {
		beta=as.vector( t(Gamma) %*% b )
	}else {
    lambda=par
    beta=as.vector( t(Gamma)%*% b + t(Lambda) %*% lambda )
	}

	J.n=fast.J.n(beta=beta,y=y,x=x,delta=delta, regularize=regularize)

  return(J.n)
}


#' @keywords internal
#' @noRd
grad.J.n = function(par, y, delta, x, Gamma, Lambda, n.gamma, n.lambda, b, regularize) {

	n.p=length(par)
	W.n=rep(0,n.p)

  if (is.null(Gamma)) {
		beta=par
		W.n=fast.W.n(beta=beta,y=y,x=x,delta=delta, regularize=regularize)$W.n
	}else if (is.null(Lambda)) {
		beta=as.vector( t(Gamma) %*% b )
		W.n.beta=fast.W.n(beta=beta,y=y,x=x,delta=delta, regularize=regularize)$W.n
		W.n=Gamma %*% W.n.beta
	}else {
    lambda=par
    beta=as.vector( t(Gamma)%*% b + t(Lambda) %*% lambda )
    W.n.beta=fast.W.n(beta=beta,y=y,x=x,delta=delta, regularize=regularize)$W.n
		W.n=Lambda %*% W.n.beta
	}

  return(W.n)
}








#test.rank.aft = function(beta, y, delta, x, Gamma=Gamma, Lambda=Lambda, n.gamma=n.gamma, n.lambda=n.lambda, b=b, score='Cox') {

#' @importFrom expm sqrtm
#' @keywords internal
#' @noRd
test.rank.aft = function(est.rank.res, score='rank') {

  if (!(score %in% c('Cox','rank')) ) {
		print('Error - score must either be Cox or rank')
		return()
	}

  y=est.rank.res$y
	x=est.rank.res$x
	delta=est.rank.res$delta
	Gamma=est.rank.res$Gamma
	Lambda=est.rank.res$Lambda
	n.gamma=est.rank.res$n.gamma
	n.lambda=est.rank.res$n.lambda
	beta=est.rank.res$beta
	if (is.null(beta)) beta=est.rank.res$beta.r
	beta.2=est.rank.res$beta.2
	if (is.null(beta.2)) beta.2=est.rank.res$beta.r2
	b=est.rank.res$b
  n.data=length(y)
	regularize=est.rank.res$regularize

#
#   calculate score and variance-covariance of score function
#


  c.matrix=rbind(Gamma,Lambda)
  W.n.beta.res=fast.W.n(beta=beta, y=y, delta=delta, x=x, regularize=regularize)
	W.n.beta=W.n.beta.res$W.n
	if (score=='Cox') {
		s.i.res=mySi.no.surv(beta=beta,y=y,delta=delta,x=x)$s
		s.i.res=s.i.res %*% t(c.matrix)
		v=cov(s.i.res)
	}else if (score=='rank') {
	  s.i=W.n.beta.res$s
		cov.x=cov(x)
		v=var(s.i)*c.matrix %*% cov.x %*% t(c.matrix)
    #	var.s=var(s.i)
	}

#
#   calculate test statistics
#

    if (is.null(Gamma)) {
	    sigma.half=sqrtm(v)
		  z.score=solve(sigma.half, W.n.beta)/sqrt(n.data)
		  test=sum( z.score^2 )
      #	sigma.inv=solve(v,W.n.beta)
      #	test=W.n.beta %*% sigma.inv/n.data
		  df=n.lambda
		  W.n=W.n.beta
      #	test.normal=NULL
		}else{
		  W.n.gamma=Gamma %*% W.n.beta
		  sigma.inv=as.matrix(solve(v)[1:n.gamma,1:n.gamma])
		  sigma.half=expm::sqrtm(sigma.inv)
		  z.score= sigma.half %*% W.n.gamma/sqrt(n.data)
		  test=sum(z.score^2)
      #	test=t(W.n.gamma) %*% sigma.inv %*% W.n.gamma/n.data
      #	test.normal=t(W.n.gamma) %*% sqrt(sigma.inv/n.data)
		  df=n.gamma
		  W.n=W.n.gamma
		}

	p.value=pchisq(test, df=df, lower.tail=FALSE)
  #	res=list(test=test, df=df, p.value=p.value, v=v, W.n=W.n, var.s=var.s)
	res=list(test=test, z.score=z.score, df=df, p.value=p.value, v=v, W.n=W.n)

	return(res)

}














#' @keywords internal
#' @noRd
simulate.aft = function(n.obs, F, F.params, G, G.params, beta, x, y.int=0, c.int=0) {
    if (length(F.params)==1) {eps=F(n.obs,F.params[1])}
	else {eps=F(n.obs,F.params[1], F.params[2])}
	eps=eps-mean(eps)
	if ('matrix' %in% class(x)) {
		y.star= as.vector(x %*% beta) + eps + y.int
		}
	else {
		y.star= x*beta + eps + y.int
		}
	t.star=exp(y.star)
	if (length(G.params)==1) {c.star=G(n.obs, G.params[1]) + c.int}
	else {c.star=G(n.obs, G.params[1], G.params[2]) + c.int}
    delta=ifelse(y.star<=c.star, 1, 0)
    y=pmin(c.star, y.star)
	t=exp(y)
	res=list( t=t, y=y, delta=delta, t.star=t.star, y.star=y.star, x=x )
	return(res)
	}









#' @keywords internal
#' @noRd
mySi.no.surv = function(beta, y, x, delta) {

	if (!('matrix' %in% class(x))) x=matrix(x,ncol=1)
	n.data=nrow(x)
	n.var=ncol(x)

	x.beta=colSums( beta*t(x) )
#   x.beta=as.vector( x %*% beta )
    e=y-x.beta
    ord=order(e)
	e=e[ord]
	delta=delta[ord]
	x=x[ord,]
	rnk=rep(0,n.data)
	rnk[ord]=1:n.data

	tab=table(delta, e)
	n.times=ncol(tab)
	events=colSums(tab)
	censored=tab[1,]
	failures=tab[2,]
    gamma0=rev( cumsum( rev(events) ) )
	cumhaz=cumsum( failures/gamma0 )
	d.Lambda=diff( c(0,cumhaz) )
    id=rep(1:n.times, times=events)
	X=rowsum(x,id)
	gamma1=apply(X, MARGIN=2, FUN=function(x) rev(cumsum(rev(x))) )

	omega0=cumsum(d.Lambda*gamma0)
	omega1=apply(d.Lambda*gamma1,MARGIN=2,FUN=cumsum)
#
#   restore ties
#
    gamma0=gamma0[id]
	gamma1=gamma1[id,,drop=FALSE]
	omega0=omega0[id]
	omega1=omega1[id,,drop=FALSE]
	cumhaz=cumhaz[id]
#
#   final calculation
#
	s = delta*(gamma0*x - gamma1) - (omega0*x - omega1)

    s=s[rnk,,drop=FALSE]
	cumhaz=cumhaz[rnk]


	res=list( s=s, gamma0=gamma0, gamma1=gamma1, omega0=omega0, omega1=omega1, cumhaz=cumhaz, n.times=n.times)
	return(res)
    }


#' @keywords internal
#' @noRd
km = function(t, delta) {
	n.obs=length(t)
    t.order=order(t)
	t.rank=rep(NA,n.obs)
	t.rank[t.order]=1:n.obs
    times=t[t.order]
    del=delta[ t.order ]
    tab=table( del,times )
	km.times=as.numeric( colnames(tab) )
    if (all( range(delta)==c(0,1) ) ) {
		fail.counts=tab[2,]
		cens.counts=tab[1,]
		n.unique.times=dim(tab)[2]
		}
	else if (all(delta==1)) {
		fail.counts=tab
		cens.counts=rep(0,length(fail.counts))
		n.unique.times=length(tab)
		}
	else if (all(delta==0)) {
		cens.counts=tab
		fail.counts=rep(0,length(cens.counts))
		n.unique.times=length(tab)
		}
    tot.counts=fail.counts+cens.counts
    risk.set=n.obs - c(0, cumsum(tot.counts[-n.unique.times]) )
	names(risk.set)=names(tot.counts)
    km=cumprod(1 - fail.counts/risk.set)
    km.all=rep(km, times=tot.counts)
	km.all.times=rep(km.times, times=tot.counts)
    km.all=km.all[t.rank]
#   res=list(km=km.all)
#   return(res)
    return(km.all)
   }
