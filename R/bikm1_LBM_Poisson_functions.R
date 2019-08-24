
##' BIKM1_LBM_Poisson fitting procedure
##'
##' Produce a blockwise estimation of a contingency matrix of observations.
##'
##' @param x contingency matrix of observations.
##' @param Hmax a positive integer less than number of rows.
##' @param Lmax a positive integer less than number of columns.The bikm1 procedure stops while the numbers of rows is higher than Hmax or the number of columns is higher than Lmax.
##' @param a hyperparameter used in the VBayes algorithm for priors on the mixing proportions.
##' By default, a=4.
##' @param alpha hyperparameter used in the VBayes algorithm for prior on the Poisson parameter.
##' By default, alpha=1.
##' @param beta hyperparameter used in the VBayes algorithm for prior on the Poisson parameter.
##' By default, beta=0.01.
##' @param Hstart a positive integer to initialize the procedure with number of row clusters.
##' By default, Hstart=2.
##' @param Lstart a positive integer to initialize the procedure with number of column clusters.
##' By default, Lstart=2.
##' @param normalization logical. To use the normalized Poisson modelling in the Latent Block Model. By default normalization=FALSE.
##' @param init_choice character string corresponding to the chosen initialization strategy used for the procedure, which can be "random" or "Gibbs" (higher time computation) or "smallVBayes" or "user".
##' By default, init_choice="smallVBayes"
##' @param userparam In the case where init_choice is "user", a list containing partitions v and w.
##' @param ntry a positive integer corresponding to the number of times which is launched the small VBayes or random initialization strategy. By default ntry=50.
##' @param criterion_choice Character string corresponding to the chosen criterion used for model selection, which can be "ICL" or "BIC".
##' By default, criterion_choice="ICL".
##' @param mc.cores a positive integer corresponding to the available number of cores for parallel computing.
##' By default, mc.cores=1.
##' @param verbose logical. To display each step and the result. By default verbose=TRUE.
##' @rdname BIKM1_LBM_Poisson-proc
##' @references Keribin, Celeux and Robert, The Latent Block Model: a useful model for high dimensional data.
##' https://hal.inria.fr/hal-01658589/document
##'
##' Govaert and Nadif. Coclustering, Wyley (2013).
##'
##' Keribin, Brault and Celeux. Estimation and Selection for the Latent Block Model on Categorical Data, Statistics and Computing (2014).
##'
##' Robert. Classification crois\'ee pour l'analyse de bases de donn\'ees de grandes dimensions de pharmacovigilance. Paris Saclay (2017).
##' @usage BIKM1_LBM_Poisson(x,Hmax,Lmax,a=4,alpha=1,beta=0.01,
##' Hstart=2,Lstart=2,normalization=FALSE,init_choice='smallVBayes',
##' userparam=NULL,ntry=50,criterion_choice='ICL', mc.cores=1,verbose=TRUE)
##'
##' @return a BIKM1_LBM_Poisson object including
##'
##' \code{model_max}: the selected model by the procedure with free energy W, theta, conditional probabilities (r_jh, t_kl), iter, empty_cluster, and the selected partitions v and w.
##'
##' \code{criterion_choice}: the chosen criterion
##'
##' \code{init_choice}: the chosen init choice
##'
##' \code{criterion tab}: matrix containing the criterion values for each selected number of row and column
##'
##' \code{W_tab}: matrix containing the free energy values for each selected number of row and column
##'
##' \code{criterion_max}: maximum of the criterion values
##'
##' \code{hopt}: the selected number of rows
##'
##' \code{lopt}: the selected number of columns
##' @examples
##' require(bikm1)
##' J=200
##' K=120
##' h=3
##' l=2
##' theta=list()
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##' theta$gamma_hl=matrix(c(1, 6,4, 1, 7, 1),ncol=2)
##' data=PoissonBlocRnd(J,K,theta)
##' res=BIKM1_LBM_Poisson(data$x,3,2,Hstart=3,Lstart=2,
##' init_choice='user',userparam=list(v=data$xrow,w=data$xcol))
##'
##' @export BIKM1_LBM_Poisson
##'


BIKM1_LBM_Poisson=function(x,Hmax,Lmax,a=4,alpha=1,beta=0.01,Hstart=2,Lstart=2,normalization=FALSE,init_choice='smallVBayes',userparam=NULL,ntry=50,criterion_choice='ICL', mc.cores=1,verbose=TRUE){

  if (normalization){


    PartRnd = function(J,proba){

      h=length(proba)

      i=sample(1:J,J, replace = FALSE)
      i1=i[1:h]
      i2=i[(h+1):J]
      v=rep(0,J)
      v[i1]=1:h

      if (J > h) {
        v[i2]=(h+1)-colSums(rep(1,h)%*%t(runif(J-h)) < cumsum(proba)%*%t(rep(1,J-h)))
        v[v==(h+1)]=h
      }

      v_jh=!(v%*%t(rep(1,h))-rep(1,J)%*%t(1:h))

      list(v=v,v_jh=v_jh)
    }
  PoissonBlocInit=function(h,l,x,normalization,start,fixed_proportion=FALSE){
    eps=10^(-16)
    J=dim(x)[1]
    K=dim(x)[2]

    # Initialization by random simulations  of centers

    theta=list()
    if (start$ale==TRUE){
      if (normalization){
        mu_i=matrix(rowSums(x),ncol=1)
        nu_j=matrix(colSums(x),ncol=K)
      }else{mu_i=matrix(rep(1,J),ncol=1)
      nu_j=matrix(rep(1,K),ncol=K)
      }
      theta$rho_h=1/h*matrix(1,h,1)
      #theta$rho_h=theta$rho_h/sum(theta$rho_h)
      theta$tau_l=1/l*matrix(1,l,1)
      #theta$tau_l=theta$tau_l/sum(theta$tau_l)
      resw=PartRnd(K,theta$tau_l)
      t_kl=resw$v_jh
      w=resw$v
      resv=PartRnd(J,theta$rho_h)
      r_jh=resv$v_jh
      v=resv$v
      n_kl=t((t(mu_i)%*%r_jh))%*%(nu_j%*%t_kl)
      theta$gamma_hl=(t(r_jh)%*%x%*%t_kl)/(n_kl)
      if (any((is.nan(theta$gamma_hl)))){
        lambda0=mean(mean(x))
        theta$gamma_hl=matrix(rpois(h,lambda0),ncol=l)
      }

    }else if(!is.null(start$theta) && !is.null(start$w)){
      theta=start$theta
      t_kl=!(start$w%*%t(rep(1,l))-rep(1,K)%*%t(1:l))
      w=start$w
    } else if(!is.null(start$theta) && !is.null(start$t_kl)){
      theta=start$theta
      t_kl=start$t_kl
      w=apply(start$t_kl,1,which.max)
    } else if(!is.null(start$theta) && !is.null(start$v)){
      theta=start$theta
      #r_jh=!(init$v*matrix(1,1,h)-matrix(1,J,1)*c(1:h))
      r_jh=!(start$v%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
      r_h=t(colSums(r_jh))
      eps10=10*eps
      log_tau_l=log(theta$tau_l)
      gamma_hl=theta$gamma_hl
      log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
      t_kl=matrix(1,K,1)%*%t(log_tau_l)- (t(nu_j)%*%((t(mu_i)%*%r_jh)%*%gamma_hl))+t(x)%*%r_jh%*%log_gamma_hl
      t_klmax=apply(t_kl,1,max)
      Tp_jl = exp(t_kl-t_klmax)
      t_kl=Tp_jl/(rowSums(Tp_jl))
      w=apply(t_kl,1,which.max)
    } else if(!is.null(start$theta) && !is.null(start$r_jh)){
      r_jh=start$r_jh
      r_h=t(colSums(r_jh))
      eps10=10*eps
      log_tau_l=log(theta$tau_l)
      gamma_hl=theta$gamma_hl
      log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
      t_kl=matrix(1,K,1)%*%t(log_tau_l)- (t(nu_j)%*%((t(mu_i)%*%r_jh)%*%gamma_hl))+t(x)%*%r_jh%*%log_gamma_hl
      t_klmax=apply(t_kl,1,max)
      Tp_jl = exp(t_kl-t_klmax)
      t_kl=Tp_jl/(rowSums(Tp_jl))
      w=apply(t_kl,1,which.max)
    } else if(!is.null(start$v) && !is.null(start$w)){

      if ((max(start$v)!=h)||(max(start$w)!=l)){
        warning('(Hstart,Lstart) need to match with the characteristics of userparam v and w')
      }else{
      r_jh=!(start$v%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
      t_kl=!(start$w%*%t(rep(1,l))-rep(1,K)%*%t(1:l))

      r_h=matrix(colSums(r_jh),ncol=1)
      t_l=matrix(colSums(t_kl),ncol=1)
      if (fixed_proportion==TRUE){
        theta$rho_h=1/h *rep(1,h)
        theta$tau_l=1/l *rep(1,l)
      }else{
        theta$rho_h=r_h/J
        theta$tau_l=t_l/K
      }
      if (normalization){
        mu_i=matrix(rowSums(x),ncol=1)
        nu_j=matrix(colSums(x),ncol=K)
      }
      else{mu_i=matrix(rep(1,J),ncol=1)
      nu_j=matrix(rep(1,K),ncol=K)
      }

      n_kl=t((t(mu_i)%*%r_jh))%*%(nu_j%*%t_kl)
      theta$gamma_hl= (t(r_jh)%*%x%*%t_kl)/ n_kl
      #if (any(is.nan(theta$gamma_hl))){
      # n_kl=r_h%*%t(t_l)
      #theta$gamma_hl=(t(r_jh)%*%data$x%*%t_kl)/(n_kl*(n_kl>0)+(n_kl<=0)*eps);
      #}
      w=start$w}

    } else if(!is.null(start$v) && !is.null(start$t_kl)){
      r_h=matrix(colSums(r_jh),ncol=1)
      t_l=matrix(colSums(t_kl),ncol=1)
      if (fixed_proportion==TRUE){
        theta$rho_h=1/h *rep(1,h)
        theta$tau_l=1/l *rep(1,l)
      }else{
        theta$rho_h=r_h/J
        theta$tau_l=t_l/K
      }

      if (normalization){
        mu_i=matrix(rowSums(x),ncol=1)
        nu_j=matrix(colSums(x),ncol=K)
      }
      else{mu_i=matrix(rep(1,J),ncol=1)
      nu_j=matrix(rep(1,K),ncol=K)
      }
      n_kl=t((t(mu_i)%*%r_jh))%*%(nu_j%*%t_kl)
      theta$gamma_hl= (t(r_jh)%*%x%*%t_kl)/ n_kl

      w=apply(start$t_kl,1,which.max)
    } else if(!is.null(start$r_jh) && !is.null(start$w)){
      r_jh=start$r_jh
      t_kl=!(start$w%*%t(rep(1,l))-rep(1,K)%*%t(1:l))

      # Computation of parameters
      r_h=colSums(r_jh)
      t_l=colSums(t_kl)

      if (fixed_proportion==TRUE){
        theta$rho_h=1/h *rep(1,h)
        theta$tau_l=1/l *rep(1,l)
      }else{
        theta$rho_h=r_h/J
        theta$tau_l=t_l/K
      }
      w=start$w

    } else if(!is.null(start$r_jh) && !is.null(start$t_kl)){

      r_jh=start$r_jh
      t_kl=start$t_kl

      # Computation of parameters
      r_h=matrix(colSums(r_jh),ncol=1)
      t_l=matrix(colSums(t_kl),ncol=1)

      if (fixed_proportion==TRUE){
        theta$rho_h=1/h *rep(1,h)
        theta$tau_l=1/l *rep(1,l)
      }else{
        theta$rho_h=r_h/J
        theta$tau_l=t_l/K
      }

      if (normalization){
        mu_i=matrix(rowSums(x),ncol=1)
        nu_j=matrix(colSums(x),ncol=K)
      }
      else{mu_i=matrix(rep(1,J),ncol=1)
      nu_j=matrix(rep(1,K),ncol=K)
      }
      n_kl=t((t(mu_i)%*%r_jh))%*%(nu_j%*%t_kl)
      theta$gamma_hl= (t(r_jh)%*%x%*%t_kl)/ n_kl
      w=apply(start$t_kl,1,which.max)

    }else if(!is.null(start$t_kl) && !is.null(start$w)&& !is.null(start$theta)){
      w=start$w
      t_kl=start$t_kl
      theta=start$theta
    }else{stop("For a not random initialization, it needs at least 2 parameters")}

    list(t_kl=t_kl,theta=theta,w=w)
  }

  BlocXemalpha0_kl_LBM =function(x,theta,t_kl,a,alpha,beta,normalization,niter=100000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
    gamma_hl=theta$gamma_hl
    rho_h=theta$rho_h
    tau_l=theta$tau_l
    h=dim(gamma_hl)[1]
    l=dim(gamma_hl)[2]
    J=dim(x)[1]
    K=dim(x)[2]
    t_l=t(colSums(t_kl))
    Ht= -sum(sum(t_kl*log(t_kl+(t_kl<=0))*(t_kl>0)))
    if (normalization==TRUE){
      mu_i=matrix(rowSums(x),ncol=1)
      nu_j=matrix(colSums(x),ncol=K)
    }else{mu_i=matrix(rep(1,J),ncol=1)
    nu_j=matrix(rep(1,K),ncol=K)
    }

    log_gamma_hl=log(gamma_hl+(gamma_hl<=0)*eps)
    log_rho_h=log(rho_h*(rho_h>0)+(rho_h<=0)*eps)
    log_tau_l=log(tau_l*(tau_l>0)+(tau_l<=0)*eps)
    # loop on the iterations of the lagorithm
    W=-Inf
    for (iter in 1:niter){

      W2=-Inf

      for (iter_int in 1:niter_int){
        # Computation of rjh
        R_jh=matrix(1,J,1)%*%t(log_rho_h)- (mu_i%*% ((nu_j%*%t_kl)%*%t(gamma_hl)))+x%*%t_kl%*%t(log_gamma_hl)
        R_jhmax=apply(R_jh, 1,max)
        Rp_jh = exp(R_jh-matrix(R_jhmax,J,1)%*%matrix(1,1,h))
        r_jh=Rp_jh/(rowSums(Rp_jh))
        Hs=-sum(sum(r_jh*(r_jh>0)*log(r_jh*(r_jh>0)+(r_jh<=0)*eps)))
        r_h=t(colSums(r_jh))

        # Computation of t_kl
        T_kl=matrix(1,K,1)%*%t(log_tau_l)- (t(nu_j)%*%((t(mu_i)%*%r_jh)%*%gamma_hl))+t(x)%*%r_jh%*%log_gamma_hl
        T_klmax=apply(T_kl, 1,max)
        Tp_kl = exp(T_kl-T_klmax)
        t_kl=Tp_kl/(rowSums(Tp_kl))
        Ht=-sum(sum(t_kl*(t_kl>0)*log(t_kl*(t_kl>0)+(t_kl<=0)*eps)))
        t_l=t(colSums(t_kl))

        W2_old=W2
        W2=sum(t_kl*t_kl) + r_h%*%log_rho_h  + Hs + Ht
        if (abs((W2-W2_old)/W2) < epsi_int) break
      }
      # Computation of parameters
      rho_h=t((r_h+a-1)/(J+h*(a-1)))

      tau_l=t((t_l+a-1)/(K+h*(a-1)))


      u_kl=t(r_jh)%*%x%*%t_kl
      n_kl=t(r_h)%*%t_l
      u_kl=u_kl*(u_kl>0)
      n_kl=n_kl*(n_kl>0)


      n_kl=t(t(mu_i)%*%r_jh)%*%(nu_j%*%t_kl)
      n_kl=n_kl*(n_kl>0)
      n_klx=t(r_jh)%*%x%*%t_kl

      if (any(any(n_kl<=0))){
        gamma_hl=(alpha-1+n_klx)/(beta+n_kl+((n_kl<=0)))*(n_kl>0)
      }else{
        gamma_hl=(alpha-1+n_klx)/(beta+n_kl)
      }
      if (any(any(gamma_hl<=0))){
        gamma_hl=gamma_hl*(gamma_hl>0)
        log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
      }else {
        log_gamma_hl=log(gamma_hl)
      }

      log_rho_h=log(rho_h*(rho_h>0)+(rho_h<=0)*eps)
      log_tau_l=log(tau_l*(tau_l>0)+(tau_l<=0)*eps)

      # Criterion
      W_old=W
      W=sum(r_h*log(r_h))- J*log(J) + sum(t_l*log(t_l)) - K*log(K) -
        sum(sum(n_kl*gamma_hl))+ sum(sum(n_klx*log_gamma_hl))+ Hs + Ht +
        (a-1)*(sum(log(rho_h))+sum(log(tau_l)))-beta*sum(sum(gamma_hl))
      if (alpha!=1){
        W=W+(alpha-1)*sum(sum(log_gamma_hl))
      }
      if (is.nan(W)){
        if (any(any(gamma_hl<0))){
          gamma_hl=gamma_hl*(gamma_hl>0)
          log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
        }
        W=sum(r_h*log(r_h+(r_h<=0)*eps))- J*log(J) + sum(t_l*log(t_l+(t_l<=0)*eps)) - K*log(K) -
          sum(sum(n_kl*gamma_hl))+ sum(sum(n_klx*log_gamma_hl))+ Hs + Ht +
          (a-1)*(sum(log(rho_h))+sum(log(tau_l)))-beta*sum(sum(gamma_hl))
        if (alpha!=1){
          W=W+(alpha-1)*sum(sum(log_gamma_hl))
        }
      }
      if (abs((W-W_old)/W) < epsi){ break}

      if (is.nan(W)) {break}
    }

    theta$gamma_hl=gamma_hl
    theta$rho_h=rho_h
    theta$tau_l=tau_l

    list(W=W,theta=theta,r_jh=r_jh,t_kl=t_kl,iter=iter)
  }




  BlocGibbs=function(x,a,normalization,alpha,beta,ini=list(n_init=10,ale=FALSE),niter=1000,eps=10^(-16),classement){
    if (is.null(ini$theta)){
      if (!is.numeric(ini$n_init)){
        stop("For a random initialization n_init must be numeric")
      }
      if (ini$n_init<=0){
        stop("n_init must be positive")
      }
      W=-Inf
      Res=list()
      #   for (i in 1:ini$n_init){
      #     init=PoissonBlocInit(h,l,data,normalisation,start=ini)
      #   }
      # }else{init=ini}
    }
    theta=ini$theta
    t_kl=ini$t_kl
    eps10=10*eps
    h=dim(theta$gamma_hl)[1]
    l=dim(theta$gamma_hl)[2]
    J=dim(x)[1]
    K=dim(x)[2]

    # loop of the iterations of the algorithm
    v=matrix(0,J,1)
    w=matrix(0,K,1)
    gamma_hl_path=list()
    rho_h_path=matrix(0,h,niter)
    tau_l_path=matrix(0,l,niter)
    Part=list()
    Part$v=matrix(0,nrow=J,ncol=h)
    Part$w=matrix(0,nrow=K,ncol=l)

    #SEM is random : 1000 runs from the same initial value

    gamma_hl=theta$gamma_hl
    rho_h=theta$rho_h
    tau_l=theta$tau_l
    t_l=colSums(t_kl)
    if (normalization){
      mu_i=matrix(rowSums(x),ncol=1)
      nu_j=matrix(colSums(x),ncol=K)
    }else{mu_i=matrix(rep(1,J),ncol=1)
    nu_j=matrix(rep(1,K),ncol=K)
    }

    for (iter in 1:niter){
      log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
      log_rho_h=log(rho_h)
      log_tau_l=log(tau_l)

      # Computation of v
      R_jh=matrix(1,J,1)%*%t(log_rho_h)- (mu_i%*% ((nu_j%*%t_kl)%*%t(gamma_hl)))+x%*%t_kl%*%t(log_gamma_hl)
      R_jhmax=apply(R_jh, 1,max)
      Rp_jh = exp(R_jh-R_jhmax)
      r_jh=Rp_jh/(rowSums(Rp_jh))
      for (iter_i in 1:J){
        v[iter_i]=min(h,h+1- sum( matrix(1,1,h)*runif(1) < cumsum(r_jh[iter_i,]) ))
      }
      r_jhn=!(v%*%matrix(1,1,h)-matrix(1,J,1)%*%(1:h))
      r_hn=colSums(r_jhn)
      r_jh=r_jhn
      r_h=r_hn

      # Computation of w
      t_kl=matrix(1,K,1)%*%t(log_tau_l)- (t(nu_j)%*%((t(mu_i)%*%r_jh)%*%gamma_hl))+t(x)%*%r_jh%*%log_gamma_hl
      t_klmax=apply(t_kl, 1,max)
      Tp_jl = exp(t_kl-t_klmax)
      t_kl=Tp_jl/(rowSums(Tp_jl))
      for (iter_j in 1:K){
        w[iter_j]=min(l,l+1- sum(matrix(1,1,l)*runif(1) < cumsum(t_kl[iter_j,])))
      }
      t_kln=!(w%*%matrix(1,1,l)-matrix(1,K,1)%*%(1:l))
      t_ln=colSums(t_kln)
      t_kl=t_kln
      t_l=t_ln


      # Computation of parameters
      nv=r_h+a
      rho_h=rgamma(h,nv,1)
      rho_h=rho_h/(sum(rho_h))
      dw=t_l+a
      tau_l=rgamma(l,dw,1)
      tau_l=tau_l/(sum(tau_l))

      n_kl=t(t(mu_i)%*%r_jh)%*%(nu_j%*%t_kl)
      n_klx=t(r_jh)%*%x%*%t_kl
      gamma_hl=matrix(rgamma(h*l,shape=alpha+n_klx,scale=1/(beta+n_kl)),ncol=l)

      # trajectoires
      gamma_hl_path=c(gamma_hl_path,list(gamma_hl))
      rho_h_path[,iter]=rho_h
      tau_l_path[,iter]=tau_l
      Part$v=Part$v+r_jh
      Part$w=Part$w+t_kl
    }


    # arguments de sortie
    theta_path=list()
    theta_path$gamma_hl=gamma_hl_path
    theta_path$rho_h=rho_h_path
    theta_path$tau_l=tau_l_path

    list(theta_path=theta_path,Part=Part)

  }

  PoissonClassement=function(theta,Part=list()){
    ind_i=order(as.vector(theta$gamma_hl%*%theta$tau_l))
    ind_j=order(as.vector(t(theta$rho_h)%*%theta$gamma_hl))
    if (is.null(Part$v)||is.null(Part$w)){
      list(rho_h=theta$rho_h[ind_i],tau_l=theta$tau_l[ind_j],gamma_hl=theta$gamma_hl[ind_i,ind_j])
    }else{
      v=matrix(0,length(Part$v),1)
      maxv=max(Part$v)
      for (i in 1:maxv){
        v=v+ind_i[i]*(Part$v==i)
      }
      w=matrix(0,length(Part$w),1)
      maxw=max(Part$w)
      for (j in 1:maxw){
        w=w+ind_j[j]*(Part$w==j)
      }
      thetares=list(rho_h=theta$rho_h[ind_i],tau_l=theta$tau_l[ind_j],gamma_hl=theta$gamma_hl[ind_i,ind_j])
      list(theta=thetares,v=v,w=w,ind_i=ind_i,ind_j=ind_j)
    }
  }

  PoissonBlocVBayes=function(x,h,l,a,alpha,beta,normalization,ini=list(n_init=1,ale=TRUE),niter=10000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
    if (is.null(ini$theta)){
      if (!is.numeric(ini$n_init)){
        stop("For a random initialization, n_init must be numeric")
      }
      if (ini$n_init<=0){
        stop("n_init must be positive")
      }
    }
    W=-Inf
    Res=list()
    for (itry in 1:ini$n_init){
      initstart=PoissonBlocInit(h,l,x,normalization,start=ini)
      ResVBayesAle=BlocXemalpha0_kl_LBM(x,initstart$theta,initstart$t_kl,a=a,alpha=alpha,beta=beta,normalization=normalization,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)
      if (is.nan(ResVBayesAle$W)){

        stop('erreur avec le prior a')}

      if (ResVBayesAle$W > W){
        W=ResVBayesAle$W
        theta_max=ResVBayesAle$theta
        r_jh_max=ResVBayesAle$r_jh
        t_kl_max=ResVBayesAle$t_kl
      }
      if (is.nan(ResVBayesAle$W)&&(itry==ini$n_init)&&(W==-Inf)){
        W=ResVBayesAle$W
        theta_max=ResVBayesAle$theta
        r_jh_max=ResVBayesAle$r_jh
        t_kl_max=ResVBayesAle$t_kl
      }
    }

    theta_path=0
    theta=theta_max
    r_jh=r_jh_max
    t_kl=t_kl_max
    v=apply(r_jh, 1, which.max)
    w=apply(t_kl, 1, which.max)
    iter=ResVBayesAle$iter
    empty_cluster=(length(unique(v))!=h)||(length(unique(w))!=l)
    list(W=W,theta=theta,r_jh=r_jh,t_kl=t_kl,iter=iter,empty_cluster=empty_cluster,v=v,w=w)
  }
  PoissonBlocGibbs=function(x,h,l,a,normalization,alpha=1,beta=0.01,niter=5000,inistart=list(n_init=10,ale=TRUE),classement=TRUE,eps=10^(-16)){
    if (is.null(inistart$theta)){
      initgibbs=PoissonBlocInit(h,l,x,normalization,start=inistart)
      initgibbs$n_init=inistart$n_init
      initgibbs$ale=inistart$ale
    }else{
      initgibbs=inistart
    }
    res=BlocGibbs(x,a=a,normalization=normalization,alpha=alpha,beta=beta,ini=initgibbs,niter=niter,eps=eps,classement=classement)
    h=dim(res$theta_path$rho_h)[1]
    l=dim(res$theta_path$tau_l)[1]
    theta=list(rho_h=matrix(0,h,1),tau_l=matrix(0,l,1),gamma_hl=matrix(0,h,l))
    if (classement){
      for (i in 1:niter){
        thetabis=list(gamma_hl=res$theta_path$gamma_hl[[i]],rho_h=res$theta_path$rho_h[,i],tau_l=res$theta_path$tau_l[,i])
        thetabisorg=PoissonClassement(thetabis)
        theta$rho_h=theta$rho_h+thetabisorg$rho_h
        theta$tau_l=theta$tau_l+thetabisorg$tau_l
        theta$gamma_hl=theta$gamma_hl+thetabisorg$gamma_hl
      }
    }else{
      for (i in 1:niter){
        theta$gamma_hl=theta$gamma_hl+res$theta_path$gamma_hl[[i]]
      }
      theta$rho_h=rowSums(res$theta_path$rho_h)
      theta$tau_l=rowSums(res$theta_path$tau_l)
    }
    theta$rho_h=theta$rho_h/niter
    theta$tau_l=theta$tau_l/niter
    theta$gamma_hl=theta$gamma_hl/niter
    v=apply(res$Part$v, 1, which.max)
    w=apply(res$Part$w, 1, which.max)
    list(theta=theta,v=v,w=w,Part=res$Part,theta_path=res$theta_path)
  }

  PoissonBlocICL =function (a,alpha,beta,x,v1,w1,normalization){


    J=dim(x)[1]
    K=dim(x)[2]
    if (!is.matrix(v1)){
      h=max(v1)
      v=!(v1%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
      #v=!(v1%*%matrix(1,1,h)-matrix(1,J,1)%*%c(1:h))
    }else{
      v=v1
    }
    if (!is.matrix(w1)){
      l=max(w1)
      w=!(w1%*%t(rep(1,l))-rep(1,K)%*%t(1:l))
      #w=!(w1%*%matrix(1,1,l)-matrix(1,K,1)%*%c(1:l))
    }else{
      w=w1
    }

    if (normalization){
      mu_i=matrix(rowSums(x),ncol=1)
      nu_j=matrix(colSums(x),ncol=K)
    }
    else{mu_i=matrix(rep(1,J),ncol=1)
    nu_j=matrix(rep(1,K),ncol=K)
    }
    nh=sum(v)
    nl=sum(w)
    nhi=matrix(colSums(v*(mu_i%*%matrix(1,1,h))),ncol=h)
    nlj=matrix(colSums(w*(t(nu_j)%*%matrix(1,1,l))),ncol=l)
    nhl=t(nhi)%*%nlj
    nhlx=t(v)%*%x%*%w
    critere=lgamma(a*h)+lgamma(a*l)-(l+h)*lgamma(a)-h*l*lgamma(alpha)-lgamma(J+h*a)-lgamma(K+l*a)+
      sum(lgamma(nh+a))+sum(lgamma(nl+a))+sum(sum(-(alpha+nhlx)*log(beta+nhl)+lgamma(alpha+nhlx)))

    if (beta>0){
      critere=critere+h*l*alpha*log(beta)
    }
    critere}

  PoissonBlocBIC =function (a,alpha,beta,x,res,normalization){
    v1=res$v
    w1=res$w
    theta=res$theta

    J=dim(x)[1]
    K=dim(x)[2]
    if (!is.matrix(v1)){
      h=max(v1)
      v=!(v1%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
    }else{
      v=v1
    }
    if (!is.matrix(w1)){
      l=max(w1)
      w=!(w1%*%t(rep(1,l))-rep(1,K)%*%t(1:l))
    }else{
      w=w1
    }

    if (normalization){
      mu_i=matrix(rowSums(x),ncol=1)
      nu_j=matrix(colSums(x),ncol=K)
    }
    else{mu_i=matrix(rep(1,J),ncol=1)
    nu_j=matrix(rep(1,K),ncol=K)
    }
    nh=sum(v)
    nl=sum(w)
    nhi=sum(v*(mu_i%*%matrix(1,1,h)))
    nlj=sum(w*(t(nu_j)%*%matrix(1,1,l)))
    nhl=t(nhi)%*%nlj
    nhlx=t(v)%*%x%*%w
    if (any(any(theta$gamma_hl<=0))){
      gamma_hl=theta$gamma_hl*(theta$gamma_hl>0)
      log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
    }else{
      gamma_hl=theta$gamma_hl
      log_gamma_hl=log(theta$gamma_hl)
    }
    if ((a==1)){
      W=res$W-((alpha-1)*sum(sum(log_gamma_hl))-beta*sum(sum(gamma_hl)))}
    else{
      W=res$W-((a-1)*(sum(log(theta$rho_h))+sum(log(theta$tau_l)))+(alpha-1)*sum(sum(log_gamma_hl))-beta*sum(sum(gamma_hl)))
      if (is.nan(W)){
        W=res$W-((a-1)*(sum(log(theta$rho_h+(theta$rho_h<=0))*(theta$rho_h>0))+sum(log(theta$tau_l+(theta$tau_l<=0))*(theta$tau_l>0)))+(alpha-1)*sum(sum(log_gamma_hl))-beta*sum(sum(gamma_hl)))}
    }
    critere=W-((h-1)*log(J)+(l-1)*log(K)+l*h*log(J*K))/2
    critere}










  ### Verif
  if (!is.numeric(mc.cores)){
    stop("mc.cores must be an integer greater than 1")
  } else if ((mc.cores<=0)||(length(mc.cores)!=1)||(floor(mc.cores)!=mc.cores)){
    stop("mc.cores must be an integer greater than 1")
  }

  ## =============================================================
  ## INITIALIZATION & PARAMETERS RECOVERY
  if ((Sys.info()[['sysname']] == "Windows")&(mc.cores>1)) {
    warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
    mc.cores <- 1
  }


  #source("Chargement.R")

  criterion_tab=matrix(-Inf,Hmax+1,Lmax+1)
  W_tab=matrix(-Inf,Hmax+1,Lmax+1)
  W=-Inf
  if (init_choice=='smallVBayes'){
    # if (mc.cores>1){
    #   Temp=parallel::mclapply(1:ntry,function(i){PoissonBlocVBayes(data,Hstart,Lstart,a,alpha,beta,
    #                                                    normalization,ini=list(n_init=1,ale=TRUE),
    #                                                    niter=50)}, mc.cores=mc.cores)
    # }else{
    #   Temp=lapply(1:ntry,function(i){PoissonBlocVBayes(data,Hstart,Lstart,a,alpha,beta,
    #                                                    normalization,ini=list(n_init=1,ale=TRUE),
    #                                                    niter=50)})
    # }

    for (i in 1:ntry){
      #res=BlocGibbs(data$x,a,b,ini=list(n_init=1,ale=TRUE),1000)
      res=PoissonBlocVBayes(x,Hstart,Lstart,a,alpha,beta,
                            normalization,ini=list(n_init=1,ale=TRUE),
                            niter=50)
      #print(W)
      if (res$W>W && res$empty_cluster==FALSE){
        W=res$W
        resopt=res
      }
      if (!exists("resopt")){
        stop(paste('The algorithm does not manage to find a configuration with  (Hstart=',as.character(Hstart),',Lstart=',as.character(Lstart),'), try normalization=FALSE or a higher ntry') )
      }
    }
  }else if (init_choice=='user'){
    if (is.null(userparam)){
      stop(paste('If you choose the user init choice, please fill the useparam field with partitions v and w' ))
    }else{
      initVBayes=list(v=userparam$v,w=userparam$w,ale=FALSE,n_init=1)
      resopt=PoissonBlocVBayes(x,Hstart,Lstart,a,alpha,beta,normalization,ini=initVBayes)

    }


  }else if (init_choice=='Gibbs'){

    resinit=PoissonBlocGibbs(x,Hstart,Lstart,a,normalization,alpha,beta,niter=5000,inistart=list(n_init=10,ale=TRUE))
    initVBayes=list(v=resinit$v,w=resinit$w,ale=FALSE,n_init=1)
    resopt=PoissonBlocVBayes(x,Hstart,Lstart,a,alpha,beta,normalization,ini=initVBayes)
    if (resopt$empty_cluster==TRUE){
      stop(paste('The algorithm does not manage to find a configuration with  (Hstart=',as.character(Hstart),',Lstart=',as.character(Lstart),'),try normalization=FALSE') )
    }

  }else if (init_choice=='random'){

    eps=10^(-16)
    J=dim(x)[1]
    K=dim(x)[2]
    # Random_Temp=function(Temp){
    #   theta=list()
    #   theta$rho_h=runif(Hstart)
    #   theta$rho_h=theta$rho_h/sum(theta$rho_h)
    #   theta$tau_l=runif(Lstart)
    #   theta$tau_l=theta$tau_l/sum(theta$tau_l)
    #   resw=PartRnd(K,theta$tau_l)
    #   t_kl=resw$v_jh
    #   w=resw$v
    #   resv=PartRnd(J,theta$rho_h)
    #   r_jh=resv$v_jh
    #   v=resv$v
    #   if (normalization){
    #     mu_i=matrix(rowSums(data$x),ncol=1)
    #     nu_j=matrix(colSums(data$x),ncol=K)
    #   }
    #   else{mu_i=matrix(rep(1,J),ncol=1)
    #   nu_j=matrix(rep(1,K),ncol=K)
    #   }
    #   n_kl=t((t(mu_i)%*%r_jh))%*%(nu_j%*%t_kl)
    #   theta$gamma_hl=(t(r_jh)%*%data$x%*%t_kl)/(n_kl)
    #   if (any(is.nan(theta$gamma_hl))){
    #     lambda0=mean(mean(data$x))
    #     theta$gamma_hl=theta$gamma_hl=matrix(rpois(h,lambda0),ncol=l)
    #   }
    #
    #   initVBayes=list(t_kl=t_kl,w=w,theta=theta,ale=FALSE,n_init=1)
    #   #res=BlocGibbs(data$x,a,b,ini=list(n_init=1,ale=TRUE),1000)
    #   PoissonBlocVBayes(data,Hstart,Lstart,a,alpha,beta,normalization,ini=initVBayes)
    # }
    #
    # if (mc.cores>1){
    #   Temp=parallel::mclapply(1:ntry,Random_Temp, mc.cores=mc.cores)
    # }else{
    #   Temp=lapply(1:ntry,Random_Temp)
    # }


    for (i in 1:ntry){
      theta=list()
      theta$rho_h=runif(Hstart)
      theta$rho_h=theta$rho_h/sum(theta$rho_h)
      theta$tau_l=runif(Lstart)
      theta$tau_l=theta$tau_l/sum(theta$tau_l)
      resw=PartRnd(K,theta$tau_l)
      t_kl=resw$v_jh
      w=resw$v
      resv=PartRnd(J,theta$rho_h)
      r_jh=resv$v_jh
      v=resv$v
      if (normalization){
        mu_i=matrix(rowSums(x),ncol=1)
        nu_j=matrix(colSums(x),ncol=K)
      }
      else{mu_i=matrix(rep(1,J),ncol=1)
      nu_j=matrix(rep(1,K),ncol=K)
      }
      n_kl=t((t(mu_i)%*%r_jh))%*%(nu_j%*%t_kl)
      theta$gamma_hl=(t(r_jh)%*%x%*%t_kl)/(n_kl)
      if (any(is.nan(theta$gamma_hl))){
        lambda0=mean(mean(x))
        theta$gamma_hl=theta$gamma_hl=matrix(rpois(h,lambda0),ncol=l)
      }

      initVBayes=list(t_kl=t_kl,w=w,theta=theta,ale=FALSE,n_init=1)
      #res=BlocGibbs(data$x,a,b,ini=list(n_init=1,ale=TRUE),1000)
      res=PoissonBlocVBayes(x,Hstart,Lstart,a,alpha,beta,normalization,ini=initVBayes)
      if (res$W>W && res$empty_cluster==FALSE){
        resopt=res}



    }
    if (!exists("resopt")){
      stop(paste('The algorithm does not manage to find a configuration with  (Hstart=',as.character(Hstart),',Lstart=',as.character(Lstart),'), try normalization=FALSE') )
    }
  } else {stop("This is not a initialization choice proposed by the procedure")}

  if (criterion_choice=='ICL'){
    Critere= PoissonBlocICL(a,alpha,beta,x,resopt$v,resopt$w,normalization)  }
  else if (criterion_choice=='BIC'){
    Critere= PoissonBlocBIC(a,alpha,beta,x,resopt,normalization)
  }
  else {stop('This is not a criterion choice proposed by the procedure')}
  criterion_tab[Hstart,Lstart]=Critere
  W_tab[Hstart,Lstart]=resopt$W

  modele=resopt
  criterion_max=Critere
  h=Hstart
  l=Lstart
  hopt=Hstart
  lopt=Lstart
  model_max=resopt

  while (h+1<=Hmax && l+1<=Lmax){
    if ((length(unique(modele$v))!=h)||(length(unique(modele$w))!=l)){
      stop(paste('The algorithm does not manage to find a configuration with (h=',as.character(h),',l=',as.character(l),'), maybe reduce Hmax and Lmax') )
    }
    W_colonne=-Inf
    W_ligne=-Inf

    f1=function(Kc){
      w=modele$w
      Xc=which(w==Kc)
      lxc=length(Xc)
      if (lxc>1){
        idxc=sample(lxc,floor(lxc/2),replace=FALSE)
        indc=Xc[idxc]
        w[indc]=l+1
        initVBayes=list(v=modele$v,w=w,ale=FALSE,n_init=1)
        resultsVBayes=PoissonBlocVBayes(x,h,l+1,a,alpha,beta,normalization,ini=initVBayes)
        iterverif=0
        while (resultsVBayes$empty_cluster==TRUE && iterverif<20){
          iterverif=iterverif+1
          w=modele$w
          idxc=sample(lxc,floor(lxc/2),replace=FALSE)
          indc=Xc[idxc]
          w[indc]=l+1
          initVBayes=list(v=modele$v,w=w,ale=FALSE,n_init=1)
          resultsVBayes=PoissonBlocVBayes(x,h,l+1,a,alpha,beta,normalization,ini=initVBayes)
        }
        if (iterverif<20){
          return(list(W=resultsVBayes$W,modele_colonne=resultsVBayes))
        }else{
          return(list(W=-Inf))
        }
      }else{
        return(list(W=-Inf))
      }
    }
    if (mc.cores>1){
      l_colonne=parallel::mclapply(1:l,f1, mc.cores=mc.cores)
    }else{
      l_colonne=lapply(1:l,f1)
    }

    for (Kc in 1:l){
      if (l_colonne[[Kc]]$W>W_colonne){
        W_colonne=l_colonne[[Kc]]$W
        modele_colonne=l_colonne[[Kc]]$modele_colonne
      }
    }
    if (W_colonne<=-Inf){
      Empty_m=TRUE
    }else{
      Empty_m=FALSE
    }
    #### Class h
    f2=function(Kl){
      v=modele$v
      Xl=which(v==Kl)
      ll=length(Xl)
      if (ll>1){
        idxl=sample(ll,floor(ll/2),replace=FALSE)
        indl=Xl[idxl]
        v[indl]=h+1
        initVBayes=list(v=v,w=modele$w,ale=FALSE,n_init=1)
        resultsVBayes=PoissonBlocVBayes(x,h+1,l,a,alpha,beta,normalization,ini=initVBayes)
        iterverif=0
        while (resultsVBayes$empty_cluster==TRUE && iterverif<20){
          iterverif=iterverif+1
          v=modele$v
          idxl=sample(ll,floor(ll/2),replace=FALSE)
          indl=Xl[idxl]
          v[indl]=h+1
          initVBayes=list(v=v,w=modele$w,ale=FALSE,n_init=1)
          resultsVBayes=PoissonBlocVBayes(x,h+1,l,a,alpha,beta,normalization,ini=initVBayes)
        }
        if (iterverif<20){
          return(list(W=resultsVBayes$W,modele_ligne=resultsVBayes))
        }else{
          return(list(W=-Inf))
        }
      }else{
        return(list(W=-Inf))
      }


    }

    if (mc.cores>1){
      l_ligne=parallel::mclapply(1:h,f2, mc.cores=mc.cores)
    }else{
      l_ligne=lapply(1:h,f2)
    }
    for (Kl in 1:h){
      if (l_ligne[[Kl]]$W>W_ligne){
        W_ligne=l_ligne[[Kl]]$W
        modele_ligne=l_ligne[[Kl]]$modele_ligne
      }
    }

    if ((W_ligne<=-Inf)&(Empty_m)){
      warning("The algorithm does not manage to find a configuration  (",as.character(h+1),',',as.character(l), ") nor (",as.character(h),',',as.character(l+1), ")")
      new(Class="BIKM1_LBM_Poisson",init_choice=init_choice,criterion_choice=criterion_choice,criterion_tab=criterion_tab,criterion_max=criterion_max,model_max=model_max,W_tab=W_tab,hopt=hopt,lopt=lopt)
    }


    if (W_ligne>W_colonne){
      modele=modele_ligne
      if (criterion_choice=='ICL'){
        Critere= PoissonBlocICL(a,alpha,beta,x,modele$v,modele$w,normalization)
      }else if (criterion_choice=='BIC'){
        Critere= PoissonBlocBIC(a,alpha,beta,x,modele,normalization)
      }
      criterion_tab[h+1,l]=Critere
      W_tab[h+1,l]=modele$W
      if (Critere>criterion_max){
        model_max=modele
        criterion_max=Critere
        hopt=h+1
        lopt=l
      }
      if (verbose){
        cat('(h,l)=(',as.character(h+1),',',as.character(l),')\n',sep = "")
      }
      h=h+1
    }else {
      modele=modele_colonne
      if (criterion_choice=='ICL'){
        Critere= PoissonBlocICL(a,alpha,beta,x,modele$v,modele$w,normalization)
      }else if (criterion_choice=='BIC'){
        Critere= PoissonBlocBIC(a,alpha,beta,x,modele,normalization)
      }
      criterion_tab[h,l+1]=Critere
      W_tab[h,l+1]=modele$W
      if (Critere>criterion_max){
        model_max=modele
        criterion_max=Critere
        hopt=h
        lopt=l+1

      }
      if (verbose){
        cat('(h,l)=(',as.character(h),',',as.character(l+1),')\n',sep = "")
      }

      l=l+1
    }
  }
  if (verbose){
    cat('The selected row (h) and column (l) clusters are (h,l)=(',as.character(hopt),',',as.character(lopt),')\n',sep = "")
  }

   new(Class="BIKM1_LBM_Poisson",init_choice=init_choice,criterion_choice=criterion_choice,criterion_tab=criterion_tab,criterion_max=criterion_max,model_max=model_max,W_tab=W_tab,hopt=hopt,lopt=lopt)
  } else {

    PartRnd = function(J,proba){

      h=length(proba)

      i=sample(1:J,J, replace = FALSE)
      i1=i[1:h]
      i2=i[(h+1):J]
      v=rep(0,J)
      v[i1]=1:h

      if (J > h) {
        v[i2]=(h+1)-colSums(rep(1,h)%*%t(runif(J-h)) < cumsum(proba)%*%t(rep(1,J-h)))
        v[v==(h+1)]=h
      }

      v_jh=!(v%*%t(rep(1,h))-rep(1,J)%*%t(1:h))

      list(v=v,v_jh=v_jh)
    }



    PoissonBlocInit=function(h,l,x,normalization,start,fixed_proportion=FALSE){
      eps=10^(-16)
      J=dim(x)[1]
      K=dim(x)[2]

      # Initialization by random simulations  of centers

      theta=list()
      if (start$ale==TRUE){
        if (normalization){
        }
        theta$rho_h=1/h*matrix(1,h,1)
        #theta$rho_h=theta$rho_h/sum(theta$rho_h)
        theta$tau_l=1/l*matrix(1,l,1)
        #theta$tau_l=theta$tau_l/sum(theta$tau_l)
        resw=PartRnd(K,theta$tau_l)
        t_kl=resw$v_jh
        w=resw$v
        resv=PartRnd(J,theta$rho_h)
        r_jh=resv$v_jh
        v=resv$v
        r_h=matrix(colSums(r_jh),ncol=1)
        t_l=matrix(colSums(t_kl),ncol=1)
        n_kl=r_h%*%t(t_l)
        theta$gamma_hl=(t(r_jh)%*%x%*%t_kl)/(n_kl)
        if (any((is.nan(theta$gamma_hl)))){
          lambda0=mean(mean(x))
          theta$gamma_hl=matrix(rpois(h,lambda0),ncol=l)
        }

      }else if(!is.null(start$theta) && !is.null(start$w)){
        theta=start$theta
        t_kl=!(start$w%*%t(rep(1,l))-rep(1,K)%*%t(1:l))
        w=start$w
      } else if(!is.null(start$theta) && !is.null(start$t_kl)){
        theta=start$theta
        t_kl=start$t_kl
        w=apply(start$t_kl,1,which.max)
      } else if(!is.null(start$theta) && !is.null(start$v)){
        theta=start$theta
        #r_jh=!(init$v*matrix(1,1,h)-matrix(1,J,1)*c(1:h))
        r_jh=!(start$v%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
        r_h=matrix(colSums(r_jh),ncol=1)
        eps10=10*eps
        log_tau_l=log(theta$tau_l)
        gamma_hl=theta$gamma_hl
        B_kl=log(1-gamma_hl+eps10)
        #t_kl=(x_ij'*r_jh)*A_kl + ones(K,1)*r_h'*B_kl+  ones(K,1)*log_tau_l'
        log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
        t_kl=matrix(1,K,1)%*%t(log_tau_l)- matrix(1,K,1)%*%t(r_h)%*%gamma_hl+t(x)%*%r_jh%*%log_gamma_hl
        t_klmax=apply(t_kl,1,max)
        Tp_jl = exp(t_kl-t_klmax)
        t_kl=Tp_jl/(rowSums(Tp_jl))
        w=apply(t_kl,1,which.max)
      } else if(!is.null(start$theta) && !is.null(start$r_jh)){
        r_jh=start$r_jh
        r_h=matrix(colSums(r_jh),ncol=1)
        eps10=10*eps
        log_tau_l=log(theta$tau_l)
        gamma_hl=theta$gamma_hl
        log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
        t_kl=matrix(1,K,1)%*%t(log_tau_l)- matrix(1,K,1)%*%t(r_h)%*%gamma_hl+t(x)%*%r_jh%*%log_gamma_hl
        t_klmax=apply(t_kl,1,max)
        Tp_jl = exp(t_kl-t_klmax)
        t_kl=Tp_jl/(rowSums(Tp_jl))
        w=apply(t_kl,1,which.max)
      } else if(!is.null(start$v) && !is.null(start$w)){

        if ((max(start$v)!=h)||(max(start$w)!=l)){
          warning('(Hstart,Lstart) need to match with the characteristics of userparam v and w')
        }else{
        r_jh=!(start$v%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
        t_kl=!(start$w%*%t(rep(1,l))-rep(1,K)%*%t(1:l))

        r_h=matrix(colSums(r_jh),ncol=1)
        t_l=matrix(colSums(t_kl),ncol=1)
        if (fixed_proportion==TRUE){
          theta$rho_h=1/h*matrix(1,h,1)
          theta$tau_l=1/l*matrix(1,l,1)
        }else{
          theta$rho_h=r_h/J
          theta$tau_l=t_l/K
        }


        n_kl=r_h%*%t(t_l)
        theta$gamma_hl= (t(r_jh)%*%x%*%t_kl)/ n_kl
        #if (any(is.nan(theta$gamma_hl))){
        # n_kl=r_h%*%t(t_l)
        #theta$gamma_hl=(t(r_jh)%*%data$x%*%t_kl)/(n_kl*(n_kl>0)+(n_kl<=0)*eps);
        #}
        w=start$w}

      } else if(!is.null(start$v) && !is.null(start$t_kl)){
        r_h=matrix(colSums(r_jh),ncol=1)
        t_l=matrix(colSums(t_kl),ncol=1)
        if (fixed_proportion==TRUE){
          theta$rho_h=1/h*matrix(1,h,1)
          theta$tau_l=1/l*matrix(1,l,1)
        }else{
          theta$rho_h=r_h/J
          theta$tau_l=t_l/K
        }


        n_kl=r_h%*%t(t_l)
        theta$gamma_hl= (t(r_jh)%*%x%*%t_kl)/ n_kl

        w=apply(start$t_kl,1,which.max)
      } else if(!is.null(start$r_jh) && !is.null(start$w)){
        r_jh=start$r_jh
        t_kl=!(start$w%*%t(rep(1,l))-rep(1,K)%*%t(1:l))

        # Computation of parameters
        r_h=matrix(colSums(r_jh),ncol=1)
        t_l=matrix(colSums(t_kl),ncol=1)

        if (fixed_proportion==TRUE){
          theta$rho_h=1/h*matrix(1,h,1)
          theta$tau_l=1/l*matrix(1,l,1)
        }else{
          theta$rho_h=r_h/J
          theta$tau_l=t_l/K
        }
        w=start$w
      } else if(!is.null(start$r_jh) && !is.null(start$t_kl)){

        r_jh=start$r_jh
        t_kl=start$t_kl

        # Computation of parameters
        r_h=matrix(colSums(r_jh),ncol=1)
        t_l=matrix(colSums(t_kl),ncol=1)

        if (fixed_proportion==TRUE){
          theta$rho_h=1/h*matrix(1,h,1)
          theta$tau_l=1/l*matrix(1,l,1)
        }else{
          theta$rho_h=r_h/J
          theta$tau_l=t_l/K
        }

        n_kl=r_h%*%t(t_l)
        theta$gamma_hl= (t(r_jh)%*%x%*%t_kl)/ n_kl
        w=apply(start$t_kl,1,which.max)
      }else if(!is.null(start$t_kl) && !is.null(start$w)&& !is.null(start$theta)){
        w=start$w
        t_kl=start$t_kl
        theta=start$theta
      }else{stop("For an initialization not random, it needs at least 2 parameters")}

      list(t_kl=t_kl,theta=theta,w=w)
    }

    BlocXemalpha1_kl_LBM=function(x,theta,t_kl,a,alpha,beta,normalization,niter=100000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
      gamma_hl=theta$gamma_hl
      rho_h=theta$rho_h
      tau_l=theta$tau_l
      h=dim(gamma_hl)[1]
      l=dim(gamma_hl)[2]
      J=dim(x)[1]
      K=dim(x)[2]
      t_l=matrix(colSums(t_kl),ncol=1)
      Ht= -sum(sum(t_kl*log(t_kl+(t_kl<=0))*(t_kl>0)))
      if (normalization==TRUE){

      }

      log_gamma_hl=log(gamma_hl+(gamma_hl<=0)*eps)
      log_rho_h=log(rho_h*(rho_h>0)+(rho_h<=0)*eps)
      log_tau_l=log(tau_l*(tau_l>0)+(tau_l<=0)*eps)
      # loop on the iterations of the lagorithm
      W=-Inf
      for (iter in 1:niter){

        W2=-Inf

        for (iter_int in 1:niter_int){
          # Computation of sik
          R_jh=matrix(1,J,1)%*%(t(log_rho_h-gamma_hl%*%t_l))+x%*%t_kl%*%t(log_gamma_hl)
          #r_jh=%*%t(log_rho_h)- matrix(1,J,1)%*%t_l%*%t(gamma_hl)+
          R_jhmax=apply(R_jh, 1,max)
          Rp_jh = exp(R_jh-matrix(R_jhmax,J,1)%*%matrix(1,1,h))
          r_jh=Rp_jh/(rowSums(Rp_jh))
          Hs=-sum(sum(r_jh*(r_jh>0)*log(r_jh*(r_jh>0)+(r_jh<=0)*eps)))
          r_h=matrix(colSums(r_jh),ncol=1)

          # Computation of t_kl
          T_kl=matrix(1,K,1)%*%t(log_tau_l)- matrix(1,K,1)%*%t(r_h)%*%gamma_hl+t(x)%*%r_jh%*%log_gamma_hl
          T_klmax=apply(T_kl, 1,max)
          Tp_kl = exp(T_kl-T_klmax)
          t_kl=Tp_kl/(rowSums(Tp_kl))
          Ht=-sum(sum(t_kl*(t_kl>0)*log(t_kl*(t_kl>0)+(t_kl<=0)*eps)))
          t_l=matrix(colSums(t_kl),ncol=1)

          W2_old=W2
          W2=sum(t_kl*t_kl) + t(r_h)%*%log_rho_h  + Hs + Ht
          if (abs((W2-W2_old)/W2) < epsi_int) break
        }
        # Computation of parameters
        rho_h=(r_h+a-1)/(J+h*(a-1))

        tau_l=(t_l+a-1)/(K+h*(a-1))


        u_kl=t(r_jh)%*%x%*%t_kl
        n_kl=r_h%*%t(t_l)
        u_kl=u_kl*(u_kl>0)
        n_kl=n_kl*(n_kl>0)



        n_klx=t(r_jh)%*%x%*%t_kl

        if (any(any(n_kl<=0))){
          gamma_hl=(alpha-1+n_klx)/(beta+n_kl+((n_kl<=0)))*(n_kl>0)
        }else{
          gamma_hl=(alpha-1+n_klx)/(beta+n_kl)
        }
        if (any(any(gamma_hl<=0))){
          gamma_hl=gamma_hl*(gamma_hl>0)
          log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
        }else {
          log_gamma_hl=log(gamma_hl)
        }

        log_rho_h=log(rho_h*(rho_h>0)+(rho_h<=0)*eps)
        log_tau_l=log(tau_l*(tau_l>0)+(tau_l<=0)*eps)

        # Criterion
        W_old=W
        W=sum(r_h*log(r_h))- J*log(J) + sum(t_l*log(t_l)) - K*log(K) -
          sum(sum(n_kl*gamma_hl))+ sum(sum(n_klx*log_gamma_hl))+ Hs + Ht +
          (a-1)*(sum(log(rho_h))+sum(log(tau_l)))-beta*sum(sum(gamma_hl))
        if (alpha!=1){
          W=W+(alpha-1)*sum(sum(log_gamma_hl))
        }
        if (is.nan(W)){
          if (any(any(gamma_hl<0))){
            gamma_hl=gamma_hl*(gamma_hl>0)
            log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
          }
          W=sum(r_h*log(r_h+(r_h<=0)*eps))- J*log(J) + sum(t_l*log(t_l+(t_l<=0)*eps)) - K*log(K) -
            sum(sum(n_kl*gamma_hl))+ sum(sum(n_klx*log_gamma_hl))+ Hs + Ht +
            (a-1)*(sum(log(rho_h))+sum(log(tau_l)))-beta*sum(sum(gamma_hl))
          if (alpha!=1){
            W=W+(alpha-1)*sum(sum(log_gamma_hl))
          }
        }
        if (abs((W-W_old)/W) < epsi){ break}

        if (is.nan(W)) {break}
      }

      theta$gamma_hl=gamma_hl
      theta$rho_h=rho_h
      theta$tau_l=tau_l

      list(W=W,theta=theta,r_jh=r_jh,t_kl=t_kl,iter=iter)
    }



    BlocGibbs=function(x,a,normalization,alpha,beta,ini=list(n_init=10,ale=FALSE),niter=1000,eps=10^(-16),classement){
      if (is.null(ini$theta)){
        if (!is.numeric(ini$n_init)){
          stop("For a random initialization n_init must be numeric")
        }
        if (ini$n_init<=0){
          stop("n_init must be positive")
        }
        W=-Inf
        Res=list()
        #   for (i in 1:ini$n_init){
        #     init=PoissonBlocInit(h,l,data,normalisation,start=ini)
        #   }
        # }else{init=ini}
      }
      theta=ini$theta
      t_kl=ini$t_kl
      eps10=10*eps
      h=dim(theta$gamma_hl)[1]
      l=dim(theta$gamma_hl)[2]
      J=dim(x)[1]
      K=dim(x)[2]

      # loop of the iterations of the algorithm
      v=matrix(0,J,1)
      w=matrix(0,K,1)
      gamma_hl_path=list()
      rho_h_path=matrix(0,h,niter)
      tau_l_path=matrix(0,l,niter)
      Part=list()
      Part$v=matrix(0,nrow=J,ncol=h)
      Part$w=matrix(0,nrow=K,ncol=l)

      #SEM is random : 1000 runs from the same initial value

      gamma_hl=theta$gamma_hl
      rho_h=theta$rho_h
      tau_l=theta$tau_l
      t_l=matrix(colSums(t_kl),ncol=1)
      if (normalization){

      }

      for (iter in 1:niter){
        log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
        log_rho_h=log(rho_h)
        log_tau_l=log(tau_l)

        # Computation of v
        R_jh=matrix(1,J,1)%*%(t(log_rho_h)- t(t_l)%*%t(gamma_hl))+x%*%t_kl%*%t(log_gamma_hl)
        R_jhmax=apply(R_jh, 1,max)
        Rp_jh = exp(R_jh-R_jhmax)
        r_jh=Rp_jh/(rowSums(Rp_jh))
        for (iter_i in 1:J){
          v[iter_i]=min(h,h+1- sum( matrix(1,1,h)*runif(1) < cumsum(r_jh[iter_i,]) ))
        }
        r_jhn=!(v%*%matrix(1,1,h)-matrix(1,J,1)%*%(1:h))
        r_hn=matrix(colSums(r_jhn),ncol=1)
        r_jh=r_jhn
        r_h=r_hn

        # Computation of w
        T_kl=matrix(1,K,1)%*%(t(log_tau_l)-t(r_h)%*%gamma_hl)+t(x)%*%r_jh%*%log_gamma_hl
        T_klmax=apply(T_kl, 1,max)
        Tp_kl = exp(T_kl-T_klmax)
        t_kl=Tp_kl/(rowSums(Tp_kl))
        for (iter_j in 1:K){
          w[iter_j]=min(l,l+1- sum(matrix(1,1,l)*runif(1) < cumsum(t_kl[iter_j,])))
        }
        t_kln=!(w%*%matrix(1,1,l)-matrix(1,K,1)%*%(1:l))
        t_ln=matrix(colSums(t_kln),ncol=1)
        t_kl=t_kln
        t_l=t_ln


        # Computation of parameters
        nv=r_h+a
        rho_h=rgamma(h,nv,1)
        rho_h=rho_h/(sum(rho_h))
        dw=t_l+a
        tau_l=rgamma(l,dw,1)
        tau_l=tau_l/(sum(tau_l))

        n_kl=t(r_h)%*%t_l
        n_klx=t(r_jh)%*%x%*%t_kl
        gamma_hl=matrix(rgamma(h*l,shape=alpha+n_klx,scale=1/(beta+n_kl)),ncol=l)

        # trajectoires
        gamma_hl_path=c(gamma_hl_path,list(gamma_hl))
        rho_h_path[,iter]=rho_h
        tau_l_path[,iter]=tau_l
        Part$v=Part$v+r_jh
        Part$w=Part$w+t_kl
      }


      # arguments de sortie
      theta_path=list()
      theta_path$gamma_hl=gamma_hl_path
      theta_path$rho_h=rho_h_path
      theta_path$tau_l=tau_l_path

      list(theta_path=theta_path,Part=Part)

    }

    PoissonClassement=function(theta,Part=list()){
      ind_i=order(as.vector(theta$gamma_hl%*%theta$tau_l))
      ind_j=order(as.vector(t(theta$rho_h)%*%theta$gamma_hl))
      if (is.null(Part$v)||is.null(Part$w)){
        list(rho_h=theta$rho_h[ind_i],tau_l=theta$tau_l[ind_j],gamma_hl=theta$gamma_hl[ind_i,ind_j])
      }else{
        v=matrix(0,length(Part$v),1)
        maxv=max(Part$v)
        for (i in 1:maxv){
          v=v+ind_i[i]*(Part$v==i)
        }
        w=matrix(0,length(Part$w),1)
        maxw=max(Part$w)
        for (j in 1:maxw){
          w=w+ind_j[j]*(Part$w==j)
        }
        thetares=list(rho_h=theta$rho_h[ind_i],tau_l=theta$tau_l[ind_j],gamma_hl=theta$gamma_hl[ind_i,ind_j])
        list(theta=thetares,v=v,w=w,ind_i=ind_i,ind_j=ind_j)
      }
    }

    PoissonBlocVBayes=function(x,h,l,a,alpha,beta,normalization,ini=list(n_init=1,ale=TRUE),niter=10000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
      if (is.null(ini$theta)){
        if (!is.numeric(ini$n_init)){
          stop("For a random initialization, n_init must be numeric")
        }
        if (ini$n_init<=0){
          stop("n_init must be positive")
        }
      }
      W=-Inf
      Res=list()
      for (itry in 1:ini$n_init){
        initstart=PoissonBlocInit(h,l,x,normalization,start=ini)
        ResVBayesAle=BlocXemalpha1_kl_LBM(x,initstart$theta,initstart$t_kl,a=a,alpha=alpha,beta=beta,normalization=normalization,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)
        if (is.nan(ResVBayesAle$W)){

          stop('erreur avec le prior a')}

        if (ResVBayesAle$W > W){
          W=ResVBayesAle$W
          theta_max=ResVBayesAle$theta
          r_jh_max=ResVBayesAle$r_jh
          t_kl_max=ResVBayesAle$t_kl
        }
        if (is.nan(ResVBayesAle$W)&&(itry==ini$n_init)&&(W==-Inf)){
          W=ResVBayesAle$W
          theta_max=ResVBayesAle$theta
          r_jh_max=ResVBayesAle$r_jh
          t_kl_max=ResVBayesAle$t_kl
        }
      }

      theta_path=0
      theta=theta_max
      r_jh=r_jh_max
      t_kl=t_kl_max
      v=apply(r_jh, 1, which.max)
      w=apply(t_kl, 1, which.max)
      iter=ResVBayesAle$iter
      empty_cluster=(length(unique(v))!=h)||(length(unique(w))!=l)
      list(W=W,theta=theta,r_jh=r_jh,t_kl=t_kl,iter=iter,empty_cluster=empty_cluster,v=v,w=w)
    }
    PoissonBlocGibbs=function(x,h,l,a,normalization,alpha=1,beta=0.01,niter=5000,inistart=list(n_init=10,ale=TRUE),classement=TRUE,eps=10^(-16)){
      if (is.null(inistart$theta)){
        initgibbs=PoissonBlocInit(h,l,x,normalization,start=inistart)
        initgibbs$n_init=inistart$n_init
        initgibbs$ale=inistart$ale
      }else{
        initgibbs=inistart
      }
      res=BlocGibbs(x,a=a,normalization=normalization,alpha=alpha,beta=beta,ini=initgibbs,niter=niter,eps=eps,classement=classement)
      h=dim(res$theta_path$rho_h)[1]
      l=dim(res$theta_path$tau_l)[1]
      theta=list(rho_h=matrix(0,h,1),tau_l=matrix(0,l,1),gamma_hl=matrix(0,h,l))
      if (classement){
        for (i in 1:niter){
          thetabis=list(gamma_hl=res$theta_path$gamma_hl[[i]],rho_h=res$theta_path$rho_h[,i],tau_l=res$theta_path$tau_l[,i])
          thetabisorg=PoissonClassement(thetabis)
          theta$rho_h=theta$rho_h+thetabisorg$rho_h
          theta$tau_l=theta$tau_l+thetabisorg$tau_l
          theta$gamma_hl=theta$gamma_hl+thetabisorg$gamma_hl
        }
      }else{
        for (i in 1:niter){
          theta$gamma_hl=theta$gamma_hl+res$theta_path$gamma_hl[[i]]
        }
        theta$rho_h=rowSums(res$theta_path$rho_h)
        theta$tau_l=rowSums(res$theta_path$tau_l)
      }
      theta$rho_h=theta$rho_h/niter
      theta$tau_l=theta$tau_l/niter
      theta$gamma_hl=theta$gamma_hl/niter
      v=apply(res$Part$v, 1, which.max)
      w=apply(res$Part$w, 1, which.max)
      list(theta=theta,v=v,w=w,Part=res$Part,theta_path=res$theta_path)
    }

    PoissonBlocICL =function (a,alpha,beta,x,v1,w1,normalization){


      J=dim(x)[1]
      K=dim(x)[2]
      if (!is.matrix(v1)){
        h=max(v1)
        v=!(v1%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
        #v=!(v1%*%matrix(1,1,h)-matrix(1,J,1)%*%c(1:h))
      }else{
        v=v1
      }
      if (!is.matrix(w1)){
        l=max(w1)
        w=!(w1%*%t(rep(1,l))-rep(1,K)%*%t(1:l))
        #w=!(w1%*%matrix(1,1,l)-matrix(1,K,1)%*%c(1:l))
      }else{
        w=w1
      }

      if (normalization){

      }

      nh=matrix(colSums(v),ncol=h)
      nl=matrix(colSums(w),ncol=l)
      nhl=t(nh)%*%nl
      nhlx=t(v)%*%x%*%w
      critere=lgamma(a*h)+lgamma(a*l)-(l+h)*lgamma(a)-h*l*lgamma(alpha)-lgamma(J+h*a)-lgamma(K+l*a)+
        sum(lgamma(nh+a))+sum(lgamma(nl+a))+sum(sum(-(alpha+nhlx)*log(beta+nhl)+lgamma(alpha+nhlx)))

      if (beta>0){
        critere=critere+h*l*alpha*log(beta)
      }
      critere}

    PoissonBlocBIC =function (a,alpha,beta,x,res,normalization){
      v1=res$v
      w1=res$w
      theta=res$theta

      J=dim(x)[1]
      K=dim(x)[2]
      if (!is.matrix(v1)){
        h=max(v1)
        v=!(v1%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
      }else{
        v=v1
      }
      if (!is.matrix(w1)){
        l=max(w1)
        w=!(w1%*%t(rep(1,l))-rep(1,K)%*%t(1:l))
      }else{
        w=w1
      }

      if (normalization){

      }

      nh=matrix(colSums(v),ncol=h)
      nl=matrix(colSums(w),ncol=l)
      nhl=t(nh)%*%nl
      nhlx=t(v)%*%x%*%w
      if (any(any(theta$gamma_hl<=0))){
        gamma_hl=theta$gamma_hl*(theta$gamma_hl>0)
        log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
      }else{
        gamma_hl=theta$gamma_hl
        log_gamma_hl=log(theta$gamma_hl)
      }
      if ((a==1)){
        W=res$W-((alpha-1)*sum(sum(log_gamma_hl))-beta*sum(sum(gamma_hl)))}
      else{
        W=res$W-((a-1)*(sum(log(theta$rho_h))+sum(log(theta$tau_l)))+(alpha-1)*sum(sum(log_gamma_hl))-beta*sum(sum(gamma_hl)))
        if (is.nan(W)){
          W=res$W-((a-1)*(sum(log(theta$rho_h+(theta$rho_h<=0))*(theta$rho_h>0))+sum(log(theta$tau_l+(theta$tau_l<=0))*(theta$tau_l>0)))+(alpha-1)*sum(sum(log_gamma_hl))-beta*sum(sum(gamma_hl)))}
      }
      critere=W-((h-1)*log(J)+(l-1)*log(K)+l*h*log(J*K))/2
      critere}










    ### Verif
    if (!is.numeric(mc.cores)){
      stop("mc.cores must be an integer greater than 1")
    } else if ((mc.cores<=0)||(length(mc.cores)!=1)||(floor(mc.cores)!=mc.cores)){
      stop("mc.cores must be an integer greater than 1")
    }

    ## =============================================================
    ## INITIALIZATION & PARAMETERS RECOVERY
    if ((Sys.info()[['sysname']] == "Windows")&(mc.cores>1)) {
      warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
      mc.cores <- 1
    }


    #source("Chargement.R")

    criterion_tab=matrix(-Inf,Hmax+1,Lmax+1)
    W_tab=matrix(-Inf,Hmax+1,Lmax+1)
    W=-Inf
    if (init_choice=='smallVBayes'){
      # if (mc.cores>1){
      #   Temp=parallel::mclapply(1:ntry,function(i){PoissonBlocVBayes(data,Hstart,Lstart,a,alpha,beta,
      #                                                    normalization,ini=list(n_init=1,ale=TRUE),
      #                                                    niter=50)}, mc.cores=mc.cores)
      # }else{
      #   Temp=lapply(1:ntry,function(i){PoissonBlocVBayes(data,Hstart,Lstart,a,alpha,beta,
      #                                                    normalization,ini=list(n_init=1,ale=TRUE),
      #                                                    niter=50)})
      # }

      for (i in 1:ntry){
        #res=BlocGibbs(data$x,a,b,ini=list(n_init=1,ale=TRUE),1000)
        res=PoissonBlocVBayes(x,Hstart,Lstart,a,alpha,beta,
                              normalization,ini=list(n_init=1,ale=TRUE),
                              niter=50)
        #print(W)
        if (res$W>W && res$empty_cluster==FALSE){
          W=res$W
          resopt=res
        }
        if (!exists("resopt")){
          stop(paste('The algorithm does not manage to find a configuration with  (Hstart=',as.character(Hstart),',Lstart=',as.character(Lstart),'), try normalization=FALSE or a higher ntry') )
        }
      }

    }else if (init_choice=='user'){
      if (is.null(userparam)){
        stop(paste('If you choose the user init choice, please fill the useparam field with partitions v and w' ))
      }else{
        initVBayes=list(v=userparam$v,w=userparam$w,ale=FALSE,n_init=1)
        resopt=PoissonBlocVBayes(x,Hstart,Lstart,a,alpha,beta,normalization,ini=initVBayes)

      }
    }else if (init_choice=='Gibbs'){

      resinit=PoissonBlocGibbs(x,Hstart,Lstart,a,normalization,alpha,beta,niter=5000,inistart=list(n_init=10,ale=TRUE))
      initVBayes=list(v=resinit$v,w=resinit$w,ale=FALSE,n_init=1)
      if (max(initVBayes$v)!=Hstart|| (max(initVBayes$w)!=Lstart)){
        stop(paste('The algorithm does not manage to find a configuration with  (Hstart=',as.character(Hstart),',Lstart=',as.character(Lstart),'),try another init_choice') )
      }
      resopt=PoissonBlocVBayes(x,Hstart,Lstart,a,alpha,beta,normalization,ini=initVBayes)
      if (resopt$empty_cluster==TRUE){
        stop(paste('The algorithm does not manage to find a configuration with  (Hstart=',as.character(Hstart),',Lstart=',as.character(Lstart),'),try another init_choice') )
      }

    }else if (init_choice=='random'){

      eps=10^(-16)
      n=dim(x)[1]
      K=dim(x)[2]


      for (i in 1:ntry){
        theta=list()
        theta$rho_h=runif(Hstart)
        theta$rho_h=theta$rho_h/sum(theta$rho_h)
        theta$tau_l=runif(Lstart)
        theta$tau_l=theta$tau_l/sum(theta$tau_l)
        resw=PartRnd(K,theta$tau_l)
        t_kl=resw$v_jh
        w=resw$v
        resv=PartRnd(J,theta$rho_h)
        r_jh=resv$v_jh
        v=resv$v
        if (normalization){

        }
        r_h=matrix(colSums(r_jh),ncol=1)
        t_l=matrix(colSums(t_kl),ncol=1)
        n_kl=r_h%*%t(t_l)
        theta$gamma_hl=(t(r_jh)%*%x%*%t_kl)/(n_kl)
        if (any(is.nan(theta$gamma_hl))){
          lambda0=mean(mean(x))
          theta$gamma_hl=theta$gamma_hl=matrix(rpois(h,lambda0),ncol=l)
        }

        initVBayes=list(t_kl=t_kl,w=w,theta=theta,ale=FALSE,n_init=1)
        #res=BlocGibbs(data$x,a,b,ini=list(n_init=1,ale=TRUE),1000)
        res=PoissonBlocVBayes(x,Hstart,Lstart,a,alpha,beta,normalization,ini=initVBayes)
        if (res$W>W && res$empty_cluster==FALSE){
          resopt=res}



      }
      if (!exists("resopt")){
        stop(paste('The algorithm does not manage to find a configuration with  (Hstart=',as.character(Hstart),',Lstart=',as.character(Lstart),'), try normalization=FALSE') )
      }
    } else {stop("This is not a initialization choice proposed by the procedure")}

    if (criterion_choice=='ICL'){
      Critere= PoissonBlocICL(a,alpha,beta,x,resopt$v,resopt$w,normalization)  }
    else if (criterion_choice=='BIC'){
      Critere= PoissonBlocBIC(a,alpha,beta,x,resopt,normalization)
    }
    else {stop('This is not a criterion choice proposed by the procedure')}
    criterion_tab[Hstart,Lstart]=Critere
    W_tab[Hstart,Lstart]=resopt$W

    modele=resopt
    criterion_max=Critere
    h=Hstart
    l=Lstart
    hopt=Hstart
    lopt=Lstart
    model_max=resopt

    while (h+1<=Hmax && l+1<=Lmax){
      if ((length(unique(modele$v))!=h)||(length(unique(modele$w))!=l)){
        stop(paste('The algorithm does not manage to find a configuration with (h=',as.character(h),',l=',as.character(l),'), maybe reduce Hmax and Lmax') )
      }
      W_colonne=-Inf
      W_ligne=-Inf

      f1=function(Kc){
        w=modele$w
        Xc=which(w==Kc)
        lxc=length(Xc)
        if (lxc>1){
          idxc=sample(lxc,floor(lxc/2),replace=FALSE)
          indc=Xc[idxc]
          w[indc]=l+1
          initVBayes=list(v=modele$v,w=w,ale=FALSE,n_init=1)
          resultsVBayes=PoissonBlocVBayes(x,h,l+1,a,alpha,beta,normalization,ini=initVBayes)
          iterverif=0
          while (resultsVBayes$empty_cluster==TRUE && iterverif<20){
            iterverif=iterverif+1
            w=modele$w
            idxc=sample(lxc,floor(lxc/2),replace=FALSE)
            indc=Xc[idxc]
            w[indc]=l+1
            initVBayes=list(v=modele$v,w=w,ale=FALSE,n_init=1)
            resultsVBayes=PoissonBlocVBayes(x,h,l+1,a,alpha,beta,normalization,ini=initVBayes)
          }
          if (iterverif<20){
            return(list(W=resultsVBayes$W,modele_colonne=resultsVBayes))
          }else{
            return(list(W=-Inf))
          }
        }else{
          return(list(W=-Inf))
        }
      }
      if (mc.cores>1){
        l_colonne=parallel::mclapply(1:l,f1, mc.cores=mc.cores)
      }else{
        l_colonne=lapply(1:l,f1)
      }

      for (Kc in 1:l){
        if (l_colonne[[Kc]]$W>W_colonne){
          W_colonne=l_colonne[[Kc]]$W
          modele_colonne=l_colonne[[Kc]]$modele_colonne
        }
      }
      if (W_colonne<=-Inf){
        Empty_m=TRUE
      }else{
        Empty_m=FALSE
      }
      #### Class h
      f2=function(Kl){
        v=modele$v
        Xl=which(v==Kl)
        ll=length(Xl)
        if (ll>1){
          idxl=sample(ll,floor(ll/2),replace=FALSE)
          indl=Xl[idxl]
          v[indl]=h+1
          initVBayes=list(v=v,w=modele$w,ale=FALSE,n_init=1)
          resultsVBayes=PoissonBlocVBayes(x,h+1,l,a,alpha,beta,normalization,ini=initVBayes)
          iterverif=0
          while (resultsVBayes$empty_cluster==TRUE && iterverif<20){
            iterverif=iterverif+1
            v=modele$v
            idxl=sample(ll,floor(ll/2),replace=FALSE)
            indl=Xl[idxl]
            v[indl]=h+1
            initVBayes=list(v=v,w=modele$w,ale=FALSE,n_init=1)
            resultsVBayes=PoissonBlocVBayes(x,h+1,l,a,alpha,beta,normalization,ini=initVBayes)
          }
          if (iterverif<20){
            return(list(W=resultsVBayes$W,modele_ligne=resultsVBayes))
          }else{
            return(list(W=-Inf))
          }
        }else{
          return(list(W=-Inf))
        }


      }

      if (mc.cores>1){
        l_ligne=parallel::mclapply(1:h,f2, mc.cores=mc.cores)
      }else{
        l_ligne=lapply(1:h,f2)
      }
      for (Kl in 1:h){
        if (l_ligne[[Kl]]$W>W_ligne){
          W_ligne=l_ligne[[Kl]]$W
          modele_ligne=l_ligne[[Kl]]$modele_ligne
        }
      }

      if ((W_ligne<=-Inf)&(Empty_m)){
        warning("The algorithm does not manage to find a configuration  (",as.character(h+1),',',as.character(l), ") nor (",as.character(h),',',as.character(l+1), ")")
        new(Class="BIKM1_LBM_Poisson",init_choice=init_choice,criterion_choice=criterion_choice,criterion_tab=criterion_tab,criterion_max=criterion_max,model_max=model_max,W_tab=W_tab,hopt=hopt,lopt=lopt)
      }


      if (W_ligne>W_colonne){
        modele=modele_ligne
        if (criterion_choice=='ICL'){
          Critere= PoissonBlocICL(a,alpha,beta,x,modele$v,modele$w,normalization)
        }else if (criterion_choice=='BIC'){
          Critere= PoissonBlocBIC(a,alpha,beta,x,modele,normalization)
        }
        criterion_tab[h+1,l]=Critere
        W_tab[h+1,l]=modele$W
        if (Critere>criterion_max){
          model_max=modele
          criterion_max=Critere
          hopt=h+1
          lopt=l
        }
        if (verbose){
          cat('(h,l)=(',as.character(h+1),',',as.character(l),')\n',sep = "")
        }
        h=h+1
      }else {
        modele=modele_colonne
        if (criterion_choice=='ICL'){
          Critere= PoissonBlocICL(a,alpha,beta,x,modele$v,modele$w,normalization)
        }else if (criterion_choice=='BIC'){
          Critere= PoissonBlocBIC(a,alpha,beta,x,modele,normalization)
        }
        criterion_tab[h,l+1]=Critere
        W_tab[h,l+1]=modele$W
        if (Critere>criterion_max){
          model_max=modele
          criterion_max=Critere
          hopt=h
          lopt=l+1

        }
        if (verbose){
         cat('(h,l)=(',as.character(h),',',as.character(l+1),')\n',sep = "")
        }

        l=l+1
      }
    }
    if (verbose){
      cat('The selected row (h) and column (l) clusters are (h,l)=(',as.character(hopt),',',as.character(lopt),')\n',sep = "")
    }

    new(Class="BIKM1_LBM_Poisson",init_choice=init_choice,criterion_choice=criterion_choice,criterion_tab=criterion_tab,criterion_max=criterion_max,model_max=model_max,W_tab=W_tab,hopt=hopt,lopt=lopt)






}

   }

##' PoissonBlocRnd function for contingency data simulation
##'
##' Produce a simulated data matrix generated under the Poisson Latent Block Model.
##'
##' @param J a positive integer specifying the number of expected rows.
##' @param K a positive integer specifying the number of expected columns.
##' @param theta a list specifying the model parameters:
##'
##' \code{rho_h}: a vector specifying the row mixing proportions.
##'
##' \code{tau_l}: a vector specifying the column mixing proportions.
##'
##' \code{gamma_hl}: a matrix specifying the distribution parameter.
##' @usage PoissonBlocRnd(J,K,theta)
##' @return a list including the arguments:
##'
##' \code{x}: simulated contingency data matrix.
##'
##' \code{xrow}: numeric vector specifying row partition.
##'
##' \code{xcol}: numeric vector specifying column partition.

##' @rdname PoissonBlocRnd-proc
##'
##' @examples
##' require(bikm1)
##' J=200
##' K=120
##' h=3
##' l=2
##' theta=list()
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##' theta$gamma_hl=matrix(c(1, 6,4, 1, 7, 1),ncol=2)
##' data=PoissonBlocRnd(J,K,theta)
##'
##' @export PoissonBlocRnd


PoissonBlocRnd =function (J,K,theta){
  PartRnd = function(J,proba){

    h=length(proba)

    i=sample(1:J,J, replace = FALSE)
    i1=i[1:h]
    i2=i[(h+1):J]
    v=rep(0,J)
    v[i1]=1:h

    if (J > h) {
      v[i2]=(h+1)-colSums(rep(1,h)%*%t(runif(J-h)) < cumsum(proba)%*%t(rep(1,J-h)))
      v[v==(h+1)]=h
    }

    v_jh=!(v%*%t(rep(1,h))-rep(1,J)%*%t(1:h))

    list(v=v,v_jh=v_jh)
  }
  resv=PartRnd(J,theta$rho_h)
  v=resv$v
  v_jh=resv$v_jh
  resw=PartRnd(K,theta$tau_l)
  w=resw$v
  w_jl=resw$v_jh
  h=dim(theta$gamma_hl)[1]
  l=dim(theta$gamma_hl)[2]
  x=matrix(rpois(J*K,v_jh%*%theta$gamma_hl%*%t(w_jl)),ncol=K)
  xrow=as.numeric(v_jh%*%c(1:h))
  xcol=as.numeric(w_jl%*%c(1:l))
  list(x=x,xrow=xrow,xcol=xcol)


}

##' PoissonBlocVisu function for visualization of contingency datasets
##'
##' Produce a plot object representing the coclustered data-sets.
##'
##' @param x contingency matrix of observations.
##' @param v a numeric vector specifying the class of rows.
##' @param w a numeric vector specifying the class of columns.
##' @usage PoissonBlocVisu(x,v,w)
##' @rdname PoissonBlocVisu-proc
##' @return a \pkg{plot} object
##' @examples
##'
##' require(bikm1)
##' J=200
##' K=120
##' h=3
##' l=2
##' theta=list()
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##' theta$gamma_hl=matrix(c(1, 6,4, 1, 7, 1),ncol=2)
##' data=PoissonBlocRnd(J,K,theta)
##' PoissonBlocVisu(data$x, data$xrow,data$xcol)
##' @export PoissonBlocVisu

PoissonBlocVisu=function(x,v,w){
  ## Initialization

  h=max(v)
  l=max(w)
  J=dim(x)[1]
  K=dim(x)[2]
  ii=c()
  jj=c()
  for (i in 1:h){
    ii=c(ii,which(v==i))
  }
  for (j in 1:l){
    jj=c(jj,which(w==j))
  }
  v_jh=((v%*%t(rep(1,h)))==(rep(1,J)%*%t(1:h)))
  w_jl=((w%*%t(rep(1,l)))==(rep(1,K)%*%t(1:l)))

  ## Affichage
  n_k=J-sort(cumsum(colSums(v_jh))-0.5,decreasing=TRUE)
  n_l=cumsum(colSums(w_jl))+0.5

  table.paint(x[ii,jj],clegend=0)

  for (i in n_k[2:h]){
    lines(c(0.5,(K+1)),c(i,i),col='blue',lwd=2)
  }
  for (i in n_l[1:(l-1)]){
    lines(c(i,i),c(0.5,(J+1)),col='blue',lwd=2)
  }

}

##' PoissonBlocVisuResum function  for visualization of contingency datasets
##'
##' Produce a plot object representing the resumed co-clustered data-sets.
##'
##' @param x contingency matrix of observations.
##' @param v a numeric vector specifying the class of rows.
##' @param w a numeric vector specifying the class of columns.
##' @rdname PoissonBlocVisuResum-proc
##' @usage PoissonBlocVisuResum(x,v,w)
##' @return a \pkg{plot} object.
##' @examples
##' require(bikm1)
##' J=200
##' K=120
##' h=3
##' l=2
##' theta=list()
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##' theta$gamma_hl=matrix(c(1, 6,4, 1, 7, 1),ncol=2)
##' data=PoissonBlocRnd(J,K,theta)
##' PoissonBlocVisuResum(data$x, data$xrow,data$xcol)
##' @export PoissonBlocVisuResum


PoissonBlocVisuResum=function(x,v,w){
  ## Initialisation

  h=max(v)
  l=max(w)
  J=dim(x)[1]
  K=dim(x)[2]
  v_jh=((v%*%t(rep(1,h)))==(rep(1,J)%*%t(1:h)))
  w_jl=((w%*%t(rep(1,l)))==(rep(1,K)%*%t(1:l)))

  ## Affichage
  nh=colSums(v_jh)
  nl=colSums(w_jl)
  nhl=nh%*%t(nl)
  nhl=t(v_jh)%*%x%*%w_jl

  table.paint(nhl/nhl)

}








##' ARI function for agreement between two partitions
##'
##' Produce a measure of agreement between two partitions.
##'
##' @param v numeric vector  specifying the class of observations.
##' @param vprime numeric vector specifying another partitions of observations.
##' @return a list including the arguments:
##'
##'   \code{ari}: value of the index.
##'
##'   \code{nv}: contingency table which the index is based on.
##'
##' @rdname ARI-proc
##' @usage ARI(v,vprime)
##' @references Hubert and Arabie. Comparing partitions. Journal of
##'   classification (1985).
##' @examples
##' \donttest{require(bikm1)
##' J=200
##' K=120
##' h=3
##' l=2
##' theta=list()
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##' theta$gamma_hl=matrix(c(1, 6,4, 1, 7, 1),ncol=2)
##' data=PoissonBlocRnd(J,K,theta)
##' res=BIKM1_LBM_Poisson(data$x,4,4,4,init_choice='random')
##' mv=ARI(res@model_max$v, data$xrow)
##' mv$ari
##' mv$nv
##' mw=ARI(res@model_max$w, data$xcol)}
##'
##' @export ARI


ARI=function(v,vprime){
  L_v=length(v)
  if (L_v!=length(vprime)){
    warning('Both partitions must contain the same number of points.')
  }
  N_v=max(v)
  N_vprime=max(vprime)
  nv=matrix(0,N_v,N_vprime)
  for (i in 1:N_v){
    for (j in 1:N_vprime){
      G1 = which(v==i)
      G2 = which(vprime==j)
      nv[i,j] = length(intersect(G1,G2))
    }
  }

  ssm = 0
  sm1 = 0
  sm2 = 0
  for (i in 1:N_v){
    for (j in 1:N_vprime){
      ssm = ssm + choose(nv[i,j],2)
    }
  }
  temp = rowSums(nv)
  for (i in 1:N_v){
    sm1 = sm1 + choose(temp[i],2)
  }
  temp = colSums(nv)
  for (i in 1:N_vprime){
    sm2 = sm2 + choose(temp[i],2)
  }
  NN = ssm - (sm1*sm2)/choose(L_v,2)
  DD = (sm1 + sm2)/2 - (sm1*sm2)/choose(L_v,2)
  ari = NN/DD
  return(list(ari=ari, nv=nv))

}








##' CARI function for agreement between coclustering partitions
##'
##' Produce a measure of agreement between two pairs of partitions for coclustering. A value of 1 means a perfect match.
##'
##' @param v numeric vector  specifying the class of  rows.
##' @param w numeric vector specifying the class of columns.
##' @param vprime numeric vector  specifying another partitions of rows.
##' @param wprime numeric vector specifying another partition of columns.
##' @return a list including the arguments:
##'
##' \code{cari}: value of the index.
##'
##' \code{nvw}: contingency table which the index is based on.
##'
##' @usage CARI(v,w,vprime,wprime)
##' @rdname CARI-proc
##' @references Robert and Vasseur. Comparing high dimensional partitions with the Coclustering Adjusted Rand Index, Preprint (2017).
##' @examples
##' \donttest{require(bikm1)
##' J=200
##' K=120
##' h=3
##' l=2
##' theta=list()
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##' theta$gamma_hl=matrix(c(1, 6,4, 1, 7, 1),ncol=2)
##' data=PoissonBlocRnd(J,K,theta)
##' res=BIKM1_LBM_Poisson(data$x,4,4,4,init_choice='smallVBayes')
##' me=CARI(res@model_max$v,res@model_max$w, data$xrow,data$xcol)
##' me$cari
##' me$nvw}
##' @export CARI

CARI=function(v,w,vprime,wprime){
  L_v=length(v)
  L_w=length(w)
  if (L_v!=length(vprime)){
    warning('Both partitions must contain the same number of points.')
  }
  if (L_w!=length(wprime)){
    warning('Both partitions must contain the same number of points.')
  }

  N_v=max(v)
  N_w=max(w)
  N_vprime=max(vprime)
  N_wprime=max(wprime)
  ari_v=ARI(v,vprime)
  ari_w=ARI(w,wprime)
  nv=ari_v[[2]]
  nw=ari_w[[2]]
  nvw=kronecker(nv,nw)
  ssm = 0
  sm1 = 0
  sm2 = 0
  for (i in 1:(N_v*N_w)){
    for (j in 1:(N_vprime*N_wprime)){
      ssm = ssm + choose(nvw[i,j],2)
    }
  }
  temp1=rowSums(nv)
  temp2=rowSums(nw)
  temp=kronecker(temp1,temp2)
  for (i in 1:(N_v*N_w)){
    sm1 = sm1 + choose(temp[i],2)
  }
  temp1=colSums(nv)
  temp2=colSums(nw)
  temp = kronecker(temp1,temp2)
  for (i in 1:(N_vprime*N_wprime)){
    sm2 = sm2 + choose(temp[i],2)
  }
  NN = ssm - sm1*sm2/choose(L_v*L_w,2)
  DD = (sm1 + sm2)/2 - sm1*sm2/choose(L_v*L_w,2)
  cari = NN/DD
  return(list(cari=cari, nvw=nvw))
}


##' PoissonBlocBIC function for the computation of the BIC criterion in the Poisson LBM
##'
##' Produce a value of the BIC criterion for coclustering partitions
##' @param a hyperparameter used in the VBayes algorithm for priors on the mixing proportions. By default, a=4.
##' @param alpha hyperparameter used in the VBayes algorithm for prior on the Poisson parameter. By default, alpha=1.
##' @param beta hyperparameter used in the VBayes algorithm for prior on the Poisson parameter. By default, beta=0.01.
##' @param v1 a numeric vector of row partitions
##' @param w1 a numeric vector of column partitions
##' @param x  contingency matrix of observations.
##' @param res a BIKM1_LBM_Poisson object
##' \code{rho_h} mixing row proportions
##' \code{tau_l} mixing column proportions
##' \code{gamma_hl} Bernoulli parameters
##'
##' @param normalization logical. To use the normalized Poisson modelling in the Latent Block Model. By default normalization=FALSE.
##' @return a value of the BIC criterion
##'
##' @usage PoissonBlocBIC(a,alpha,beta,v1,w1,x,res,normalization)
##' @rdname PoissonBlocBIC-proc
##' @examples
##' \donttest{require(bikm1)
##' J=200
##' K=120
##' h=3
##' l=2
##' theta=list()
##' theta$rho_h=1/h*matrix(1,h,1)
##' theta$tau_l=1/l*matrix(1,l,1)
##' theta$gamma_hl=matrix(c(1, 6,4, 1, 7, 1),ncol=2)
##' data=PoissonBlocRnd(J,K,theta)
##' res=BIKM1_LBM_Poisson(data$x,3,3,4,init_choice='smallVBayes')
##' bic=PoissonBlocBIC(res@model_max$v,res@model_max$w,data$x,res,normalization=TRUE)}
##' @export PoissonBlocBIC


PoissonBlocBIC =function (a=4,alpha=1,beta=0.01,v1,w1,x,res,normalization=TRUE){
  # v1=res@model_max$v
  # w1=res@model_max$w
   theta=res@model_max$theta
  #data=list()
  #x=data$x
  J=dim(x)[1]
  K=dim(x)[2]
  if (!is.matrix(v1)){
    h=max(v1)
    v=!(v1%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
  }else{
    v=v1
  }
  if (!is.matrix(w1)){
    l=max(w1)
    w=!(w1%*%t(rep(1,l))-rep(1,K)%*%t(1:l))
  }else{
    w=w1
  }

  if (normalization){

  }

  nh=matrix(colSums(v),ncol=h)
  nl=matrix(colSums(w),ncol=l)
  nhl=t(nh)%*%nl
  nhlx=t(v)%*%x%*%w
  if (any(any(theta$gamma_hl<=0))){
    gamma_hl=theta$gamma_hl*(theta$gamma_hl>0)
    log_gamma_hl=log(gamma_hl+(gamma_hl<=0))*(gamma_hl>0)
  }else{
    gamma_hl=theta$gamma_hl
    log_gamma_hl=log(theta$gamma_hl)
  }
  if ((a==1)){
    W=res@model_max$W-((alpha-1)*sum(sum(log_gamma_hl))-beta*sum(sum(gamma_hl)))}
  else{
    W=res@model_max$W-((a-1)*(sum(log(theta$rho_h))+sum(log(theta$tau_l)))+(alpha-1)*sum(sum(log_gamma_hl))-beta*sum(sum(gamma_hl)))
    if (is.nan(W)){
      W=res@model_max$W-((a-1)*(sum(log(theta$rho_h+(theta$rho_h<=0))*(theta$rho_h>0))+sum(log(theta$tau_l+(theta$tau_l<=0))*(theta$tau_l>0)))+(alpha-1)*sum(sum(log_gamma_hl))-beta*sum(sum(gamma_hl)))}
  }
  critere=W-((h-1)*log(J)+(l-1)*log(K)+l*h*log(J*K))/2
  critere}


##' PoissonBlocICL function for the computation of the ICL criterion in the Poisson LBM
##'
##' Produce a value of the ICL criterion for coclustering partitions
##' @param a hyperparameter used in the VBayes algorithm for priors on the mixing proportions. By default, a=4.
##' @param alpha hyperparameter used in the VBayes algorithm for prior on the Poisson parameter. By default, alpha=1.
##' @param beta hyperparameter used in the VBayes algorithm for prior on the Poisson parameter. By default, beta=0.01.
##' @param x  contingency matrix of observations.
##' @param v1 a numeric vector specifying the class of rows.
##' @param w1 a numeric vector specifying the class of columns.
##' @param normalization logical. To use the normalized Poisson modelling in the Latent Block Model. By default normalization=FALSE.
##' @return a value of the ICL criterion
##' @usage PoissonBlocICL(a,alpha,beta,x,v1,w1,normalization)
##' @rdname PoissonBlocICL-proc
##' @examples
##' \donttest{require(bikm1)
##' J=200
##' K=120
##' h=3
##' l=2
##' theta=list()
##' theta$rho_h=(1/h)*matrix(1,h,1)
##' theta$tau_l=(1/l)*matrix(1,l,1)
##' theta$gamma_hl=matrix(c(1, 6,4, 1, 7, 1),ncol=2)
##' data=PoissonBlocRnd(J,K,theta)
##' res=BIKM1_LBM_Poisson(data$x,4,4,4,init_choice='smallVBayes')
##' icl=PoissonBlocICL(4,1,0.01,data$x,res@model_max$v,res@model_max$w, normalization=FALSE)}
##' @export PoissonBlocICL


PoissonBlocICL =function (a=4,alpha=1,beta=0.01,x,v1,w1,normalization=TRUE){

  #x=data$x
  J=dim(x)[1]
  K=dim(x)[2]
  if (!is.matrix(v1)){
    h=max(v1)
    v=!(v1%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
    #v=!(v1%*%matrix(1,1,h)-matrix(1,J,1)%*%c(1:h))
  }else{
    v=v1
  }
  if (!is.matrix(w1)){
    l=max(w1)
    w=!(w1%*%t(rep(1,l))-rep(1,K)%*%t(1:l))
    #w=!(w1%*%matrix(1,1,l)-matrix(1,K,1)%*%c(1:l))
  }else{
    w=w1
  }

  if (normalization){

  }

  nh=matrix(colSums(v),ncol=h)
  nl=matrix(colSums(w),ncol=l)
  nhl=t(nh)%*%nl
  nhlx=t(v)%*%x%*%w
  critere=lgamma(a*h)+lgamma(a*l)-(l+h)*lgamma(a)-h*l*lgamma(alpha)-lgamma(J+h*a)-lgamma(K+l*a)+
    sum(lgamma(nh+a))+sum(lgamma(nl+a))+sum(sum(-(alpha+nhlx)*log(beta+nhl)+lgamma(alpha+nhlx)))

  if (beta>0){
    critere=critere+h*l*alpha*log(beta)
  }
  critere}


