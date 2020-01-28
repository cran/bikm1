
##' BIKM1_MLBM_Binary fitting procedure
##'
##' Produce a blockwise estimation of  double matrices of observations.
##'
##' @param x  matrix of observations (1rst matrix).
##' @param y matrix of observations (2nd matrix).
##' @param Gmax a positive integer less than number of rows.
##' @param Hmax a positive integer less than number of columns of the 1st matrix.
##' @param Lmax a positive integer less than number of columns of the 2nd matrix. The bikm1 procedure stops while the numbers of rows is higher than Gmax or the number of columns is higher than Hmax or the numbers of columns(2nd matrix) is higher than Lmax.
##' @param a hyperparameter used in the VBayes algorithm for priors on the mixing proportions.
##' By default, a=4.
##' @param b hyperparameter used in the VBayes algorithm for prior on the Bernoulli parameter.
##' By default, b=1.
##' @param Gstart a positive integer to initialize the procedure with number of row clusters.
##' By default, Gstart=2.
##' @param Hstart a positive integer to initialize the procedure with number of column clusters.
##' By default, Hstart=2.
##' @param Lstart a positive integer to initialize the procedure with number of column clusters.
##' By default, Lstart=2.
##' @param init_choice character string corresponding to the chosen initialization strategy used for the procedure, which can be "random" or "smallVBayes" or "user".
##' By default, init_choice="smallVBayes".
##' @param userparam In the case where init_choice is "user", a list containing partitions z,v and w.
##' @param ntry a positive integer corresponding to the number of times which is launched the small VBayes initialization strategy. By default ntry=100.
##' @param criterion_choice Character string corresponding to the chosen criterion used for model selection, which can be "ICL" as for now.
##' By default, criterion_choice="ICL".
##' @param mc.cores a positive integer corresponding to the available number of cores for parallel computing.
##' By default, mc.cores=1.
##' @param verbose logical. To display each step and the result. By default verbose=TRUE.
##' @rdname BIKM1_MLBM_Binary-proc
##' @references Govaert and Nadif. Co-clustering, Wyley (2013).
##'
##' Keribin, Brault and Celeux. Estimation and Selection for the Latent Block Model on Categorical Data, Statistics and Computing (2014).
##'
##' Robert. Classification crois\'ee pour l'analyse de bases de donn\'ees de grandes dimensions de pharmacovigilance. Paris Saclay (2017).
##' @usage BIKM1_MLBM_Binary(x,y,Gmax,Hmax,Lmax,a=4,b=1,
##' Gstart=2,Hstart=2,Lstart=2,init_choice='smallVBayes',userparam=NULL,
##' ntry=50,criterion_choice='ICL', mc.cores=1,verbose=TRUE)
##'
##' @return a BIKM1_MLBM_Binary object including
##'
##' \code{model_max}: the selected model by the procedure including free energy W, theta, conditional probabilities (s_ig, r_jh,t_kl), iter, empty_cluster, and the selected partitions z,v and w.
##'
##' \code{criterion_choice}: the chosen criterion
##'
##' \code{init_choice}: the chosen init_choice
##'
##' \code{criterion_tab}: matrix containing the criterion values for each selected number of row and column
##'
##' \code{W_tab}: matrix containing the free energy values for each selected number of row and column
##'
##' \code{criterion_max}: maximum of the criterion values
##'
##' \code{gopt}: the selected number of rows
##'
##' \code{hopt}: the selected number of columns (1rst matrix)
##'
##' \code{lopt}: the selected number of columns (2nd matrix)
##'
##' @examples
##'
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' K=120
##' g=3
##' h=2
##' l=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##' theta$alpha_gh=matrix(runif(6),ncol=h)
##' theta$beta_gl=matrix(runif(6),ncol=l)
##' data=BinBlocRnd_MLBM(n,J,K,theta)
##' res=BIKM1_MLBM_Binary(data$x,data$y,3,2,2,Gstart=3,Hstart=2,Lstart=2,init_choice='user',
##' userparam=list(z=data$xrow,v=data$xcolx,w=data$xcoly))
##'
##'
##' @export BIKM1_MLBM_Binary
##'


BIKM1_MLBM_Binary=function(x,y,Gmax,Hmax,Lmax,a=4,b=1,Gstart=2,Hstart=2,Lstart=2,init_choice='smallVBayes',userparam=NULL,ntry=50,criterion_choice='ICL', mc.cores=1,verbose=TRUE){

  BinBlocInit_MLBM=function(g,h,l,x,y,start,fixed_proportion=FALSE){

    eps=10^(-16)
    n=dim(x)[1]
    J=dim(x)[2]
    K=dim(y)[2]
    # Initialisation par tirage aleatoire des centres


    if(start$ale==TRUE){
      theta=list()
      #theta$pi_g=runif(g)
      #theta$pi_g=theta$pi_g/sum(theta$pi_g)
      #theta$rho_h=runif(m)
      #theta$rho_h=theta$rho_h/sum(theta$rho_h)
      #theta$tau_l=runif(s)
      #theta$tau_l=theta$tau_l/sum(theta$tau_l)

      theta$pi_g=1/g *rep(1,g)
      theta$rho_h=1/h *rep(1,h)
      theta$tau_l=1/l *rep(1,l)

      resx=PartRnd(J,theta$rho_h)
      resy=PartRnd(K,theta$tau_l)
      u_ilx=x%*%resx$z_ik
      u_ily=y%*%resy$z_ik
      alpha_gh=u_ilx[sample(1:n,g),]/(matrix(1,g,1)%*%colSums(resx$z_ik))
      #alpha_gh=matrix(runif(g*m),ncol=m)
      alpha_gh[alpha_gh<eps]=eps
      alpha_gh[alpha_gh>1-eps]=1-eps
      theta$alpha_gh=alpha_gh
      beta_gl=u_ily[sample(1:n,g),]/(matrix(1,g,1)%*%colSums(resy$z_ik))
      #beta_gl=matrix(runif(g*s),ncol=s)
      beta_gl[beta_gl<eps]=eps
      beta_gl[beta_gl>1-eps]=1-eps
      theta$beta_gl=beta_gl
      list(r_jh=resx$z_ik,t_kl=resy$z_ik,theta=theta,v=resx$z,w=resy$z)
    }else if (any(names(start)=='theta')&& any(names(start)=='v') && any(names(start)=='w'))

        {theta=start$theta
        r_jh=matrix(mapply(as.numeric,(!(start$v%*%matrix(1,1,h)-matrix(1,J,1)%*%(1:h)))),ncol=h)
        t_kl=matrix(mapply(as.numeric,(!(start$w%*%matrix(1,1,l)-matrix(1,K,1)%*%(1:l)))),ncol=l)
        list(r_jh=r_jh,t_kl=t_kl,theta=theta)

     }else if (any(names(start)=='z')&& any(names(start)=='theta')){

            theta=start$theta
            z=start$z
            # Computation of t_jl
            s_ig=matrix(mapply(as.numeric,(!(z%*%matrix(1,1,g)-matrix(1,n,1)%*%(1:g)))),ncol=g)
            s_g=t(colSums(s_ig))
            eps10=10*eps
            log_rho_h=log(theta$rho_h)
            log_tau_l=log(theta$tau_l)
            alpha_gh=theta$alpha_gh
            beta_gl=theta$beta_gl
            A_gh=log((alpha_gh+eps10)/(1-alpha_gh+eps10))
            A_gl=log((beta_gl+eps10)/(1-beta_gl+eps10))
            B_gh=log(1-alpha_gh+eps10)
            B_gl=log(1-beta_gl+eps10)
            R_jh=(t(x)%*%s_ig)%*%A_gh+rep(1,J)%*%s_g%*%B_gh+matrix(rep(1,J),ncol=1)%*%t(log_rho_h)
            R_jhmax=apply(R_jh, 1,max)
            Tp_jh = exp(R_jh-R_jhmax)
            r_jh=Tp_jh/(rowSums(Tp_jh))

            T_kl=(t(y)%*%s_ig)%*%A_gl+rep(1,K)%*%s_g%*%B_gl+matrix(rep(1,K),ncol=1)%*%t(log_tau_l)
            T_klmax=apply(T_kl, 1,max)
            Tp_kl = exp(T_kl-T_klmax)
            t_kl=Tp_kl/(rowSums(Tp_kl))
            list(r_jh=r_jh,t_kl=t_kl,theta=theta)
          }else if  (any(names(start)=='z')&& any(names(start)=='v') && any(names(start)=='w')){
            if ((max(start$z)!=g)||(max(start$v)!=h)||(max(start$w)!=l)){
              warning('(Gstart,Hstart,start) need to match with the characteristics of userparam z, v and w')
            }else{


            s_ig=matrix(mapply(as.numeric,(!(start$z%*%matrix(1,1,g)-matrix(1,n,1)%*%(1:g)))),ncol=g)
            r_jh=matrix(mapply(as.numeric,(!(start$v%*%matrix(1,1,h)-matrix(1,J,1)%*%(1:h)))),ncol=h)
            t_kl=matrix(mapply(as.numeric,(!(start$w%*%matrix(1,1,l)-matrix(1,K,1)%*%(1:l)))),ncol=l)
            # Computation of parameters
            s_g=colSums(s_ig)
            r_h=colSums(r_jh)
            t_l=colSums(t_kl)
            if (fixed_proportion==TRUE){
              theta$pi_g=1/g * rep(1,g)
              theta$rho_h=1/h * rep(1,h)
              theta$tau_l=1/l * rep(1,l)
            }else{
              theta$pi_g=s_g/n
              theta$rho_h=r_h/J
              theta$tau_l=t_l/K}

            theta$alpha_gh = (t(s_ig)%*%x%*%r_jh)/(s_g%*%t(r_h))
            theta$beta_gl = (t(s_ig)%*%y%*%t_kl)/(s_g%*%t(t_l))
            list(r_jh=r_jh,t_kl=t_kl,theta=theta)}
          }else if (any(names(start)=='theta')&& any(names(start)=='r_jh') && any(names(start)=='t_kl'))
          {list(r_jh=r_jh,t_kl=t_kl,theta=theta) }
          else{warning('error, you need to provide a correct initialization with at least 2 parameters')}
  }
  maxeps=function(i){
    eps=10^(-16)
    max(eps,min(i,1-eps))
  }
  BlocXemAlpha_kl_MLBM=function(x,y,theta,r_jh,t_kl,a=a,b=b,niter=10000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){

    alpha_gh=theta$alpha_gh
    alpha_gh=matrix(mapply(maxeps,alpha_gh),ncol=dim(alpha_gh)[2])
    beta_gl=theta$beta_gl
    beta_gl=matrix(mapply(maxeps,beta_gl),ncol=dim(beta_gl)[2])
    pi_g=theta$pi_g
    rho_h=theta$rho_h
    tau_l=theta$tau_l
    g=dim(alpha_gh)[1]
    h=dim(alpha_gh)[2]
    l=dim(beta_gl)[2]
    n=dim(x)[1]
    J=dim(x)[2]
    K=dim(y)[2]
    r_h=t(colSums(r_jh))
    t_l=t(colSums(t_kl))


    # boucle sur les it?rations de l'algorithme
    W=-Inf
    oldtheta=theta
    for(iter in 1:niter){
      A_gh=log(alpha_gh/(1-alpha_gh))
      A_gl=log(beta_gl/(1-beta_gl))
      B_gh=log(1-alpha_gh)
      B_gl=log(1-beta_gl)
      log_pi_g=log(pi_g+(pi_g==0)*eps)
      log_rho_h=log(rho_h+(rho_h==0)*eps)
      log_tau_l=log(tau_l+(tau_l==0)*eps)

      W2=-Inf

      ############
      # Etape E
      #############

      for (iter_int in 1:niter_int){
        # Computation of sik
        S_ig=(x%*%r_jh)%*%t(A_gh)+rep(1,n)%*%r_h%*%t(B_gh)+matrix(1,n,1)%*%t(log_pi_g)+(y%*%t_kl)%*%t(A_gl)+rep(1,n)%*%t_l%*%t(B_gl)
        S_igmax=apply(S_ig, 1,max)
        Sp_ig = exp(S_ig-S_igmax)
        s_ig=Sp_ig/(rowSums(Sp_ig))

        #for ( iter_i in 1:n){
        #z[iter_i]=min(g,g+1- sum(rep(1,g)*runif(1) < cumsum(s_ig[iter_i,]) ))
        #}

        #s_ign=matrix(mapply(as.numeric,(!(z%*%matrix(1,1,g)-matrix(1,n,1)%*%(1:g)))),ncol=g)
        #s_gn=t(colSums(s_ign))
        #s_ig=s_ign
        #s_g=s_gn


        Hs=-sum(sum(s_ig*(s_ig>0)*log(s_ig*(s_ig>0)+(s_ig<=0)*eps)))
        s_g=t(colSums(s_ig))

        # Computation of t_jl
        R_jh=(t(x)%*%s_ig)%*%A_gh+rep(1,J)%*%s_g%*%B_gh+matrix(rep(1,J),ncol=1)%*%t(log_rho_h)
        R_jhmax=apply(R_jh, 1,max)
        Tp_jh = exp(R_jh-R_jhmax)
        r_jh=Tp_jh/(rowSums(Tp_jh))
        Htx=-sum(sum(r_jh*(r_jh>0)*log(r_jh+(r_jh<=0)*eps)))
        r_h=t(colSums(r_jh))

        T_kl=(t(y)%*%s_ig)%*%A_gl+rep(1,K)%*%s_g%*%B_gl+matrix(rep(1,K),ncol=1)%*%t(log_tau_l)
        T_klmax=apply(T_kl, 1,max)
        Tp_kl = exp(T_kl-T_klmax)
        t_kl=Tp_kl/(rowSums(Tp_kl))
        Hty=-sum(sum(t_kl*(t_kl>0)*log(t_kl+(t_kl<=0)*eps)))
        t_l=t(colSums(t_kl))

        #for (iter_j in 1:p){
        #v[iter_j]=min(m,m+1- sum(rep(1,m)*runif(1) < cumsum(r_jh[iter_j,]) ))
        #}

        #for (iter_k in 1:q){
        #w[iter_k]=min(s,s+1- sum(rep(1,s)*runif(1) < cumsum(t_kl[iter_k,]) ))
        #}

        #t_jlnx=matrix(mapply(as.numeric,(!(v%*%matrix(1,1,m)-matrix(1,p,1)%*%(1:m)))),ncol=m)
        #t_lnx=t(colSums(t_jlnx))
        #r_jh=t_jlnx
        #r_h=t_lnx

        #t_jlny=matrix(mapply(as.numeric,(!(w%*%matrix(1,1,s)-matrix(1,q,1)%*%(1:s)))),ncol=s)
        #t_lny=t(colSums(t_jlny))
        #t_kl=t_jlny
        #t_l=t_lny

        #empty_cluster=any(t_l<10.e-5)||any(s_g<10.e-5)

        W2_old=W2
        W2=sum(r_jh*R_jh)+sum(t_kl*T_kl)+s_g%*%log_pi_g  + Hs + Htx+Hty
        W2=W2[1,]
        if (abs((W2-W2_old)/W2) < epsi_int) break
      }


      ##############
      #  Etape M
      ##############"
      # Computation of parameters
      pi_g=t((s_g+a-1)/(n+g*(a-1)))
      rho_h=t((r_h+a-1)/(J+h*(a-1)))
      tau_l=t((t_l+a-1)/(K+l*(a-1)))
      u_klx=t(s_ig)%*%x%*%r_jh
      n_klx=t(s_g)%*%r_h
      u_klx=u_klx*(u_klx>0)
      n_klx=n_klx*(n_klx>0)
      alpha_gh=(u_klx+b-1)/(n_klx+2*(b-1))

      u_kly=t(s_ig)%*%y%*%t_kl
      n_kly=t(s_g)%*%t_l
      u_kly=u_kly*(u_kly>0)
      n_kly=n_kly*(n_kly>0)
      beta_gl=(u_kly+b-1)/(n_kly+2*(b-1))


      if (any(any((u_klx+b-1)<=0))){
        alpha_gh=(u_klx+b-1)/(n_klx+2*(b-1)+eps*(n_klx<=0))*((u_klx+b-1)>0)
      }
      alpha_gh=matrix(mapply(maxeps,alpha_gh),ncol=dim(alpha_gh)[2])

      if (any(any((u_kly+b-1)<=0))){
        beta_gl=(u_kly+b-1)/(n_kly+2*(b-1)+eps*(n_kly<=0))*((u_kly+b-1)>0)
      }
      beta_gl=matrix(mapply(maxeps,beta_gl),ncol=dim(beta_gl)[2])
      # Criterion
      W_old=W
      W=sum(s_g*log(s_g))- n*log(n) + sum(r_h*log(r_h))+sum(t_l*log(t_l))-J*log(J)-K*log(K)+
        sum(u_klx*log(u_klx)+(n_klx-u_klx)*log(n_klx-u_klx)-n_klx*log(n_klx))+
        sum(u_kly*log(u_kly)+(n_kly-u_kly)*log(n_kly-u_kly)-n_kly*log(n_kly))+
        Hs+Htx+Hty+
        (a-1)*(sum(log(pi_g))+sum(log(rho_h))+sum(log(tau_l)))+
        (b-1)*sum(log(alpha_gh)+log(1-alpha_gh))+(b-1)*sum(log(beta_gl)+log(1-beta_gl))

      if (is.nan(W)||(W=Inf)){
        if (b==1){
          W=sum(s_g*log(s_g+(s_g<=0)*eps))-n*log(n)+sum(r_h*log(r_h+(r_h<=0)*eps))+sum(t_l*log(t_l+(t_l<=0)*eps))-J*log(J)-K*log(K)+
            sum(u_klx*log(u_klx+(u_klx<=0)*eps)+(n_klx-u_klx)*((n_klx-u_klx)>0)*log((n_klx-u_klx)*((n_klx-u_klx)>0)+(((n_klx-u_klx)<=0)*eps))-n_klx*(n_klx>0)*log(n_klx+(n_klx==0)*eps))+
            sum(u_kly*log(u_kly+(u_kly<=0)*eps)+(n_kly-u_kly)*((n_kly-u_kly)>0)*log((n_kly-u_kly)*((n_kly-u_kly)>0)+(((n_kly-u_kly)<=0)*eps))-n_kly*(n_kly>0)*log(n_kly+(n_kly==0)*eps))+
            Hs + Htx+Hty+
            (a-1)*(sum(log(pi_g+(pi_g<=0))*(pi_g>0))+sum(log(rho_h+(rho_h<=0))*(rho_h>0))+sum(log(tau_l+(tau_l<=0))*(tau_l>0)))
        }else{
          W=sum(s_g*log(s_g+(s_g<=0)*eps))-n*log(n)+sum(r_h*log(r_h+(r_h<=0)*eps))+sum(t_l*log(t_l+(t_l<=0)*eps))-J*log(J)-K*log(K)+
            sum(u_klx*log(u_klx+(u_klx<=0)*eps)+(n_klx-u_klx)*((n_klx-u_klx)>0)*log((n_klx-u_klx)*((n_klx-u_klx)>0)+(((n_klx-u_klx)<=0)*eps))-n_klx*log(n_klx+(n_klx==0)*eps))+
            sum(u_kly*log(u_kly+(u_kly<=0)*eps)+(n_kly-u_kly)*((n_kly-u_kly)>0)*log((n_kly-u_kly)*((n_kly-u_kly)>0)+(((n_kly-u_kly)<=0)*eps))-n_kly*log(n_kly+(n_kly==0)*eps))+
            Hs + Htx+Hty+
            (a-1)*(sum(log(pi_g+(pi_g<=0))*(pi_g>0))+sum(log(rho_h+(rho_h<=0))*(rho_h>0))+sum(log(tau_l+(tau_l<=0))*(tau_l>0)))+
            (b-1)*sum(log(alpha_gh+(alpha_gh<=0))*(alpha_gh>0)+log(1-alpha_gh+(alpha_gh>=1))*(alpha_gh<1))+(b-1)*sum(log(beta_gl+(beta_gl<=0))*(beta_gl>0)+log(1-beta_gl+(beta_gl>=1))*(beta_gl<1))
        }
      }
      normtheta=max(max(abs(oldtheta$pi_g-theta$pi_g)),max(abs(oldtheta$rho_h-theta$rho_h)),max(abs(oldtheta$tau_l-theta$tau_l)),max((max(abs(oldtheta$alpha_gh-theta$alpha_gh)))),max((max(abs(oldtheta$beta_gl-theta$beta_gl)))))
      if ((abs((W-W_old)/W) < epsi) && (normtheta<10^(-7)))
      {
        break
      }
      else
      {
        oldtheta=theta
      }


      if (is.nan(W)) break
    }
    theta$alpha_gh=alpha_gh
    theta$beta_gl=beta_gl
    theta$pi_g=pi_g
    theta$rho_h=rho_h
    theta$tau_l=tau_l
    #z=apply(s_ig, 1, which.max)
    #w=apply(T_jl, 1, which.max)
    #print(c("W= ", W))
    ResVBayesAle=list(W=W,theta=theta,s_ig=s_ig,r_jh=r_jh,t_kl=t_kl,iter=iter,normtheta=normtheta)
  }

  BinBlocVBayes_MLBM=function(x,y,g,h,l,a,b,ini=list(n_init=1,ale=TRUE,user=FALSE),niter=100000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
    if (is.null(ini$theta)){
      if (!is.numeric(ini$n_init)){
        stop("For a random initialization, n_init must be numeric")
      }
      if (ini$n_init<=0){
        stop("n_init must be positive")
      }
    }
      W2=-Inf
      Res=list()
      for (i in 1:ini$n_init){
        #niter=1
        #niter_int=1
        #eps=10^(-16)
        #epsi=10^(-16)
        #epsi_int=10^(-16)
        initVBayes1=BinBlocInit_MLBM(g,h,l,x,y,start=ini)
        ResVBayesAle=BlocXemAlpha_kl_MLBM(x,y,initVBayes1$theta,initVBayes1$r_jh,initVBayes1$t_kl,a=a,b=b,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)

        ### ? enlever
        #print(c("ResVBayesAle$W= ", ResVBayesAle$W))
        W2=ResVBayesAle$W
        theta_max=ResVBayesAle$theta
        s_ig_max=ResVBayesAle$s_ig
        r_jh_max=ResVBayesAle$r_jh
        t_kl_max=ResVBayesAle$t_kl
        ### fin ? enlever
        if ((ResVBayesAle$W>W2)||(is.nan(ResVBayesAle$W))||(i=ini$n_init)||(W2=-Inf)){
          W2=ResVBayesAle$W
          theta_max=ResVBayesAle$theta
          s_ig_max=ResVBayesAle$s_ig
          r_jh_max=ResVBayesAle$r_jh
          t_kl_max=ResVBayesAle$t_kl

        }
      }# fin iteration

    theta_path=0
    theta=theta_max
    s_ig=s_ig_max
    r_jh=r_jh_max
    t_kl=t_kl_max

    #matrice de classification
    z=apply(s_ig, 1, which.max)
    v=apply(r_jh, 1, which.max)
    w=apply(t_kl, 1, which.max)

    #calcul des marges lignes et colonnes
    Repart=list()
    Repart$z=matrix(0,g,1)

    for (ig in 1:g){
      Repart$z[ig]=sum(z==ig)
    }
    Empty=list()
    Empty$z=sum(Repart$z==0)
    Repart$v=matrix(0,h,1)
    Repart$w=matrix(0,l,1)
    for (im in 1:h){
      Repart$v[im]=sum(v==im)
    }
    for (is in 1:l){
      Repart$w[is]=sum(w==is)
    }
    Empty$v=sum(Repart$v==0)
    Empty$w=sum(Repart$w==0)

    list(W=W2,theta=theta,s_ig=s_ig,r_jh=r_jh,t_kl=t_kl,z=z,v=v,w=w,Empty=Empty,Repart=Repart,iter=ResVBayesAle$iter,normtheta=ResVBayesAle$normtheta)
  }

  PartRnd = function(n,proba){

    g=length(proba)

    i=sample(1:n,n, replace = FALSE)
    i1=i[1:g]
    i2=i[(g+1):n]
    z=rep(0,n)
    z[i1]=1:g

    if (n > g) {
      z[i2]=(g+1)-colSums(rep(1,g)%*%t(runif(n-g)) < cumsum(proba)%*%t(rep(1,n-g)))
      z[z==(g+1)]=g
    }

    z_ik=!(z%*%t(rep(1,g))-rep(1,n)%*%t(1:g))

    list(z=z,z_ik=z_ik)
  }



  BinBlocICL =function (a,b,x,y,z1,v1,w1){

    b1x=b
    b1y=b
    b2x=b
    b2y=b
    n=dim(x)[1]
    p=dim(x)[2]
    q=dim(y)[2]
    if (!is.matrix(z1)){
      g=max(z1)
      z=!(z1%*%t(rep(1,g))-rep(1,n)%*%t(1:g))
    }else{
      z=z1
    }
    if (!is.matrix(v1)){
      h=max(v1)
      v=!(v1%*%t(rep(1,h))-rep(1,p)%*%t(1:h))
    }else{
      v=v1
    }

    if (!is.matrix(w1)){
      l=max(w1)
      w=!(w1%*%t(rep(1,l))-rep(1,q)%*%t(1:l))
    }else{
      w=w1
    }

    nk=colSums(z)
    nlx=colSums(v)
    Nklx=nk%*%t(nlx)
    nklx=t(z)%*%x%*%v

    nly=colSums(w)
    Nkly=nk%*%t(nly)
    nkly=t(z)%*%y%*%w
    #critere=lgamma(g*a)+lgamma(m*a)-(g+m)*lgamma(a)+m*g*(lgamma(2*b)-2*lgamma(b))-lgamma(n+g*a)-lgamma(d+m*a)+
    #sum(lgamma(a+nk))+sum(lgamma(a+nl))+
    #sum(lgamma(b+nkl)+lgamma(b+Nkl-nkl)-lgamma(2*b+Nkl))

    critere=lgamma(g*a)+lgamma(h*a)+lgamma(l*a)-(g+h+l)*lgamma(a)+g*h*(lgamma(b1x+b2x)-
                                                                         lgamma(b1x)*lgamma(b2x))+g*l*(lgamma(b1y+b2y)-lgamma(b1y)*lgamma(b2y))-lgamma(n+g*a)-lgamma(p+h*a)-lgamma(q+l*a)+
      sum(lgamma(a+nk))+sum(lgamma(a+nlx))+sum(lgamma(a+nly))+
      sum(sum(lgamma(b1x+nklx)+lgamma(b2x+Nklx-nklx)-lgamma(b1x+b2x+Nklx)))+
      sum(sum(lgamma(b1y+nkly)+lgamma(b1y+Nkly-nkly)-lgamma(b1y+b2y+Nkly)))
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

  criterion_tab=array(-Inf,dim=c(Gmax+1,Hmax+1,Lmax+1))
  W_tab=array(-Inf,dim=c(Gmax+1,Hmax+1,Lmax+1))
  W=-Inf
  if (init_choice=='smallVBayes'){
    # if (mc.cores>1){
    #   Temp=parallel::mclapply(1:ntry,function(i){PoissonBlocVBayes(data,Gstart,Hstart,a,alpha,beta,
    #                                                    normalization,ini=list(n_init=1,ale=TRUE),
    #                                                    niter=50)}, mc.cores=mc.cores)
    # }else{
    #   Temp=lapply(1:ntry,function(i){PoissonBlocVBayes(data,Gstart,Hstart,a,alpha,beta,
    #                                                    normalization,ini=list(n_init=1,ale=TRUE),
    #                                                    niter=50)})
    # }

    for (i in 1:ntry){
      #res=BlocGibbs(data$x,a,b,ini=list(n_init=1,ale=TRUE),1000)
      res=BinBlocVBayes_MLBM(x,y,Gstart,Hstart,Lstart,a,b
                            ,ini=list(n_init=1,ale=TRUE),
                            niter=50)
      if (res$W>W && res$Empty$v==FALSE && res$Empty$w==FALSE && res$Empty$z==FALSE){
        W=res$W
        relopt=res
      }
    }
      if (!exists("relopt")){
        stop(paste('The algorithm does not manage to find a configuration with  (Gstart=',as.character(Gstart),',Hstart=',as.character(Hstart),',Lstart=',as.character(Lstart),')' ))
      }

  # }else if (init_choice=='Gibbs'){
  #
  #   resinit=PoissonBlocGibbs(data,Gstart,Hstart,a,normalization,alpha,beta,niter=5000,inistart=list(n_init=10,ale=TRUE))
  #   initVBayes=list(z=resinit$z,w=resinit$w,ale=FALSE,n_init=1)
  #   relopt=PoissonBlocVBayes(data,Gstart,Hstart,a,alpha,beta,normalization,ini=initVBayes)
  #   if (relopt$empty_cluster==TRUE){
  #     stop(paste('The algorithm does not manage to find a configuration with  (Gstart=',as.character(Gstart),',Hstart=',as.character(Hstart),'),try normalization=FALSE') )
  #   }
  }else if (init_choice=='user'){
    if (is.null(userparam)){
      stop(paste('If you choose the user init choice, please fill the useparam field with partitions z, v and w' ))
    }else{
      initVBayes=list(z=userparam$z,v=userparam$v,w=userparam$w,ale=FALSE,n_init=1)
      relopt=BinBlocVBayes_MLBM(x,y,Gstart,Hstart,Lstart,a,b,ini=initVBayes)

    }
  }else if (init_choice=='random'){

    eps=10^(-16)
    n=dim(x)[1]
    J=dim(x)[2]
    K=dim(y)[2]
    # Random_Temp=function(Temp){
    #   theta=list()
    #   theta$pi_g=runif(Gstart)
    #   theta$pi_g=theta$pi_g/sum(theta$pi_g)
    #   theta$rho_l=runif(Hstart)
    #   theta$rho_l=theta$rho_l/sum(theta$rho_l)
    #   resw=PartRnd(d,theta$rho_l)
    #   t_jl=resw$z_ik
    #   w=resw$z
    #   resz=PartRnd(n,theta$pi_g)
    #   s_ig=resz$z_ik
    #   z=resz$z
    #   if (normalization){
    #     mu_i=matrix(rowSums(data$x),ncol=1)
    #     nu_j=matrix(colSums(data$x),ncol=d)
    #   }
    #   else{mu_i=matrix(rep(1,n),ncol=1)
    #   nu_j=matrix(rep(1,d),ncol=d)
    #   }
    #   n_kl=t((t(mu_i)%*%s_ig))%*%(nu_j%*%t_jl)
    #   theta$lambda_kl=(t(s_ig)%*%data$x%*%t_jl)/(n_kl)
    #   if (any(is.nan(theta$lambda_kl))){
    #     lambda0=mean(mean(data$x))
    #     theta$lambda_kl=theta$lambda_kl=matrix(rpois(g,lambda0),ncol=m)
    #   }
    #
    #   initVBayes=list(t_jl=t_jl,w=w,theta=theta,ale=FALSE,n_init=1)
    #   #res=BlocGibbs(data$x,a,b,ini=list(n_init=1,ale=TRUE),1000)
    #   PoissonBlocVBayes(data,Gstart,Hstart,a,alpha,beta,normalization,ini=initVBayes)
    # }
    #
    # if (mc.cores>1){
    #   Temp=parallel::mclapply(1:ntry,Random_Temp, mc.cores=mc.cores)
    # }else{
    #   Temp=lapply(1:ntry,Random_Temp)
    # }


    for (i in 1:ntry){
      theta=list()
      #theta$pi_g=runif(g)
      #theta$pi_g=theta$pi_g/sum(theta$pi_g)
      #theta$rho_h=runif(m)
      #theta$rho_h=theta$rho_h/sum(theta$rho_h)
      #theta$tau_l=runif(s)
      #theta$tau_l=theta$tau_l/sum(theta$tau_l)

      theta$pi_g=1/g *rep(1,g)
      theta$rho_h=1/h *rep(1,h)
      theta$tau_l=1/l *rep(1,l)

      resx=PartRnd(J,theta$rho_h)
      resy=PartRnd(K,theta$tau_l)
      u_ilx=x%*%resx$z_ik
      u_ily=y%*%resy$z_ik
      alpha_gh=u_ilx[sample(1:n,g),]/(matrix(1,g,1)%*%colSums(resx$z_ik))
      #alpha_gh=matrix(runif(g*m),ncol=m)
      alpha_gh[alpha_gh<eps]=eps
      alpha_gh[alpha_gh>1-eps]=1-eps
      theta$alpha_gh=alpha_gh
      beta_gl=u_ily[sample(1:n,g),]/(matrix(1,g,1)%*%colSums(resy$z_ik))
      #beta_gl=matrix(runif(g*s),ncol=s)
      beta_gl[beta_gl<eps]=eps
      beta_gl[beta_gl>1-eps]=1-eps
      theta$beta_gl=beta_gl
      #list(r_jh=resx$z_ik,t_kl=resy$z_ik,theta=theta,v=resx$z,w=resy$z)

      initVBayes=list(r_jh=resx$z_ik,t_kl=resy$z_ik,theta=theta,v=resx$z,w=resy$z,ale=FALSE,n_init=1)
      #res=BlocGibbs(data$x,a,b,ini=list(n_init=1,ale=TRUE),1000)
      res=BinBlocVBayes_MLBM(x,y,Gstart,Hstart,Lstart,a,b,ini=initVBayes)
      if (res$W>W && res$Empty$z==FALSE && res$Empty$v==FALSE && res$Empty$w==FALSE){
        relopt=res}



    }
    if (!exists("relopt")){
      stop(paste('The algorithm does not manage to find a configuration with  (Gstart=',as.character(Gstart),',Hstart=',as.character(Hstart),',','Lstart=',as.character(Lstart),')') )
    }

   }else {stop("This is not a initialization choice proposed by the procedure")}

  if (criterion_choice=='ICL'){
     Critere= BinBlocICL_MLBM(a,b,x,y,relopt$z,relopt$v,relopt$w)
  # else if (criterion_choice=='BIC'){
  #   Critere= PoissonBlocBIC(a,alpha,beta,data,relopt,normalization)
  # }
  }else {stop('This is not a criterion choice proposed by the procedure')}
  criterion_tab[Gstart,Hstart,Lstart]=Critere
  W_tab[Gstart,Hstart,Lstart]=relopt$W

  modele=relopt
  criterion_max=Critere
  g=Gstart
  h=Hstart
  l=Lstart
  gopt=Gstart
  hopt=Hstart
  lopt=Lstart
  model_max=relopt

  while (g+1<=Gmax && h+1<=Hmax && l+1<=Lmax){
    if ((length(unique(modele$z))!=g)||(length(unique(modele$v))!=h)||(length(unique(modele$w))!=l)){
      stop(paste('The algorithm does not manage to find a configuration with (g=',as.character(g),',h=',as.character(h),',s=',as.character(l),'), maybe reduce Gmax, Hmax or Lmax') )
    }

    W_colonnex=-Inf
    W_colonney=-Inf
    W_ligne=-Inf

    f1=function(Kc){
      v=modele$v
      Xc=which(v==Kc)
      lxc=length(Xc)
      if (lxc>1){
        idxc=sample(lxc,floor(lxc/2),replace=FALSE)
        indc=Xc[idxc]
        v[indc]=h+1
        initVBayes=list(z=modele$z,v=v,w=modele$w,ale=FALSE,n_init=1)
        resultsVBayes=BinBlocVBayes_MLBM(x,y,g,h+1,l,a,b,ini=initVBayes)
        iterverif=0
        while ((resultsVBayes$Empty$z==TRUE || resultsVBayes$Empty$v==TRUE || resultsVBayes$Empty$w==TRUE ) && iterverif<20){
          iterverif=iterverif+1
          v=modele$v
          idxc=sample(lxc,floor(lxc/2),replace=FALSE)
          indc=Xc[idxc]
          v[indc]=h+1
          initVBayes=list(z=modele$z,v=v,w=modele$w,ale=FALSE,n_init=1)
          resultsVBayes=BinBlocVBayes_MLBM(x,y,g,h+1,l,a,b,ini=initVBayes)
        }
        if (iterverif<20){
          return(list(W=resultsVBayes$W,modele_colonnex=resultsVBayes))
        }else{
          return(list(W=-Inf))
        }
      }else{
        return(list(W=-Inf))
      }
    }
    if (mc.cores>1){
      l_colonnex=parallel::mclapply(1:h,f1, mc.cores=mc.cores)
    }else{
      l_colonnex=lapply(1:h,f1)
    }

    for (Kc in 1:h){
      if (l_colonnex[[Kc]]$W>W_colonnex){
        W_colonnex=l_colonnex[[Kc]]$W
        modele_colonnex=l_colonnex[[Kc]]$modele_colonnex
      }
    }
    if (W_colonnex<=-Inf){
      Empty_m=TRUE
    }else{
      Empty_m=FALSE
    }
    #### Class g
    f2=function(Kl){
      z=modele$z
      Xl=which(z==Kl)
      ll=length(Xl)
      if (ll>1){
        idxl=sample(ll,floor(ll/2),replace=FALSE)
        indl=Xl[idxl]
        z[indl]=g+1
        initVBayes=list(z=z,v=modele$v,w=modele$w,ale=FALSE,n_init=1)
        resultsVBayes=BinBlocVBayes_MLBM(x,y,g+1,h,l,a,b,ini=initVBayes)
        iterverif=0
        while ((resultsVBayes$Empty$w==TRUE || resultsVBayes$Empty$v==TRUE || resultsVBayes$Empty$z==TRUE)  && iterverif<20){
          iterverif=iterverif+1
          z=modele$z
          idxl=sample(ll,floor(ll/2),replace=FALSE)
          indl=Xl[idxl]
          z[indl]=g+1
          initVBayes=list(z=z,v=modele$v,w=modele$w,ale=FALSE,n_init=1)
          resultsVBayes=BinBlocVBayes_MLBM(x,y,g+1,h,l,a,b,ini=initVBayes)
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
      l_ligne=parallel::mclapply(1:g,f2, mc.cores=mc.cores)
    }else{
      l_ligne=lapply(1:g,f2)
    }
    for (Kl in 1:g){
      if (l_ligne[[Kl]]$W>W_ligne){
        W_ligne=l_ligne[[Kl]]$W
        modele_ligne=l_ligne[[Kl]]$modele_ligne
      }
    }
##Class l
    f3=function(Kcy){
      w=modele$w
      Xcy=which(w==Kcy)
      lxcy=length(Xcy)
      if (lxcy>1){
        idxcy=sample(lxcy,floor(lxcy/2),replace=FALSE)
        indcy=Xcy[idxcy]
        w[indcy]=l+1
        initVBayes=list(z=modele$z,v=modele$v,w=w,ale=FALSE,n_init=1)
        resultsVBayes=BinBlocVBayes_MLBM(x,y,g,h,l+1,a,b,ini=initVBayes)
        iterverif=0
        while ((resultsVBayes$Empty$w==TRUE || resultsVBayes$Empty$v==TRUE || resultsVBayes$Empty$z==TRUE)  && iterverif<20){
          iterverif=iterverif+1
          w=modele$w
          idxcy=sample(lxcy,floor(lxcy/2),replace=FALSE)
          indcy=Xcy[idxcy]
          w[indcy]=l+1
          initVBayes=list(z=modele$z,v=modele$v,w=w,ale=FALSE,n_init=1)
          resultsVBayes=BinBlocVBayes_MLBM(x,y,g,h,l+1,a,b,ini=initVBayes)
        }
        if (iterverif<20){
          return(list(W=resultsVBayes$W,modele_colonney=resultsVBayes))
        }else{
          return(list(W=-Inf))
        }
      }else{
        return(list(W=-Inf))
      }
    }
    if (mc.cores>1){
      l_colonney=parallel::mclapply(1:l,f3, mc.cores=mc.cores)
    }else{
      l_colonney=lapply(1:l,f3)
    }

    for (Kcy in 1:l){
      if (l_colonney[[Kcy]]$W>W_colonney){
        W_colonney=l_colonney[[Kcy]]$W
        modele_colonney=l_colonney[[Kcy]]$modele_colonney
      }
    }
    if (W_colonney<=-Inf){
      Empty_s=TRUE
    }else{
      Empty_s=FALSE
    }




    if ((W_ligne<=-Inf)&(Empty_m) &(Empty_s)){
      warning("The algorithm does not manage to find a configuration  (",as.character(g+1),',',as.character(h),',',as.character(l), ") nor (",as.character(g),',',as.character(h+1),',',as.character(l),'nor (', as.character(g),',',as.character(h),',',as.character(l+1),")")
      new(Class="BIKM1_MLBM_Binary",init_choice=init_choice,criterion_choice=criterion_choice,criterion_tab=criterion_tab,criterion_max=criterion_max,model_max=model_max,W_tab=W_tab,gopt=gopt,hopt=hopt,lopt=lopt)
    }



    if (W_ligne>W_colonnex && W_ligne>W_colonney){
      modele=modele_ligne
      if (criterion_choice=='ICL'){
        Critere= BinBlocICL_MLBM(a,b,x,y,modele$z,modele$v,modele$w)
      # }else if (criterion_choice=='BIC'){
      #   Critere= PoissonBlocBIC(a,alpha,beta,data,modele,normalization)
      }
      criterion_tab[g+1,h,l]=Critere
      W_tab[g+1,h,l]=modele$W
      if (Critere>criterion_max){
        model_max=modele
        criterion_max=Critere
        gopt=g+1
        hopt=h
        lopt=l
      }
      if (verbose){
        cat('(g,h,l)=(',as.character(g+1),',',as.character(h),',',as.character(l),')\n',sep = "")
      }
      g=g+1
    }else if (W_colonnex>W_ligne && W_colonnex>W_colonney) {
      modele=modele_colonnex
      if (criterion_choice=='ICL'){
        Critere= BinBlocICL_MLBM(a,b,x,y,modele$z,modele$v,modele$w)
      # }else if (criterion_choice=='BIC'){
      #   Critere= PoissonBlocBIC(a,alpha,beta,data,modele,normalization)
      }
      criterion_tab[g,h+1,l]=Critere
      W_tab[g,h+1,l]=modele$W
      if (Critere>criterion_max){
        model_max=modele
        criterion_max=Critere
        gopt=g
        hopt=h+1
        lopt=l

      }
      if (verbose){
        cat('(g,h,l)=(',as.character(g),',',as.character(h+1),',',as.character(l),')\n',sep = "")
      }

      h=h+1
    }
    else {modele=modele_colonney
    if (criterion_choice=='ICL'){
      Critere= BinBlocICL_MLBM(a,b,x,y,modele$z,modele$v,modele$w)
      # }else if (criterion_choice=='BIC'){
      #   Critere= PoissonBlocBIC(a,alpha,beta,data,modele,normalization)
    }
    criterion_tab[g,h,l+1]=Critere
    W_tab[g,h,l+1]=modele$W
    if (Critere>criterion_max){
      model_max=modele
      criterion_max=Critere
      gopt=g
      hopt=h
      lopt=l+1

    }
    if (verbose){
      cat('(g,h,l)=(',as.character(g),',',as.character(h),',',as.character(l+1),')\n',sep = "")
    }

    l=l+1

    }
}
  if (verbose){
    cat('The selected row (g), 1st matrix column (h)  and 2nd matrix column (l) clusters are (g,h,l)=(',as.character(gopt),',',as.character(hopt),',',as.character(lopt),')\n',sep = "")
  }
  new(Class="BIKM1_MLBM_Binary",init_choice=init_choice,criterion_choice=criterion_choice,criterion_tab=criterion_tab,criterion_max=criterion_max,model_max=model_max,W_tab=W_tab,gopt=gopt,hopt=hopt,lopt=lopt)
}

##' BinBlocRnd_MLBM function for binary double data matrix simulation
##'
##' Produce two simulated data matrices generated under the Binary Multiple Latent Block Model.
##'
##' @param n a positive integer specifying the number of expected rows.
##' @param J a positive integer specifying the number of expected columns of the first matrix.
##' @param K a positive integer specifying the number of expected columns of the second matrix.
##' @param theta a list specifying the model parameters:
##'
##' \code{pi_g}: a vector specifying the row mixing proportions.
##'
##' \code{rho_h}: a vector specifying the first matrix column mixing proportions.
##'
##' \code{tau_l}: a vector specifying the second matrix column mixing proportions.
##'
##' \code{alpha_gh}: a matrix specifying the distribution parameter of the first matrix.
##'
##'  \code{beta_gl}: a matrix specifying the distribution parameter of the second matrix.
##' @usage BinBlocRnd_MLBM(n,J,K,theta)
##' @return a list including the arguments:
##'
##' \code{x}: simulated first  data matrix.
##' \code{y}: simulated second data matrix.
##'
##' \code{xrow}: numeric vector specifying row partition.
##'
##' \code{xcolx}: numeric vector specifying first matrix column partition.
##'
##'  \code{xcoly}: numeric vector specifying second matrix column partition.

##' @rdname BinBlocRnd_MLBM-proc
##'
##' @examples
##' require(bikm1)
##'  set.seed(42)
  ##' n=200
  ##' J=120
  ##' K=120
  ##' g=3
  ##' h=2
  ##' l=2
  ##' theta=list()
  ##' theta$pi_g=1/g *matrix(1,g,1)
  ##' theta$rho_h=1/h *matrix(1,h,1)
  ##' theta$tau_l=1/l *matrix(1,l,1)
##'theta$alpha_gh=matrix(runif(6),ncol=h)
##'theta$beta_gl=matrix(runif(6),ncol=l)
##' data=BinBlocRnd_MLBM(n,J,K,theta)
##'
##' @export BinBlocRnd_MLBM


BinBlocRnd_MLBM =function (n,J,K,theta){
  PartRnd = function(n,proba){

    g=length(proba)

    i=sample(1:n,n, replace = FALSE)
    i1=i[1:g]
    i2=i[(g+1):n]
    z=rep(0,n)
    z[i1]=1:g

    if (n > g) {
      z[i2]=(g+1)-colSums(rep(1,g)%*%t(runif(n-g)) < cumsum(proba)%*%t(rep(1,n-g)))
      z[z==(g+1)]=g
    }

    z_ik=!(z%*%t(rep(1,g))-rep(1,n)%*%t(1:g))

    list(z=z,z_ik=z_ik)
  }
  resz=PartRnd(n,theta$pi_g)
  z=resz$z
  z_ik=resz$z_ik

  resv=PartRnd(J,theta$rho_h)
  v=resv$z
  w_jlx=resv$z_ik

  resw=PartRnd(K,theta$tau_l)
  w=resw$z
  w_jly=resw$z_ik
  g=dim(theta$alpha_gh)[1]
  h=dim(theta$alpha_gh)[2]
  l=dim(theta$beta_gl)[2]
  x=matrix(rbinom(n*J,1,z_ik%*%theta$alpha_gh%*%t(w_jlx)),ncol=J)
  y=matrix(rbinom(n*K,1,z_ik%*%theta$beta_gl%*%t(w_jly)),ncol=K)
  xrow=as.numeric(z_ik%*%c(1:g))
  xcolx=as.numeric(w_jlx%*%c(1:h))
  xcoly=as.numeric(w_jly%*%c(1:l))
  list(x=x,y=y,xrow=xrow,xcolx=xcolx,xcoly=xcoly)


}

##' BinBlocVisu_MLBM function for visualization of double matrix datasets
##'
##' Produce a plot object representing the co-clustered data-sets.
##'
##' @param x first data matrix of observations.
##' @param y second data matrix of observations.
##' @param z a numeric vector specifying the class of rows.
##' @param v a numeric vector specifying the class of columns (1rst matrix).
##' @param w a numeric vector specifying the class of columns (2nd matrix).
##' @usage BinBlocVisu_MLBM(x,y,z,v,w)
##' @rdname BinBlocVisu_MLBM-proc
##' @return a \pkg{plot} object
##' @examples
##'
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' K=120
##' g=3
##' h=2
##' l=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##'theta$alpha_gh=matrix(runif(6),ncol=h)
##'theta$beta_gl=matrix(runif(6),ncol=l)
##' data=BinBlocRnd_MLBM(n,J,K,theta)
##' BinBlocVisu_MLBM(data$x,data$y, data$xrow,data$xcolx,data$xcoly)
##' @export BinBlocVisu_MLBM

BinBlocVisu_MLBM=function(x,y,z,v,w){
  ## Initialisation
  xbase=x
  ybase=y
  g=max(z)
  hbase=max(v)
  l=max(w)
  n=dim(xbase)[1]
  p=dim(xbase)[2]
  q=dim(ybase)[2]
  vw=c(v,w+max(v))
  x=cbind(xbase,ybase)
  d=p+q
  h=hbase+l
  ii=c()
  jj=c()
  for (i in 1:g){
    ii=c(ii,which(z==i))
  }
  for (j in 1:h){
    jj=c(jj,which(vw==j))
  }
  z_ik=((z%*%t(rep(1,g)))==(rep(1,n)%*%t(1:g)))
  #z_ik=z%*%matrix(1,1,g)
  w_jl=((vw%*%t(rep(1,h)))==(rep(1,d)%*%t(1:h)))

  ## Affichage
  n_k=n-sort(cumsum(colSums(z_ik))-0.5,decreasing=TRUE)
  n_l=cumsum(colSums(w_jl))+0.5

  table.paint(x[ii,jj],clegend=0)

  for (i in n_k[2:g]){
    lines(c(0.5,(d+1)),c(i,i),col='blue',lwd=2)
  }
  for (i in n_l[1:(h-1)]){
    lines(c(i,i),c(0.5,(n+1)),col='blue',lwd=2)
  }
  for (i in n_l[hbase]){
    lines(c(i,i),c(0.5,(n+1)),col='red',lwd=2)
  }
}

##' BinBlocVisuResum_MLBM function  for visualization of double matrix datasets
##'
##' Produce a plot object representing the resumed co-clustered data-sets.
##'
##' @param x binary matrix of observations.
##' @param y binary second matrix of observations.
##' @param z a numeric vector specifying the class of rows.
##' @param v a numeric vector specifying the class of columns (1rst matrix).
##' @param w a numeric vector specifying the class of columns (2nd matrix).
##' @rdname BinBlocVisuResum_MLBM-proc
##' @usage BinBlocVisuResum_MLBM(x,y,z,v,w)
##' @return a \pkg{plot} object.
##' @examples
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' K=120
##' g=3
##' h=2
##' l=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##'theta$alpha_gh=matrix(runif(6),ncol=h)
##'theta$beta_gl=matrix(runif(6),ncol=l)
##' data=BinBlocRnd_MLBM(n,J,K,theta)
##' BinBlocVisuResum_MLBM(data$x,data$y, data$xrow,data$xcolx,data$xcoly)
##' @export BinBlocVisuResum_MLBM


BinBlocVisuResum_MLBM=function(x,y,z,v,w){
  ## Initialisation
  xbase=x
  ybase=y
  g=max(z)
  hbase=max(v)
  l=max(w)
  n=dim(xbase)[1]
  p=dim(xbase)[2]
  q=dim(ybase)[2]
  w=c(v,w+max(v))
  x=cbind(xbase,ybase)
  d=p+q
  h=hbase+l
  #m=max(w)
  z_ik=((z%*%t(rep(1,g)))==(rep(1,n)%*%t(1:g)))
  w_jl=((w%*%t(rep(1,h)))==(rep(1,d)%*%t(1:h)))


  ## Affichage
  nk=colSums(z_ik)
  nlx=colSums(w_jl)
  Nklx=nk%*%t(nlx)
  nklx=t(z_ik)%*%x%*%w_jl
  n_k=n-sort(cumsum(colSums(z_ik))-0.5,decreasing=TRUE)
  n_l=cumsum(colSums(w_jl))+0.5

  table.paint(nklx/Nklx)
  for (i in n_k[2:g]){
    lines(c(0.5,(d+1)),c(i,i),col='blue',lwd=2)
  }
  for (i in n_l[1:(h-1)]){
    lines(c(i,i),c(0.5,(n+1)),col='blue',lwd=2)
  }
  for (i in n_l[hbase]){
    lines(c(i,i),c(0.5,(n+1)),col='red',lwd=2)
  }


}

##' CE_simple function for agreement between clustering partitions
##'
##' Produce a measure of agreement between two  partitions for clustering. A value of 1 means a perfect match.
##'
##' @param z numeric vector  specifying the class of  rows.
##' @param zprime numeric vector  specifying the class of  rows.
##'
##' @return the value of the index.(between 0 and 1). A value of 0 corresponds to a perfect match.
##'
##' @usage CE_simple(z,zprime)
##' @rdname CE_simple-proc
##' @examples
##' \donttest{
##' require(bikm1)
##' set.seed(42)
##' z=floor(runif(4)*3)
##' zprime=floor(runif(4)*3)
##' error=CE_simple(z,zprime)
##' }
##' @export CE_simple



CE_simple=function(z,zprime){
  #PartDist returns the distance between the two partitions.
  #  This distance is computed with the best arrangement.
  #
  #Inputs:
  # z  : 1st partition
  # zprime  : 2nd partition
  #Outputs:
  # d  :  distance
  # p1 :  best permutation
  # p2 :  inverse best permutation

  #Functions:
  #31/05/2010, Gerard Govaert, UTC
  #####################################################"

  #t='Adlab:PartDist';
  #if nargin ~= 2, error(t,'2 arguments.'); end

  # n1=length(z)
  # n2=length(zprime)
  # if (n1!=n2) warning('Size of partitions incorrect.')
  # g1=max(z)
  # permutation=permutations(factorial(g1),g1,1:g1)
  # nperm=dim(permutation)[1]
  # dd=matrix(0,1,nperm)
  # for (i in 1:nperm){
  #   p=t(permutation[i,])
  #   dd[1,i]=sum(p[1,z]!=zprime)
  # }
  # d=min(dd)
  # i=which.min(dd)
  # d=d/n1
  # p1=permutation[i,]
  # tmp=sort(p1,index.return=TRUE)
  # p2=tmp$ix
  # list(d=d,p1=p1,p2=p2)

   g=max(length(unique(z)),length(unique(zprime)))
    n=length(z)
    if (max(z)>g){
      val=sort(unique(z),decreasing = FALSE)
      for (it in 1:length(val)){
        if (val[it]!=it){
          z[z==val[it]]=it
        }
      }
    }
    zsik=!(zprime%*%t(rep(1,g))-rep(1,n)%*%t(1:g));
    zik=!(z%*%t(rep(1,g))-rep(1,n)%*%t(1:g));
    if (g<7){
      #P<-permute::allPerms(g)
      P<-permutations(g,g,1:g)
    }else{
      #P1<-permute::allPerms(6)
      P1<-permutations(6,6,1:6)
      if (g>=7){
        P=cbind(rep(7,nrow(P1)),P1)
        for (it in 2:7){
          P2=matrix(7,nrow=nrow(P1),ncol=7)
          P2[,-it]=P1
          P<-rbind(P,P2)
        }
        if (g>=8){
          P1=P
          P=cbind(rep(8,nrow(P1)),P1)
          for (it in 2:8){
            P2=matrix(8,nrow=nrow(P1),ncol=8)
            P2[,-it]=P1
            P<-rbind(P,P2)
          }
          if (g==9){
            P1=P
            P=cbind(rep(9,nrow(P1)),P1)
            for (it in 2:9){
              P2=matrix(9,nrow=nrow(P1),ncol=9)
              P2[,-it]=P1
              P<-rbind(P,P2)
            }
          }
        }
      }
    }
    d=Inf
    for (iter in 1:nrow(P)){
      Val=sum(abs(zik-zsik[,P[iter,]]))/(2*n)
      if (Val<d){
        d=Val
        if (Val==0){
          break
        }
      }
    }
    d








  }


##' CE_LBM function for agreement between co-clustering partitions
##'
##' Produce a measure of agreement between two pairs of partitions for co-clustering. A value of 1 means a perfect match.
##'
##' @param z numeric vector  specifying the class of  rows.
##' @param w numeric vector  specifying the class of  columns.
##' @param zprime numeric vector  specifying another partition of  rows.
##' @param wprime numeric vector  specifying another partition of  columns.
##'
##' @return the value of the index. (between 0 and 1). A value of 0 corresponds to a perfect match.
##'
##' @usage CE_LBM(z,w,zprime,wprime)
##' @rdname CE_LBM-proc
##' @examples
##' \donttest{
##' require(bikm1)
##' set.seed(42)
##' z=floor(runif(4)*2)
##' zprime=floor(runif(4)*2)
##' w=floor(runif(4)*3)
##' wprime=floor(runif(4)*3)
##' error=CE_LBM(z,w,zprime,wprime)
##' }
##' @export CE_LBM

CE_LBM=function(z,w,zprime,wprime){


  L_z=length(z)
  L_w=length(w)
  if (L_z!=length(zprime)){
    warning('Both partitions z and zprime must contain the same number of points.')
  }
  if (L_w!=length(wprime)){
    warning('Both partitions w and wprime must contain the same number of points.')
  }


  e1=CE_simple(z,zprime)
  e2=CE_simple(w,wprime)
  CE=e1+e2-e1*e2
  CE
}




##' CE_MLBM function for agreement between co-clustering partitions in the MBLM
##'
##' Produce a measure of agreement between two triplets of partitions for co-clustering. A value of 1 means a perfect match.
##'
##' @param z numeric vector  specifying the class of  rows.
##' @param v numeric vector specifying the class of column partitions for the first matrix.
##' @param w numeric vector specifying the class of column partitions for the second matrix.
##' @param zprime numeric vector  specifying another partitions of rows.
##' @param vprime numeric vector specifying another partition of columns for the first matrix.
##' @param wprime numeric vector specifying another partition of columns for the second matrix.
##' @return the value of the index (between 0 and 1). A value of 0 corresponds to a perfect match.
##'
##' @usage CE_MLBM(z,v,w,zprime,vprime,wprime)
##' @rdname CE_MLBM-proc
##' @examples
##' \donttest{
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' K=120
##' g=2
##' h=2
##' l=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##'theta$alpha_gh=matrix(runif(4),ncol=h)
##'theta$beta_gl=matrix(runif(4),ncol=l)
##' data=BinBlocRnd_MLBM(n,J,K,theta)
##' res=BIKM1_MLBM_Binary(data$x,data$y,2,2,2,4,init_choice='smallVBayes')
##' error=CE_MLBM(res@model_max$z,res@model_max$v,res@model_max$w,data$xrow,data$xcolx,data$xcoly)
##' }
##' @export CE_MLBM




CE_MLBM=function(z,v,w,zprime,vprime,wprime){
  L_z=length(z)
  L_v=length(v)
  L_w=length(w)
  if (L_z!=length(zprime)){
    warning('Both partitions z and zprime must contain the same number of points.')
  }
  if (L_v!=length(vprime)){
    warning('Both partitions v and vprime must contain the same number of points.')
  }
  if (L_w!=length(wprime)){
    warning('Both partitions w and wprime must contain the same number of points.')
  }


  w=c(v,w+max(v))
  wprime=c(vprime,wprime+max(vprime))

  CE=CE_LBM(z,w,zprime,wprime)
  CE
}


##' BinBlocICL_MLBM function  for computation of the ICL criterion in the MLBM
##'
##' Produce a plot object representing the resumed co-clustered data-sets.
##'
##' @param a an hyperparameter for priors on the mixing proportions. By default, a=4.
##' @param b an hyperparameter for prior on the Bernoulli parameter. By default, b=1.
##' @param x binary matrix of observations (1rst matrix).
##' @param y binary matrix of observations (2nd matrix).
##'

##' @param z1 a numeric vector specifying the class of rows.
##' @param v1 a numeric vector specifying the class of columns (1rst matrix).
##' @param w1 a numeric vector specifying the class of columns (2nd matrix).
##' @rdname BinBlocICL_MLBM-proc
##' @usage BinBlocICL_MLBM(a,b,x,y,z1,v1,w1)
##' @return a value of the ICL criterion.
##' @examples
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' K=120
##' g=2
##' h=2
##' l=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##'theta$alpha_gh=matrix(runif(4),ncol=h)
##'theta$beta_gl=matrix(runif(4),ncol=l)
##' data=BinBlocRnd_MLBM(n,J,K,theta)
##' res=BIKM1_MLBM_Binary(data$x,data$y,2,2,2,4,init_choice='smallVBayes')
##' BinBlocICL_MLBM(a=4,b=1,data$x,data$y, data$xrow,data$xcolx,data$xcoly)
##' @export BinBlocICL_MLBM





BinBlocICL_MLBM =function (a,b,x,y,z1,v1,w1){
  b1x=b
  b1y=b
  b2x=b
  b2y=b

  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(y)[2]
  if (!is.matrix(z1)){
    g=max(z1)
    z=!(z1%*%t(rep(1,g))-rep(1,n)%*%t(1:g))
  }else{
    z=z1
  }
  if (!is.matrix(v1)){
    h=max(v1)
    v=!(v1%*%t(rep(1,h))-rep(1,p)%*%t(1:h))
  }else{
    v=v1
  }

  if (!is.matrix(w1)){
    l=max(w1)
    w=!(w1%*%t(rep(1,l))-rep(1,q)%*%t(1:l))
  }else{
    w=w1
  }

  nk=colSums(z)
  nlx=colSums(v)
  Nklx=nk%*%t(nlx)
  nklx=t(z)%*%x%*%v

  nly=colSums(w)
  Nkly=nk%*%t(nly)
  nkly=t(z)%*%y%*%w
  #critere=lgamma(g*a)+lgamma(m*a)-(g+m)*lgamma(a)+m*g*(lgamma(2*b)-2*lgamma(b))-lgamma(n+g*a)-lgamma(d+m*a)+
  #sum(lgamma(a+nk))+sum(lgamma(a+nl))+
  #sum(lgamma(b+nkl)+lgamma(b+Nkl-nkl)-lgamma(2*b+Nkl))

  critere=lgamma(g*a)+lgamma(h*a)+lgamma(l*a)-(g+h+l)*lgamma(a)+g*h*(lgamma(b1x+b2x)-
                                                                       lgamma(b1x)*lgamma(b2x))+g*l*(lgamma(b1y+b2y)-lgamma(b1y)*lgamma(b2y))-lgamma(n+g*a)-lgamma(p+h*a)-lgamma(q+l*a)+
    sum(lgamma(a+nk))+sum(lgamma(a+nlx))+sum(lgamma(a+nly))+
    sum(sum(lgamma(b1x+nklx)+lgamma(b2x+Nklx-nklx)-lgamma(b1x+b2x+Nklx)))+
    sum(sum(lgamma(b1y+nkly)+lgamma(b1y+Nkly-nkly)-lgamma(b1y+b2y+Nkly)))
  critere}








