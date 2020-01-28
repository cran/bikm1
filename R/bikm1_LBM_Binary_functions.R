
##' BIKM1_LBM_Binary fitting procedure
##'
##' Produce a blockwise estimation of a contingency matrix of observations.
##'
##' @param x binary matrix of observations.
##' @param Gmax a positive integer less than number of rows.
##' @param Hmax a positive integer less than number of columns.The bikm1 procedure stops while the numbers of rows is higher than Gmax or the number of columns is higher than Hmax.
##' @param a hyperparameter used in the VBayes algorithm for priors on the mixing proportions.
##' By default, a=4.
##' @param b hyperparameter used in the VBayes algorithm for prior on the Bernoulli parameter.
##' By default, b=1.
##' @param Gstart a positive integer to initialize the procedure with number of row clusters.
##' By default, Gstart=2.
##' @param Hstart a positive integer to initialize the procedure with number of column clusters.
##' By default, Hstart=2.
##' @param init_choice a character string corresponding to the chosen initialization strategy used for the procedure, which can be "random" or "smallVBayes" or "user".
##' By default, init_choice="smallVBayes".
##' @param userparam in the case where init_choice is "user", a list containing partitions z and w. By default userparam=NULL.
##' @param ntry a positive integer corresponding to the number of times which is launched the small VBayes or random initialization strategy. By default ntry=100.
##' @param criterion_choice a character string corresponding to the chosen criterion used for model selection, which can be "ICL" as for now.
##' By default, criterion_choice="ICL".
##' @param mc.cores a positive integer corresponding to the available number of cores for parallel computing.
##' By default, mc.cores=1.
##' @param verbose logical. To display each step and the result. By default verbose=TRUE.
##' @rdname BIKM1_LBM_Binary-proc
##' @references Govaert and Nadif. Co-clustering, Wyley (2013).
##'
##' Keribin, Brault and Celeux. Estimation and Selection for the Latent Block Model on Categorical Data, Statistics and Computing (2014).
##'
##' Robert. Classification crois\'ee pour l'analyse de bases de donn\'ees de grandes dimensions de pharmacovigilance. Paris Saclay (2017).
##' @usage BIKM1_LBM_Binary(x,Gmax,Hmax,a=4,b=1,
##' Gstart=2,Hstart=2,init_choice='smallVBayes',userparam=NULL,
##' ntry=50,criterion_choice='ICL', mc.cores=1,verbose=TRUE)
##'
##' @return a BIKM1_LBM_Binary object including
##'
##' \code{model_max}: the selected model by the procedure with free energy W, theta, conditional probabilities (s_ig, r_jh), iter, empty_cluster, and the selected partitions z and w.
##'
##' \code{criterion_choice}: the chosen criterion
##'
##' \code{init_choice}: the chosen init choice
##'
##' \code{criterion tab}:  the matrix containing the criterion values for each selected number of row and column
##'
##' \code{W_tab}: the matrix containing the free energy values for each selected number of row and column
##'
##' \code{criterion_max}: the maximum of the criterion values
##'
##' \code{gopt}: the selected number of rows
##'
##' \code{hopt}: the selected number of columns
##' @examples
##'
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' g=3
##' h=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##'theta$alpha_gh=matrix(runif(6),ncol=h)
##' data=BinBlocRnd_LBM(n,J,theta)
##' res=BIKM1_LBM_Binary(data$x,3,2,Gstart=3,Hstart=2,
##' init_choice='user',userparam=list(z=data$xrow,v=data$xcol))
##'
##' @export BIKM1_LBM_Binary
##'


BIKM1_LBM_Binary=function(x,Gmax,Hmax,a=4,b=1,Gstart=2,Hstart=2,init_choice='smallVBayes',userparam=NULL,ntry=50,criterion_choice='ICL', mc.cores=1,verbose=TRUE){
  #data=list()
  #data$x=x
  #x_ij=data$x
  BinBlocInit_LBM=function(g,h,x,start,fixed_proportion=FALSE){
    #x_ij=data$x
    eps=10^(-16)
    n=dim(x)[1]
    J=dim(x)[2]
    # Initialisation par tirage al?atoire des centres


    if(start$ale==TRUE){
      theta=list()
      #theta$pi_g=runif(g)
      #theta$pi_g=theta$pi_g/sum(theta$pi_g)
      #theta$rho_h=runif(m)
      #theta$rho_h=theta$rho_h/sum(theta$rho_h)
      #theta$rho_hy=runif(s)
      #theta$rho_hy=theta$rho_hy/sum(theta$rho_hy)

      theta$pi_g=1/g *rep(1,g)
      theta$rho_h=1/h *rep(1,h)

      resx=PartRnd(J,theta$rho_h)
      u_ilx=x%*%resx$z_ik
      alpha_gh=u_ilx[sample(1:n,g),]/(matrix(1,g,1)%*%colSums(resx$z_ik))
      #alpha_gh=matrix(runif(g*m),ncol=m)
      alpha_gh[alpha_gh<eps]=eps
      alpha_gh[alpha_gh>1-eps]=1-eps
      theta$alpha_gh=alpha_gh


      list(r_jh=resx$z_ik,theta=theta,v=resx$z)
    }else if (any(names(start)=='theta')&& any(names(start)=='v') )

        {theta=start$theta
        r_jh=matrix(mapply(as.numeric,(!(start$v%*%matrix(1,1,h)-matrix(1,J,1)%*%(1:h)))),ncol=h)

        list(r_jh=r_jh,theta=theta)

     }else if (any(names(start)=='z')&& any(names(start)=='theta')){

            theta=start$theta
            z=start$z
            # Computation of r_jh
            s_ig=matrix(mapply(as.numeric,(!(z%*%matrix(1,1,g)-matrix(1,n,1)%*%(1:g)))),ncol=g)
            s_k=t(colSums(s_ig))
            eps10=10*eps
            log_rho_h=log(theta$rho_h)
            alpha_gh=theta$alpha_gh

            A_klx=log((alpha_gh+eps10)/(1-alpha_gh+eps10))

            B_klx=log(1-alpha_gh+eps10)

            r_jh=(t(x)%*%s_ig)%*%A_klx+rep(1,J)%*%s_k%*%B_klx+matrix(rep(1,J),ncol=1)%*%t(log_rho_h)
            t_kHmax=apply(r_jh, 1,max)
            Tp_jlx = exp(r_jh-t_kHmax)
            r_jh=Tp_jlx/(rowSums(Tp_jlx))


            list(r_jh=r_jh,theta=theta)
          }else if  (any(names(start)=='z')&& any(names(start)=='v') ){
            if ((max(start$z)!=g)||(max(start$v)!=h)){
              warning('(Gstart,Hstart) need to match with the characteristics of userparam z and v')
            }else{
            s_ig=matrix(mapply(as.numeric,(!(start$z%*%matrix(1,1,g)-matrix(1,n,1)%*%(1:g)))),ncol=g)
            r_jh=matrix(mapply(as.numeric,(!(start$v%*%matrix(1,1,h)-matrix(1,J,1)%*%(1:h)))),ncol=h)

            # Computation of parameters
            s_k=colSums(s_ig)
            t_lx=colSums(r_jh)

            if (fixed_proportion==TRUE){
              theta$pi_g=1/g * rep(1,g)
              theta$rho_h=1/h * rep(1,h)

            }else{
              theta$pi_g=s_k/n
              theta$rho_h=t_lx/J
              }

            theta$alpha_gh = (t(s_ig)%*%x%*%r_jh)/(s_k%*%t(t_lx))

            list(r_jh=r_jh,theta=theta)}
          }else if (any(names(start)=='theta')&& any(names(start)=='r_jh') )
          {list(r_jh=r_jh,theta=theta) }
          else{warning('For a not random initialization, it needs at least 2 parameters')}
  }
  maxeps=function(i){
    eps=10^(-16)
    max(eps,min(i,1-eps))
  }
  BlocXemalpha_gh_LBM=function(x,theta,r_jh,a,b,niter=100000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
    alpha_gh=theta$alpha_gh
    pi_g=theta$pi_g
    rho_h=theta$rho_h
    g=dim(alpha_gh)[1]
    h=dim(alpha_gh)[2]
    n=dim(x)[1]
    J=dim(x)[2]
    t_l=t(colSums(r_jh))
    Ht= -sum(r_jh*log(r_jh+eps*(r_jh<=0)))

    # boucle sur les it?rations de l'algorithme
    W=-Inf
    for (iter in 1:niter){
      A_kl=log((alpha_gh*(alpha_gh>0)+(alpha_gh<=0)*eps)/((1-alpha_gh)*(alpha_gh<1)+(alpha_gh>=1)*eps));
      B_kl=log((1-alpha_gh)*(alpha_gh<1)+(alpha_gh>=0)*eps)
      log_pi_g=log(pi_g+(pi_g<=0)*eps)
      log_rho_h=log(rho_h+(rho_h<=0)*eps)

      W2=-Inf
      for (iter_int in 1:niter_int){
        # Computation of sig
        S_ig=(x%*%r_jh)%*%t(A_kl)+rep(1,n)%*%t_l%*%(t(B_kl))+matrix(rep(1,n),ncol=1)%*%t(log_pi_g);
        r_jGmax=apply(S_ig, 1,max)
        Sp_ik = exp(S_ig-r_jGmax)
        s_ig=Sp_ik/(rowSums(Sp_ik))
        Hs=-sum(sum(s_ig*(s_ig>0)*log(s_ig*(s_ig>0)+(s_ig<=0)*eps)))
        s_k=t(colSums(s_ig))

        # Computation of r_jh
        R_jh=(t(x)%*%s_ig)%*%A_kl+rep(1,J)%*%s_k%*%B_kl+matrix(rep(1,J),ncol=1)%*%t(log_rho_h)
        t_kHmax=apply(R_jh, 1,max)
        Tp_jl = exp(R_jh-t_kHmax)
        r_jh=Tp_jl/(rowSums(Tp_jl))
        Ht=-sum(sum(r_jh*(r_jh>0)*log(r_jh*(r_jh>0)+(r_jh<=0)*eps)))
        t_l=t(colSums(r_jh))
        empty_cluster=any(t_l<10.e-5)||any(s_k<10.e-5)
        W2_old=W2
        W2=sum(r_jh*r_jh) + s_k%*%log_pi_g  + Hs + Ht
        if (abs((W2-W2_old)/W2) < epsi_int) break
      }
      # Computation of parameters
      pi_g=t((s_k+a-1)/(n+g*(a-1)))
      rho_h=t((t_l+a-1)/(J+g*(a-1)))
      u_kl=t(s_ig)%*%x%*%r_jh
      n_kl=t(s_k)%*%t_l
      u_kl=u_kl*(u_kl>0)
      n_kl=n_kl*(n_kl>0)
      alpha_gh=(u_kl+b-1)/(n_kl+2*(b-1))
      if (any(is.nan(alpha_gh))){
        alpha_gh=(u_kl+b-1)/(n_kl*(n_kl>0)+2*(b-1)+(n_kl<=0)*eps)
      }
      # Criterion
      W_old=W
      W=sum(s_k*log(s_k))- n*log(n) + sum(t_l*log(t_l))-J*log(J)+
        sum(u_kl*log(u_kl)+(n_kl-u_kl)*log(n_kl-u_kl)-n_kl*log(n_kl))+
        Hs+Ht+
        (a-1)*(sum(log(pi_g))+sum(log(rho_h)))+
        (b-1)*sum(log(alpha_gh)+log(1-alpha_gh))

      if (is.nan(W)){
        if (b==1){
          W=sum(s_k*log(s_k+(s_k<=0)*eps))-n*log(n)+sum(t_l*log(t_l+(t_l<=0)*eps))-J*log(J)+
            sum(u_kl*log(u_kl+(u_kl<=0)*eps)+(n_kl-u_kl)*((n_kl-u_kl)>0)*log((n_kl-u_kl)*((n_kl-u_kl)>0)+(((n_kl-u_kl)<=0)*eps))-n_kl*log(n_kl+(n_kl<=0)*eps))+
            Hs + Ht+
            (a-1)*(sum(log(pi_g*(pi_g>0)+(pi_g<=0)*eps))+
                     sum(log(rho_h*(rho_h>0)+(rho_h<=0)*eps)))
        }else{
          W=sum(s_k*log(s_k+(s_k<=0)*eps))-n*log(n)+sum(t_l*log(t_l+(t_l<=0)*eps))-J*log(J)+
            sum(u_kl*log(u_kl+(u_kl<=0)*eps)+(n_kl-u_kl)*((n_kl-u_kl)>0)*log((n_kl-u_kl)*((n_kl-u_kl)>0)+(((n_kl-u_kl)<=0)*eps))-n_kl*log(n_kl+(n_kl<=0)*eps))+
            Hs + Ht+
            (a-1)*(sum(log(pi_g))+sum(log(rho_h)))+
            (b-1)*sum(log(alpha_gh*(alpha_gh>0)+(alpha_gh<=0)*eps)+
                        log((1-alpha_gh)*(alpha_gh<1)+(alpha_gh>=1)*eps))
        }
      }
      if (abs((W-W_old)/W) < epsi) break
    }
    theta$alpha_gh=alpha_gh
    theta$pi_g=pi_g
    theta$rho_h=rho_h
    z=apply(s_ig, 1, which.max)
    v=apply(r_jh, 1, which.max)
    list(W=W,theta=theta,s_ig=s_ig,r_jh=r_jh,iter=iter,empty_cluster=empty_cluster,z=z,v=v)
  }
  BinBlocVBayes_LBM=function(x,g,h,a,b,ini=list(n_init=1,ale=TRUE,user=FALSE),niter=100000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
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
        initVBayes=BinBlocInit_LBM(g,h,x,start=ini)
        ResVBayesAle=BlocXemalpha_gh_LBM(x,initVBayes$theta,initVBayes$r_jh,a=a,b=b,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)

        ### ? enlever
        #print(c("ResVBayesAle$W= ", ResVBayesAle$W))
        W2=ResVBayesAle$W
        theta_max=ResVBayesAle$theta
        s_ig_max=ResVBayesAle$s_ig
        r_jh_max=ResVBayesAle$r_jh
        ### fin ? enlever
        if ((ResVBayesAle$W>W2)||(is.nan(ResVBayesAle$W))||(i=ini$n_init)||(W2=-Inf)){
          W2=ResVBayesAle$W
          theta_max=ResVBayesAle$theta
          s_ig_max=ResVBayesAle$s_ig
          r_jh_max=ResVBayesAle$r_jh
        }
      }# fin iteration

    theta_path=0
    theta=theta_max
    s_ig=s_ig_max
    r_jh=r_jh_max

    #matrice de classification
    z=apply(s_ig, 1, which.max)
    v=apply(r_jh, 1, which.max)

    #calcul des marges lignes et colonnes
    Repart=list()
    Repart$z=matrix(0,dim(s_ig)[2],1)

    for (ig in 1:dim(s_ig)[2]){
      Repart$z[ig]=sum(z==ig)
    }
    Empty=list()
    Empty$z=sum(Repart$z==0)
    Repart$v=matrix(0,dim(r_jh)[2],1)
    for (im in 1:dim(r_jh)[2]){
      Repart$v[im]=sum(v==im)
    }

    Empty$v=sum(Repart$v==0)


    list(W=W2,theta=theta,s_ig=s_ig,r_jh=r_jh,z=z,v=v,Empty=Empty,Repart=Repart,iter=ResVBayesAle$iter,normtheta=ResVBayesAle$normtheta)
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



  BinBlocICL =function (a,b,x,z1,v1){


    n=dim(x)[1]
    J=dim(x)[2]
    if (!is.matrix(z1)){
      g=max(z1)
      z=!(z1%*%t(rep(1,g))-rep(1,n)%*%t(1:g))
    }else{
      z=z1
    }
    if (!is.matrix(v1)){
      h=max(v1)
      v=!(v1%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
    }else{
      v=v1
    }

    nk=colSums(z)
    nl=colSums(v)
    Nkl=nk%*%t(nl)
    nkl=t(z)%*%x%*%v
    critere=lgamma(g*a)+lgamma(h*a)-(g+h)*lgamma(a)+h*g*(lgamma(2*b)-2*lgamma(b))-lgamma(n+g*a)-lgamma(J+h*a)+
      sum(lgamma(a+nk))+sum(lgamma(a+nl))+
      sum(lgamma(b+nkl)+lgamma(b+Nkl-nkl)-lgamma(2*b+Nkl))
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

  criterion_tab=matrix(-Inf,Gmax+1,Hmax+1)
  W_tab=matrix(-Inf,Gmax+1,Hmax+1)
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
      res=BinBlocVBayes_LBM(x,Gstart,Hstart,a,b
                            ,ini=list(n_init=1,ale=TRUE),
                            niter=50)
      if (res$W>W && res$Empty$v==FALSE && res$Empty$z==FALSE){
        W=res$W
        resopt=res
      }
    }
      if (!exists("resopt")){
        stop(paste('The algorithm does not manage to find a configuration with  (Gstart=',as.character(Gstart),',Hstart=',as.character(Hstart),')' ))
      }

  # }else if (init_choice=='Gibbs'){
  #
  #   resinit=PoissonBlocGibbs(data,Gstart,Hstart,a,normalization,alpha,beta,niter=5000,inistart=list(n_init=10,ale=TRUE))
  #   initVBayes=list(z=resinit$z,w=resinit$w,ale=FALSE,n_init=1)
  #   resopt=PoissonBlocVBayes(data,Gstart,Hstart,a,alpha,beta,normalization,ini=initVBayes)
  #   if (resopt$empty_cluster==TRUE){
  #     stop(paste('The algorithm does not manage to find a configuration with  (Gstart=',as.character(Gstart),',Hstart=',as.character(Hstart),'),try normalization=FALSE') )
  #   }
  }else if (init_choice=='user'){
    if (is.null(userparam)){
      stop(paste('If you choose the user init choice, please fill the useparam field with partitions z and v' ))
    }else{
      initVBayes=list(z=userparam$z,v=userparam$v,ale=FALSE,n_init=1)
      resopt=BinBlocVBayes_LBM(x,Gstart,Hstart,a,b,ini=initVBayes)

    }

  }else if (init_choice=='random'){

    eps=10^(-16)
    n=dim(x)[1]
    J=dim(x)[2]
    # Random_Temp=function(Temp){
    #   theta=list()
    #   theta$pi_g=runif(Gstart)
    #   theta$pi_g=theta$pi_g/sum(theta$pi_g)
    #   theta$rho_h=runif(Hstart)
    #   theta$rho_h=theta$rho_h/sum(theta$rho_h)
    #   resw=PartRnd(d,theta$rho_h)
    #   r_jh=resw$z_ik
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
    #   n_kl=t((t(mu_i)%*%s_ig))%*%(nu_j%*%r_jh)
    #   theta$lambda_kl=(t(s_ig)%*%data$x%*%r_jh)/(n_kl)
    #   if (any(is.nan(theta$lambda_kl))){
    #     lambda0=mean(mean(data$x))
    #     theta$lambda_kl=theta$lambda_kl=matrix(rpois(g,lambda0),ncol=m)
    #   }
    #
    #   initVBayes=list(r_jh=r_jh,w=w,theta=theta,ale=FALSE,n_init=1)
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
      #theta$rho_hy=runif(s)
      #theta$rho_hy=theta$rho_hy/sum(theta$rho_hy)

      theta$pi_g=1/Gstart *rep(1,Gstart)
      theta$rho_h=1/Hstart *rep(1,Hstart)


      resx=PartRnd(J,theta$rho_h)

      u_ilx=x%*%resx$z_ik

      alpha_gh=u_ilx[sample(1:n,Gstart),]/(matrix(1,Gstart,1)%*%colSums(resx$z_ik))
      #alpha_gh=matrix(runif(g*m),ncol=m)
      alpha_gh[alpha_gh<eps]=eps
      alpha_gh[alpha_gh>1-eps]=1-eps
      theta$alpha_gh=alpha_gh

      #list(r_jh=resx$z_ik,r_jhy=resy$z_ik,theta=theta,wx=resx$z,wy=resy$z)

      initVBayes=list(r_jh=resx$z_ik,theta=theta,v=resx$z,ale=FALSE,n_init=1)
      #res=BlocGibbs(data$x,a,b,ini=list(n_init=1,ale=TRUE),1000)
      res=BinBlocVBayes_LBM(x,Gstart,Hstart,a,b,ini=initVBayes)
      if (res$W>W && res$Empty$z==FALSE && res$Empty$v==FALSE ){
        resopt=res}



    }
    if (!exists("resopt")){
      stop(paste('The algorithm does not manage to find a configuration with  (Gstart=',as.character(Gstart),',Hstart=',as.character(Hstart),')') )
    }

   }else {stop("This is not a initialization choice proposed by the procedure")}

  if (criterion_choice=='ICL'){
     Critere= BinBlocICL(a,b,x,resopt$z,resopt$v)
  # else if (criterion_choice=='BIC'){
  #   Critere= PoissonBlocBIC(a,alpha,beta,data,resopt,normalization)
  # }
  }else {stop('This is not a criterion choice proposed by the procedure')}
  criterion_tab[Gstart,Hstart]=Critere
  W_tab[Gstart,Hstart]=resopt$W

  modele=resopt
  criterion_max=Critere
  g=Gstart
  h=Hstart
  gopt=Gstart
  hopt=Hstart
  model_max=resopt

  while (g+1<=Gmax && h+1<=Hmax ){
    if ((length(unique(modele$z))!=g)||(length(unique(modele$v))!=h)){
      stop(paste('The algorithm does not manage to find a configuration with (g=',as.character(g),',h=',as.character(h),'), maybe reduce Gmax, Hmax or change init_choice ') )
    }

    W_colonnex=-Inf
    W_ligne=-Inf

    f1=function(Kc){
      v=modele$v
      Xc=which(v==Kc)
      lxc=length(Xc)
      if (lxc>1){
        idxc=sample(lxc,floor(lxc/2),replace=FALSE)
        indc=Xc[idxc]
        v[indc]=h+1
        initVBayes=list(z=modele$z,v=v,ale=FALSE,n_init=1)
        resultsVBayes=BinBlocVBayes_LBM(x,g,h+1,a,b,ini=initVBayes)
        iterverif=0
        while ((resultsVBayes$Empty$z==TRUE || resultsVBayes$Empty$v==TRUE ) && iterverif<20){
          iterverif=iterverif+1
          v=modele$v
          idxc=sample(lxc,floor(lxc/2),replace=FALSE)
          indc=Xc[idxc]
          v[indc]=h+1
          initVBayes=list(z=modele$z,v=v,ale=FALSE,n_init=1)
          resultsVBayes=BinBlocVBayes_LBM(x,g,h+1,a,b,ini=initVBayes)
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
        initVBayes=list(z=z,v=modele$v,ale=FALSE,n_init=1)
        resultsVBayes=BinBlocVBayes_LBM(x,g+1,h,a,b,ini=initVBayes)
        iterverif=0
        while ((resultsVBayes$Empty$v==TRUE || resultsVBayes$Empty$z==TRUE)  && iterverif<20){
          iterverif=iterverif+1
          z=modele$z
          idxl=sample(ll,floor(ll/2),replace=FALSE)
          indl=Xl[idxl]
          z[indl]=g+1
          initVBayes=list(z=z,v=modele$v,ale=FALSE,n_init=1)
          resultsVBayes=BinBlocVBayes_LBM(x,g+1,h,a,b,ini=initVBayes)
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



    if ((W_ligne<=-Inf)&(Empty_m) ){
      warning("The algorithm does not manage to find a configuration  (",as.character(g+1),',',as.character(h), ") nor (",as.character(g),',',as.character(h+1),')')
      new(Class="BIKM1_LBM_Binary",init_choice=init_choice,criterion_choice=criterion_choice,criterion_tab=criterion_tab,criterion_max=criterion_max,model_max=model_max,W_tab=W_tab,gopt=gopt,hopt=hopt)
    }



    if (W_ligne>W_colonnex ){
      modele=modele_ligne
      if (criterion_choice=='ICL'){
        Critere= BinBlocICL(a,b,x,modele$z,modele$v)
      # }else if (criterion_choice=='BIC'){
      #   Critere= PoissonBlocBIC(a,alpha,beta,data,modele,normalization)
      }
      criterion_tab[g+1,h]=Critere
      W_tab[g+1,h]=modele$W
      if (Critere>criterion_max){
        model_max=modele
        criterion_max=Critere
        gopt=g+1
        hopt=h

      }
      if (verbose){
        cat('(g,h)=(',as.character(g+1),',',as.character(h),')\n',sep = "")
      }
      g=g+1
    }else if (W_colonnex>W_ligne ) {
      modele=modele_colonnex
      if (criterion_choice=='ICL'){
        Critere= BinBlocICL(a,b,x,modele$z,modele$v)
      # }else if (criterion_choice=='BIC'){
      #   Critere= PoissonBlocBIC(a,alpha,beta,data,modele,normalization)
      }
      criterion_tab[g,h+1]=Critere
      W_tab[g,h+1]=modele$W
      if (Critere>criterion_max){
        model_max=modele
        criterion_max=Critere
        gopt=g
        hopt=h+1


      }
      if (verbose){
        cat('(g,h)=(',as.character(g),',',as.character(h+1),')\n',sep = "")
      }

      h=h+1
    }

}
  if (verbose){
    cat('The selected row (g) and column (h) clusters are (g,h)=(',as.character(gopt),',',as.character(hopt),')\n',sep = "")
  }
  new(Class="BIKM1_LBM_Binary",init_choice=init_choice,criterion_choice=criterion_choice,criterion_tab=criterion_tab,criterion_max=criterion_max,model_max=model_max,W_tab=W_tab,gopt=gopt,hopt=hopt)
}

##' BinBlocRnd_LBM function for binary data matrix simulation
##'
##' Produce a data matrix generated under the Binary Latent Block Model.
##'
##' @param n a positive integer specifying the number of expected rows.
##' @param J a positive integer specifying the number of expected columns.
##' @param theta a list specifying the model parameters:
##'
##' \code{pi_g}: a vector specifying the row mixing proportions.
##'
##' \code{rho_h}: a vector specifying the matrix column mixing proportions.
##'
##' \code{alpha_gh}: a matrix specifying the distribution parameter of the matrix.

##' @usage BinBlocRnd_LBM(n,J,theta)
##' @return a list including the arguments:
##'
##' \code{x}: simulated  data matrix.
##'
##' \code{xrow}: numeric vector specifying row partition.
##'
##' \code{xcol}: numeric vector specifying  column partition.
##'

##' @rdname BinBlocRnd_LBM-proc
##'
##' @examples
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' g=3
##' h=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##'theta$alpha_gh=matrix(runif(6),ncol=h)
##' data=BinBlocRnd_LBM(n,J,theta)
##'
##' @export BinBlocRnd_LBM


BinBlocRnd_LBM =function (n,J,theta){
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

  reswx=PartRnd(J,theta$rho_h)
  w=reswx$z
  w_jlx=reswx$z_ik



  g=dim(theta$alpha_gh)[1]
  h=dim(theta$alpha_gh)[2]

  x=matrix(rbinom(n*J,1,z_ik%*%theta$alpha_gh%*%t(w_jlx)),ncol=J)

  xrow=as.numeric(z_ik%*%c(1:g))
  xcol=as.numeric(w_jlx%*%c(1:h))

  list(x=x,xrow=xrow,xcol=xcol)


}

##' BinBlocVisu_LBM function for visualization of binary matrix datasets
##'
##' Produce a plot object representing the co-clustered data-sets.
##'
##' @param x data matrix of observations.
##'
##'
##' @param z a numeric vector specifying the class of rows.
##' @param v a numeric vector specifying the class of columns.
##' @usage BinBlocVisu_LBM(x,z,v)
##' @rdname BinBlocVisu_LBM-proc
##' @return a \pkg{plot} object
##' @examples
##'
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' g=3
##' h=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$alpha_gh=matrix(runif(6),ncol=h)
##' data=BinBlocRnd_LBM(n,J,theta)
##' BinBlocVisu_LBM(data$x,data$xrow,data$xcol)
##' @export BinBlocVisu_LBM

BinBlocVisu_LBM=function(x,z,v){
  ## Initialisation

  g=max(z)
  h=max(v)
  n=dim(x)[1]
  J=dim(x)[2]
  ii=c()
  jj=c()
  for (i in 1:g){
    ii=c(ii,which(z==i))
  }
  for (j in 1:h){
    jj=c(jj,which(v==j))
  }
  z_ik=((z%*%t(rep(1,g)))==(rep(1,n)%*%t(1:g)))
  w_jl=((v%*%t(rep(1,h)))==(rep(1,J)%*%t(1:h)))

  ## Affichage
  n_k=n-sort(cumsum(colSums(z_ik))-0.5,decreasing=TRUE)
  n_l=cumsum(colSums(w_jl))+0.5

  table.paint(x[ii,jj],clegend=0)

  for (i in n_k[2:g]){
    lines(c(0.5,(J+1)),c(i,i),col='blue',lwd=2)
  }
  for (i in n_l[1:(h-1)]){
    lines(c(i,i),c(0.5,(n+1)),col='blue',lwd=2)
  }

}

##' BinBlocVisuResum_LBM function  for visualization of binary matrix data-sets
##'
##' Produce a plot object representing the resumed co-clustered data-sets.
##'
##' @param x binary matrix of observations.
##'
##'
##' @param z a numeric vector specifying the class of rows.
##' @param v a numeric vector specifying the class of columns.
##' @rdname BinBlocVisuResum_LBM-proc
##' @usage BinBlocVisuResum_LBM(x,z,v)
##' @return a \pkg{plot} object.
##' @examples
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' g=3
##' h=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##'theta$alpha_gh=matrix(runif(6),ncol=h)
##' data=BinBlocRnd_LBM(n,J,theta)
##' BinBlocVisuResum_LBM(data$x,data$xrow,data$xcol)
##' @export BinBlocVisuResum_LBM


BinBlocVisuResum_LBM=function(x,z,v){
  ## Initialisation

  g=max(z)
  h=max(v)
  n=dim(x)[1]
  J=dim(x)[2]
  z_ik=((z%*%t(rep(1,g)))==(rep(1,n)%*%t(1:g)))
  w_jl=((v%*%t(rep(1,h)))==(rep(1,J)%*%t(1:h)))

  ## Affichage
  nk=colSums(z_ik)
  nl=colSums(w_jl)
  Nkl=nk%*%t(nl)
  nkl=t(z_ik)%*%x%*%w_jl

  table.paint(nkl/Nkl)

}


##' BinBlocICL_LBM function  for computation of the ICL criterion in the Binary LBM
##'
##' Produce a value of the ICL criterion in the Binary LBM.
##'
##' @param a an hyperparameter for priors on the mixing proportions. By default, a=4.
##' @param b an hyperparameter for prior on the Bernoulli parameter. By default, b=1.
##' @param x  contingency matrix of observations.
##' @param z1 a numeric vector specifying the class of rows.
##' @param v1 a numeric vector specifying the class of columns.
##' @rdname BinBlocICL_LBM-proc
##' @usage BinBlocICL_LBM(a,b,x,z1,v1)
##' @return a value of the ICL criterion.
##' @examples
##' require(bikm1)
##' set.seed(42)
##' n=200
##' J=120
##' g=3
##' h=2
##' theta=list()
##' theta$pi_g=1/g *matrix(1,g,1)
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$alpha_gh=matrix(runif(6),ncol=h)
##' data=BinBlocRnd_LBM(n,J,theta)
##' BinBlocICL_LBM(a=4,b=1,data$x, data$xrow,data$xcol)
##' @export BinBlocICL_LBM






BinBlocICL_LBM =function (a=4,b=1,x,z1,v1){

  #x=data$x
  n=dim(x)[1]
  J=dim(x)[2]
  if (!is.matrix(z1)){
    g=max(z1)
    z=!(z1%*%t(rep(1,g))-rep(1,n)%*%t(1:g))
  }else{
    z=z1
  }
  if (!is.matrix(v1)){
    h=max(v1)
    v=!(v1%*%t(rep(1,h))-rep(1,J)%*%t(1:h))
  }else{
    v=v1
  }

  nk=colSums(z)
  nl=colSums(v)
  Nkl=nk%*%t(nl)
  nkl=t(z)%*%x%*%v
  critere=lgamma(g*a)+lgamma(h*a)-(g+h)*lgamma(a)+h*g*(lgamma(2*b)-2*lgamma(b))-lgamma(n+g*a)-lgamma(J+h*a)+
    sum(lgamma(a+nk))+sum(lgamma(a+nl))+
    sum(lgamma(b+nkl)+lgamma(b+Nkl-nkl)-lgamma(2*b+Nkl))
  critere}












