##' Class "BIKM1_MLBM_Binary"
##'
##' Class of object returned by the \code{\link{BIKM1_MLBM_Binary}} function.
##'
##' @section Slots: \describe{
##'
##' \item{\code{model_max}: }{The selected model by the procedure with free energy W, theta, conditional probabilities (s_ig, r_jh, t_kl), iter, empty_cluster, and the selected partitions z, v and w. }
##'
##' \item{\code{criterion_choice}: }{A character string corresponding to the chosen criterion used for model selection, which can be "ICL" or "BIC".}
##'
##' \item{\code{init_choice}: }{A character string corresponding to the chosen initialization strategy used for the procedure, which can be "random" or "Gibbs" or "smallVBayes".}
##'
##' \item{\code{criterion_tab}: }{The matrix corresponding to the values of the chosen criterion for pairs of numbers of clusters visited by the BIKM1_MLBM_Binary function. The matrix rows design the numbers of row clusters. If a pair is not visited, by default, the value is -Inf.}
##'
##'
##' \item{\code{W_tab}: }{The matrix corresponding to the values of the free energy (minimizer of the loglikelihood in the algorithm) for pairs of numbers of clusters visited by the procedure. The matrix rows design the numbers of row clusters. If a pair is not visited, by default, the value is -Inf.}
##'
##'
##' \item{\code{criterion_max}: }{Numeric indicating the maximum of the criterion values, calculated on the pairs of numbers of clusters visited by the BIKM1_MLBM_Binary function.}
##'
##'
##' \item{\code{gopt}: }{An integer value indicating the number of row clusters selected by the BIKM1_MLBM_Binary function.}
##'
##'
##'  \item{\code{hopt}: }{An integer value indicating the number of column clusters for the first matrix selected by the BIKM1_MLBM_Binary function.}
##'
##'   \item{\code{lopt}: }{An integer value indicating the number of row clusters for the second matrix selected by the BIKM1_MLBM_Binary function.}
##' }

##'
##'
##'
##' @aliases BIKM1_MLBM_Binary-class
##'
##' @docType class
##'
##' @keywords class
##'
##'
##'
##' @rdname BIKM1_MLBM_Binary-class
##'
##' @exportClass BIKM1_MLBM_Binary
##'
##' @examples
##' \donttest{
##' require(bikm1)
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
##' res=BIKM1_MLBM_Binary(data$x,data$y,3,3,3,4,init_choice='smallVBayes')}

setClass(
  Class="BIKM1_MLBM_Binary",
  representation=representation(
    model_max="list",
    criterion_choice="character",
    init_choice="character",
    criterion_tab="array",
    W_tab="array",
    criterion_max="numeric",
    gopt="numeric",
    hopt="numeric",
    lopt="numeric")
)


##'
##' Print method for a BIKM1_MLBM_Binary object
##'
##' Print method for a \code{\linkS4class{BIKM1_MLBM_Binary}} object
##'
##'
##' @param x in the print method, a BIKM1_MLBM_Binary object
##'
##' @param ... in the print method, additional parameters (ignored)
##'
##' @export print
##' @aliases print,BIKM1_MLBM_Binary-method
##' @examples
##' \donttest{require(bikm1)
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
##' res=BIKM1_MLBM_Binary(data$x,data$y,3,3,3,4)
##' print(res)}
##'


setMethod("print","BIKM1_MLBM_Binary",
          function(x,...){
            result=list(x@criterion_choice,x@criterion_max,x@gopt,x@hopt,x@lopt)
            names(result)=c("criterion_choice","criterion_max","Selected number of row clusters (g)","Selected number of column clusters 1st matrix (h)","Selected number of column clusters 2nd matrix (l)")
            print(result)
          }
)


##'  Show method for a BIKM1_MLBM_Binary object
##'
##'  show method for a \code{\linkS4class{BIKM1_MLBM_Binary}} object

##'
##'
##'
##' @param object a BIKM1_MLBM_Binary object
##' @export show
##' @aliases show,BIKM1_MLBM_Binary-method
##' @examples
##' \donttest{require(bikm1)
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
##' res=BIKM1_MLBM_Binary(data$x,data$y,3,3,3,4)
##' show(res)}
##'


setMethod("show","BIKM1_MLBM_Binary",
          function(object){
            result=list(object@criterion_choice,object@criterion_max,object@gopt,object@hopt,object@lopt)
            names(result)=c("criterion_choice","criterion_max","Selected number of row clusters (g)","Selected number of column clusters 1st matrix (h)","Selected number of column clusters 2nd matrix (l)")
            print(result)
          }
)



##'
##' Summary method for a BIKM1_MLBM_Binary object
##'
##' Produce a summary of informations of a \code{BIKM1_MLBM_Binary} object
##'
##'
##'
##' @param object in the summary method, a BIKM1_MLBM_Binary object
##' @param ... in the summary method, additional parameters (ignored)
##' @export summary
##'
##' @aliases summary,BIKM1_MLBM_Binary-method
##'
##' @examples
##' \donttest{require(bikm1)
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
##' res=BIKM1_MLBM_Binary(data$x,data$y,3,3,3,4)
##' summary(res)}
##'
##'

setMethod("summary","BIKM1_MLBM_Binary",
          function(object,...){
            cat("\nCriterion Choice : ",as.character(object@criterion_choice),"\n",sep="")
            cat("\nInitialization Choice : ",as.character(object@init_choice),"\n",sep="")
            cat("\nMaximum Criterion value: ",as.character(object@criterion_max),"\n",sep="")
            cat("\nEstimated number of row clusters (g): ",as.character(object@gopt),"\n",sep="")
            cat("\nEstimated number of column clusters  1st matrix (h): ",as.character(object@hopt),"\n",sep="")
            cat("\nEstimated number of column clusters  2nd matrix (l): ",as.character(object@lopt),"\n",sep="")
            cat("\nEstimated row proportion:",sep="")
            print(object@model_max$theta$pi_g)
            cat("\nEstimated column proportion 1st matrix:",sep="")
            print(object@model_max$theta$rho_h)
            cat("\nEstimated column proportion 2nd matrix:",sep="")
            print(object@model_max$theta$tau_l)
            cat("\nEmpty cluster z : ",as.character(object@model_max$Empty$z),"\n",sep="")
            cat("\nEmpty cluster v : ",as.character(object@model_max$Empty$v),"\n",sep="")
            cat("\nEmpty cluster w : ",as.character(object@model_max$Empty$w),"\n",sep="")


          }
)
##'
##'
##'
##'
##' Plot method for a \code{\linkS4class{BIKM1_MLBM_Binary}} object
##'
##'
##' Produce respectively a plot of two-dimensional segmentation of a \code{BIKM1_MLBM_Binary} fit, and a boxplot of conditional posteriors for each row and column cluster.
##'
##' @param x an object of class \code{BIKM1_MLBM_Binary}.
##' @param y a list specifying :
##'
##' \code{x}: the first matrix of observations
##'
##' \code{y}: the second matrix of observations.
##'
##' @param ... in the plot method, additional parameters (ignored)
##'
##' @return Two \pkg{plot} and on \pkg{ggplot2} object.
##'
##' @export plot
##'
##' @aliases plot,BIKM1_MLBM_Binary-method
##' @examples
##' \donttest{require(bikm1)
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
##' res=BIKM1_MLBM_Binary(data$x,data$y,3,3,3,4)
##' plot(res,data)}




setMethod(
  f="plot",
  signature="BIKM1_MLBM_Binary",
  definition=function(x,y,...){

    oldpar=par(mfrow=c(1,2),oma=c(0,0,3,0))
    on.exit(par(oldpar))
    #par(layout(c(1,2,rep(c(3,4),5))))
    BinBlocVisu_MLBM(y$x,y$y,rep(1,dim(y$x)[1]),rep(1,dim(y$x)[2]),rep(1,dim(y$y)[2]))
    BinBlocVisu_MLBM(y$x,y$y,x@model_max$z,x@model_max$v,x@model_max$w)

    mtext('Initial partitions',side = 3,outer = TRUE,adj = 0.1)

    mtext('Estimated partitions',side = 3,outer = TRUE,adj = 0.9)
    #refgraph=list(ask=par()$ask,
                  # mgp=par()$mgp,
                  # oma=par()$oma,
                  # xaxs=par()$xaxs,
                  # mfrow=par()$mfrow,
                  # cex.axis=par()$cex.axis)

    # Niter=sum(sum(x@criterion_tab>-Inf))
    # pos=matrix(0,Niter,3)
    # g=2
    # m=2
    # s=2
    # pos[1,1]=x@criterion_tab[g,m]
    # pos[1,2]=g
    # pos[1,3]=m
    # for (iter in 2:Niter){
    #   if (x@criterion_tab[g+1,m]>-Inf){
    #     g=g+1
    #   }else{
    #     m=m+1
    #   }
    #   pos[iter,1]=x@criterion_tab[g,m]
    #   pos[iter,2]=g
    #   pos[iter,3]=m
    # }
    # dev.new(width=20)
    #
    #
    #
    # # dev.new(width=14)
    # # par(mfrow=c(1,2),oma=c(0,0,3,0))
    # # #par(layout(c(1,2,rep(c(3,4),5))))
    # # BinBlocVisuResum_MLBM(y,rep(1,dim(y$x)[1]),rep(1,dim(y$x)[2]),rep(1,dim(y$y)[2]))
    # # mtext('Reorganized data matrix with initial partitions',side = 3,outer = TRUE,adj = 0.1)
    # # BinBlocVisuResum_MLBM(y,x@model_max$z,x@model_max$v,x@model_max$w)
    # # mtext('Reorganized data matrix with estimated partitions',side = 3,outer = TRUE,adj = 0.9)
    # # #
    #
    # #ev.new(width=14)
    # #library(grid)
    # grid.newpage()
    # pushViewport(viewport(layout = grid.layout(1,2)))
    # dat=data.frame(pos)
    # names(dat)=c("Criterion","Row","Column")
    # Criterion='Criterion'
    # Row='Row'
    # Column='Column'
    # p1<-ggplot(dat,aes(x=Row,y=Criterion))+geom_line(linetype = "dashed")+
    #   geom_point()+
    #   labs(x="Number of row clusters",y="Criterion values")
    #
    # print(p1, vp = viewport(layout.pos.row = 1,
    #                         layout.pos.col = 1))
    # p2<-ggplot(dat,aes(x=Column,y=Criterion))+geom_line(linetype = "dashed")+
    #   geom_point()+
    #   labs(x="Number of column clusters",y="Criterion values")
    #
    # print(p2, vp = viewport(layout.pos.row = 1,
    #                         layout.pos.col = 2))


    # par(mfrow = c(1,2))
    # plot(pos[,2],pos[,1],type='b', xlab='Number of row clusters',ylab='Criterion values')
    #
    # plot(pos[,3],pos[,1],type='b', xlab='Number of column clusters',ylab='Criterion values')
    # title(main="Evolution of Criterion values",outer=TRUE,line=-2)


    dev.new(width=20)
    oldpar=par(mfrow=c(1,2),oma=c(3,3,3,3))
    on.exit(par(oldpar))
    di=dim(x@model_max$s_ig)[1]
    dit=dim(x@model_max$r_jh)[1]
    dis=dim(x@model_max$t_kl)[1]
    c=matrix(0,1,di)
    s=x@model_max$s_ig
    t=x@model_max$r_jh
    ty=x@model_max$t_kl
    a=matrix(apply(s,1,max),1,di)
    for (i in 1:di){
      c[1,i]=(which(s[i,]==max(s[i,])))

    }
    b=sort(c[1,],decreasing=TRUE,index.return=TRUE)
    ordonne= a[1,b$ix]
    uX = unique(b$x)

    n = histc(b$x,sort(uX))
    p=n$cnt



    p=cumsum(p)
    matbox_ligne=list()
    matbox_ligne[[1]]=ordonne[1:p[1]]
    g=dim(x@model_max$s_ig)[2]
    for (i in 1:g){
      if (i!=g){
        matbox_ligne[[i+1]]= ordonne[(p[i]+1):p[i+1]]
      }
    }

    at=matrix(apply(t,1,max),1,dit)
    ct=matrix(0,1,dit)
    for (i in 1:dit){
      ct[1,i]=(which(t[i,]==max(t[i,])))

    }
    bt=sort(ct[1,],decreasing=TRUE,index.return=TRUE)
    ordonnet= at[1,bt$ix]
    uXt = unique(bt$x)

    nt = histc(bt$x,sort(uXt))
    pt=nt$cnt



    pt=cumsum(pt)
    matbox_colonne=list()
    matbox_colonne[[1]]=ordonnet[1:pt[1]]
    m=dim(x@model_max$r_jh)[2]
    for (i in 1:m){
      if (i!=m){
        matbox_colonne[[i+1]]= ordonnet[(pt[i]+1):pt[i+1]]
      }
    }


    aty=matrix(apply(ty,1,max),1,dis)
    cty=matrix(0,1,dis)
    for (i in 1:dis){
      cty[1,i]=(which(ty[i,]==max(ty[i,])))

    }
    bty=sort(cty[1,],decreasing=TRUE,index.return=TRUE)
    ordonnety= aty[1,bty$ix]
    uXty = unique(bty$x)

    nty = histc(bty$x,sort(uXty))
    pty=nty$cnt



    pty=cumsum(pty)
    matbox_colonney=list()
    matbox_colonney[[1]]=ordonnety[1:pty[1]]
    ss=dim(x@model_max$t_kl)[2]
    for (i in 1:ss){
      if (i!=ss){
        matbox_colonney[[i+1]]= ordonnety[(pty[i]+1):pty[i+1]]
      }
    }



    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1,3)))

    dat=data.frame(melt(matbox_ligne))

    names(dat)=c("y","x")
    dat$x=as.factor(dat$x)
    p1<-ggplot(dat,aes(x=x,y=y,fill=x)) + geom_boxplot()+
      labs(x="row cluster",y="Conditionnal posterior",
           title="")+
      theme(legend.position='none')+ylim(0,1)
    print(p1, vp = viewport(layout.pos.row = 1,
                            layout.pos.col = 1))
    dat=data.frame(melt(matbox_colonne))
    names(dat)=c("y","x")
    dat$x=as.factor(dat$x)
    p2<-ggplot(dat,aes(x=x,y=y,fill=x)) + geom_boxplot()+
      labs(x="column cluster 1rst matrix",y="Conditionnal posterior",
           title=" ")+
      theme(legend.position='none')+ylim(0,1)
    print(p2, vp = viewport(layout.pos.row = 1,
                            layout.pos.col = 2))
    daty=data.frame(melt(matbox_colonney))
    names(daty)=c("y","x")
    daty$x=as.factor(daty$x)
    p3<-ggplot(daty,aes(x=x,y=y,fill=x)) + geom_boxplot()+
      labs(x="column cluster 2nd matrix",y="",
           title="Conditional posterior for each column cluster ")+
      theme(legend.position='none')+ylim(0,1)
    print(p3, vp = viewport(layout.pos.row = 1,
                            layout.pos.col = 3))



    #par(refgraph)
  }
)



