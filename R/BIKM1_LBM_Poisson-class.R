##' Class "BIKM1_LBM_Poisson"
##'
##' Class of object returned by the \code{\link{BIKM1_LBM_Poisson}} function.
##'
##' @section Slots: \describe{
##'
##' \item{\code{model_max}: }{The selected model by the procedure with free energy W, theta, conditional probabilities (r_jh, t_kl), iter, empty_cluster, and the selected partitions v and w. }
##'
##' \item{\code{criterion_choice}: }{A character string corresponding to the chosen criterion used for model selection, which can be "ICL" or "BIC".}
##'
##' \item{\code{init_choice}: }{A character string corresponding to the chosen initialization strategy used for the procedure, which can be "random" or "Gibbs" or "smallVBayes".}
##'
##' \item{\code{criterion_tab}: }{The matrix corresponding to the values of the chosen criterion for pairs of numbers of clusters visited by the BIKM1_LBM_Poisson function. The matrix rows design the numbers of row clusters. If a pair is not visited, by default, the value is -Inf.}
##'
##'
##' \item{\code{W_tab}: }{The matrix corresponding to the values of the free energy (minimizer of the loglikelihood in the algorithm) for pairs of numbers of clusters visited by the procedure. The matrix rows design the numbers of row clusters. If a pair is not visited, by default, the value is -Inf.}
##'
##'
##' \item{\code{criterion_max}: }{Numeric indicating the maximum of the criterion values, calculated on the pairs of numbers of clusters visited by the BIKM1_LBM_Poisson function.}
##'
##'
##' \item{\code{lopt}: }{An Integer value indicating the number of row clusters selected by the BIKM1_LBM_Poisson function.}
##'
##'
##'  \item{\code{hopt}: }{An integer value indicating the number of column clusters selected by the BIKM1_LBM_Poisson function.}
##'
##'
##' }

##'
##'
##'
##' @aliases BIKM1_LBM_Poisson-class
##'
##' @docType class
##'
##' @keywords class
##'
##'
##'
##' @rdname BIKM1_LBM_Poisson-class
##'
##' @exportClass BIKM1_LBM_Poisson
##'


##' @examples
##'
##' require(bikm1)
##' set.seed(42)
##' J=200
##' K=120
##' h=3
##' l=2
##' theta=list()
##' theta$rho_h=1/h *matrix(1,h,1)
##' theta$tau_l=1/l *matrix(1,l,1)
##' theta$gamma_hl=matrix(floor(runif(h*l)*20+1),ncol=l)
##' data=PoissonBlocRnd(J,K,theta)
##' res=BIKM1_LBM_Poisson(data$x,3,3,4,init_choice='smallVBayes')

setClass(
  Class="BIKM1_LBM_Poisson",
  representation=representation(
    model_max="list",
    criterion_choice="character",
    init_choice="character",
    criterion_tab="matrix",
    W_tab="matrix",
    criterion_max="numeric",
    hopt="numeric",
    lopt="numeric")
)


##'
##' Print method for a BIKM1_LBM_Poisson object
##'
##' Print method for a \code{\linkS4class{BIKM1_LBM_Poisson}} object
##'
##'
##' @param x in the print method, a BIKM1_LBM_Poisson object
##'
##' @param ... in the print method, additional parameters (ignored)
##'
##' @export print
##' @aliases print,BIKM1_LBM_Poisson-method
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
##' res=BIKM1_LBM_Poisson(data$x,3,2,4,init_choice='random')
##' print(res)}
##'


setMethod("print","BIKM1_LBM_Poisson",
          function(x,...){
            result=list(x@criterion_choice,x@criterion_max,x@hopt,x@lopt)
            names(result)=c("criterion_choice","criterion_max","Selected number of row clusters (h)","Selected number of column clusters (l)")
            print(result)
          }
)


##'  Show method for a BIKM1_LBM_Poisson object
##'
##'  show method for a \code{\linkS4class{BIKM1_LBM_Poisson}} object

##'
##'
##'
##' @param object a BIKM1_LBM_Poisson object
##' @export show
##' @aliases show,BIKM1_LBM_Poisson-method
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
##' show(res)}
##'


setMethod("show","BIKM1_LBM_Poisson",
          function(object){
            result=list(object@criterion_choice,object@criterion_max,object@hopt,object@lopt)
            names(result)=c("criterion_choice","criterion_max","Selected number of row clusters (h)","Selected number of column clusters (l)")
            print(result)
          }
)



##'
##' Summary method for a BIKM1_LBM_Poisson object
##'
##' Produce a summary of informations of a \code{BIKM1_LBM_Poisson} object
##'
##'
##'
##' @param object in the summary method, a BIKM1_LBM_Poisson object
##' @param ... in the summary method, additional parameters (ignored)
##' @export summary
##'
##' @aliases summary,BIKM1_LBM_Poisson-method
##'
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
##' summary(res)}
##'
##'

setMethod("summary","BIKM1_LBM_Poisson",
          function(object,...){
            cat("\nCriterion Choice : ",as.character(object@criterion_choice),"\n",sep="")
            cat("\nInitialization Choice : ",as.character(object@init_choice),"\n",sep="")
            cat("\nMaximum Criterion value: ",as.character(object@criterion_max),"\n",sep="")
            cat("\nEstimated number of row clusters (h): ",as.character(object@hopt),"\n",sep="")
            cat("\nEstimated number of column clusters (l): ",as.character(object@lopt),"\n",sep="")
            cat("\nEstimated row proportion:",sep="")
            print(object@model_max$theta$rho_h)
            cat("\nEstimated column proportions:",sep="")
            print(object@model_max$theta$tau_l)
            cat("\nEstimated Poisson parameters:\n",sep="")
            print(object@model_max$theta$gamma_hl)
            cat("\nEmpty cluster: ",as.character(object@model_max$empty_cluster),"\n",sep="")


          }
)
##'
##'
##'
##'
##' Plot method for a \code{\linkS4class{BIKM1_LBM_Poisson}} object
##'
##'
##' Produce respectively one plot of two-dimensional segmentation of a \code{BIKM1_LBM_Poisson} fit, an evolution of the criterion as a function of the numbers of rows and columns, and a boxplot of conditional posteriors for each row and column cluster.
##'
##' @param x an object of class \code{BIKM1_LBM_Poisson}.
##' @param y a list specifying
##'
##' \code{x} : contingency matrix of observations.
##'
##' @param ... in the plot method, additional parameters (ignored)
##'
##' @return Two \pkg{plot} and two \pkg{ggplot2} object.
##'
##' @export plot
##'
##' @aliases plot,BIKM1_LBM_Poisson-method
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
##' res=BIKM1_LBM_Poisson(data$x,3,3,4,init_choice='random')
##' plot(res,data)}




setMethod(
  f="plot",
  signature="BIKM1_LBM_Poisson",
  definition=function(x,y,...){

    oldpar=par(mfrow=c(1,2),oma=c(1,1,1,1))
    on.exit(par(oldpar))
    PoissonBlocVisu(y$x,rep(1,dim(y$x)[1]),rep(1,dim(y$x)[2]))
    PoissonBlocVisu(y$x,x@model_max$v,x@model_max$w)
    mtext('Initial partitions',side = 3,outer = TRUE,adj = 0.1,cex.lab=0.3,
          cex.main=0.3)

    mtext('Estimated partitions',side = 3,outer = TRUE,adj = 0.9,cex.lab=0.3,
          cex.main=0.3)
    #refgraph=list(ask=par()$ask,
                  # mgp=par()$mgp,
                  # oma=par()$oma,
                  # xaxs=par()$xaxs,
                  # mfrow=par()$mfrow,
                  # cex.axis=par()$cex.axis)

    Niter=sum(sum(x@criterion_tab>-Inf))
    pos=matrix(0,Niter,3)
    h=2
    l=2
    pos[1,1]=x@criterion_tab[h,l]
    pos[1,2]=h
    pos[1,3]=l
    if (Niter>1){
    for (iter in 2:Niter){
      if (x@criterion_tab[h+1,l]>-Inf){
        h=h+1
      }else{
        l=l+1
      }
      pos[iter,1]=x@criterion_tab[h,l]
      pos[iter,2]=h
      pos[iter,3]=l
    }
    }
    #dev.new(width=30)
    #oldpar<-par(mfrow=c(1,2),oma=c(2,2,2,2))
    #on.exit(par(oldpar))
    #par(layout(c(1,2,rep(c(3,4),5))))






    # par(mfrow=c(1,2),oma=c(0,0,3,0))
    # #par(layout(c(1,2,rep(c(3,4),5))))
    # PoissonBlocVisuResum(y,rep(1,dim(y$x)[1]),rep(1,dim(y$x)[2]))
    # mtext('Reorganized data matrix with initial partitions',side = 3,outer = TRUE,adj = 0.1,cex.lab=0.3,
    #       cex.main=0.3)
    # PoissonBlocVisuResum(y,x@model_max$z,x@model_max$w)
    # mtext('Reorganized data matrix with estimated partitions',side = 3,outer = TRUE,adj = 0.9,cex.lab=0.3,
    #       cex.main=0.3)


    dev.new(width=14)
    if (Niter>1){
    #library(grid)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1,2)))
    dat=data.frame(pos)
    names(dat)=c("Criterion","Row","Column")
    Criterion='Criterion'
    Row='Row'
    Column='Column'
    p1<-ggplot(dat,aes(x=Row,y=Criterion))+geom_line(linetype = "dashed")+
      geom_point()+
      labs(x="Number of row clusters",y="Criterion values")+
      theme(plot.title = element_text(size=12),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12))

    print(p1, vp = viewport(layout.pos.row = 1,
                                    layout.pos.col = 1))
    p2<-ggplot(dat,aes(x=Column,y=Criterion))+geom_line(linetype = "dashed")+
      geom_point()+
      labs(x="Number of column clusters",y="Criterion values")+
      theme(plot.title = element_text(size=12),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12))

    print(p2, vp = viewport(layout.pos.row = 1,
                            layout.pos.col = 2))

    }


    # par(mfrow = c(1,2))
    # plot(pos[,2],pos[,1],type='b', xlab='Number of row clusters',ylab='Criterion values')
    # #
    # plot(pos[,3],pos[,1],type='b', xlab='Number of column clusters',ylab='Criterion values')
    # title(main="Evolution of Criterion values",outer=TRUE,line=-2)



    di=dim(x@model_max$r_jh)[1]
    dit=dim(x@model_max$t_kl)[1]
    c=matrix(0,1,di)
    s=x@model_max$r_jh
    t=x@model_max$t_kl
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
    g=dim(x@model_max$r_jh)[2]
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
    m=dim(x@model_max$t_kl)[2]
    for (i in 1:m){
      if (i!=m){
        matbox_colonne[[i+1]]= ordonnet[(pt[i]+1):pt[i+1]]
      }
    }

    dev.new(width=14)
    grid.newpage()

    pushViewport(viewport(layout = grid.layout(1,2)))

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
      labs(x="column cluster",y="Conditionnal posterior",
           title=" ")+
      theme(legend.position='none')+ylim(0,1)
    #grid.arrange(p1,p2, ncol=2, nrow = 1)
    print(p2, vp = viewport(layout.pos.row = 1,
                            layout.pos.col = 2))


    #par(refgraph)

  }
)



