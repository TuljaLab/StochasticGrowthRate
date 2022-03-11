### Before running on Mac, you need to have Xcode installed, with command-line tools
# From Xcode preferences go to Downloads, and click on Install by command-line tools
# Then you need to download gfortran from http://cran.r-project.org/bin/macosx/tools/


#k=# matrices
# X=list of matrices
# d=dimension of matrices
# n=order of term to compute



library(lpSolve)
library(expm)
library(inline)
library(Rcpp)
library(RcppArmadillo)
#library(rbenchmark)
#setwd('new_stoch_growth')

sourceCpp('senstoch.cpp')

unitvec=function(i,n){
    c(rep(0,i-1),1,rep(0,n-i))
}

### Input nonneg matrix m, output with rows normed to sum to 1
rownormal=function(m){
    msum=apply(m,1,function(v) notzero(sum(v)))
    m/matrix(msum,dim(m)[1],dim(m)[2])
}

### Input nonneg matrix m, output with columns normed to sum to 1
colnormal=function(m){
    msum=apply(m,2,function(v) notzero(sum(v)))
    m/matrix(msum,dim(m)[1],dim(m)[2],byrow=TRUE)
}

# Eliminate 0/0 problem
notzero=function(x){
    sapply(x,function(xx) xx+(xx==0))
}




### Start with a vector giving an MC realisation in reverse order, preceded
###   

## Generate iid sequence
iidseq=function(p,n){
    vs=cumsum(p)
    r=runif(n)
    vapply(r,function(rr) sum(rr>vs)+1,1)
}





# source('c++ code/mcovorsrc.R')

# PP=list(matrix(c(.5,.3,.5,.7),2,2),matrix(c(.9,.1,.1,.9),2,2))

# I=inhomseq(PP,c(.5,.5),20)
# Generating sequence of matrices

# W=matrix(rnorm(10),5,2)
# 
# arraycovariance(W)
# 
# Prev=reverseMC(P1)

#benchmark(out <- SimulateU(unlist(X),Prev,7,8,10,20,6,.01),out <- SimulateU1(unlist(X),Prev,7,8,10,20,6,.01), replications = 5,columns = c("test", "elapsed", "replications"))

##############################################################
# GSens: Growth and sensitivity object
# 
#       Input all data
#       Output complete collection of calls to SGmean, sensitivity1, sensitivity 2
# 
# Inputs:
#   X=list of projection matrices
#   P=transition matrix or vector of iid probabilities
#   N=vector of number of samples to take, or single total number of samples
#   T=initial number of steps for simulating U or V
#   maxerror= tolerance for systematic error in U or V
#   maxbias= tolerance for systematic error in computing zeta
#   whichside = 0 if using V estimates, anything else for U estimates
#
##############################################################

check_inputs=function(mats,P){
    
}

GSens=function(mats,P,N,NU,NV=NU,TV=50,basestate=0,maxerror=1e-6,maxbias=1e-6,whichside=1){
    SGm=SGmean(mats,P,NU,T,maxerror,whichside)
    s1=sensitivity1(mats,P,N,NU,NV,T,maxerror)
    s2=sensitivity2(mats,P,N,NU,NV,T,basestate,maxerror,maxbias)
    list(Growth=SGm,Sensitivity1=s1,Sensitivity2=s2)
}

check_matrices <- function(proj_matrices, transitions){
    if (is.vector(transitions)){
        K <- length(transitions)
        type <- 'iid'
        nu <- transitions
        transition_matrix <- matrix(transitions,K,K,byrow = TRUE)
    }
    else {
        transition_matrix <- transitions
        K <- dim(transition_matrix)[1]
        if ( !dim(transition_matrix)[2]==K ){stop("Transition matrix not square")}
        nu <- statdist(transition_matrix)
        if(!verifystochastic(transition_matrix)){stop("Transition matrix not stochastic")}
        type <- 'Markov'
    }
    if (sum(sapply(proj_matrices,function(x) abs(dim(x)[1]-dim(x)[2])))>0){stop("Projection matrices not all square")}
    if (length(unique(c(sapply(proj_matrices,dim))))>1) {stop("Projection matrices not all the same size")}
    list(K = K, nu = nu, TM = transition_matrix, type =type)
}

##############################################################
# SGmean: Stochastic growth estimation
# 
#       First estimate the variances using 100 samples for each starting state.
#       Use those estimates to decide an approximately optimal distribution of the remaining
#           samples for minimising the variance of the estimate.
# 
# Inputs:
#   proj_matrices=list of projection matrices
#   transitions=transition matrix or vector of iid probabilities
#   N=vector of number of samples to take, or single total number of samples
#   T=initial number of steps for simulating U or V
#   maxerror= tolerance for systematic error in U or V
#   whichside = 0 if using V estimates, anything else for U estimates
#
##############################################################
SGmean=function(proj_matrices,transitions,N,T=50,maxerror=1e-6,whichside=1){
    checkM <- check_matrices(proj_matrices = proj_matrices , transitions = transitions )
        # Check dimensions of matrices, etc., return appropriate transition matrix, K, etc
    K <- checkM$K
    transition_matrix <- checkM$TM
    nu <- checkM$nu
    N2Nvec <- N_to_Nvec(N,K,proj_matrices, transition_matrix, nu, T,maxerror,whichside)
    # Calculate reasonable number of samples for each state;
    #   Return vector of Nsamples, and the samples themselves (so they can be folded into the estimate) 
    Nvec <- N2Nvec[[1]]
    firstsg <- N2Nvec[[2]]
    newsg=stochgrowthall(unlist(proj_matrices),transition_matrix,dim(proj_matrices[[1]])[1],Nvec,T,maxerror,whichside)
    if (length(N)==1){#Might as well not throw away those extra samples
        NV0=pmax(ceiling(N/K/10),100)
        newsg[1,]=(NV0*firstsg[1,]+Nvec*newsg[1,])/(NV0+Nvec)
        newsg[2,]=(NV0*NV0*firstsg[2,]+Nvec*Nvec*newsg[2,])/(NV0+Nvec)^2
        newsg[3,]=(NV0*firstsg[3,]+Nvec*newsg[3,])/(NV0+Nvec)
    }
    nu=matrix(nu,1,K)
    results=apply(newsg*rbind(nu,nu*nu,nu),1,sum)
    names(results)=c('Estimate','Variance','Bias bound')
    results
}

# Function to break up the stated number of U samples into a sample size for each state
N_to_Nvec <- function(N,K,proj_matrices,transition_matrix, nu, T,maxerror,whichside=1){
    if (!length(N)==K & length(N)>1){stop("Wrong number of sample sizes")}
    if (length(N)==1){
        NV0=pmax(ceiling(N/K/10),100) # First step, sample equal (small) numbers for each state to get a crude estimate of SE
        firstsg=stochgrowthall(unlist(proj_matrices),transition_matrix,dim(proj_matrices[[1]])[1],rep(NV0,K),T,maxerror,whichside)
        SEestimates=firstsg[2,]
        Nvec=pmax(ceiling(.9*nu*sqrt(SEestimates)*N/sum(nu*sqrt(SEestimates))),25)
    }
    else {Nvec <- N}
    list(Nvec,firstsg)
}
##########################################################################################
# SGr: Wrapper function for taking inputs and submitting to C++ function
#     stochgrowthall to obtain stochastic growth rates.
# 
# Input 1) list of matrices X
#       2) Transition matrix P
#       3) Number of samples to take; Could be either a single number (to be later divided up)
#         or a vector of length K, saying how many samples to be taken for each starting state.
# 
# Output: Estimate of stochastic growth rate based on samples from U or V respectively.
#       
#
##########################################################################################

SGr=function(X,P,N,T=50,maxerror=1e-6,whichside=1){
    K=dim(P)[1]
    if (!length(N)==K & length(N)>1){stop("Wrong number of sample sizes")}
    if (!dim(P)[1]==dim(P)[2]){stop("Transition matrix not square")}
    if (sum(sapply(X,function(x) abs(dim(x)[1]-dim(x)[2])))>0){stop("Projection matrices not all square")}
    if (length(N)==1){
        nu=statdist(P)
        Nvec=pmax(ceiling(nu*N),25)
    }
    else {Nvec=N}
    stochgrowthall(unlist(X),P,dim(X[[1]])[1],Nvec,T,maxerror,whichside)
}

################################################################################################
#   sensitivity1
#
#       Compute sensitivity of stochastic growth wrt changes in matrix elements
#
#    Input: 
#             mats=list of projection matrices
#             P= transition matrix
#             NU,NV= number of samples in computing samples of U,V
#             TV= depth of first try at computing V or U
#             maxerror= limit on systematic error used in sampling U and V
#
#       Variables: K=# matrices
#                   d= dimensions of matrices
#
#        Output: List contatining K lists of 3 dxd matrices (of beta mean, variance, bias)
#               (i,j) element corresponds to an increase the probability of transition i->j
#               The column corresponding to the basestate is by definition 0
#
################################################################################################

sensitivity1=function(proj_matrices,transitions,NU,NV=NU,TV=50,maxerror=1e-6){
    checkM <- check_matrices(proj_matrices, transitions)
    K <- checkM$K
    transition_matrix <- checkM$TM
    nu <- checkM$nu
    matsu=unlist(proj_matrices)
    d=dim(proj_matrices[[1]])[1]
    as1=allsens1(matsu,transition_matrix,d,NU,NV,TV,maxerror)
    lapply(seq(K),function(i) list('Mean'=nu[i]*as1[,seq(d),i],'Variance'=nu[i]*nu[i]*as1[,seq(d)+d,i],'Bias'=nu[i]*as1[,seq(d)+2*d,i]))
}

# Input a matrix, output TRUE or FALSE depending on whether it is a stochastic matrix
verifystochastic=function(p){
    psum=apply(p,1,sum)
    all(psum==1)&all(p>=0)
}


#############################################################################
#
#           sensitivity2: Compute the frame for sensitivity to changes in driving process
#
#       Input:
#             mats=list of projection matrices
#             P= transition matrix
#             N= number of samples in computing zeta (irrelevant in iid case)
#             NU,NV= number of samples in computing xi
#             TV= depth of first try at computing V or U
#             basestate= state used as pivot for zetas
#             maxerror= limit on systematic error used in sampling U and V
#             maxbias= limit on systematic error used in choosing depth for zeta
#
#        Output: List contatining 3 KxK matrices (of beta mean, variance, bias) 
#               (i,j) element corresponds to an increase the probability of transition i->j
#               The column corresponding to the basestate is by definition 0
#############################################################################

sensitivity2=function(proj_matrices,transitions,N,NU,NV=NU,TV=50,basestate=0,maxerror=.001,maxbias=.001){
    checkM <- check_matrices(proj_matrices , transitions)
    K <- checkM$K
    transition_matrix <- checkM$TM
    nu <- checkM$nu
    matsu <- unlist(proj_matrices)
    d <- dim(proj_matrices[[1]])[1]
    zetareturn <- allzetas(proj_matrices,transitions,N,TV,basestate,maxerror,maxbias)
        # Note that zetareturn will be all 0s if transitions is a vector
    if (checkM$type == 'Markov') {
        xis <- allxi(matsu,transition_matrix,d,NU,NV,TV,maxerror)
        numat <- matrix(nu,K,K)
        betamean <- numat*(matrix(zetareturn[,1],K,K,byrow=TRUE)+xis[,seq(K)])
        betavar <- numat*numat*(matrix(zetareturn[,2],K,K,byrow=TRUE)+xis[,K+seq(K)])
        betabias <- numat*(matrix(zetareturn[,3],K,K,byrow=TRUE)+xis[,2*K+seq(K)])
        # Note: These coefficients are defined only up to translating a whole row
        #       by a constant. So we choose to make one element in each row 0,
        #       and the rest positive.
        for (i in seq(K)){
            wm=which.min(betamean[i,])
            wv=betavar[i,wm]
            wb=betabias[i,wm]
            betamean[i,]=betamean[i,]-min(betamean[i,])
            betavar[i,]=betavar[i,]+wv
            betavar[i,wm]=0
            betabias[i,]=betabias[i,]+wb
            betabias[i,wm]=0
        }
    }
    else{ # iid case. 
        xis <- allxi_iid(unlist(proj_matrices),transitions,dim(proj_matrices[[1]])[1],NU,NV,TV,maxerror)
        betamean <- xis[seq(K)]
        betavar <- xis[seq(K)+K]
        betabias <- xis[2*K+1]
        wm=which.min(betamean)
        wv=betavar[wm]
        wb=betabias[wm]
        betamean=betamean-betamean[wm]
        betavar=betavar+wv
        betavar[wm]=0
        betabias=betabias+wb
        betabias[wm]=0
    }
    list("mean"=betamean,"variance"=betavar,"bias"=betabias)
}

################################################################################################
#               zeta
# 
# Wrapper to partition the number of samples and call zetaintern twice.
# 
#         N=# samples
#         starts=pair of states being compared
#         Note: Zeta takes a transition matrix P. The function allzetas
#               splits off the iid case.
# 
################################################################################################
zeta=function(mats,P,N,starts,TV=50,maxerr=.001,maxb=.001,maxstep=1){
    K=length(mats)
    if (!dim(P)[1]==dim(P)[2]){stop("Transition matrix not square")}
    if (sum(sapply(mats,function(x) abs(dim(x)[1]-dim(x)[2])))>0){stop("Projection matrices not all square")}
    matsu=unlist(mats)
    d=dim(mats[[1]])[1]
    if (dim(P)[1]!=K || dim(P)[2]!=K || dim(mats[[1]])[2]!=d ){stop("Dimensions don't match.")}
    biasbound=maxb/maxstep*.9  # Conservative bound
    TY=QTchooseT(P,biasbound) # Choose depth of calculation to keep the bias below the desired level.
    transitionsQ=coupletransitions(P,TY,starts) # Vector that can be folded into two sequences of KxK matrices
    # of reverse transitions
    q <- matrix(transitionsQ[seq(1,K*TY)+2*K*K*(TY-1)],K,TY)
    ### q is the matrix of coupling probabilities
    ###     (which comes out tacked onto the end of transitionsQ)
    qResidue=max(1-sum(q),0) # Weight by stationary distribution, though in principle
    # the bias bound is independent of the distribution in the final state.
    q[,TY]=q[,TY]+qResidue*statdist(P)
    ## First step: Call zetaintern with all N's equal to 25, to get a baseline for estimating the error.
    allN=matrix(25,K,TY)
    firstzeta=zetaintern(matsu,transitionsQ[seq(1,K*K*(TY-1))],transitionsQ[seq(1,K*K*(TY-1))+K*K*(TY-1)],P,d,allN,TV,starts,maxerr,maxb,maxstep)
    qsum=sum((q*sqrt(firstzeta$Variance))) 
    allN=pmax(ceiling(q*sqrt(firstzeta$Variance)*N/qsum)-25,0)
    allN=ifelse(allN>=2,allN,0)  # Don't want N=1, because the variance is then undefined for that component.
    ## Second step: Call zetaintern with N's that are the remaining size above 25 based on the preliminary estimate of variance.
    secondzeta=zetaintern(matsu,transitionsQ[seq(1,K*K*(TY-1))],transitionsQ[seq(1,K*K*(TY-1))+K*K*(TY-1)],P,d,allN,TV,starts,maxerr,maxb,maxstep)
    jointzetamean=sum(q*(25*firstzeta$Mean+allN*secondzeta$Mean)/(25+allN))
    jointzetabias=sum(q*(25*firstzeta$Bias+allN*secondzeta$Bias)/(25+allN)) + sum(firstzeta$Cutoff2)*qResidue
    # Note: There's another component of bias that we only can compute after we have all zetas together.
    jointzetavar=sum(q*q*(625*firstzeta$Variance+allN*allN*secondzeta$Variance)/(25+allN)/(25+allN))
    c(jointzetamean,jointzetavar,jointzetabias,qResidue) # Note, qResidue always the same, except for different
    # sample lengths, but we need to
    # return it somewhere, so that it can be used to compute last piece of bias
}

allzetas=function(proj_matrices,transitions,N,TV=50,basestate=0,maxerror=.001,maxbias=.001){
    checkM <- check_matrices(proj_matrices , transitions)
    K <- checkM$K
    transition_matrix <- checkM$TM
    nu <- checkM$nu
    matsu=unlist(proj_matrices)
    d=dim(proj_matrices[[1]])[1]
    if (checkM$type == 'Markov'){
        outputs=matrix(0,K,4)
        mstep=rep(0,K)
        # We need to estimate the maximum possible
        for (i in seq(K)){
            u=SimulateU(matsu,transitions,d,ceiling(N/K),30,i-1,maxerror); # //NV/K and 30 are guesses
            mstep[i]=maxU(c(proj_matrices[[i]]),u[seq(d),]); # We only send one
            #  							// matrix at a time to maxU to save on storing u.
        }
    stepmax=max(mstep)  # This should be the largest possible
    for (i in seq(0,K-1)){
        if (i==basestate){outputs[i+1,]=rep(0,4)}
        else{outputs[i+1,]=zeta(proj_matrices,transitions,N,c(i,basestate),TV,maxerror,maxbias,stepmax)}
            #submit to zeta for every pair (basestate, something else); (basestate,basestate) returns 0
        }
        outputs[,3]=outputs[,3]+outputs[,4]*(max(outputs[,1])-min(outputs[,1]))
        return(outputs[,1:3])
    }
    else{
        matrix(0,K,3)
    }
}

#break

#####################################################################################
#
#                               Examples
#
#####################################################################################


testmats=list(matrix(c(0,.8,2,.3),2,2),matrix(c(.5,.7,3,.6),2,2),matrix(c(0,.5,4,.5),2,2))
testu=unlist(testmats)
testP=matrix(c(.4,.4,.2,.3,.4,.3,.5,0,.5),3,3,byrow=TRUE)
testP2=matrix(c(.4,.35,.25,.3,.4,.3,.5,0,.5),3,3,byrow=TRUE)
testP4=matrix(c(.4,.4,.2,.3,.4,.3,.48,0,.52),3,3,byrow=TRUE)
testP5=matrix(c(.4,.4,.2,.3,.4,.3,.5,.01,.49),3,3,byrow=TRUE)
testP6=matrix(c(.4,.4,.2,.3,.4,.3,.4,0,.6),3,3,byrow=TRUE)
testP7=matrix(c(.35,.45,.2,.3,.4,.3,.5,0,.5),3,3,byrow=TRUE)
testP8=matrix(c(.35,.4,.25,.3,.4,.3,.5,0,.5),3,3,byrow=TRUE)

testmats2=list(matrix(c(0,.6,2,.3),2,2),matrix(c(.5,.7,3,.6),2,2),matrix(c(0,.5,4,.5),2,2))
testmats3=list(matrix(c(0,.8,2.5,.3),2,2),matrix(c(.5,.7,3,.6),2,2),matrix(c(0,.5,4,.5),2,2))

SGmean(proj_matrices = testmats , transitions = testP , N = 10000)

# Horvitz matrices
# Sources http://www.jstor.org/stable/suppl/10.1086/378648/suppl_file/020306appxA.txt
# http://www.jstor.org/stable/suppl/10.1086/378648/suppl_file/020306appxB.txt


P0=matrix(c(.9525,.81,.6202,.4304,.2937,.157,.0785,.0475,.1425,.1898,.1898,.1367,.1367,.0785,0,.0475,.095,.1898,.1898,.1367,.1367,0,0,.095,.0475,.1898,.1898,.1367,0,0,0,.1425,.0475,.1898,.1898,0,0,0,0,.1425,.0475,.1898,0,0,0,0,0,.1425,.19),7,7)

#P001=matrix(c(7703,2298,0,0,0,0,0,))

P1=t(matrix(c(
    9190,190,190,137,137,78,78,6893,2298,190,190,136,136,157,0,6893,2298,189,189,137,294,0,0,6893,2297,190,190,430,0,0,0,4595,4595,190,620,0,0,0,0,2298,6892,810,0,0,0,0,0,2298,7702)/10000,7,7))

P1mod=t(matrix(c(
    9090,190,190,137,137,78,178,6893,2298,190,190,136,136,157,0,6893,2298,189,189,137,294,0,0,6893,2297,190,190,430,0,0,0,4595,4595,190,620,0,0,0,0,2298,6892,810,0,0,0,0,0,2298,7702)/10000,7,7))

#Moved .01 probability from (1,1) to (1,7)

Y=list(c(0,.1,0,0,0,0,0,0,0,0,.7,0,0,0,0,0,0,0,.95,.01,0,0,0,0,0,0,.17,.66,.17,0,0,0,.9,0,0,0,.96,.04,0,0,2.3,0,0,0,0,.88,0,0,2.6,0,0,0,0,.04,.96,0,.8,0,0,0,0,0,0,1),
       c(0,.1,0,0,0,0,0,0,1.8,0,.7,0,0,0,0,0,2.3,0,.76,.09,0,0,0,0,12.6,0,.05,.66,.29,0,0,0,73.4,0,.04,.06,.56,.29,.04,0,153.4,0,.02,0,.04,.67,.27,0,568.9,0,0,0,.07,0,.6,.33,1431.7,0,0,0,0,0,0,1),
       c(0,.1,0,0,0,0,0,0,.07,0,.07,0,0,0,0,0,5,0,.66,.25,.02,0,0,0,65.3,0,.04,.5,.33,.08,0,0,189.1,0,0,0,.52,.32,.03,.02,330.7,0,0,0,0,.76,.17,.03,679.7,0,0,0,0,.17,.5,.33,1142.8,0,0,0,0,0,.17,.83),
       c(0,.1,0,0,0,0,0,0,0,0,.7,0,0,0,0,0,1.8,0,.52,.19,.03,.03,0,0,14.3,0,.04,.32,.16,.36,0,0,18.1,0,0,.06,.42,.38,.02,0,61.2,0,0,0,.05,.43,.3,.16,125.2,0,0,0,.02,.1,.24,.56,179.4,0,.08,0,0,0,0,.92),
       c(0,.1,0,0,0,0,0,0,0.7,0,.7,0,0,0,0,0,1.3,0,.43,.22,0,0,0,0,62.1,0,0,.55,.27,.05,0,0,306.9,0,0,0,.43,.45,.06,0,579.7,0,0,0,0,.64,.31,.03,890.6,0,0,0,.1,.15,.5,.2,1843,0,0,0,0,.13,.13,.75),
       c(0,.1,0,0,0,0,0,0,0,0,.7,0,0,0,0,0,1.1,0,.7,.08,.01,0,0,0,58.6,0,.1,.49,.27,.12,0,0,190.3,0,0,.07,.5,.23,.07,.03,481.5,0,0,0,.03,.51,.33,.1,702,0,0,0,.03,.21,.38,.28,1508.1,0,0,0,0,0,0,.94),
       c(0,.1,0,0,0,0,0,0,0,0,.7,0,0,0,0,0,11.4,0,.57,.12,.01,0,0,0,110,.0,.2,.4,.4,0,0,0,790.7,0,.01,.11,.1,.5,.13,.03,1450.6,0,0,0,.1,.48,.29,.14,3216.2,0,0,0,0,0,.33,.67,4066.9,0,0,0,0,0,0,1)
)

Ymod=list(c(0,.1,0,.2,0,0,0,0,0,0,.7,0,0,0,0,0,0,0,.95,.01,0,0,0,0,0,0,.17,.66,.17,0,0,0,.9,0,0,0,.96,.04,0,0,2.3,0,0,0,0,.88,0,0,2.6,0,0,0,0,.04,.96,0,.8,0,0,0,0,0,0,1),
          c(0,.1,0,0,0,0,0,0,1.8,0,.7,0,0,0,0,0,2.3,0,.76,.09,0,0,0,0,12.6,0,.05,.66,.29,0,0,0,73.4,0,.04,.06,.56,.29,.04,0,153.4,0,.02,0,.04,.67,.27,0,568.9,0,0,0,.07,0,.6,.33,1431.7,0,0,0,0,0,0,1),
          c(0,.1,0,0,0,0,0,0,.07,0,.07,0,0,0,0,0,5,0,.66,.25,.02,0,0,0,65.3,0,.04,.5,.33,.08,0,0,189.1,0,0,0,.52,.32,.03,.02,330.7,0,0,0,0,.76,.17,.03,679.7,0,0,0,0,.17,.5,.33,1142.8,0,0,0,0,0,.17,.83),
          c(0,.1,0,0,0,0,0,0,0,0,.7,0,0,0,0,0,1.8,0,.52,.19,.03,.03,0,0,14.3,0,.04,.32,.16,.36,0,0,18.1,0,0,.06,.42,.38,.02,0,61.2,0,0,0,.05,.43,.3,.16,125.2,0,0,0,.02,.1,.24,.56,179.4,0,.08,0,0,0,0,.92),
          c(0,.1,0,0,0,0,0,0,0.7,0,.7,0,0,0,0,0,1.3,0,.43,.22,0,0,0,0,62.1,0,0,.55,.27,.05,0,0,306.9,0,0,0,.43,.45,.06,0,579.7,0,0,0,0,.64,.31,.03,890.6,0,0,0,.1,.15,.5,.2,1843,0,0,0,0,.13,.13,.75),
          c(0,.1,0,0,0,0,0,0,0,0,.7,0,0,0,0,0,1.1,0,.7,.08,.01,0,0,0,58.6,0,.1,.49,.27,.12,0,0,190.3,0,0,.07,.5,.23,.07,.03,481.5,0,0,0,.03,.51,.33,.1,702,0,0,0,.03,.21,.38,.28,1508.1,0,0,0,0,0,0,.94),
          c(0,.1,0,0,0,0,0,0,0,0,.7,0,0,0,0,0,11.4,0,.57,.12,.01,0,0,0,110,.0,.2,.4,.4,0,0,0,790.7,0,.01,.11,.1,.5,.13,.03,1450.6,0,0,0,.1,.48,.29,.14,3216.2,0,0,0,0,0,.33,.67,4066.9,0,0,0,0,0,0,1)
)

X=lapply(Y,function(y) matrix(y,8,8))
Xmod=lapply(Ymod,function(y) matrix(y,8,8))
Xu=unlist(X)

# SGmean(proj_matrices=X,transitions=P1,N=100000)

# [1] 2.335509e-01 4.292320e-06 8.910421e-08

# s1=sensitivity1(X,P1,1000)

# > s1[[1]]
# $Mean
# [,1]        [,2]        [,3]         [,4]        [,5]         [,6]         [,7]         [,8]
# [1,]  0.01286207 0.002256441  0.01574393 0.0005279406 0.001019373 0.0003047901 0.0000622153 5.184919e-05
# [2,]  0.14064219 0.022888573  0.16444406 0.0054741292 0.010717643 0.0032718196 0.0006679588 5.319845e-04
# [3,]  0.18831595 0.033803079  0.24330890 0.0082395017 0.015368043 0.0048232204 0.0009855367 8.174750e-04
# [4,]  1.15357182 0.149339259  1.26770207 0.0422536210 0.078577662 0.0240152588 0.0052253436 4.287337e-03
# [5,]  2.24892283 0.277106806  2.32294001 0.0769392329 0.150073215 0.0457695903 0.0098023008 7.796048e-03
# [6,]  1.76974024 0.276839779  1.92516745 0.0681001318 0.120153764 0.0377936601 0.0076599439 6.279537e-03
# [7,]  5.22142004 0.981367580  6.56555920 0.2157563031 0.415088919 0.1233589255 0.0261018442 2.112696e-02
# [8,] 12.40239970 2.007800503 13.55003090 0.4583612420 0.865597067 0.2615662622 0.0546087822 4.522816e-02
# 
# $Variance
# [,1]         [,2]         [,3]         [,4]         [,5]         [,6]         [,7]         [,8]
# [1,] 1.501975e-06 2.595676e-08 5.527168e-08 6.019776e-11 3.904930e-10 3.294057e-11 1.050571e-12 5.246970e-13
# [2,] 1.708122e-04 2.743976e-06 4.867115e-06 5.233458e-09 3.682473e-08 3.381954e-09 9.562475e-11 4.129887e-11
# [3,] 3.388803e-04 5.914306e-06 8.841156e-06 8.999232e-09 6.823657e-08 6.629221e-09 1.647969e-10 7.857292e-11
# [4,] 1.083721e-02 1.428650e-04 2.415816e-04 2.192641e-07 1.604881e-06 1.435564e-07 4.460662e-09 1.759043e-09
# [5,] 4.090334e-02 4.777548e-04 8.757892e-04 8.455223e-07 5.629069e-06 5.167065e-07 1.773359e-08 7.144160e-09
# [6,] 2.937817e-02 4.299231e-04 2.113432e-03 2.615945e-06 8.631443e-06 8.116895e-07 3.388985e-08 1.826573e-08
# [7,] 2.600379e-01 4.759981e-03 1.069499e-02 1.079454e-05 5.573842e-05 4.945337e-06 1.882488e-07 8.489873e-08
# [8,] 1.288559e+00 2.118539e-02 7.020504e-02 7.385954e-05 3.902477e-04 3.280149e-05 1.302347e-06 6.707575e-07
# 
# $Bias
# [,1]         [,2]         [,3]         [,4]         [,5]         [,6]         [,7]         [,8]
# [1,] 1.191623e-07 1.116182e-07 1.112687e-07 9.322099e-08 1.046598e-07 9.600034e-08 1.009637e-07 9.706462e-08
# [2,] 2.753293e-07 1.173082e-07 1.516048e-07 9.444026e-08 1.019803e-07 9.913167e-08 1.013882e-07 1.062936e-07
# [3,] 3.787768e-07 1.344776e-07 1.620501e-07 1.076077e-07 1.050968e-07 1.075543e-07 9.739424e-08 1.004807e-07
# [4,] 1.872713e-06 2.770385e-07 4.555632e-07 1.254116e-07 1.157819e-07 1.074300e-07 1.060349e-07 1.047849e-07
# [5,] 2.801043e-06 4.083374e-07 8.116983e-07 1.185785e-07 1.436021e-07 1.256604e-07 9.130448e-08 9.785114e-08
# [6,] 2.622181e-06 3.996164e-07 1.043347e-06 1.431800e-07 1.439688e-07 1.114560e-07 1.063526e-07 9.984157e-08
# [7,] 8.392116e-06 9.890373e-07 2.568575e-06 1.922041e-07 2.544592e-07 1.489533e-07 1.133315e-07 1.110829e-07
# [8,] 1.811955e-05 2.279736e-06 6.101397e-06 3.446453e-07 6.228306e-07 2.083362e-07 1.302004e-07 1.215632e-07

# > SGmean(Xmod,P1,100000)
# [1] 3.452544e-01 4.125799e-06 7.950817e-08

# s2=sensitivity2(X,P1,10000,500)
# 
# > s2
# $mean
# [,1]       [,2]       [,3]      [,4]      [,5]      [,6]      [,7]
# [1,]    0 0.39819106 1.06440920 1.5455268 2.2138087 3.2228336 4.4775992
# [2,]    0 0.04138872 0.09680221 0.1559391 0.2445112 0.3597096 0.5081297
# [3,]    0 0.03713765 0.08087182 0.1190507 0.2104891 0.3149138 0.4445068
# [4,]    0 0.04017335 0.08115001 0.1241095 0.1807485 0.2791917 0.3823390
# [5,]    0 0.02553069 0.06972826 0.1299994 0.2054757 0.3118477 0.4336626
# [6,]    0 0.05504393 0.12053052 0.2144009 0.3301829 0.4886798 0.6775851
# [7,]    0 0.03852828 0.08686189 0.1794192 0.2953095 0.4424452 0.6160436
# 
# $variance
# [,1]         [,2]         [,3]         [,4]         [,5]         [,6]         [,7]
# [1,]    0 6.294626e-04 5.779874e-04 6.108461e-04 6.370242e-04 6.124284e-04 6.677312e-04
# [2,]    0 1.655878e-05 1.984939e-05 1.887580e-05 1.853706e-05 1.730501e-05 1.944866e-05
# [3,]    0 1.331835e-05 1.506914e-05 1.483163e-05 1.388791e-05 1.382952e-05 1.398683e-05
# [4,]    0 1.069623e-05 1.137104e-05 1.033353e-05 1.037369e-05 1.097223e-05 1.174663e-05
# [5,]    0 1.328988e-05 1.394504e-05 1.410879e-05 1.439337e-05 1.407199e-05 1.556999e-05
# [6,]    0 1.846955e-05 2.117830e-05 1.851712e-05 1.877111e-05 1.911898e-05 1.775998e-05
# [7,]    0 2.943529e-05 3.140280e-05 3.047464e-05 2.813832e-05 2.812648e-05 2.712696e-05
# 
# $bias
# [,1]         [,2]         [,3]         [,4]         [,5]         [,6]         [,7]
# [1,]    0 3.058810e-04 2.965343e-04 3.064996e-04 3.266349e-04 3.350635e-04 3.771617e-04
# [2,]    0 3.421539e-05 3.356554e-05 3.488219e-05 3.806640e-05 4.089013e-05 4.289307e-05
# [3,]    0 2.981940e-05 3.069929e-05 3.142247e-05 3.347204e-05 3.495817e-05 3.807631e-05
# [4,]    0 2.431428e-05 2.397463e-05 2.623044e-05 2.710588e-05 2.883903e-05 3.202155e-05
# [5,]    0 2.950036e-05 3.019744e-05 3.056974e-05 3.102973e-05 3.373306e-05 3.917898e-05
# [6,]    0 4.583045e-05 4.595810e-05 4.778320e-05 4.914067e-05 5.378836e-05 5.829750e-05
# [7,]    0 4.444509e-05 4.453496e-05 4.706111e-05 4.807976e-05 5.138786e-05 5.719820e-05

# SGmean(X,P1mod,100000)  Should increase by about .01*4.478=.0448. Was .233. Becomes .279
#[1] 2.795000e-01 4.633401e-06 8.344288e-08
# So that's about right
