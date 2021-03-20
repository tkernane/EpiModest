#' Estimation of the parameters of SIR Epidemic Model
#' @export
#' @param I numeric vector variable for the number of infected individuals
#' @param S numeric vector variable for the number of scuceptible individuals
#' @param m numeric variable for the number of imputed data between each pair of observations
#' @param beta0 numeric variable representing the starting infection rate for the MCMC algorithm
#' @param mu0 numeric variable representing the starting recovery rate for the MCMC algorithm
#' @param cv numeric variable representing the number of iterations for the MCMC algorithm
#' @param cv numeric variable representing the number of iterations for the MCMC algorithm
#' @param bur numeric variable representing the burning period from which we begin the computation of the mean posterior values
SIRMCMC=function(I,S,m,cv,bur,beta0,mu0){
    library(GIGrvg)
    library(mvtnorm)
    ####################################
    ####################################
    ####################################

    MHAR=function(x123,x333,x223,dt1,n,b,m){
        instru=dnorm(.5*(x123[1]+x223[1]),.5*(x123[1]+x223[1]),sqrt(0.5*dt1*b*x123[1]*(1-x123[1]-x123[2])/n),log=TRUE)+
            dnorm(.5*(x123[2]+x223[2]),.5*(x123[2]+x223[2]),sqrt(0.5*dt1*m*(1-x123[1]-x123[2])/n),log=TRUE)

        sigma11=diag(c(b*.5*(x123[1]+x223[1])*(1-(.5*(x123[1]+x223[1]))-(.5*(x123[2]+x223[2])))/n,
                       m*(1-(.5*(x123[1]+x223[1]))-(.5*(x123[2]+x223[2])))/n))
        mu11=c(-b*.5*(x123[1]+x223[1])*(1-(.5*(x123[1]+x223[1]))-(.5*(x123[2]+x223[2])))*dt1,
               m*(1-(.5*(x123[1]+x223[1]))-(.5*(x123[2]+x223[2])))*dt1)

        sigma1=diag(c(b*x123[1]*(1-x123[1]-x123[2])/n,m*(1-x123[1]-x123[2])/n))
        mu1=c(-b*x123[1]*(1-x123[1]-x123[2])*dt1,m*(1-x123[1]-x123[2])*dt1)

        cible=dmvnorm(c(.5*(x123[1]+x223[1]),.5*(x123[2]+x223[2]))-x123,mu1,dt1*sigma1,log=TRUE)+
            dmvnorm(x223-c(.5*(x123[1]+x223[1]),.5*(x123[2]+x223[2])),mu11,dt1*sigma11,log=TRUE)

        c=exp(cible-instru)

        repeat{
            z=vector("numeric",2)
            z[1]=rnorm(1,.5*(x123[1]+x223[1]),sqrt(0.5*dt1*b*x123[1]*(1-x123[1]-x123[2])/n))
            z[2]=rnorm(1,.5*(x123[2]+x223[2]),sqrt(0.5*dt1*m*(1-x123[1]-x123[2])/n))
            sigma2=diag(c(b*z[1]*(1-z[1]-z[2])/n,m*(1-z[1]-z[2])/n))
            mu2=c(-b*z[1]*(1-z[1]-z[2])*dt1,m*(1-z[1]-z[2])*dt1)

            cibl=dmvnorm(z-x123,mu1,dt1*sigma1)*
                dmvnorm(x223-z,mu2,dt1*sigma2)

            inst=dnorm(z[1],.5*(x123[1]+x223[1]),sqrt(0.5*dt1*b*x123[1]*(1-x123[1]-x123[2])/n))*
                dnorm(z[2],.5*(x123[2]+x223[2]),sqrt(0.5*dt1*m*(1-x123[1]-x123[2])/n))

            if(z[1]<0 | z[2]<0 | z[1]>1 | z[2]>1) {alpha=0} else {alpha=min(1,cibl/(c*inst))}

            u1=runif(1)

            if(u1<alpha) break
        }
        sigma222=diag(c(b*x333[1]*(1-x333[1]-x333[2])/n,m*(1-x333[1]-x333[2])/n))
        mu222=c(-b*x333[1]*(1-x333[1]-x333[2])*dt1,m*(1-x333[1]-x333[2])*dt1)

        cibl3=dmvnorm(x333-x123,mu1,dt1*sigma1)*
            dmvnorm(x223-x333,mu222,dt1*sigma222)

        inst3=dnorm(x333[1],.5*(x123[1]+x223[1]),sqrt(sigma1[1,1]*0.5*dt1))*
            dnorm(x333[2],.5*(x123[2]+x223[2]),sqrt(sigma1[2,2]*0.5*dt1))


        if ( cibl3<c*inst3 ){
            alpha1=1
        } else if ( (cibl3>c*inst3) & (cibl<c*inst) ) {
            alpha1=(c*inst3)/cibl3
        } else alpha1=min(1,exp(log(cibl)-log(cibl3)+log(inst3)-log(inst)))


        u2=runif(1)
        if(u2<alpha1) {x=z} else {x=x333}

        return(x)

    }

    ####################################################################
    ####################################################################
    ####   data initialization     #####
    ######################################


    ########### I: vector of infected
    ########### S: vector of susceptibles


    T=length(I)       ### the simulation horizon
    a=I[1]            ### number of infected at t=0
    n=S[1]            ### number of susceptibles at t=0


    ############# m: number of missing data
    m1=m+1
    dt1=1/m1
    n1=n+a


    x111=S/n1
    x222=I/n1
    x333=1-x111-x222


    para1=1
    para2=0.2


    h1=matrix(c(x111,x333),ncol=length(x111),byrow=T)


    ############## cv: number of iterations to perform

    param=matrix(rep(0,cv*2),ncol=2)   #### matrix for storing parameters beta and mu (2 columns), the insertion is done by line

    param[1,1]=beta0                   #### initial value of the parameter beta

    param[1,2]=mu0                     #### initial value of the parameter mu

    d1=length(h1[1, ])
    N=((d1-1)*m1)+1
    y1=vector("numeric",N)
    y2=vector("numeric",N)

    y1[N]=h1[1,d1]
    y2[N]=h1[2,d1]




    for(i in 1:(d1-1)){

        y1[((i-1)*m1)+1]=h1[1,i]
        y2[((i-1)*m1)+1]=h1[2,i]

        for(j in 1:(m1-1)){

            y1[(m1*(i-1))+1+j]=(-(abs(h1[1,i]-h1[1,i+1])/m1)*j)+h1[1,i]
        }

    }

    for(i in 1:(d1-1)){
        if(y2[(i*m1)+1]<1+(a/n)-y1[((i-1)*m1)+1]){x1=runif(m1-1,y2[((i-1)*m1)+1],y2[(i*m1)+1])
        x1=sort(x1,decreasing = TRUE)
        for(j in 1:(m1-1)){
            y2[(m1*(i-1))+1+j]=x1[m1-j]}}
        else{x1=runif(m1-1,y2[((i-1)*m1)+1],1-y1[((i-1)*m1)+1])
        x1=sort(x1,decreasing = TRUE)
        for(j in 1:(m1-1)){
            y2[(m1*(i-1))+1+j]=x1[m1-j]}}

    }


    yn1=matrix(rep(0,(cv*N)),ncol=N)    #################### matrix for storing data for susceptible individuals, the data is inserted by line.
    ####################
    yn1[1, ]=y1

    yn2=matrix(rep(0,(cv*N)),ncol=N)    #################### matrix for storing data for removed individuals, the data is inserted by line.
    ####################
    yn2[1, ]=y2



    ################################################################
    ####         simulation of parameters and missing data           #####
    ###############################################################
    t=1
    while(t<cv){
        t=t+1

        for(i in 1:(N)){
            if((i-1)%%m1==0){yn1[t,i]=yn1[t-1,i];yn2[t,i]=yn2[t-1,i]}
            else{
                x=c(yn1[t,i-1],yn2[t,i-1])
                y=c(yn1[t-1,i],yn2[t-1,i])
                z=c(yn1[t-1,i+1],yn2[t-1,i+1])

                sim=MHAR(x,y,z,dt1,n1,param[t-1,1],param[t-1,2])
                yn1[t,i]=sim[1]
                yn2[t,i]=sim[2]
            }

        }
        lambda1=-(N/2)+para1
        lambda2=-(N/2)+para1
        sigma1=(n1/dt1)*(sum(((yn1[t,-1]-yn1[t,-N])^2)/
                                 (yn1[t,-N]*(1-yn1[t,-N]-yn2[t,-N]))))
        sigma2=(n1/dt1)*(sum(((yn2[t,-1]-yn2[t,-N])^2)/
                                 (1-yn1[t,-N]-yn2[t,-N])))
        gamma1=(n1*sum(yn1[t,-N]*(1-yn1[t,-N]-yn2[t,-N])*dt1))+(2*para2)
        gamma2=(n1*sum((1-yn1[t,-N]-yn2[t,-N])*dt1))+(2*para2)

        a1=rgig(1,lambda1,sigma1,gamma1)
        a2=rgig(1,lambda2,sigma2,gamma2)

        param[t,1]=a1
        param[t,2]=a2


        if(t%%1000==0) print(t)


    }
    l=list()
    l[[1]]=param
    l[[2]]=yn1
    l[[3]]=yn2

    return(l)
beta<-mean(a[[1]][bur:cv,1])
}
