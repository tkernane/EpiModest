S_change=rep(0,2)
I_change=rep(0,2)
E_change=rep(0,2)
R_change=rep(0,2)
T_change=rep(0,2)
N=150
I=5
E=3
R=0
k=1
S=N-(I+E)
t=0
a<-2.32  #taux d'infection
b<-2      #période d'incubation
c<-1.47     #taux de rétablisement
T_change[1]=t
S_change[1]=S
E_change[1]=E
I_change[1]=I
R_change[1]=R
while (t<T_fin && E>0 && I>0){
    k=k+1
    p1=(a*S*I)/N  #probabilité d'infection 
    p2=b*E        #probabilité d'expo
    p3=c*I        #probabilité de rétablissement  
    p<-c(p1,p2,p3)   #matrice des probabilité
    r<-runif(2)   #nombres aléatoire 
    tt<-(-log(r[1])/sum(p))    #temps du premier occurence 
    valeur_comp<-cumsum(p)/sum(p);
    select_cass=which.max((valeur_comp > r[2]))   #l'événement qui se produise a l'instant tt  
    switch(select_cass,(S<-(S-1)):(E<-(E+1)),(E<-(E-1)):(I<-(I+1)),(I<-(I-1)):(R<-(R+1)))
   
    t<-(t+tt);
    S_change[k]<-S;
    E_change[k]<-E;
    I_change[k]<-I;
    R_change[k]<-R;
    T_change[k]<-t;
}
MA<-list(N,T_change,S_change,E_change,I_change,R_change)
#save('MA')
par(mfrow=c(2,2))
plot(T_change,S_change,type="l",col="blue",xlab="Temps",ylab="Suceptible",panel.first = grid(15,15))
plot(T_change,E_change,type="l",col="green",xlab="Temps",ylab="Exposed",panel.first = grid(15,15))
plot(T_change,I_change,type="l",col="red",xlab="Temps",ylab="Infected",panel.first = grid(15,15))
plot(T_change,R_change,type="l",col="black",xlab="Temps",ylab="Recovered",panel.first = grid(15,15))


##################superposé les graphes################################


plot(T_change,S_change,type="l",col="blue",xlab="Temps",ylab="Individus",ylim=c(0,N))
lines(T_change,E_change,type="l",col="green")
lines(T_change,I_change,type="l",col="red")
lines(T_change,R_change,type="l",col="black")
op <- par(bg="antiquewhite1")
legend("topright",legend=c("Suceptibles","Exposed","Infected","Recovred"),col = c("blue","green","red","black"),pch=c(16,16,16),cex = 0.6)
