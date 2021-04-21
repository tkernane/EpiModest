#SIR<-function(T_fin,N,I,a,b){
S_change=rep(0,2)
I_change=rep(0,2)
R_change=rep(0,2)
T_change=rep(0,2)
a<-2.32  #taux d'infection
b<-2      #période d'incubation 
T_fin<-20 #temps
N<-260    #taille de la population
I<-5      #individus infecté
R=0
k=1
S=N-I
t=0
T_change[1]=t
S_change[1]=S
I_change[1]=I
R_change[1]=R
while (t<T_fin && I>0){
    k=k+1
    p1=(a*S*I)/N  #probabilité d'infection 
    p2=b*I        #probabilité de récupération 
    p<-c(p1,p2)   #matrice des probabilité
    r<-runif(2)   #nombres aléatoire 
    tt<-(-log(r[1])/sum(p))    #temps du premier occurence 
    valeur_comp<-cumsum(p)/sum(p);
    select_cass=which.max((valeur_comp > r[2]))   #l'événement qui se produise a l'instant tt  
    switch(select_cass,(S<-(S-1)):(I<-(I+1)),(I<-(I-1)):(R<-(R+1)))
   
    t<-(t+tt);
    S_change[k]<-S;
    I_change[k]<-I;
    R_change[k]<-R;
    T_change[k]<-t;
}
MA<-list(N,T_change,S_change,I_change,R_change)
#save('MA')
par(mfrow=c(2,2))
plot(T_change,S_change,type="l",col="blue",xlab="Temps",ylab="Suceptible",panel.first = grid(15,15))
plot(T_change,I_change,type="l",col="red",xlab="Temps",ylab="Infected",panel.first = grid(15,15))
plot(T_change,R_change,type="l",col="green",xlab="Temps",ylab="Recovered",panel.first = grid(15,15))




plot(T_change,S_change,type="l",col="blue",xlab="Temps",ylab="Individus",ylim=c(0,N))
lines(T_change,I_change,type="l",col="red")
lines(T_change,R_change,type="l",col="green")
op <- par(bg="antiquewhite1")
legend("topright",legend=c("Suceptibles","Infected","Recovred"),col = c("blue","red","green"),pch=c(16,16,16),cex = 0.6)
