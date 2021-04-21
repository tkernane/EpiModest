S_change<-rep(0,2)
E_change<-rep(0,2)
I_change<-rep(0,2)
T_change<-rep(0,2)
a<-2.32  #taux d'infection
b<-2      #période d'incubation
c<-1.47     #taux de rétablisement 
T_fin<-20 #temps
N<-260    #taille de la population
I<-5      #individus infecté
E<-3      #individus exposé
R<-0
S<-N-(I+E+R)    #individus suscéptible
k<-1
t<-0
dt<-T_fin/N   #pas de descritisation
S_change[1]<-S
E_change[1]<-E
I_change[1]<-I
R_change[1]<-R
TT<-seq(0,T_fin,by=dt)
TTT<-TT[-1]
T_change<-TTT
Z1<-runif(N)
Z2<-runif(N)
Z3<-runif(N)
i<-1

while((E_change[i]>0) && ((I_change[i]>0) && (TTT[i]<T_fin))){
i<-i+1
S_change[i]<-S_change[i-1]-(a*S_change[i-1]*I_change[i-1]*dt/N)-(sqrt(a*S_change[i-1]*I_change[i-1]*dt/N)*Z1[i-1])
E_change[i]<-E_change[i-1]+(a*S_change[i-1]*I_change[i-1]*dt/N)-(b*E_change[i-1]*dt)+(sqrt(a*S_change[i-1]*I_change[i-1]*dt/N)*Z1[i-1])-(sqrt(b*E_change[i-1]*dt)*Z2[i-1])
I_change[i]<-I_change[i-1]+(b*E_change[i-1]*dt)-(c*I_change[i-1]*dt)+(sqrt(b*E_change[i-1]*dt)*Z2[i-1])-(sqrt(c*I_change[i-1]*dt))*Z3[i-1]
R_change[i]<-N-(S_change[i]+E_change[i]+I_change[i])
}

SS<-S
EE<-E
II<-I
if (E_change[length(E_change)]<0){
I_change<-I_change[1:length(I_change)-1]
E_change<-E_change[1:length(E_change)-1]
S_change<-S_change[1:length(S_change)-1]
R_change<-R_change[1:length(R_change)-1]
}
if (I_change[length(I_change)]<0){
I_change<-I_change[1:length(I_change)-1]
E_change<-E_change[1:length(E_change)-1]
S_change<-S_change[1:length(S_change)-1]
R_change<-R_change[1:length(R_change)-1]
}
A<-TTT[1:length(S_change)]
R_change<-R_change[1:length(S_change)]
MA<-list(A,S_change,E_change,I_change,R_change)

plot(A,S_change,type="l",col="blue",xlab="Temps",ylab="Individus",ylim=c(0,N))
lines(A,E_change,type="l",col="green")
lines(A,I_change,type="l",col="red")
lines(A,R_change,type="l",col="black")
op <- par(bg="antiquewhite1")
legend("topright",legend=c("Suceptibles","Exposed","Infected","Recovred"),col = c("blue","green","red","black"),pch=c(16,16,16),cex = 0.6)
