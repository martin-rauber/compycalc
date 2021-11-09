####################################################################################
#corr_100EC.R: correction to 100% EC-yield
####################################################################################
#This part of the code was written by Gary Salazar: gary.salazar@dcb.unibe.ch
####################################################################################
 
  library(MASS)
alfa<-function(b) { #fixed T for 3 sunset steps
  alfa1<-exp(-t1*Kref*exp(Ea/R*((1/Tref)-(b/temp1))))
  alfa2<-alfa1*exp(-t2*Kref*exp(Ea/R*((1/Tref)-(b/temp2))))
  alfa3<-alfa2*exp(-t3*Kref*exp(Ea/R*((1/Tref)-(b/temp3))))
  alfa<-alfa3
  return(alfa)
}

SS_for_a<-function(par) {
  a<-par[1:n];b<-par[(n+1):length(par)]
  alfa.nv<-alfa(b) # non volatile
  alfa.v<-alfa(rep(1,length(b)))  # volatile
  Ys<-(alfa.nv+(a*alfa.v))/(1+a)
  Fs<-((a*alfa.v*Fv)+(alfa.nv*Fnv))/(Ys*(1+a))
  sum_sqrs<-1/((F_meas_sig^2+Y_meas_sig^2))*sum((Fs-F_meas)^2+(Ys-Y_meas)^2)
  return(sum_sqrs)
}

####### Fit all data with Fixed T, individuals as and bs
data<-list(c(0.42,0.56,0.62,0.76),c(0.12,0.40),c(0.61,0.70,0.91),c(0.5,0.62,0.78),c(0.3,0.56,0.77),c(0.74,0.8),c(0.48,0.92),c(0.74,0.85),
           c(0.66,0.74,0.77),c(0.73,0.92),c(0.7,0.77,0.8,0.83),c(0.39,0.65),c(.75,.81,.84),c(.75,.85,.93))
lin_par<-rbind(c(0.29,-.068),c(.31,.48),c(.32,.15),c(.3,.12),c(.32,.08),c(.23,.5),c(.53,.02),c(.40,.19),c(.43,.07),c(.26,.18),c(.48,-.14),
               c(.21,.03),c(.09,.11),c(.05,.001))
locations<-c("Gote wint", "Pay wint", "Pay sum", "Duben fall", "Duben fall2", "Pay Jan", "Rov Jan", "Zur Nov", "Chi Jan", "Zur Apr", "Chi Jan", "Ber Jan", "Zur Jul", "Ber Jun")

data_n<-8
F_meas<-c(); Y_meas<-c()
for (counter in c(1:14)) {  #data collection from paper
  data_n<-counter
  F_m<-(data[[data_n]]*lin_par[data_n,1])+lin_par[data_n,2]; Y_m<-data[[data_n]]
  if (counter==1) {F_meas<-F_m; Y_meas<-Y_m}
  if (counter>1) {F_meas<-c(F_meas,F_m); Y_meas<-c(Y_meas,Y_m)} 
}

#####

#get EC-yield data
Y_meas = EC_yield_mean_summary

#get EC F14C raw data
F_meas = F14C_EC_raw_data

#####

b0<-1.5 #----> b IS NOT FIXED <-----------
a<-c(); alfa.nv.all<-c(); alfa.v.all<-c(); Ys.all<-c(); Fs.all<-c(); linmodel_slope<-c(); F0.all<-c()
n<-length(Y_meas); Fnv<-0; Fv<-1.11; Fnv_sig<-0.03; Fv_sig<-0.03*Fv; F_meas_sig<-0.03;Y_meas_sig<-0.01
Kref<-0.00003; Tref<-673; Ea<-90; R<-8.31E-3
#three sunset steps temperature and times
temp1<-648; temp2<-698; temp3<-973
t1<-240; t2<-120; t3<-120

#calc a, b using whole data together
initial_par<-c(rep(1,n),rep(b0,n))          #vector of initial a's and b's
low_par<-c(rep(0.005,n),rep(0.90,n))       #vector of lower bound a's and b's
up_par<-c(rep(4,n),rep(4*b0,n))             #vector of upper bound a's and b's
model3<-optim(initial_par,SS_for_a,lower=low_par,upper=up_par,method="L-BFGS-B",hessian=TRUE) #error MINIMIZATION changing a & T on non linear model
a<-model3$par[1:n]; b<-model3$par[(n+1):(2*n)]

#H<-model3$hessian/2
#C<-ginv(H);variance<-diag(C)                   #diagonal of inverse hessian
#a_sig<-sqrt(abs(variance[1:n])); b_sig<-sqrt(abs(variance[(n+1):(2*n)])) #vector of uncertainties of a's and b's from hessian

alfa.nv.all<-alfa(b) # non volatile with fixed 3 thermal steps
alfa.v.all<-alfa(rep(1,length(b)))  # volatile with fixed 3 thermal steps
Ys.all<-(alfa.nv.all+(a*alfa.v.all))/(1+a) #calculated yield
Fs.all<-((a*alfa.v.all*Fv)+(alfa.nv.all*Fnv))/(Ys.all*(1+a)) #calculated Fm
linmodel_slope<-((a*Fv/(1+a))-F_meas)/(1-Y_meas)  #slope of linear model
F0.all<-((a*Fv)+Fnv)/(1+a)  #Fm correction by extrapolation to yield=1 non linear model

#error calculation for the above section by propagation
#F_meas_sig<-(sum((F_meas-Fs.all)^2))/n;Y_meas_sig<-(sum((Y_meas-Ys.all)^2))/n 

plot(F_meas,Fs.all,xlim=c(0,1),ylim=c(0,1))
points(Y_meas,Ys.all,col="blue")
rsqr2<-(cor(c(Ys.all,Fs.all),c(Y_meas,F_meas)))^2


###end Gary###
F0.all


