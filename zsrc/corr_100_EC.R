####################################################################################
#correction to 100% EC-yield
####################################################################################

  ###Gary###
  
  library(MASS)
SS_for_a<-function(par) {
  a<-par[1:n_taken];b<-b_fixed;T<-par[(n_taken+1):length(par)]
  alfa.nv<-exp(log(alfa0)*exp(Ea_RT0-(1/T))) # non volatile
  alfa.v<-exp(log(alfa0)*exp(Ea_RT0-(b/T)))  # volatile
  Ys<-(alfa.nv+(a*alfa.v))/(1+a)
  Fs<-((a*alfa.v*Fv)+(alfa.nv*Fnv))/(Ys*(1+a))
  sum_sqrs<-1/((F_meas_sig^2+Y_meas_sig^2))*sum((Fs-Fm_taken)^2+(Ys-Y_taken)^2)
  return(sum_sqrs)
}

alfas<-function(par) {
  a<-par[1];b<-par[2];T<-par[3:length(par)]
  alfa.nv<-exp(log(alfa0)*exp(Ea_RT0-(1/T))) # non volatile
  alfa.v<-exp(log(alfa0)*exp(Ea_RT0-(b/T)))  # volatile
  return(c(alfa.v))
}

SS_forT<-function(par) {
  a<-par[1];b<-par[2];T<-par[3:length(par)]
  alfa.nv<-exp(log(alfa0)*exp(Ea_RT0-(1/T))) # non volatile
  alfa.v<-exp(log(alfa0)*exp(Ea_RT0-(b/T)))  # volatile
  Ys<-(alfa.nv+(a*alfa.v))/(1+a)
  Fs<-((a*alfa.v*Fv)+(alfa.nv*Fnv))/(Ys*(1+a))
  sum_sqrs<-sum((Fs-F_meas)^2+(Ys-Y_meas)^2)
  return(sum_sqrs)
}

Fs_forT<-function(par) {
  a<-par[1];b<-par[2];T<-par[3:length(par)]
  alfa.nv<-exp(log(alfa0)*exp(Ea_RT0-(1/T))) # non volatile
  alfa.v<-exp(log(alfa0)*exp(Ea_RT0-(b/T)))  # volatile
  Ys<-(alfa.nv+(a*alfa.v))/(1+a)
  Fs<-((a*alfa.v*Fv)+(alfa.nv*Fnv))/(Ys*(1+a))
  return(c(Ys,Fs))
}

####### Fit all data with Fixed b, individuals a?s and T?s
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

b_fixed<-0.86 #----> b IS FIXED <-----------
Ea_RT0<-7.5   # ----> max Ea/RT0 value <---
a<-c(); T<-c(); alfa.nv.all<-c(); alfa.v.all<-c(); Ys.all<-c(); Fs.all<-c(); linmodel_slope<-c(); F0.all<-c(); a_sig<-c(); T_sig<-c()
n<-length(Y_meas); Fnv<-0; Fv<-1.11  #main parameter values
Fnv_sig<-0.03; Fv_sig<-0.03*Fv; F_meas_sig<-0.03;Y_meas_sig<-0.01 #main parameters precision for error propagation

#calc a, a_sigma,T, T_sigma using whole data together
alfa0<-0.045; T0<-0.10
initial_par<-c(rep(1,n),rep(T0,n))          #vector of initial a's and T's
low_par<-c(rep(0.01,n),rep(0.1*T0,n))       #vector of lower bound a's and T's
up_par<-c(rep(Ea_RT0,n),rep(0.25,n))             #vector of upper bound a's and T's
n_taken<-n; Fm_taken<-F_meas; Y_taken<-Y_meas
model3<-optim(initial_par,SS_for_a,lower=low_par,upper=up_par,method="L-BFGS-B",hessian=TRUE) #error MINIMIZATION changing a & T on non linear model
a<-model3$par[1:n]; T<-model3$par[(n+1):(2*n)]
H<-model3$hessian/2
C<-ginv(H);variance<-diag(C)                   #diagonal of inverse hessian
a_sig<-sqrt(abs(variance[1:n])); T_sig<-sqrt(abs(variance[(n+1):(2*n)])) #vector of uncertainties of a's and T's from hessian

#calc a, a_sigma,T, T_sigma each sample individually
for (samples in c(1:n)) {
  alfa0<-0.045; T0<-0.10;n_taken<-1;H<-c()
  initial_par<-c(rep(1,n_taken),rep(T0,n_taken))          #vector of initial a's and T's
  low_par<-c(rep(0.01,n_taken),rep(0.1*T0,n_taken))       #vector of lower bound a's and T's
  up_par<-c(rep(Ea_RT0,n_taken),rep(0.25,n_taken))             #vector of upper bound a's and T's
  Fm_taken<-F_meas[samples]; Y_taken<-Y_meas[samples]
  model3<-optim(initial_par,SS_for_a,lower=low_par,upper=up_par,method="L-BFGS-B",hessian=TRUE) #error MINIMIZATION changing a & T on non linear model
  a[samples]<-model3$par[1]; T[samples]<-model3$par[2]
  H<-model3$hessian/2
  C<-ginv(H);variance<-diag(C)
  a_sig[samples]<-sqrt(abs(variance[1])); T_sig[samples]<-sqrt(abs(variance[2]))
}

alfa.nv.all<-alfa0^(exp(Ea_RT0-(1/T))) # non volatile
alfa.v.all<-alfa0^(exp(Ea_RT0-(b_fixed/T)))  # volatile
Ys.all<-(alfa.nv.all+(a*alfa.v.all))/(1+a) #calculated yield
Fs.all<-((a*alfa.v.all*Fv)+(alfa.nv.all*Fnv))/(Ys.all*(1+a)) #calculated Fm
linmodel_slope<-((a*Fv/(1+a))-F_meas)/(1-Y_meas)  #slope of linear model
F0.all<-((a*Fv)+Fnv)/(1+a)  #Fm correction by extrapolation to yield=1 non linear model

#error calculation for the above section by propagation
#F_meas_sig<-(sum((F_meas-Fs.all)^2))/n;Y_meas_sig<-(sum((Y_meas-Ys.all)^2))/n 

F0_err_a<-(Fv/(1+a))-(Fnv/((1+a)^2))-(a*Fv/((1+a)^2))
F0_err_Fv<-a/(1+a); F0_err_Fnv<-1/(1+a)
alfav_err_T<-b_fixed*log(alfa0)*exp(Ea_RT0-(b_fixed/T))*(1/T^2)*alfa0^(exp(Ea_RT0-(b_fixed/T)))
alfanv_err_T<-log(alfa0)*exp(Ea_RT0-(1/T))*(1/T^2)*alfa0^(exp(Ea_RT0-(1/T)))
F0_err_alfav<-(a*Fv/(1+a))-(((a^2*Fv)+(a*Fnv))/((1+a)^2))
F0_err_alfanv<-(Fnv/(1+a))-(((a*Fv)+Fnv)/((1+a)^2))
F0_sig_total<-sqrt((F0_err_a*a_sig)^2 +(F0_err_Fv*Fv_sig)^2 + (F0_err_Fnv*Fnv_sig)^2 + (F0_err_alfav*alfav_err_T*T_sig)^2 + (F0_err_alfanv*alfanv_err_T*T_sig)^2)

plot(F_meas,Fs.all)
points(Y_meas,Ys.all,col="blue")
rsqr2<-(cor(c(Ys.all,Fs.all),c(Y_meas,F_meas)))^2
lm(c(Ys.all,Fs.all)~c(Y_meas,F_meas))
write.table(cbind(linmodel_slope,F0.all,F0_sig_total),"clipboard",sep=",", row.names=FALSE)
write.table(cbind(a,a_sig,T, T_sig),"clipboard",sep=",", row.names=FALSE)

###end Gary###
F0.all


