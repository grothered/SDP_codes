##R code to read the SDP file outputs.
n_best=array(scan('n_best'),dim=c(55,55,10))
n_best[,,1:10]=n_best[,,10:1]
#So n_best is the optimal monitoring effort 
#-- e.g. n_best[20,12,5] is the optimal monitoring effort associated with Means(20), Varis(12), and time step 5 -- so if we have 10 time steps, then n_best[,,10] is the single time step optimization, and n_best[,,1] is the 10 time step optimization
TR_best=array(scan('TR_best'),dim=c(55,55,10))
TR_best[,,1:10]=TR_best[,,10:1]
#So TR_best is the expected return, indexed in the same way as n_best.


#A plot
Mbb=seq(-3,2,len=56)
Means=0.5*(Mbb[2:56]+Mbb[1:55])
Vbb=seq(0.0001,1,len=56)
Varis=0.5*(Vbb[2:56]+Vbb[1:55])

library(fields)
png('optimalmon_10.png',width=5.4,height=5,units='in',res=600)
drape.plot(Means, Varis, n_best[,,10],xlab='Prior Mean',ylab='Prior Variance',main='Optimal monitoring effort \n optimized over a single time step',theta=-50)
dev.off()

png('optimalmon_1.png',width=5.4,height=5,units='in',res=600)
drape.plot(Means, Varis, n_best[,,1],xlab='Prior Mean',ylab='Prior Variance',main='Optimal monitoring effort \n optimized over 10 time steps', theta=-50)
dev.off()

png('Expected_return_1.png',width=5.4,height=5,units='in',res=600)
drape.plot(Means, Varis, TR_best[,,1],xlab='Prior Mean',ylab='Prior Variance',main='Expected return \n optimized over 10 time steps',theta=-50)
dev.off()

png('Expected_return_10.png',width=5.4,height=5,units='in',res=600)
drape.plot(Means, Varis, TR_best[,,10],xlab='Prior Mean',ylab='Prior Variance',main='Expected return \n optimized over a single time step',theta=-50)
dev.off()
