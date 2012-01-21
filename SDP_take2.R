#################
#Useful functions
#################


#Function to calculate the probability of a given value of the data mean, given the prior mean and variance, and the monitoring effort (the variance in the data is also used, though is not a function argument)
pr_mu_dta<-function(Mpr,Vpr,nj,mu_dta){
	#Mean of belief pdf
	m1= Mpr
	#Variance of belief pdf
	if(nj>0){
	sig1= Vpr + Vdta/nj
	}else{
	sig1=Vpr
	}
	
	dnorm(mu_dta,mean=m1, sd=sqrt(sig1))
}

#Function to call a fortran version of pr_mu_dta. Assumes we have already used dyn.load
pr_mu_dta_f<-function(Mpr,Vpr,nj,mu_dta){
	#Predefine output vector
	pr_out=mu_dta*0
	#Call fortran routine

	output<-.Fortran("pr_mu_dta_f", Mpr=as.double(Mpr), Vpr=as.double(Vpr),Vdta=Vdta, nj=as.integer(nj), mu_dta=as.double(mu_dta), len_mu_dta=as.integer(length(mu_dta)), pr_out=as.double(pr_out)) 
	#return(output$pr_out)
	output$pr_out
	}

#This function takes the prior mean, prior variance, monitoring effort, data mean, critical value of the data mean, and logging decision in the j time-step, and outputs the prior mean and variance in the j+1 time step
mu_sigma_jp1<-function(Mpr, Vpr,nj, mu_dta, mu_crit, LDj){
	
	
	#Vectorised form
	Vpo=(1/Vpr+nj/Vdta)^(-1)
	Mpo=(Mpr+mu_dta*Vpr*nj/Vdta)/(1+Vpr*nj/Vdta) #Note - this is the usual formula divided by Vdta/nj -- done to prevent division by zero.
	

	#If logging is intented, apply logging effects. Note that in the case with nj=0, mu_crit is manipulated to ensure that this still gives the correct answer
		if((LDj==1)){
			Mpo= Mpo+lsrm*(mu_dta>mu_crit)
			Vpo= Vpo+lsrvar*(mu_dta>mu_crit)
		}
	
	#Apply population growth
	Mpo=Mpo+grm
	Vpo=Vpo+grvar
	
	#Output Mean and variance at j+1 time step
	cbind(Mpo,Vpo)		
} 

#This calls fortran code to replicate the function mu_sigma_jp1
mu_sigma_jp1_f<-function(Mpr, Vpr, nj, mu_dta, mu_crit, LDj){

	#predefine some local variables	
	Mpo=mu_dta*1
	Vpo=Mpo

#SUBROUTINE mu_sigma_jp1_f(Mpr, Vpr,Vdta,nj, mu_dta,len_mu_dta, mu_crit, LDj, Mpo, Vpo, lsrm, lsrvar, grm, grvar)
	output<- .Fortran('mu_sigma_jp1_f', Mpr=as.double(Mpr), Vpr=as.double(Vpr), Vdta=as.double(Vdta), nj=as.integer(nj),
	mu_dta=as.double(mu_dta), len_mu_dta=as.integer(length(mu_dta)), mu_crit=as.double(mu_crit), LDj=as.integer(LDj),
	Mpo=as.double(Mpo), Vpo=as.double(Vpo), lsrm=as.double(lsrm), lsrvar=as.double(lsrvar), grm=as.double(grm), 
	grvar=as.double(grvar) )

	cbind(output$Mpo, output$Vpo)	
	#output
	}

#This calculates the integrand in a big_nasty integral
useful<-function(mu_dta, Mpr, Vpr, nj, mu_crit,LDj, TR_jp1){

	#Next time step value of Mpr and Vpr
	#m_v=mu_sigma_jp1(Mpr,Vpr,nj,mu_dta,mu_crit,LDj) #Native R version
	m_v=mu_sigma_jp1_f(Mpr,Vpr,nj,mu_dta,mu_crit,LDj) #Fortran version

	#Now, we can do some crude rounding of these values, and associate them with a value of TR_jp1. Might need to think about a better way to do this sometime
	#Mu_jp1=findInterval(m_v[,1],Mbb, all.inside=T) #An index
	#Var_jp1=findInterval(m_v[,2],Varbb, all.inside=T) #An index
	Mu_jp1=findInterval_f(m_v[,1],Mbb) #An index, fortran version
	Var_jp1=findInterval_f(m_v[,2],Varbb) #An index, fortran version

	#Probability density of a given value of mu_dta
	#f1=pr_mu_dta(Mpr,Vpr,nj,mu_dta) #Native R version -- slow
	f1=pr_mu_dta_f(Mpr,Vpr,nj,mu_dta) #Fortran version -- fast
	TR_jp1[cbind(Mu_jp1,Var_jp1)]*f1
	}

#This calls a fortran version of useful (the function above)
useful_f<-function(mu_dta, Mpr, Vpr, nj, mu_crit, LDj, TR_jp1){
	#SUBROUTINE useful_f(mu_dta,len_mu_dta, Mpr, Vpr,Vdta, nj, mu_crit,LDj, TR_jp1, lsrm, &
	# lsrvar, grm, grvar, vals_mu_dta, Mbb, len_Mbb, Varbb, len_Varbb)

	vals_mu_dta=mu_dta*0 #Eventual output values of this function

	output<-.Fortran('useful_f',
		mu_dta=as.double(mu_dta), len_mu_dta=as.integer(length(mu_dta)), 
		Mpr=as.double(Mpr), Vpr=as.double(Vpr), Vdta=as.double(Vdta), nj=as.integer(nj), 
		mu_crit=as.double(mu_crit), LDj=as.integer(LDj), TR_jp1=as.double(TR_jp1),
		lsrm=as.double(lsrm), lsrvar=as.double(lsrvar), grm=as.double(grm), grvar=as.double(grvar),
		vals_mu_dta=as.double(vals_mu_dta), Mbb=as.double(Mbb), len_Mbb=as.integer(length(Mbb)), 
		Varbb=as.double(Varbb), len_Varbb=as.integer(length(Varbb)) 
		)
	
	output$vals_mu_dta
	}

#This calls a fortran version of useful (the function above)
useful_f2<-function(mu_dta, Mpr, Vpr, nj, mu_crit, LDj, TR_jp1){
	#SUBROUTINE useful_f(mu_dta,len_mu_dta, Mpr, Vpr,Vdta, nj, mu_crit,LDj, TR_jp1, lsrm, &
	# lsrvar, grm, grvar, vals_mu_dta, Mbb, len_Mbb, Varbb, len_Varbb)

	vals_mu_dta=mu_dta*0 #Eventual output values of this function
	Meanz=Means
	Variz=Varis

	output<-.Fortran('useful_f2',
		mu_dta=as.double(mu_dta), len_mu_dta=as.integer(length(mu_dta)), 
		Mpr=as.double(Mpr), Vpr=as.double(Vpr), Vdta=as.double(Vdta), nj=as.integer(nj), 
		mu_crit=as.double(mu_crit), LDj=as.integer(LDj), TR_jp1=as.double(TR_jp1),
		lsrm=as.double(lsrm), lsrvar=as.double(lsrvar), grm=as.double(grm), grvar=as.double(grvar),
		vals_mu_dta=as.double(vals_mu_dta), Mbb=as.double(Mbb), len_Mbb=as.integer(length(Mbb)), 
		Varbb=as.double(Varbb), len_Varbb=as.integer(length(Varbb)),Meanz=as.double(Means),Variz=as.double(Varis),1) 
		
	
	output$vals_mu_dta
	}
#This calls the fortran version of findInterval
findInterval_f<-function(Vec1, Vec2){
	indout=Vec1*0
	output<-.Fortran('findInterval_f', Vec1=as.double(Vec1), len_Vec1=as.integer(length(Vec1)),
	Vec2=as.double(Vec2), len_Vec2=as.integer(length(Vec2)), Index_vec1=as.integer(indout))

	output$Index_vec1
}


trapz_integral_f<-function(x,y){
	len_x=length(x)
	outval=0
	output<-.Fortran('trapz_integral_f',x=as.double(x),y=as.double(y), len_x=as.integer(length(x)), outval=as.double(outval))
	
	output$outval
	}

		
simps_integral_f<-function(x,y){
	len_x=length(x)
	outval=0
	if(len_x%%2==0){
		stop('Error in Simpsons Rule Integration: Need an odd number of points (even number of intervals)')
		}
	output<-.Fortran('simps_integral_f',x=as.double(x),y=as.double(y), len_x=as.integer(length(x)), outval=as.double(outval))
	
	output$outval
	}

#TR_n_LD=TR_n_LD_loop(sizen, tt, NTS, Mpr, Vpr, nj, 400, Pr_logging)
TR_n_LD_loop<-function(tt,NTS,Mpr,Vpr,nj,int_nos,Pr_logging,mu_crit, TR_jp1){

		TR_n_LD=matrix(0,ncol=2,nrow=length(nj))
		len_Mbb=length(Mbb)
		len_Varbb=length(Varbb)
		sizen=length(nj)
		
		output<- .Fortran('TR_n_LD_loop', Mpr=as.double(Mpr), Vpr=as.double(Vpr), Vdta=as.double(Vdta), nj=as.integer(nj), 
			mu_crit=as.double(mu_crit), Pr_logging= as.double(Pr_logging), sizen= as.integer(sizen), 
			tt= as.integer(tt), NTS= as.integer(NTS), int_nos= as.integer(int_nos), lsrm= as.double(lsrm), 
			lsrvar= as.double(lsrvar), grm= as.double(grm), grvar= as.double(grvar),Lr= as.double(Lr), Mc= as.double(Mc), 
			TR_n_LD= as.double(TR_n_LD),Mbb= as.double(Mbb),len_Mbb = as.integer(len_Mbb),
			Varbb=as.double(Varbb), len_Varbb=as.integer(len_Varbb),Means= as.double(Means), Varis=as.double(Varis)
			, TR_jp1=as.double(TR_jp1))		

#		(Mpr, Vpr, Vdta, nj, mu_crit,Pr_logging, sizen, tt, NTS, & 
#                                int_nos, lsrm, lsrvar, grm, grvar, Lr, Mc, TR_n_LD, & 
#                                Mbb, len_Mbb, Varbb, len_Varbb, Meanz, Variz)        
		cbind(output$TR_n_LD[1:length(nj)],output$TR_n_LD[(1:length(nj))+length(nj)])

		}

Starttime=Sys.time();
###############################
# 1. INITIAL PARAMETERS - user defined
############################### 

#Time stepping
NTS = 10; #  the number of time steps that we run the model for 

#Ecological(-ish) parameters
grm= 0.1 # the log(growth rate) mean
grvar= 0.04; # the log(growth rate) variance 
lsrm= -0.1;  # The log(survival rate) mean (i.e. survival following logging)   
lsrvar= 0.04 ;#The log(survival rate) variance 

Vdta= 0.1; #variance in the data (constant)- if we divide by the monitoring effort n, we get the (standard error)^2 of the data mean. 

#Money
Lr= 300000; # logging return ($) , not considering monitoring costs
Mc= 1000; # monitoring cost for one unit of monitoring ($)
dr= 0.8 # future discounting rate -- the value of (X dollars) at 1 time step in the future is (dr*X dollars) 

#Decision-theory
T= -0.5; # acceptable log(number of gliders) for logging to be permitted
p= 0.9; # if Pr( glider_density>T )>p , then logging is permitted

############################
#2. NUMERICAL PARAMETERS - so the solution should not be sensitive to these things, as long as they take 'reasonable' values
############################
n_max = 55; # the upper bound on the sampling effort (so there will be at MOST n_max units of monitoring done)

numMeans= 210; #the number of possible log(glider_density) means in the state-space 
numVaris= 55; #the number of population log(glider_density) variances in the state space
sizen=n_max+1 ; # the number of elements in the monitoring vector - so if this is n+1, the the values are 0,1,2,3,4,...n - or if it is (n+1)/2, then the values are 0,2,4,6,8,10,12,...

lmbb=-3 #Lower limit of log(glider_density) mean
umbb=2 #Upper limit of log(glider_density) mean

lvbb=0.0001 #Lower limit of log(glider_density) variance
uvbb=1 #Upper limit of log(glider_density) variance

####################################
#3. Setting up discrete state-space
####################################

#Mean bin boundaries
Mbb= seq(lmbb,umbb,len=numMeans+1);
#Centre of mean bins
Means=0.5*(Mbb[1:numMeans]+ Mbb[2:(numMeans+1)]);
if(min(diff(Means)) > grm){
	print('Some Mean bins will not register growth')	
}else{
	print('These Mean bins seem okay')
}

#NoMBins = length(Means) # number of mean bins 
#num_bbT = length(Mbb[Mbb<T]) # the number of mean bins below the threshold
#num_bbTCC = length(Mbb[Mbb>=T]) #the number of mean bins between the threshold and the carrying capacity.


#Variance bin boundaries - space these geometrically for now
#Varbb=lvbb*((uvbb/lvbb)^(1/(numVaris+1)) )^(seq(1,numVaris+1))
Varbb=seq(lvbb,uvbb,len=numVaris+1)
#Geometric centre of the variance bins
#Varis=sqrt(Varbb[1:numVaris]*Varbb[2:(numVaris+1)]);
Varis=0.5*(Varbb[1:numVaris]+ Varbb[2:(numVaris+1)]);


	##Some plotting
if(FALSE){
	Real_means=matrix(NA,nrow=numMeans,ncol=numVaris)
	Real_varis=Real_means
	for(i in 1:numMeans){
		for(j in 1:numVaris){
		Real_means[i,j]=exp(Means[i]+Varis[j]/2)
		Real_varis[i,j]=exp(2*Means[i]+Varis[j])*(exp(Varis[j])-1)
		}
	}
	par(mfrow=c(2,1))
	library(fields)
	image.plot(Means, Varis, log(Real_means,10),xlab='Transformed means',ylab='Transformed varis')
	image.plot(Means, Varis, log(Real_varis,10),xlab='Transformed means',ylab='Transformed varis')
}

##Load library containing fortran subroutines, which cover the most computationally intensive parts of the code.
##To make the libary (*nix), take the file SDP_take2.f90, open a terminal in its directory, and type:  R CMD SHLIB SDP_take2.f90
dyn.load("SDP_take2.so")

nj= ceiling( seq(0,n_max, len=sizen)) ; #Monitoring vector. Note that if n_max is too small, the code will warn (or at least try!)

######
###Predefine some storage arrays
######

#Total return at step j+1
TR_jp1 = matrix(0,nrow=numMeans,ncol=numVaris)

#Arrays which store the return, monitoring decision and logging decision for every time step
TR_best=array(0,dim=c(numMeans,numVaris,NTS))
n_best=TR_best
LD_best=TR_best

#Total return for a given nj and logging decision, where Mpr and Vpr are assumed (this appears in an inner loop)
TR_n_LD=matrix(0,nrow=sizen,ncol=2)



#The looping begins
for(tt in NTS:1){
	#print(paste('time-step',tt))

	for (jj in 1:length(Means)){ 
		
		print(paste("mean ",jj, 'time-step', tt))
		
		for (k in 1:length(Varis)){
			#print(paste('varis',k))
		
			Mpr= Means[jj];#Prior Mean
			Vpr=Varis[k];#Prior variance
			
			#At this stage, we can already calculate the posterior mean and variance assocaited with a given amount of monitoring
			#Posterior variance
			Vpo=(1/Vpr+nj/Vdta)^(-1)
			
			#mu_crit - critical data value for logging to be permitted
			mu_crit=rep(NA,len=sizen)
			mu_crit[2:sizen]=1/Vpr*( T*(Vdta/nj[2:sizen] +Vpr) - Mpr*Vdta/nj[2:sizen] - qnorm(1-p,0,1)*sqrt(Vpo[2:sizen])*(Vdta/nj[2:sizen] +Vpr) ) #Note that this was totally wrong in previous versions
			#mu_crit[1]= 9999999 #In this case nj=0, and there is no monitoring
			
			#The probability that logging is permitted
			Pr_logging=rep(NA,len=sizen)
			Pr_logging[2:sizen]= 1 - pnorm(mu_crit[2:sizen],Mpr,sqrt(Vpr+Vdta/nj[2:sizen]))
			if(qnorm(1-p,Mpr,sqrt(Vpr)) > T){
			Pr_logging[1]=1
			mu_crit[1]=-999999999 #Note that this is a 'trick' to get the code mu_sigma_jp1 to work in the case without monitoring
					}else{
			Pr_logging[1]=0
			mu_crit[1]= 999999999 #Note that this is a 'trick' to get the code mu_sigma_jp1 to work in the case without monitoring
				}
		if(max(Pr_logging>1)){
			print('MAX_pr_logging>1')
				}

			#Here, we calculate the total return for each amount of monitoring and each logging decision. Thus, we can find the optimal values 
		int_nos=4800*3+1
		if(int_nos%%2==0){
			stop('int_nos must be odd for simpson integration')
			}
		TR_n_LD=TR_n_LD_loop(tt, NTS, Mpr, Vpr, nj, int_nos, Pr_logging, mu_crit, TR_jp1)
if(FALSE){
		fzzs=0
		for(LDj in c(0,1)){
			for(nn in 1:sizen){

			#Here we calculate the integral of the total return function. Of course, there is not point doing this on the NTS time step. 
			if(tt!=NTS){
				ats=seq(Mpr-6*sqrt(Vpr+Vdta/max(nj[nn],1)*(nj[nn]>0)), Mpr+6*sqrt(Vpr+Vdta/max(nj[nn],1)*(nj[nn]>0)),len=400)
				fzzs=useful_f2(ats,Mpr,Vpr,nj[nn],mu_crit[nn],LDj,TR_jp1)
				if(max(fzzs)>0){				
					#Trapezoidal rule
					bigint=trapz_integral_f(ats,fzzs)

						}else{
						bigint=0
					}

					}else{
				bigint=0
				}
			bigint=as.numeric(bigint)
			#print(bigint)
			#Total return for each monitoring effort and logging decision
		
			TR_n_LD[nn,LDj+1]= bigint+ Pr_logging[nn]*LDj*Lr - nj[nn]*Mc
			#Need to be careful about the bounds on this integration, and about the interpolation in useful
			}

		}
}		
		#print(bigint)
		#if(max(abs(TR_n_LDa-TR_n_LD))>1){
		#stop(paste("TR_DIFFERENCE", max(abs(TR_n_LDa-TR_n_LD))))
		#}
		#Find max value of TR_n_LD in each column
		v=max.col(t(TR_n_LD),ties.method='first')
		#print(TR_n_LD)
		#print(paste('ERROR (now sleep): Length v = ', length(v), TR_n_LD))
		if(TR_n_LD[v[1],1]>TR_n_LD[v[2],2]){ #So here, NOT choosing to log is optimal - this is probably an unusual case
			LD_best[jj,k,tt]=0
			n_best[jj,k,tt]=nj[v[1]]
			TR_best[jj,k,tt]=TR_n_LD[v[1],1]
			print(paste('######## NOT choosing to log is optimal.', 'No_logging :$', TR_n_LD[v[1],1],'Logging :$', TR_n_LD[v[2],2], 'nbest=',n_best[jj,k,tt]))
			}else{#So here, choosing to log is optimal - this is probably the normal case
			LD_best[jj,k,tt]=1
			n_best[jj,k,tt]=nj[v[2]]
			TR_best[jj,k,tt]=TR_n_LD[v[2],2]
			#print(paste('No Logging :$', TR_n_LD[v[1],1],'; Logging :$',TR_n_LD[v[2],2], 'n_best=', n_best[jj,k,tt] ))
			}

			
		if(n_best[jj,k,tt] == max(nj)){
			print('The upper bound on nj was chosen as optimal -- need to increase the possible amount of monitoring')
			stop()
			}
			
		
			}#End of variance loop
		}#End of means loop

		TR_jp1=TR_best[,,tt]*(dr)
	
}#End of time loop

save(TR_best, n_best, LD_best,file='Outputs.R')


