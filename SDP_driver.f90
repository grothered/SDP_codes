!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Program to run this stuff.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM SDP_run

use SDP_routines

IMPLICIT NONE

!Time stepping
INTEGER:: NTS = 2; !  the number of time steps that we run the model for 

!Ecological(-ish) parameters
REAL(dp):: grm= 0.1_dp ! the log(growth rate) mean
REAL(dp):: grvar= 0.04_dp! the log(growth rate) variance 
REAL(dp):: lsrm= -0.1_dp  ! The log(survival rate) mean (i.e. survival following logging)   
REAL(dp):: lsrvar= 0.04_dp !The log(survival rate) variance 

REAL(dp):: Vdta= 0.1_dp !variance in the data (constant)- if we divide by the monitoring effort n, we get the (standard error)^2 of the data mean. 

!Money
REAL(dp):: Lr= 300000._dp ! logging return ($) , not considering monitoring costs
REAL(dp):: Mc= 1000._dp ! monitoring cost for one unit of monitoring ($)
REAL(dp):: dr= 1.0_dp ! future discounting rate -- the value of (X dollars) at 1 time step in the future is (dr*X dollars) 

!Decision-theory
REAL(dp):: T= -0.5_dp ! acceptable log(number of gliders) for logging to be permitted
REAL(dp):: p=  0.9_dp ! if Pr( glider_density>T )>p , then logging is permitted

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!2. NUMERICAL PARAMETERS - so the solution should not be sensitive to these things, as long as they take 'reasonable' values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER:: n_max = 55 ! the upper bound on the sampling effort (so there will be at MOST n_max units of monitoring done)
INTEGER:: numMeans= 55 !the number of possible log(glider_density) means in the state-space 
INTEGER:: numVaris= 55 !the number of population log(glider_density) variances in the state space
INTEGER:: sizen =56 !=n_max+1  ! the number of elements in the monitoring vector - so if this is n+1, the the values are 0,1,2,3,4,...n - or if it is (n+1)/2, then the values are 0,2,4,6,8,10,12,...
INTEGER:: int_nos=2001

REAL(dp):: lmbb=-3._dp !Lower limit of log(glider_density) mean
REAL(dp):: umbb=2._dp !Upper limit of log(glider_density) mean

REAL(dp):: lvbb=0.0001_dp !Lower limit of log(glider_density) variance
REAL(dp):: uvbb=1._dp !Upper limit of log(glider_density) variance

!Work/temp type variables
INTEGER:: i, tt,jj,k, v(2), tmp(1)

!Functions
!REAL(dp):: r8_normal_01_cdf_inverse

REAL(dp):: x, Mpr, Vpr
Logical:: writer

!Means and variances and monitoring, and various arrays to store the results
REAL(dp):: Mbb, Means, Varbb, Varis, Vpo, TR_jp1, TR_best, TR_n_LD, mu_crit, Pr_logging
INTEGER:: nj, n_best, LD_best 
Allocatable Mbb(:), Means(:), Varbb(:), Varis(:), Vpo(:), nj(:), TR_jp1(:,:), &
    TR_best(:,:,:), TR_n_LD(:,:), n_best(:,:,:), LD_best(:,:,:), mu_crit(:), Pr_logging(:)

Allocate( Mbb(numMeans+1), Means(numMeans), Varbb(numVaris+1), Varis(numVaris), Vpo(sizen), &
    nj(sizen),TR_jp1(numMeans,numVaris), TR_best(numMeans,numVaris,NTS),& 
    TR_n_LD(sizen,2), n_best(numMeans,numVaris,NTS), LD_best(numMeans,numVaris,NTS), mu_crit(sizen), & 
    Pr_logging(sizen) )

open(1,file='n_best')
open(2,file='TR_best')
open(3,file='Dumper')

!Define mean bin boundaries
x=(umbb-lmbb)/(1._dp*numMeans)
DO i=0,numMeans
Mbb(i+1)=lmbb+i*x
END DO
Means=0.5_dp*( Mbb(1:numMeans)+Mbb(2:(numMeans+1)) )!Mean bins

!Define variance bin boundaries
x=(uvbb-lvbb)/(1._dp*numVaris)
DO i=0,numVaris
Varbb(i+1)=lvbb+(i*1._dp)*x
END DO
Varis=0.5_dp*( Varbb(1:numVaris)+Varbb(2:(numVaris+1)) )! Variance bins

!Define monitoring vector
nj(1)=0
DO i=1,sizen-1
nj(i+1)=i
END DO

!Predefine
TR_jp1 = 0._dp
TR_best = 0._dp
TR_n_LD=0._dp
n_best=0._dp
LD_best=0._dp

!!!!!!!!!!!!!!!!!!
! Begin main loop
!!!!!!!!!!!!!!!!!!

DO tt = NTS, 1, -1  !Loop backwards through time-steps
!$OMP PARALLEL SHARED(LD_best, n_best, TR_best) PRIVATE(jj, k, Mpr, Vpr, Vpo, x, mu_crit, Pr_logging, writer, TR_n_LD, tmp, v)
!$OMP DO SCHEDULE(DYNAMIC, 10)
    DO jj=1, numMeans !Loop through all means

        print*, 'Mean ', jj, ' Time ', tt 

        DO k=1, numVaris !Loop through all variances
        
                Mpr=Means(jj)   !Prior mean
                Vpr=Varis(k)    !Prior variance

                Vpo = (1.0_dp/Vpr+(1._dp*nj)/Vdta)**(-1.0_dp) ! Posterior variance for a given amount of monitoring

                x=r8_normal_01_cdf_inverse(1._dp-p) !Useful variable
                
                mu_crit(2:sizen) = 1._dp/Vpr*( T*(Vdta/(1._dp*nj(2:sizen)) +Vpr) - Mpr*Vdta/(1._dp*nj(2:sizen)) & 
                        - x*(Vpo(2:sizen))**0.5_dp*(Vdta/(1._dp*nj(2:sizen)) +Vpr) ) !Data value required for logging to be permitted. 
                
                !Note -- the alnorm function gives slightly different results to
                !R's pnorm -- i.e. different in the 7-8 decimal place
                DO i=2,sizen !For all monitoring choices except n=0
                x= (mu_crit(i)-Mpr)/(Vpr+Vdta/(1._dp*nj(i)))**0.5 !Standardised normal deviate
                Pr_logging(i)=1._dp- alnorm(x, .FALSE.) !Probability of logging
                END DO
                
                !Case without monitoring, n=0 
                x=alnorm((T-Mpr)/(Vpr)**0.5, .FALSE.) !Standardized deviate
                IF(x<(1._dp-p)) THEN !Logging is permitted
                Pr_logging(1)=1._dp
                mu_crit(1)=-999999999999._dp
                ELSE
                Pr_logging(1)=0._dp
                mu_crit(1)=9999999999999_dp
                END IF 
                !write(1,*) Pr_logging(sizen/2)                

                !Up to here, I have checked I am getting the same results as the
                !R code to within about 1.0E-08  
                IF(k==1) THEN
                    writer=.TRUE.
                ELSE
                    writer=.FALSE.
                END IF              
                call TR_n_LD_loop(Mpr, Vpr, Vdta, nj, mu_crit,Pr_logging, sizen, tt, NTS, & 
                                int_nos, lsrm, lsrvar, grm, grvar, Lr, Mc, TR_n_LD, & 
                                Mbb, size(Mbb), Varbb, size(Varbb), Means, Varis, TR_jp1,writer)        
                tmp = maxloc(TR_n_LD(:,1))
                v(1)=tmp(1)
                tmp = maxloc(TR_n_LD(:,2))
                v(2)=tmp(1)
                !v = maxloc(TR_n_LD,1)
                !print*, v !, size(TR_n_LD(:,1))!, size(TR_n_LD(1,:)) !, maxloc(TR_n_LD,1)
               
                IF(TR_n_LD(v(1),1)>TR_n_LD(v(2),2)) THEN
                print*, 'NOT CHOOSING TO LOG IS OPTIMAL', 'No Logging:$ ', TR_n_LD(v(1),1), 'Logging:$ ', TR_n_LD(v(2),2)
                        LD_best(jj,k,tt)=0
                        n_best(jj,k,tt)=nj(v(1))
                        TR_best(jj,k,tt)=TR_n_LD(v(1),1)
                ELSE
                        LD_best(jj,k,tt)=1
                        n_best(jj,k,tt)=nj(v(2))
                        TR_best(jj,k,tt)=TR_n_LD(v(2),2)
                END IF 

                !Check that the maximal value of n is okay
                IF(n_best(jj,k,tt)==nj(sizen)) THEN
                        print*, 'The max value of n is selected, consider increasing n_max'
                        stop
                END IF
                !IF((k==1).AND.(tt==1)) write(3,*) TR_n_LD(:,2) !print*, Means(jj), TR_best(jj,k,tt), n_best(jj,k,tt)
        END DO !End variances loop
    END DO !End means loop
!$OMP END DO
!$OMP END PARALLEL
    
    write(1,*) n_best(:,:,tt)
    write(2,*) TR_best(:,:,tt)

    TR_jp1=TR_best(:,:,tt)*(dr)
END DO !End time loop

    !Close files for writing
    close(1)
    close(2)
    close(3)
END PROGRAM
