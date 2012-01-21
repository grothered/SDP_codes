!Fortran routines to take on some of the more computationally
!intensive stuff in the greedy forrestor algorithm.
module SDP_routines
IMPLICIT NONE
INTEGER, PARAMETER, PUBLIC:: dp = selected_real_kind(12,60)
contains

!Calculate the PRobability_density_of_MU_DTA_Fortran (i.e. Pr_mu_dta, the probability of a given data mean)
SUBROUTINE pr_mu_dta_f(Mpr,Vpr, Vdta, nj,mu_dta,len_mu_dta, pr_out)
        implicit none
        INTEGER, intent(in):: nj, len_mu_dta !monitoring effort, length of the mu_dta vector
        REAL(dp), intent(in):: Mpr, Vpr, Vdta !Prior Mean, Prior variance, data variance
        REAL(dp), intent(in):: mu_dta !The values where we need the probability density to be evaluated
        REAL(dp), intent(out):: pr_out !The probability vector associated with mu_dta 

        DIMENSION:: mu_dta(len_mu_dta), pr_out(len_mu_dta)
       
        !Local varaibles
        REAL(dp):: sig1_sq, pi=atan(1._dp)*4._dp

        !Calculate the variance of the mu_data belief pdf, accounting for the
        !case with no monitoring (nj==0) 
        IF(nj>0) THEN
                sig1_sq= Vpr+Vdta/(1._dp*nj)
        ELSE
                sig1_sq=Vpr
        END IF

        !Gaussian density
        pr_out=  1._dp/dsqrt(2._dp*pi*sig1_sq)*dexp(- ( (mu_dta - Mpr)**2/(2._dp*sig1_sq)) ) 
        return

END SUBROUTINE pr_mu_dta_f


!This function takes the prior mean, prior variance, monitoring effort, data mean, critical value of the data mean, and logging decision in the j time-step, and outputs the prior mean and variance in the j+1 time step
SUBROUTINE mu_sigma_jp1_f(Mpr, Vpr,Vdta,nj, mu_dta,len_mu_dta, mu_crit, LDj, Mpo, Vpo, lsrm, lsrvar, grm, grvar)
        implicit none
        INTEGER, intent(in):: len_mu_dta, nj,LDj !length of mu_dta, monitoring effort, Logging decision variable
        REAL(dp), intent(in):: Mpr, Vpr,Vdta, mu_crit, mu_dta !Prior mean, Prior variance, data variance, critical value of data for logging, values of the monitoring data -- a vector,
        REAL(dp), intent(in):: lsrm, lsrvar, grm, grvar !logging survival rate mean, logging survival rate variance, growth rate mean, growth rate variance
        REAL(dp), intent(out):: Mpo, Vpo ! Posterior mean, posterior variance -- both vectors
       
        DIMENSION:: mu_dta(len_mu_dta), Mpo(len_mu_dta), Vpo(len_mu_dta) 

        !Local variables
        INTEGER:: i
        
    !Vectorised form of posterior mean and variance
    Vpo= (1._dp/Vpr+nj*1._dp/Vdta)**(-1._dp)
    Mpo= (Mpr+mu_dta*Vpr*(1._dp*nj)/Vdta)/(1._dp+Vpr*(1._dp*nj)/Vdta) !Note - this is the usual formula divided by Vdta/nj -- done to prevent division by zero.

    !If logging is intented, and it occurs, apply logging effects.
        IF(LDj.eq.1) THEN
                DO i = 1, len_mu_dta
                        !Check if logging occured - and apply growth if so
                        IF(mu_dta(i).ge.mu_crit) THEN
                        Mpo(i)= Mpo(i)+lsrm        
                        Vpo(i)= Vpo(i)+lsrvar
                        END IF
                END DO
        END IF
        
        !Apply population growth
        Mpo= Mpo+grm
        Vpo= Vpo+grvar
       
        return

END SUBROUTINE mu_sigma_jp1_f

!!This subroutine outputs the integrand in the SDP equation where we calculate
!the total return
SUBROUTINE useful_f(mu_dta,len_mu_dta, Mpr, Vpr,Vdta, nj, mu_crit,LDj, TR_jp1, lsrm, &
 lsrvar, grm, grvar, vals_mu_dta, Mbb, len_Mbb, Varbb, len_Varbb)
        implicit none
        INTEGER, INTENT(in):: len_mu_dta, nj, LDj, len_Mbb, len_Varbb !length of mu_dta, Monitoring effort, Logging decision, length of Mbb, length of Varbb
        REAL(dp), INTENT(in):: Mpr, Vpr, Vdta, mu_dta, mu_crit, lsrm, lsrvar, grm, grvar, Mbb, Varbb, TR_jp1 !Prior Mean, Prior variance, data variance, data mean, critical value of data mean, logging survival rate mean, logging survival rate variance, growth rate mean, growth rate variance, Mean bin boundaries, Variance bin boundaries
        REAL(dp), INTENT(out):: vals_mu_dta !The output values -- the integrand at mu_dta
        Dimension:: mu_dta(len_mu_dta), vals_mu_dta(len_mu_dta), Mbb(len_Mbb), Varbb(len_Varbb)
        Dimension:: TR_jp1((len_Mbb-1), (len_Varbb-1))

        !Local variables
        REAL(dp):: Mpo(len_mu_dta), Vpo(len_mu_dta), pr_out(len_mu_dta) !Posterior mean, Posterior variance, probability density at mu_dta
        INTEGER:: Mu_jp1(len_mu_dta), Var_jp1(len_mu_dta) !The index associating the posterior mean and variance with a value in Means and Varis
        INTEGER:: kk, k1, k2, i
        !Calculate posterior mean and variance 
        call mu_sigma_jp1_f(Mpr, Vpr, Vdta, nj, mu_dta, len_mu_dta, mu_crit, LDj, Mpo, Vpo, lsrm, lsrvar, grm, grvar)

        !Associate the posterior mean and variance with a point in Means and Varis
        !Here we mimic the function of findInterval in R. So mu_jp1 and Var_jp1
        !will contain index's
        call findInterval_f(Mpo, len_mu_dta, Mbb, len_Mbb, Mu_jp1)
        call findInterval_f(Vpo, len_mu_dta, Varbb, len_Varbb, Var_jp1)
        !Calculate the probability density of different data
        call pr_mu_dta_f(Mpr, Vpr, Vdta, nj,mu_dta,len_mu_dta, pr_out)
       
        DO i = 1, len_mu_dta
        vals_mu_dta(i) = TR_jp1(Mu_jp1(i), Var_jp1(i))*pr_out(i)
        END DO 

        return
END SUBROUTINE useful_f

!!This subroutine outputs the integrand in the SDP equation where we calculate
!the total return -- slight variation on above
SUBROUTINE useful_f2(mu_dta,len_mu_dta, Mpr, Vpr,Vdta, nj, mu_crit,LDj, TR_jp1, lsrm, &
 lsrvar, grm, grvar, vals_mu_dta, Mbb, len_Mbb, Varbb, len_Varbb, Meanz, Variz, md)
        implicit none
        INTEGER, INTENT(in):: len_mu_dta, nj, LDj, len_Mbb, len_Varbb, md !length of mu_dta, Monitoring effort, Logging decision, length of Mbb, length of Varbb
        REAL(dp), INTENT(in):: Mpr, Vpr, Vdta, mu_dta, mu_crit, lsrm, lsrvar, grm, grvar, Mbb, Varbb, TR_jp1, Meanz, Variz !Prior Mean, Prior variance, data variance, data mean, critical value of data mean, logging survival rate mean, logging survival rate variance, growth rate mean, growth rate variance, Mean bin boundaries, Variance bin boundaries, Total_return, Means, Variances
        REAL(dp), INTENT(out):: vals_mu_dta !The output values -- the integrand at mu_dta
        Dimension:: mu_dta(len_mu_dta), vals_mu_dta(len_mu_dta), Mbb(len_Mbb), Varbb(len_Varbb)
        Dimension:: TR_jp1((len_Mbb-1), (len_Varbb-1)), Meanz(len_Mbb-1), Variz(len_Varbb-1)

        !Local variables
        REAL(dp):: Mpo(len_mu_dta), Vpo(len_mu_dta), pr_out(len_mu_dta) !Posterior mean, Posterior variance, probability density at mu_dta
        INTEGER:: Mu_jp1(len_mu_dta), Var_jp1(len_mu_dta) !The index associating the posterior mean and variance with a value in Means and Varis
        INTEGER:: kk, k1, k2, i, m1, v1, jj1, jj2,ier
        REAL(dp):: dtdx, dtdy,dt2dxdy, eps, d1, d2, d3, d4, TR_top, TR_bot, a1, a2, a3, TR_jp1_tmp((len_Mbb-1), (len_Varbb-1))
        Logical:: bi_cubic=.FALSE., bi_linear=.TRUE.

        
        !Calculate posterior mean and variance 
        call mu_sigma_jp1_f(Mpr, Vpr, Vdta, nj, mu_dta, len_mu_dta, mu_crit, LDj, Mpo, Vpo, lsrm, lsrvar, grm, grvar)

        !Associate the posterior mean and variance with a point in Means and Varis
        !Here we mimic the function of findInterval in R. So mu_jp1 and Var_jp1
        !will contain index's
        call findInterval_f(Mpo, len_mu_dta, Mbb, len_Mbb, Mu_jp1)
        call findInterval_f(Vpo, len_mu_dta, Varbb, len_Varbb, Var_jp1)
        !Calculate the probability density of different data
        call pr_mu_dta_f(Mpr, Vpr, Vdta, nj,mu_dta,len_mu_dta, pr_out)

        TR_jp1_tmp=TR_jp1 !Predefine TR_jp1_tmp, to get around the fact that it appears as intent(inout) in rgbi3p
        IF(bi_cubic) THEN 
                !Perform bicubic spline interpolation to get TR_jp1 at the values
                !Mpo,Vpo
                call rgbi3p(md,len_Mbb-1, len_Varbb-1, Meanz, Variz, TR_jp1_tmp, len_mu_dta, Mpo,Vpo, vals_mu_dta,ier)
                IF(ier.ne.0) THEN
                        print*, 'Error in bicubic interpolation, ier = ', ier 
                !        stop ! Meanz
                ELSE
                !       print*, 'Bicubic okay'
                END IF
               
                !Get the integrand values
                vals_mu_dta=vals_mu_dta*pr_out
                RETURN
        END IF


        !THE FOLLOWING CODE needs to be cleaned up. It includes treatments of
        !inverse distance weighted averages using various different points, and bilinear interpolation
 
        DO i = 1, len_mu_dta
        !vals_mu_dta(i) = TR_jp1(Mu_jp1(i), Var_jp1(i))*pr_out(i)
         
                !Interpolate the value of TR_jp1, except on the boundaries
                !Boundary case 
                IF((((Mu_jp1(i).eq.len_Mbb-1).or.(Var_jp1(i).eq. len_Varbb-1)).or. & 
                                                        ((Mu_jp1(i).eq.1).or.(Var_jp1(i).eq.1)))) THEN
                        vals_mu_dta(i) = TR_jp1(Mu_jp1(i), Var_jp1(i))*pr_out(i)
                ELSE

 
        !               !Krigging type estimate of TR_jp1 -- calculate distances, then
        !               !use inverse distance weighting to get our estimate
                        eps=10.0_dp**(-6)

                        !Interior case - this will tell us if we are less than
                        !or greater than the nearby Mean/Variance point
                        !print*, Mu_jp1(i)
                        IF(Mpo(i)>Meanz(Mu_jp1(i)))THEN
                                m1=1
                        ELSE
                                m1=-1
                        END IF
                        
                        IF(Vpo(i)>Variz(Var_jp1(i))) THEN
                                v1=1
                        ELSE
                                v1=-1
                        END IF
                        
                        !Try bi-linear interpolation
                        IF(bi_linear) THEN
                                !Calculate distances for bi-linear weighted
                                !averages
                                d1=abs(Mpo(i) - Meanz(Mu_jp1(i)))
                                d2=abs(Mpo(i) - Meanz(Mu_jp1(i)+m1))
                                d3=abs(Vpo(i) - Variz(Var_jp1(i)))
                                d4=abs(Vpo(i) - Variz(Var_jp1(i)+v1))
                                !Weighted average on Var_jp1(i) side
                                a1=( TR_jp1(Mu_jp1(i),Var_jp1(i))*d2 +  TR_jp1(Mu_jp1(i)+m1,Var_jp1(i))*d1)/(d2+d1)
                                !Weighted average on Var_jp1(i)+v1 side
                                a2=(TR_jp1(Mu_jp1(i),Var_jp1(i)+v1)*d2 +  TR_jp1(Mu_jp1(i)+m1,Var_jp1(i)+v1)*d1)/(d2+d1)
                                !Final interpolated value
                                a3=(a1*d4+a2*d3)/(d4+d3)
                                !Final value
                                vals_mu_dta(i)=a3*pr_out(i)
                        ELSE
                                !If we are only one point inside the boundary, then use a simple
                                !four point inverse distance weighted average estimate
                                IF((((Mu_jp1(i).eq.len_Mbb-2).or.(Var_jp1(i).eq. len_Varbb-2)).or. & 
                                                                        ((Mu_jp1(i).eq.2).or.(Var_jp1(i).eq.2)))) THEN
                                
                                        d1=( (Mpo(i)-Meanz(Mu_jp1(i)))**2 + (Vpo(i)-Variz(Var_jp1(i)))**2 )**0.5_dp +eps
                                        d2=( (Mpo(i)-Meanz(Mu_jp1(i)+m1))**2 + (Vpo(i)-Variz(Var_jp1(i)))**2 )**0.5_dp+eps
                                        d3=( (Mpo(i)-Meanz(Mu_jp1(i)))**2 + (Vpo(i)-Variz(Var_jp1(i)+v1))**2 )**0.5_dp+eps
                                        d4=( (Mpo(i)-Meanz(Mu_jp1(i)+m1))**2 + (Vpo(i)-Variz(Var_jp1(i)+v1))**2 )**0.5_dp+eps
                                       
                                        !Inverse weighted average value of TR_jp1
                                        vals_mu_dta(i)= pr_out(i)*( TR_jp1(Mu_jp1(i),Var_jp1(i))/d1 + &
                                                                TR_jp1(Mu_jp1(i)+m1,Var_jp1(i))/d2 + &
                                                                TR_jp1(Mu_jp1(i),Var_jp1(i)+v1)/d3 + &
                                                                TR_jp1(Mu_jp1(i)+m1,Var_jp1(i)+v1)/d4 &
                                                                )/(1._dp/d1+1._dp/d2+1._dp/d3+1._dp/d4)


                                ELSE
                                !In this case we are 2 or more points inside the boundary, so we use a nine point krigging estimate
                                        d1=0._dp
                                        TR_top=0._dp !This will be the numerator in the krigging estimate
                                        TR_bot=0._dp !This will be the denominator in the krigging estimate
                                        DO jj1= -1,1 !Loop through means in neighbourhood
                                                DO jj2=-1,1 !Loop through variances in neighbourhood
                                                d1=( (Mpo(i)-Meanz(Mu_jp1(i)+jj1))**2 + & 
                                                        (Vpo(i)-Variz(Var_jp1(i)+jj2))**2 )**0.5_dp +eps

                                                TR_top=TR_top+TR_jp1(Mu_jp1(i)+jj1, Var_jp1(i)+jj2)/d1
                                                TR_bot=TR_bot+1._dp/d1

                                                END DO
                                        END DO
                                        vals_mu_dta(i)=pr_out(i)*TR_top/TR_bot !Note that this is a nine krigging estimate point extension 

                                END IF
                        END IF !End IF(bi_linear)

                END IF
         

        END DO 


        return
END SUBROUTINE useful_f2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This turns useful_f2 into a function -- good for some adaptive integration
! routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FUNCTION useful_f2_func(mu_dta,len_mu_dta, Mpr, Vpr,Vdta, nj, mu_crit,LDj, TR_jp1, lsrm, &
! lsrvar, grm, grvar, vals_mu_dta, Mbb, len_Mbb, Varbb, len_Varbb, Meanz, Variz, md)
!        implicit none
!        INTEGER, INTENT(in):: len_mu_dta, nj, LDj, len_Mbb, len_Varbb, md !length of mu_dta, Monitoring effort, Logging decision, length of Mbb, length of Varbb
!        REAL(dp), INTENT(in):: Mpr, Vpr, Vdta, mu_dta, mu_crit, lsrm, lsrvar, grm, grvar, Mbb, Varbb, TR_jp1, Meanz, Variz !Prior Mean, Prior variance, data variance, data mean, critical value of data mean, logging survival rate mean, logging survival rate variance, growth rate mean, growth rate variance, Mean bin boundaries, Variance bin boundaries, Total_return, Means, Variances
!        REAL(dp), INTENT(out):: vals_mu_dta !The output values -- the integrand at mu_dta
!        Dimension:: mu_dta(len_mu_dta), vals_mu_dta(len_mu_dta), Mbb(len_Mbb), Varbb(len_Varbb)
!        Dimension:: TR_jp1((len_Mbb-1), (len_Varbb-1)), Meanz(len_Mbb-1), Variz(len_Varbb-1)
!
!call useful_f2(mu_dta,len_mu_dta, Mpr, Vpr,Vdta, nj, mu_crit,LDj, TR_jp1, lsrm, &
! lsrvar, grm, grvar, vals_mu_dta, Mbb, len_Mbb, Varbb, len_Varbb, Meanz, Variz, md)
!
!useful_f2_func=vals_mu_dta
!
!return 
!END FUNCTION useful_f2_func

!This mimics, in a limited way, the subroutine findInterval in R
SUBROUTINE findInterval_f(Vec1, len_Vec1, Vec2, len_Vec2, Index_vec1)
        IMPLICIT NONE
        INTEGER, INTENT(in):: len_Vec1, len_Vec2 !The length of the first and second vectors.
        REAL(dp), INTENT(in):: Vec1, Vec2 !The vector to be binned, the vector defining the bins. BOTH SHOULD BE SORTED - DIFFERENT TO THE R VERSION
        INTEGER, INTENT(out):: Index_vec1 !The indexes showing the bin of Vec1 in Vec2

        Dimension::Vec1(len_Vec1), Vec2(len_Vec2), Index_vec1(len_Vec1)
       !local variables
        INTEGER:: kk, k1, i       

 
        k1 = 0 !Predefine - this will store the index in Vec2 that the value in Vec1(i) is greater than Vec2[k1]
        DO i = 1, len_Vec1
                !NOTICE THAT THE BINS ARE SHAPED LIKE [ ), same as findInterval. 
                IF(Vec1(i) <= Vec2(1)) THEN !This mimics the all.inside=T functionality of findInterval
                        Index_Vec1(i)=1
                        goto 3456
                END IF
                
                IF(Vec1(i) >= Vec2(len_Vec2)) THEN
                        Index_Vec1(i)=len_Vec2-1
                        k1=len_Vec2
                        goto 3456
                END IF
                !We use k1 to reduce the sorting time -- however it only works
                !for Vec1 increasing. Due to the logging threshold and
                !mortality, this is not always true. Thus, if Vec1 is not increasing, reset k1
                IF(i>1) THEN
                IF(Vec1(i)<Vec1(i-1)) k1=0
                END IF

                !Scroll through Vec2, and find the index that Vec1(i) is between
                DO kk = max(k1,1), (len_Vec2-1)
                        IF((Vec1(i) >= Vec2(kk) ).and.(Vec1(i)< Vec2(kk+1))) THEN
                                Index_Vec1(i)=kk
                                k1=kk
                                goto 3456
                        END IF
                END DO !kk
                print*, 'Missing something', k1, len_Vec1, len_Vec2, Vec1((len_Vec1-1):len_Vec1), Vec2(len_Vec2-1)
                stop

                3456 continue
                        
        END DO !i

        return

END SUBROUTINE findInterval_f

!!!!!!!!!!!!Subroutine to compute the trapozoidal rule for approximating an integral
SUBROUTINE trapz_integral_f(x,y,len_x, output)
        IMPLICIT NONE
        INTEGER, INTENT(IN):: len_x
        REAL(dp), INTENT(IN):: x,y
        REAL(dp), INTENT(OUT):: output
        DIMENSION:: x(len_x), y(len_x)
       
        !Local variables
        INTEGER:: i
        output=0._dp !Predefine
        !This should be clear
        DO i=2,(len_x-1)
        output=output+y(i)*(x(i+1)-x(i-1))*0.5_dp
        END DO
        output=output+y(1)*0.5_dp*(x(2)-x(1))
        output=output+y(len_x)*0.5_dp*(x(len_x)-x(len_x-1))

END SUBROUTINE trapz_integral_f

!!!!!!!!!!!!Subroutine to compute the trapozoidal rule for approximating an integral
SUBROUTINE simps_integral_f(x,y,len_x, output)
        IMPLICIT NONE
        INTEGER, INTENT(IN):: len_x
        REAL(dp), INTENT(IN):: x,y
        REAL(dp), INTENT(OUT):: output
        DIMENSION:: x(len_x), y(len_x)
       
        !Local variables
        INTEGER:: i
        REAL(dp):: tmp
        
        !Check that the number of intervals is okay
        IF(mod(len_x,2).EQ.0) THEN
                PRINT*, 'number of INTERVALS in simpson rule integration should be even -- so & 
                       & we need an ODD number of function evaluation points'
                stop
        END IF

        !Implement Simpsons rule: Integral = dx/3*[ first + last + 2*(interior
        !odds) + 4*(interior_evens)] ;  where the indexes are 1,2,..len_x, the
        !interior_evens are 2,4,6,...,(len_x-1), and the interior_odds are
        !3,5,7,..len_x-2 -- note that my indexing is different to the wikipedia
        !page
        output=y(1)+y(len_x) ! First + last
        tmp=0._dp
        DO i= 1, (len_x-1)/2
        tmp=tmp+y(2*i) !Sum of evens
        END DO
        output=output+4._dp*tmp ! + 4* interior evens

        tmp=0._dp
        DO i= 2, (len_x-1)/2
        tmp=tmp+y(2*i-1) !Sum of odds
        END DO
        output=output+2._dp*tmp !+ 2* interior odds
       
        output=output*(x(2)-x(1))/3.0_dp 

END SUBROUTINE simps_integral_f


!!!!!!!!!THE following subroutine calculates the optimal monitoring/logging strategy and return, given Mpr, Vpr,
SUBROUTINE TR_n_LD_loop(Mpr, Vpr, Vdta, nj, mu_crit,Pr_logging, sizen, tt, NTS, & 
                                int_nos, lsrm, lsrvar, grm, grvar, Lr, Mc, TR_n_LD, & 
                                Mbb, len_Mbb, Varbb, len_Varbb, Meanz, Variz, TR_jp1, writer)        
        IMPLICIT NONE
        INTEGER, INTENT(IN):: nj, sizen, tt, NTS, int_nos, len_Mbb, len_Varbb
        REAL(dp), INTENT(IN):: Mpr, Vpr, Vdta, mu_crit, Pr_logging, Lr, Mc, lsrm, lsrvar, grm, grvar
        REAL(dp), INTENT(IN):: Mbb, Varbb, Meanz, Variz, TR_jp1
        REAL(dp), INTENT(OUT):: TR_n_LD
        LOGICAL, INTENT(IN)::writer
        DIMENSION:: Pr_logging(sizen), mu_crit(sizen), nj(sizen), TR_n_LD(sizen, 2), Mbb(len_Mbb), &
                        Varbb(len_Varbb), Meanz(len_Mbb-1), Variz(len_Varbb-1), TR_jp1(len_Mbb-1, len_Varbb-1)
        !Local Variables
        REAL(dp):: bigint, ats(int_nos), fzzs(int_nos), ats_low, ats_high, ats_range=5._dp, errr, pr_tmp(int_nos)
        REAL(dp):: tmp1, tmp2
        INTEGER:: LDj, nn, i,md, bb,bb2(1)
        fzzs=0._dp
        md=1 !This is useful in the bicubic interpolation.
        DO LDj = 0,1 !Loop over both logging strategies
                DO nn = 1,sizen !Loop over all monitoring efforts
                        IF (tt.ne.NTS) THEN !Unless we are on the last time step, we need to integrate over future time-steps to get the optimal strategy

                                !The following values define the domain over
                                !which we will do the integration.
                                IF(nj(nn)>0) THEN
                                        ats_low=Mpr-ats_range*sqrt(Vpr+Vdta/max(nj(nn)*1._dp,1._dp))
                                        ats_high=Mpr+ats_range*sqrt(Vpr+Vdta/max(nj(nn)*1._dp,1._dp))
                                ELSE
                                        ats_low=Mpr-ats_range*sqrt(Vpr)
                                        ats_high=Mpr+ats_range*sqrt(Vpr)
                                END IF
                                DO i=1,int_nos !ats are the points we integrate over. int_nos is the number of points that we integrate over
                                ats(i)= ats_low+ (ats_high-ats_low)/(1._dp*(int_nos-1))*(1._dp*(i-1))
                                END DO
                                
                                !do i=1,int_nos-1
                                !if(ats(i)>ats(i+1)) THEN
                                !print*, 'ats is not sorted', i, ats(i:i+1), ats_low, ats_high 
                                !END IF
                                !end do
                                !Now we evaluate the integrand at ats, using the
                                !useful_f2 function
                                call useful_f2(ats,int_nos, Mpr, Vpr,Vdta, nj(nn), mu_crit(nn), LDj, TR_jp1, lsrm, &
                                                 lsrvar, grm, grvar, fzzs, Mbb, len_Mbb, Varbb, len_Varbb, Meanz, Variz,md)
                                md=2
                                !fzzs == the TR*Pr values which we must integrate over
                                
                                !IF((writer).AND.(nj(nn)<20).AND.(LDj==1)) THEN
                                !        print*, 'writing'
                                !        call pr_mu_dta_f(Mpr, Vpr, Vdta, nj(nn),ats,int_nos, pr_tmp)
                                !        write(3,*) fzzs !/pr_tmp, minloc(abs(ats-mu_crit(nn))) 
                                !END IF

                                !Calculate integral, accounting for
                                !discontinuity
                                IF(maxval(fzzs)>0._dp) THEN
                                        bb2=minloc(abs(ats-mu_crit(nn))) 
                                        bb=bb2(1)
                                        if(ats(bb)>mu_crit(nn)) bb=bb-1 ! ats(bb) is where the discontinuity in the integral will be

                                                !call trapz_integral_f(ats,fzzs,int_nos, bigint) !So bigint is the value of the integral
                                                !call simps_integral_f(ats,fzzs,int_nos, bigint) !So bigint is the value of the integral
                                                IF((bb>4).AND.(bb<int_nos-4)) THEN
                                                        call cubint(bb,ats(1:bb),fzzs(1:bb),1,bb,tmp1,errr) !Integrate First continuous section
                                                        call cubint(int_nos - bb,ats(bb+1:int_nos),fzzs(bb+1:int_nos),&
                                                                           1,int_nos-bb,tmp2,errr) !Integrate second continuous section
                                                        bigint=tmp1+tmp2 +  & 
                                                          (mu_crit(nn)-ats(bb))*(fzzs(bb))+(ats(bb+1)-mu_crit(nn))*fzzs(bb+1) !Total integral
                                                ELSE
                                                        call cubint(int_nos,ats,fzzs,1,int_nos,bigint,errr) !Integrate First continuous section
                                                END IF
                                        !call wedint(int_nos,ats(2)-ats(1),fzzs,bigint)
                                        !call simp(useful_f2_func, ats_low,ats_high, 0.001, bigint)
                                ELSE
                                        bigint=0._dp
                                END IF
                        ELSE                     
                        bigint=0._dp
                        END IF 
                !Calculate the value of the total return, and store it
                TR_n_LD(nn,LDj+1)= bigint + Pr_logging(nn)*(LDj*Lr) - (nj(nn)*Mc)*1._dp                            

                END DO
        END DO
!print*, 'par_fortran', bigint
END SUBROUTINE TR_n_LD_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGINNING OF THE BICUBIC INTERPOLATION FUNCTIONS, sourced from
! toms760.f90 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rgbi3p(md, nxd, nyd, xd, yd, zd, nip, xi, yi, zi, ier)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-11  Time: 10:11:03

! Rectangular-grid bivariate interpolation
! (a master subroutine of the RGBI3P/RGSF3P subroutine package)

! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 1995/08

! This subroutine performs interpolation of a bivariate function, z(x,y), on a
! rectangular grid in the x-y plane.  It is based on the revised Akima method.

! In this subroutine, the interpolating function is a piecewise function
! composed of a set of bicubic (bivariate third-degree) polynomials, each
! applicable to a rectangle of the input grid in the x-y plane.
! Each polynomial is determined locally.

! This subroutine has the accuracy of a bicubic polynomial, i.e., it
! interpolates accurately when all data points lie on a surface of a
! bicubic polynomial.

! The grid lines can be unevenly spaced.

! The input arguments are
!   MD  = mode of computation
!       = 1 for new XD, YD, or ZD data (default)
!       = 2 for old XD, YD, and ZD data,
!   NXD = number of the input-grid data points in the x coordinate
!         (must be 2 or greater),
!   NYD = number of the input-grid data points in the y coordinate
!         (must be 2 or greater),
!   XD  = array of dimension NXD containing the x coordinates of the
!         input-grid data points (must be in a monotonic increasing order),
!   YD  = array of dimension NYD containing the y coordinates of the
!         input-grid data points (must be in a monotonic increasing order),
!   ZD  = two-dimensional array of dimension NXD*NYD
!         containing the z(x,y) values at the input-grid data points,
!   NIP = number of the output points at which interpolation
!         of the z value is desired (must be 1 or greater),
!   XI  = array of dimension NIP containing the x coordinates
!         of the output points,
!   YI  = array of dimension NIP containing the y coordinates
!         of the output points.

! The output arguments are
!   ZI  = array of dimension NIP where the interpolated z
!         values at the output points are to be stored,
!   IER = error flag
!       = 0 for no errors
!       = 1 for NXD = 1 or less
!       = 2 for NYD = 1 or less
!       = 3 for identical XD values or XD values out of sequence
!       = 4 for identical YD values or YD values out of sequence
!       = 5 for NIP = 0 or less.

! N.B. The workspace has been removed from the argument list.
!   WK  = three dimensional array of dimension 3*NXD*NYD used internally
!         as a work area.

! The very fisrt call to this subroutine and the call with a new XD, YD, and
! ZD array must be made with MD=1.  The call with MD=2 must be preceded by
! another call with the same XD, YD, and ZD arrays.  Between the call with
! MD=2 and its preceding call, the WK array must not be disturbed.

! The constant in the PARAMETER statement below is
!   NIPIMX = maximum number of output points to be processed at a time.
! The constant value has been selected empirically.

! This subroutine calls the RGPD3P, RGLCTN, and RGPLNL subroutines.


! Specification statements
!     .. Parameters ..

INTEGER, INTENT(IN)   :: md
INTEGER, INTENT(IN)   :: nxd
INTEGER, INTENT(IN)   :: nyd
REAL(dp), INTENT(IN)      :: xd(nxd)
REAL(dp), INTENT(IN)      :: yd(nyd)
REAL(dp), INTENT(IN OUT)  :: zd(nxd,nyd)
INTEGER, INTENT(IN)   :: nip
REAL(dp), INTENT(IN OUT)  :: xi(nip)
REAL(dp), INTENT(IN OUT)  :: yi(nip)
REAL(dp), INTENT(IN OUT)  :: zi(nip)
INTEGER, INTENT(OUT)  :: ier

!     ..
!     .. Local Scalars ..
INTEGER, PARAMETER  :: nipimx=51

INTEGER  :: iip, ix, iy, nipi
!     ..
!     .. Local Arrays ..
INTEGER  :: inxi(nipimx), inyi(nipimx)

! Allocate workspace
REAL(dp)  :: wk(3,nxd,nyd)
!     ..
!     .. External Subroutines ..
! EXTERNAL         rglctn, rgpd3p, rgplnl
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC        MIN
!     ..

! Preliminary processing
! Error check
IF (nxd <= 1) GO TO 40
IF (nyd <= 1) GO TO 50
DO  ix = 2,nxd
  IF (xd(ix) <= xd(ix-1)) THEN
         print*, ix,nxd,nyd, xd
         GO TO 60
  END IF
END DO
DO  iy = 2,nyd
  IF (yd(iy) <= yd(iy-1)) GO TO 70
END DO
IF (nip <= 0) GO TO 80
ier = 0

! Calculation
! Estimates partial derivatives at all input-grid data points (for MD=1).
IF (md /= 2) THEN
  CALL rgpd3p(nxd, nyd, xd, yd, zd, wk)
END IF

! DO-loop with respect to the output point
! Processes NIPIMX output points, at most, at a time.
DO  iip = 1,nip,nipimx
  nipi = MIN(nip- (iip-1),nipimx)
! Locates the output points.
  CALL rglctn(nxd, nyd, xd, yd, nipi, xi(iip), yi(iip), inxi, inyi)

! Calculates the z values at the output points.
  CALL rgplnl(nxd, nyd, xd, yd, zd, wk, nipi, xi(iip), yi(iip), inxi, inyi, &
              zi(iip))
END DO
RETURN

! Error exit
40 WRITE (*,FMT=9000)
ier = 1
GO TO 90
50 WRITE (*,FMT=9010)
ier = 2
GO TO 90
60 WRITE (*,FMT=9020) ix,xd(ix)
ier = 3
GO TO 90
70 WRITE (*,FMT=9030) iy,yd(iy)
ier = 4
GO TO 90
80 WRITE (*,FMT=9040)
ier = 5
90 WRITE (*,FMT=9050) nxd,nyd,nip
RETURN

! Format statements for error messages
9000 FORMAT (/' *** RGBI3P Error 1: NXD = 1 or less')
9010 FORMAT (/' *** RGBI3P Error 2: NYD = 1 or less')
9020 FORMAT (/' *** RGBI3P Error 3: Identical XD values or',  &
    ' XD values out of sequence'/ '    IX =', i6, ',  XD(IX) =', e11.3)
9030 FORMAT (/' *** RGBI3P Error 4: Identical YD values or',  &
    ' YD values out of sequence',/,'    IY =',i6,',  YD(IY) =', e11.3)
9040 FORMAT (/' *** RGBI3P Error 5: NIP = 0 or less')
9050 FORMAT ('    NXD =', i5,',  NYD =', i5,',  NIP =', i5/)
END SUBROUTINE rgbi3p



SUBROUTINE rgsf3p(md, nxd, nyd, xd, yd, zd, nxi, xi, nyi, yi, zi, ier)

! Rectangular-grid surface fitting
! (a master subroutine of the RGBI3P/RGSF3P subroutine package)

! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 1995/08

! This subroutine performs surface fitting by interpolating values of a
! bivariate function, z(x,y), on a rectangular grid in the x-y plane.
! It is based on the revised Akima method.

! In this subroutine, the interpolating function is a piecewise function
! composed of a set of bicubic (bivariate third-degree) polynomials, each
! applicable to a rectangle of the input grid in the x-y plane.
! Each polynomial is determined locally.

! This subroutine has the accuracy of a bicubic polynomial, i.e., it fits the
! surface accurately when all data points lie on a surface of a bicubic
! polynomial.

! The grid lines of the input and output data can be unevenly spaced.

! The input arguments are
!   MD  = mode of computation
!       = 1 for new XD, YD, or ZD data (default)
!       = 2 for old XD, YD, and ZD data,
!   NXD = number of the input-grid data points in the x
!         coordinate (must be 2 or greater),
!   NYD = number of the input-grid data points in the y
!         coordinate (must be 2 or greater),
!   XD  = array of dimension NXD containing the x coordinates
!         of the input-grid data points (must be in a
!         monotonic increasing order),
!   YD  = array of dimension NYD containing the y coordinates
!         of the input-grid data points (must be in a
!         monotonic increasing order),
!   ZD  = two-dimensional array of dimension NXD*NYD
!         containing the z(x,y) values at the input-grid data points,
!   NXI = number of output grid points in the x coordinate
!         (must be 1 or greater),
!   XI  = array of dimension NXI containing the x coordinates
!         of the output grid points,
!   NYI = number of output grid points in the y coordinate
!         (must be 1 or greater),
!   YI  = array of dimension NYI containing the y coordinates
!         of the output grid points.

! The output arguments are
!   ZI  = two-dimensional array of dimension NXI*NYI, where the interpolated
!         z values at the output grid points are to be stored,
!   IER = error flag
!       = 0 for no error
!       = 1 for NXD = 1 or less
!       = 2 for NYD = 1 or less
!       = 3 for identical XD values or XD values out of sequence
!       = 4 for identical YD values or YD values out of sequence
!       = 5 for NXI = 0 or less
!       = 6 for NYI = 0 or less.

! N.B. The workspace has been removed from the argument list.
!   WK  = three-dimensional array of dimension 3*NXD*NYD used internally
!         as a work area.

! The very first call to this subroutine and the call with a new XD, YD, or
! ZD array must be made with MD=1.  The call with MD=2 must be preceded by
! another call with the same XD, YD, and ZD arrays.  Between the call with
! MD=2 and its preceding call, the WK array must not be disturbed.

! The constant in the PARAMETER statement below is
!   NIPIMX = maximum number of output points to be processed at a time.
! The constant value has been selected empirically.

! This subroutine calls the RGPD3P, RGLCTN, and RGPLNL subroutines.


! Specification statements
!     .. Parameters ..

INTEGER, INTENT(IN)   :: md
INTEGER, INTENT(IN)   :: nxd
INTEGER, INTENT(IN)   :: nyd
REAL(dp), INTENT(IN)      :: xd(nxd)
REAL(dp), INTENT(IN)      :: yd(nyd)
REAL(dp), INTENT(IN OUT)  :: zd(nxd,nyd)
INTEGER, INTENT(IN)   :: nxi
REAL(dp), INTENT(IN OUT)  :: xi(nxi)
INTEGER, INTENT(IN)   :: nyi
REAL(dp), INTENT(IN)      :: yi(nyi)
REAL(dp), INTENT(IN OUT)  :: zi(nxi,nyi)
INTEGER, INTENT(OUT)  :: ier

!     ..
!     .. Local Scalars ..
INTEGER, PARAMETER  :: nipimx=51

INTEGER  :: ix, ixi, iy, iyi, nipi
!     ..
!     .. Local Arrays ..
REAL(dp)     :: yii(nipimx)
INTEGER  :: inxi(nipimx), inyi(nipimx)

! Allocate workspace
REAL(dp)  :: wk(3,nxd,nyd)
!     ..
!     .. External Subroutines ..
! EXTERNAL         rglctn,rgpd3p,rgplnl
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC        MIN
!     ..

! Preliminary processing
! Error check
IF (nxd <= 1) GO TO 60
IF (nyd <= 1) GO TO 70
DO  ix = 2,nxd
  IF (xd(ix) <= xd(ix-1)) GO TO 80
END DO
DO  iy = 2,nyd
  IF (yd(iy) <= yd(iy-1)) GO TO 90
END DO
IF (nxi <= 0) GO TO 100
IF (nyi <= 0) GO TO 110
ier = 0

! Calculation
! Estimates partial derivatives at all input-grid data points
! (for MD=1).
IF (md /= 2) THEN
  CALL rgpd3p(nxd, nyd, xd, yd, zd, wk)
END IF

! Outermost DO-loop with respect to the y coordinate of the output grid points
DO  iyi = 1,nyi
  DO  ixi = 1,nipimx
    yii(ixi) = yi(iyi)
  END DO

! Second DO-loop with respect to the x coordinate of the output grid points
! Processes NIPIMX output-grid points, at most, at a time.
  DO  ixi = 1,nxi,nipimx
    nipi = MIN(nxi- (ixi-1), nipimx)
! Locates the output-grid points.
    CALL rglctn(nxd, nyd, xd, yd, nipi, xi(ixi), yii, inxi, inyi)

! Calculates the z values at the output-grid points.
    CALL rgplnl(nxd, nyd, xd, yd, zd, wk, nipi, xi(ixi), yii, inxi, inyi,  &
                zi(ixi,iyi))
  END DO
END DO
RETURN

! Error exit
60 WRITE (*,FMT=9000)
ier = 1
GO TO 120
70 WRITE (*,FMT=9010)
ier = 2
GO TO 120
80 WRITE (*,FMT=9020) ix,xd(ix)
ier = 3
GO TO 120
90 WRITE (*,FMT=9030) iy,yd(iy)
ier = 4
GO TO 120
100 WRITE (*,FMT=9040)
ier = 5
GO TO 120
110 WRITE (*,FMT=9050)
ier = 6
120 WRITE (*,FMT=9060) nxd,nyd,nxi,nyi
RETURN

! Format statements for error messages
9000 FORMAT (/' *** RGSF3P Error 1: NXD = 1 or less')
9010 FORMAT (/' *** RGSF3P Error 2: NYD = 1 or less')
9020 FORMAT (/' *** RGSF3P Error 3: Identical XD values or',  &
    ' XD values out of sequence',/,'    IX =',i6,',  XD(IX) =', e11.3)
9030 FORMAT (/' *** RGSF3P Error 4: Identical YD values or',  &
    ' YD values out of sequence',/,'    IY =',i6,',  YD(IY) =', e11.3)
9040 FORMAT (/' *** RGSF3P Error 5: NXI = 0 or less')
9050 FORMAT (/' *** RGSF3P Error 6: NYI = 0 or less')
9060 FORMAT ('    NXD =', i5, ',  NYD =', i5, ',  NXI =', i5,',  NYI =', i5 /)
END SUBROUTINE rgsf3p



!     ..
! Statement Function definitions
! z2f(xx1,xx2,zz0,zz1) = (zz1-zz0)*xx2/xx1 + zz0
! z3f(xx1,xx2,xx3,zz0,zz1,zz2) = ((zz2-zz0)* (xx3-xx1)/xx2 -  &
!    (zz1-zz0)* (xx3-xx2)/xx1)* (xx3/ (xx2-xx1)) + zz0

FUNCTION z2f(xx1, xx2, zz0, zz1) RESULT(fn_val)

REAL(dp), INTENT(IN)  :: xx1, xx2, zz0, zz1
REAL(dp)              :: fn_val

fn_val = (zz1-zz0)*xx2/xx1 + zz0
RETURN
END FUNCTION z2f



FUNCTION z3f(xx1, xx2, xx3, zz0, zz1, zz2) RESULT(fn_val)

REAL(dp), INTENT(IN)  :: xx1, xx2, xx3, zz0, zz1, zz2
REAL(dp)              :: fn_val

fn_val = ((zz2-zz0)*(xx3-xx1)/xx2 - (zz1-zz0)*(xx3-xx2)/xx1) *  &
          (xx3/(xx2-xx1)) + zz0
RETURN
END FUNCTION z3f



SUBROUTINE rgpd3p(nxd, nyd, xd, yd, zd, pdd)

! Partial derivatives of a bivariate function on a rectangular grid
! (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)

! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 1995/08

! This subroutine estimates three partial derivatives, zx, zy, and
! zxy, of a bivariate function, z(x,y), on a rectangular grid in
! the x-y plane.  It is based on the revised Akima method that has
! the accuracy of a bicubic polynomial.

! The input arguments are
!   NXD = number of the input-grid data points in the x
!         coordinate (must be 2 or greater),
!   NYD = number of the input-grid data points in the y
!         coordinate (must be 2 or greater),
!   XD  = array of dimension NXD containing the x coordinates of the
!         input-grid data points (must be in a monotonic increasing order),
!   YD  = array of dimension NYD containing the y coordinates of the
!         input-grid data points (must be in a monotonic increasing order),
!   ZD  = two-dimensional array of dimension NXD*NYD
!         containing the z(x,y) values at the input-grid data points.

! The output argument is
!   PDD = three-dimensional array of dimension 3*NXD*NYD,
!         where the estimated zx, zy, and zxy values at the
!         input-grid data points are to be stored.


! Specification statements
!     .. Scalar Arguments ..

INTEGER, INTENT(IN)  :: nxd
INTEGER, INTENT(IN)  :: nyd
REAL(dp), INTENT(IN)     :: xd(nxd)
REAL(dp), INTENT(IN)     :: yd(nyd)
REAL(dp), INTENT(IN)     :: zd(nxd,nyd)
REAL(dp), INTENT(OUT)    :: pdd(3,nxd,nyd)

!     ..
!     .. Local Scalars ..
REAL(dp) :: b00, b00x, b00y, b01, b10, b11, cx1, cx2, cx3, cy1, cy2,  &
        cy3, disf, dnm, dz00, dz01, dz02, dz03, dz10, dz11, dz12,  &
        dz13, dz20, dz21, dz22, dz23, dz30, dz31, dz32, dz33,  &
        dzx10, dzx20, dzx30, dzxy11, dzxy12, dzxy13, dzxy21,  &
        dzxy22, dzxy23, dzxy31, dzxy32, dzxy33, dzy01, dzy02,  &
        dzy03, epsln, pezx, pezxy, pezy, smpef, smpei, smwtf,  &
        smwti, sx, sxx, sxxy, sxxyy, sxy, sxyy, sxyz, sxz, sy, syy,  &
        syz, sz, volf, wt, x0, x1, x2, x3, y0, y1, y2,  &
        y3, z00, z01, z02, z03, z10, z11, z12, z13, z20, z21, z22,  &
        z23, z30, z31, z32, z33, zxdi, zxydi, zydi
INTEGER :: ipex, ipey, ix0, ix1, ix2, ix3, iy0, iy1, iy2, iy3, nx0, ny0
!     ..
!     .. Local Arrays ..
REAL(dp)    :: b00xa(4), b00ya(4), b01a(4), b10a(4), cxa(3,4), cya(3,4),   &
           sxa(4), sxxa(4), sya(4), syya(4), xa(3,4), ya(3,4),   &
           z0ia(3,4), zi0a(3,4)
! INTEGER :: idlt(3,4)
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC        MAX
!     ..
!     .. Statement Functions ..
! REAL(dp) :: z2f,z3f
!     ..
! Data statements
! DATA ((idlt(jxy,jpexy),jpexy=1,4),jxy=1,3)/-3,-2,-1,1,-2,-1,1,2,-1,1,2,3/
INTEGER, SAVE  :: idlt(3,4) = RESHAPE(   &
                  (/ -3,-2,-1, 1,-2,-1, 1,2,-1, 1,2,3 /), (/ 3, 4 /) )
!     ..

! Calculation
! Initial setting of some local variables
nx0 = MAX(4,nxd)
ny0 = MAX(4,nyd)

! Double DO-loop with respect to the input grid points
DO  iy0 = 1,nyd
  DO  ix0 = 1,nxd
    x0 = xd(ix0)
    y0 = yd(iy0)
    z00 = zd(ix0,iy0)

! Part 1.  Estimation of ZXDI
! Initial setting
    smpef = 0.0
    smwtf = 0.0
    smpei = 0.0
    smwti = 0.0
! DO-loop with respect to the primary estimate
    DO  ipex = 1,4
! Selects necessary grid points in the x direction.
      ix1 = ix0 + idlt(1,ipex)
      ix2 = ix0 + idlt(2,ipex)
      ix3 = ix0 + idlt(3,ipex)
      IF ((ix1 < 1) .OR. (ix2 < 1) .OR. (ix3 < 1) .OR.  &
          (ix1 > nx0) .OR. (ix2 > nx0) .OR. (ix3 > nx0)) CYCLE
! Selects and/or supplements the x and z values.
      x1 = xd(ix1) - x0
      z10 = zd(ix1,iy0)
      IF (nxd >= 4) THEN
        x2 = xd(ix2) - x0
        x3 = xd(ix3) - x0
        z20 = zd(ix2,iy0)
        z30 = zd(ix3,iy0)
      ELSE IF (nxd == 3) THEN
        x2 = xd(ix2) - x0
        z20 = zd(ix2,iy0)
        x3 = 2*xd(3) - xd(2) - x0
        z30 = z3f(x1,x2,x3,z00,z10,z20)
      ELSE IF (nxd == 2) THEN
        x2 = 2*xd(2) - xd(1) - x0
        z20 = z2f(x1,x2,z00,z10)
        x3 = 2*xd(1) - xd(2) - x0
        z30 = z2f(x1,x3,z00,z10)
      END IF
      dzx10 = (z10-z00)/x1
      dzx20 = (z20-z00)/x2
      dzx30 = (z30-z00)/x3
! Calculates the primary estimate of partial derivative zx as
! the coefficient of the bicubic polynomial.
      cx1 = x2*x3/ ((x1-x2)* (x1-x3))
      cx2 = x3*x1/ ((x2-x3)* (x2-x1))
      cx3 = x1*x2/ ((x3-x1)* (x3-x2))
      pezx = cx1*dzx10 + cx2*dzx20 + cx3*dzx30
! Calculates the volatility factor and distance factor in the x
! direction for the primary estimate of zx.
      sx = x1 + x2 + x3
      sz = z00 + z10 + z20 + z30
      sxx = x1*x1 + x2*x2 + x3*x3
      sxz = x1*z10 + x2*z20 + x3*z30
      dnm = 4.0*sxx - sx*sx
      b00 = (sxx*sz-sx*sxz)/dnm
      b10 = (4.0*sxz-sx*sz)/dnm
      dz00 = z00 - b00
      dz10 = z10 - (b00+b10*x1)
      dz20 = z20 - (b00+b10*x2)
      dz30 = z30 - (b00+b10*x3)
      volf = dz00**2 + dz10**2 + dz20**2 + dz30**2
      disf = sxx

! Calculates the EPSLN value, which is used to decide whether or
! not the volatility factor is essentially zero.
      epsln = (z00**2+z10**2+z20**2+z30**2)*1.0E-12
! Accumulates the weighted primary estimates of zx and their weights.
      IF (volf > epsln) THEN
! - For a finite weight.
        wt = 1.0/ (volf*disf)
        smpef = smpef + wt*pezx
        smwtf = smwtf + wt
      ELSE
! - For an infinite weight.
        smpei = smpei + pezx
        smwti = smwti + 1.0
      END IF

! Saves the necessary values for estimating zxy
      xa(1,ipex) = x1
      xa(2,ipex) = x2
      xa(3,ipex) = x3
      zi0a(1,ipex) = z10
      zi0a(2,ipex) = z20
      zi0a(3,ipex) = z30
      cxa(1,ipex) = cx1
      cxa(2,ipex) = cx2
      cxa(3,ipex) = cx3
      sxa(ipex) = sx
      sxxa(ipex) = sxx
      b00xa(ipex) = b00
      b10a(ipex) = b10
    END DO

! Calculates the final estimate of zx.
    IF (smwti < 0.5) THEN
! - When no infinite weights exist.
      zxdi = smpef/smwtf
    ELSE
! - When infinite weights exist.
      zxdi = smpei/smwti
    END IF
! End of Part 1.

! Part 2.  Estimation of ZYDI
! Initial setting
    smpef = 0.0
    smwtf = 0.0
    smpei = 0.0
    smwti = 0.0
! DO-loop with respect to the primary estimate
    DO  ipey = 1,4
! Selects necessary grid points in the y direction.
      iy1 = iy0 + idlt(1,ipey)
      iy2 = iy0 + idlt(2,ipey)
      iy3 = iy0 + idlt(3,ipey)
      IF ((iy1 < 1) .OR. (iy2 < 1) .OR. (iy3 < 1) .OR.  &
          (iy1 > ny0) .OR. (iy2 > ny0) .OR. (iy3 > ny0)) CYCLE
! Selects and/or supplements the y and z values.
      y1 = yd(iy1) - y0
      z01 = zd(ix0,iy1)
      IF (nyd >= 4) THEN
        y2 = yd(iy2) - y0
        y3 = yd(iy3) - y0
        z02 = zd(ix0,iy2)
        z03 = zd(ix0,iy3)
      ELSE IF (nyd == 3) THEN
        y2 = yd(iy2) - y0
        z02 = zd(ix0,iy2)
        y3 = 2*yd(3) - yd(2) - y0
        z03 = z3f(y1,y2,y3,z00,z01,z02)
      ELSE IF (nyd == 2) THEN
        y2 = 2*yd(2) - yd(1) - y0
        z02 = z2f(y1,y2,z00,z01)
        y3 = 2*yd(1) - yd(2) - y0
        z03 = z2f(y1,y3,z00,z01)
      END IF
      dzy01 = (z01-z00)/y1
      dzy02 = (z02-z00)/y2
      dzy03 = (z03-z00)/y3
! Calculates the primary estimate of partial derivative zy as
! the coefficient of the bicubic polynomial.
      cy1 = y2*y3/ ((y1-y2)* (y1-y3))
      cy2 = y3*y1/ ((y2-y3)* (y2-y1))
      cy3 = y1*y2/ ((y3-y1)* (y3-y2))
      pezy = cy1*dzy01 + cy2*dzy02 + cy3*dzy03
! Calculates the volatility factor and distance factor in the y
! direction for the primary estimate of zy.
      sy = y1 + y2 + y3
      sz = z00 + z01 + z02 + z03
      syy = y1*y1 + y2*y2 + y3*y3
      syz = y1*z01 + y2*z02 + y3*z03
      dnm = 4.0*syy - sy*sy
      b00 = (syy*sz-sy*syz)/dnm
      b01 = (4.0*syz-sy*sz)/dnm
      dz00 = z00 - b00
      dz01 = z01 - (b00+b01*y1)
      dz02 = z02 - (b00+b01*y2)
      dz03 = z03 - (b00+b01*y3)
      volf = dz00**2 + dz01**2 + dz02**2 + dz03**2
      disf = syy

! Calculates the EPSLN value, which is used to decide whether or
! not the volatility factor is essentially zero.
      epsln = (z00**2+z01**2+z02**2+z03**2)*1.0E-12
! Accumulates the weighted primary estimates of zy and their weights.
      IF (volf > epsln) THEN
! - For a finite weight.
        wt = 1.0/ (volf*disf)
        smpef = smpef + wt*pezy
        smwtf = smwtf + wt
      ELSE
! - For an infinite weight.
        smpei = smpei + pezy
        smwti = smwti + 1.0
      END IF
! Saves the necessary values for estimating zxy
      ya(1,ipey) = y1
      ya(2,ipey) = y2
      ya(3,ipey) = y3
      z0ia(1,ipey) = z01
      z0ia(2,ipey) = z02
      z0ia(3,ipey) = z03
      cya(1,ipey) = cy1
      cya(2,ipey) = cy2
      cya(3,ipey) = cy3
      sya(ipey) = sy
      syya(ipey) = syy
      b00ya(ipey) = b00
      b01a(ipey) = b01
    END DO

! Calculates the final estimate of zy.
    IF (smwti < 0.5) THEN
! - When no infinite weights exist.
      zydi = smpef/smwtf
    ELSE
! - When infinite weights exist.
      zydi = smpei/smwti
    END IF
! End of Part 2.

! Part 3.  Estimation of ZXYDI
! Initial setting
    smpef = 0.0
    smwtf = 0.0
    smpei = 0.0
    smwti = 0.0
! Outer DO-loops with respect to the primary estimates in the x direction
    DO  ipex = 1,4
      ix1 = ix0 + idlt(1,ipex)
      ix2 = ix0 + idlt(2,ipex)
      ix3 = ix0 + idlt(3,ipex)
      IF ((ix1 < 1) .OR. (ix2 < 1) .OR. (ix3 < 1) .OR.  &
          (ix1 > nx0) .OR. (ix2 > nx0) .OR. (ix3 > nx0)) CYCLE
! Retrieves the necessary values for estimating zxy in the x direction.
      x1 = xa(1,ipex)
      x2 = xa(2,ipex)
      x3 = xa(3,ipex)
      z10 = zi0a(1,ipex)
      z20 = zi0a(2,ipex)
      z30 = zi0a(3,ipex)
      cx1 = cxa(1,ipex)
      cx2 = cxa(2,ipex)
      cx3 = cxa(3,ipex)
      sx = sxa(ipex)
      sxx = sxxa(ipex)
      b00x = b00xa(ipex)
      b10 = b10a(ipex)

! Inner DO-loops with respect to the primary estimates in the y direction
      DO  ipey = 1,4
        iy1 = iy0 + idlt(1,ipey)
        iy2 = iy0 + idlt(2,ipey)
        iy3 = iy0 + idlt(3,ipey)
        IF ((iy1 < 1) .OR. (iy2 < 1) .OR. (iy3 < 1) .OR. (iy1 > ny0) .OR.  &
            (iy2 > ny0) .OR. (iy3 > ny0)) CYCLE
! Retrieves the necessary values for estimating zxy in the y direction.
        y1 = ya(1,ipey)
        y2 = ya(2,ipey)
        y3 = ya(3,ipey)
        z01 = z0ia(1,ipey)
        z02 = z0ia(2,ipey)
        z03 = z0ia(3,ipey)
        cy1 = cya(1,ipey)
        cy2 = cya(2,ipey)
        cy3 = cya(3,ipey)
        sy = sya(ipey)
        syy = syya(ipey)
        b00y = b00ya(ipey)
        b01 = b01a(ipey)
! Selects and/or supplements the z values.
        IF (nyd >= 4) THEN
          z11 = zd(ix1,iy1)
          z12 = zd(ix1,iy2)
          z13 = zd(ix1,iy3)
          IF (nxd >= 4) THEN
            z21 = zd(ix2,iy1)
            z22 = zd(ix2,iy2)
            z23 = zd(ix2,iy3)
            z31 = zd(ix3,iy1)
            z32 = zd(ix3,iy2)
            z33 = zd(ix3,iy3)
          ELSE IF (nxd == 3) THEN
            z21 = zd(ix2,iy1)
            z22 = zd(ix2,iy2)
            z23 = zd(ix2,iy3)
            z31 = z3f(x1,x2,x3,z01,z11,z21)
            z32 = z3f(x1,x2,x3,z02,z12,z22)
            z33 = z3f(x1,x2,x3,z03,z13,z23)
          ELSE IF (nxd == 2) THEN
            z21 = z2f(x1,x2,z01,z11)
            z22 = z2f(x1,x2,z02,z12)
            z23 = z2f(x1,x2,z03,z13)
            z31 = z2f(x1,x3,z01,z11)
            z32 = z2f(x1,x3,z02,z12)
            z33 = z2f(x1,x3,z03,z13)
          END IF
        ELSE IF (nyd == 3) THEN
          z11 = zd(ix1,iy1)
          z12 = zd(ix1,iy2)
          z13 = z3f(y1,y2,y3,z10,z11,z12)
          IF (nxd >= 4) THEN
            z21 = zd(ix2,iy1)
            z22 = zd(ix2,iy2)
            z31 = zd(ix3,iy1)
            z32 = zd(ix3,iy2)
          ELSE IF (nxd == 3) THEN
            z21 = zd(ix2,iy1)
            z22 = zd(ix2,iy2)
            z31 = z3f(x1,x2,x3,z01,z11,z21)
            z32 = z3f(x1,x2,x3,z02,z12,z22)
          ELSE IF (nxd == 2) THEN
            z21 = z2f(x1,x2,z01,z11)
            z22 = z2f(x1,x2,z02,z12)
            z31 = z2f(x1,x3,z01,z11)
            z32 = z2f(x1,x3,z02,z12)
          END IF
          z23 = z3f(y1,y2,y3,z20,z21,z22)
          z33 = z3f(y1,y2,y3,z30,z31,z32)
        ELSE IF (nyd == 2) THEN
          z11 = zd(ix1,iy1)
          z12 = z2f(y1,y2,z10,z11)
          z13 = z2f(y1,y3,z10,z11)
          IF (nxd >= 4) THEN
            z21 = zd(ix2,iy1)
            z31 = zd(ix3,iy1)
          ELSE IF (nxd == 3) THEN
            z21 = zd(ix2,iy1)
            z31 = z3f(x1,x2,x3,z01,z11,z21)
          ELSE IF (nxd == 2) THEN
            z21 = z2f(x1,x2,z01,z11)
            z31 = z2f(x1,x3,z01,z11)
          END IF
          z22 = z2f(y1,y2,z20,z21)
          z23 = z2f(y1,y3,z20,z21)
          z32 = z2f(y1,y2,z30,z31)
          z33 = z2f(y1,y3,z30,z31)
        END IF
! Calculates the primary estimate of partial derivative zxy as
! the coefficient of the bicubic polynomial.
        dzxy11 = (z11-z10-z01+z00)/ (x1*y1)
        dzxy12 = (z12-z10-z02+z00)/ (x1*y2)
        dzxy13 = (z13-z10-z03+z00)/ (x1*y3)
        dzxy21 = (z21-z20-z01+z00)/ (x2*y1)
        dzxy22 = (z22-z20-z02+z00)/ (x2*y2)
        dzxy23 = (z23-z20-z03+z00)/ (x2*y3)
        dzxy31 = (z31-z30-z01+z00)/ (x3*y1)
        dzxy32 = (z32-z30-z02+z00)/ (x3*y2)
        dzxy33 = (z33-z30-z03+z00)/ (x3*y3)
        pezxy = cx1* (cy1*dzxy11+cy2*dzxy12+cy3*dzxy13) +  &
                cx2* (cy1*dzxy21+cy2*dzxy22+cy3*dzxy23) +  &
                cx3* (cy1*dzxy31+cy2*dzxy32+cy3*dzxy33)
! Calculates the volatility factor and distance factor in the x
! and y directions for the primary estimate of zxy.
        b00 = (b00x+b00y)/2.0
        sxy = sx*sy
        sxxy = sxx*sy
        sxyy = sx*syy
        sxxyy = sxx*syy
        sxyz = x1* (y1*z11+y2*z12+y3*z13) + x2* (y1*z21+y2*z22+y3*z23) +  &
               x3* (y1*z31+y2*z32+y3*z33)
        b11 = (sxyz-b00*sxy-b10*sxxy-b01*sxyy)/sxxyy
        dz00 = z00 - b00
        dz01 = z01 - (b00+b01*y1)
        dz02 = z02 - (b00+b01*y2)
        dz03 = z03 - (b00+b01*y3)
        dz10 = z10 - (b00+b10*x1)
        dz11 = z11 - (b00+b01*y1+x1* (b10+b11*y1))
        dz12 = z12 - (b00+b01*y2+x1* (b10+b11*y2))
        dz13 = z13 - (b00+b01*y3+x1* (b10+b11*y3))
        dz20 = z20 - (b00+b10*x2)
        dz21 = z21 - (b00+b01*y1+x2* (b10+b11*y1))
        dz22 = z22 - (b00+b01*y2+x2* (b10+b11*y2))
        dz23 = z23 - (b00+b01*y3+x2* (b10+b11*y3))
        dz30 = z30 - (b00+b10*x3)
        dz31 = z31 - (b00+b01*y1+x3* (b10+b11*y1))
        dz32 = z32 - (b00+b01*y2+x3* (b10+b11*y2))
        dz33 = z33 - (b00+b01*y3+x3* (b10+b11*y3))
        volf = dz00**2 + dz01**2 + dz02**2 + dz03**2 +  &
            dz10**2 + dz11**2 + dz12**2 + dz13**2 +  &
            dz20**2 + dz21**2 + dz22**2 + dz23**2 +  &
            dz30**2 + dz31**2 + dz32**2 + dz33**2
        disf = sxx*syy
! Calculates EPSLN.
        epsln = (z00**2 + z01**2 + z02**2 + z03**2 + z10**2 +   &
                z11**2 + z12**2 + z13**2 + z20**2 + z21**2 + z22**2 +   &
                z23**2 + z30**2 + z31**2 + z32**2 + z33**2)* 1.0E-12
! Accumulates the weighted primary estimates of zxy and their weights.
        IF (volf > epsln) THEN
! - For a finite weight.
          wt = 1.0/ (volf*disf)
          smpef = smpef + wt*pezxy
          smwtf = smwtf + wt
        ELSE
! - For an infinite weight.
          smpei = smpei + pezxy
          smwti = smwti + 1.0
        END IF
      END DO
    END DO

! Calculates the final estimate of zxy.
    IF (smwti < 0.5) THEN
! - When no infinite weights exist.
      zxydi = smpef/smwtf
    ELSE
! - When infinite weights exist.
      zxydi = smpei/smwti
    END IF
! End of Part 3

    pdd(1,ix0,iy0) = zxdi
    pdd(2,ix0,iy0) = zydi
    pdd(3,ix0,iy0) = zxydi
  END DO
END DO
RETURN
END SUBROUTINE rgpd3p



SUBROUTINE rglctn(nxd, nyd, xd, yd, nip, xi, yi, inxi, inyi)

! Location of the desired points in a rectangular grid
! (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)

! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 1995/08

! This subroutine locates the desired points in a rectangular grid
! in the x-y plane.

! The grid lines can be unevenly spaced.

! The input arguments are
!   NXD  = number of the input-grid data points in the x
!          coordinate (must be 2 or greater),
!   NYD  = number of the input-grid data points in the y
!          coordinate (must be 2 or greater),
!   XD   = array of dimension NXD containing the x coordinates of the
!          input-grid data points (must be in a monotonic increasing order),
!   YD   = array of dimension NYD containing the y coordinates of the
!          input-grid data points (must be in a monotonic increasing order),
!   NIP  = number of the output points to be located (must be 1 or greater),
!   XI   = array of dimension NIP containing the x coordinates
!          of the output points to be located,
!   YI   = array of dimension NIP containing the y coordinates
!          of the output points to be located.

! The output arguments are
!   INXI = integer array of dimension NIP where the interval
!          numbers of the XI array elements are to be stored,
!   INYI = integer array of dimension NIP where the interval
!          numbers of the YI array elements are to be stored.
! The interval numbers are between 0 and NXD and between 0 and NYD,
! respectively.


! Specification statements
!     .. Scalar Arguments ..

INTEGER, INTENT(IN)   :: nxd
INTEGER, INTENT(IN)   :: nyd
REAL(dp), INTENT(IN)      :: xd(nxd)
REAL(dp), INTENT(IN)      :: yd(nyd)
INTEGER, INTENT(IN)   :: nip
REAL(dp), INTENT(IN)      :: xi(nip)
REAL(dp), INTENT(IN)      :: yi(nip)
INTEGER, INTENT(OUT)  :: inxi(nip)
INTEGER, INTENT(OUT)  :: inyi(nip)

!     ..
!     .. Local Scalars ..
REAL(dp)     :: xii, yii
INTEGER  :: iip, imd, imn, imx, ixd, iyd, nintx, ninty

!     ..
! DO-loop with respect to IIP, which is the point number of the output point
DO  iip = 1,nip
  xii = xi(iip)
  yii = yi(iip)
! Checks if the x coordinate of the IIPth output point, XII, is
! in a new interval.  (NINTX is the new-interval flag.)
  IF (iip == 1) THEN
    nintx = 1
  ELSE
    nintx = 0
    IF (ixd == 0) THEN
      IF (xii > xd(1)) nintx = 1
    ELSE IF (ixd < nxd) THEN
      IF ((xii < xd(ixd)) .OR. (xii > xd(ixd+1))) nintx = 1
    ELSE
      IF (xii < xd(nxd)) nintx = 1
    END IF
  END IF

! Locates the output point by binary search if XII is in a new interval.
! Determines IXD for which XII lies between XD(IXD) and XD(IXD+1).
  IF (nintx == 1) THEN
    IF (xii <= xd(1)) THEN
      ixd = 0
    ELSE IF (xii < xd(nxd)) THEN
      imn = 1
      imx = nxd
      imd = (imn+imx)/2
      10 IF (xii >= xd(imd)) THEN
        imn = imd
      ELSE
        imx = imd
      END IF
      imd = (imn+imx)/2
      IF (imd > imn) GO TO 10
      ixd = imd
    ELSE
      ixd = nxd
    END IF
  END IF
  inxi(iip) = ixd

! Checks if the y coordinate of the IIPth output point, YII, is
! in a new interval.  (NINTY is the new-interval flag.)
  IF (iip == 1) THEN
    ninty = 1
  ELSE
    ninty = 0
    IF (iyd == 0) THEN
      IF (yii > yd(1)) ninty = 1
    ELSE IF (iyd < nyd) THEN
      IF ((yii < yd(iyd)) .OR. (yii > yd(iyd+1))) ninty = 1
    ELSE
      IF (yii < yd(nyd)) ninty = 1
    END IF
  END IF

! Locates the output point by binary search if YII is in a new interval.
! Determines IYD for which YII lies between YD(IYD) and YD(IYD+1).
  IF (ninty == 1) THEN
    IF (yii <= yd(1)) THEN
      iyd = 0
    ELSE IF (yii < yd(nyd)) THEN
      imn = 1
      imx = nyd
      imd = (imn+imx)/2
      20 IF (yii >= yd(imd)) THEN
        imn = imd
      ELSE
        imx = imd
      END IF
      imd = (imn+imx)/2
      IF (imd > imn) GO TO 20
      iyd = imd
    ELSE
      iyd = nyd
    END IF
  END IF
  inyi(iip) = iyd
END DO
RETURN
END SUBROUTINE rglctn



SUBROUTINE rgplnl(nxd, nyd, xd, yd, zd, pdd, nip, xi, yi, inxi, inyi, zi)

! Polynomials for rectangular-grid bivariate interpolation and surface fitting
! (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)

! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 1995/08

! This subroutine determines a polynomial in x and y for a rectangle of the
! input grid in the x-y plane and calculates the z value for the desired
! points by evaluating the polynomial for rectangular-grid bivariate
! interpolation and surface fitting.

! The input arguments are
!   NXD  = number of the input-grid data points in the x
!          coordinate (must be 2 or greater),
!   NYD  = number of the input-grid data points in the y
!          coordinate (must be 2 or greater),
!   XD   = array of dimension NXD containing the x coordinates of the
!          input-grid data points (must be in a monotonic increasing order),
!   YD   = array of dimension NYD containing the y coordinates of the
!          input-grid data points (must be in a monotonic increasing order),
!   ZD   = two-dimensional array of dimension NXD*NYD
!          containing the z(x,y) values at the input-grid data points,
!   PDD  = three-dimensional array of dimension 3*NXD*NYD
!          containing the estimated zx, zy, and zxy values
!          at the input-grid data points,
!   NIP  = number of the output points at which interpolation
!          is to be performed,
!   XI   = array of dimension NIP containing the x coordinates
!          of the output points,
!   YI   = array of dimension NIP containing the y coordinates
!          of the output points,
!   INXI = integer array of dimension NIP containing the
!          interval numbers of the input grid intervals in the
!          x direction where the x coordinates of the output points lie,
!   INYI = integer array of dimension NIP containing the
!          interval numbers of the input grid intervals in the
!          y direction where the y coordinates of the output points lie.

! The output argument is
!   ZI   = array of dimension NIP, where the interpolated z
!          values at the output points are to be stored.


! Specification statements
!     .. Scalar Arguments ..

INTEGER, INTENT(IN)  :: nxd
INTEGER, INTENT(IN)  :: nyd
REAL(dp), INTENT(IN)     :: xd(nxd)
REAL(dp), INTENT(IN)     :: yd(nyd)
REAL(dp), INTENT(IN)     :: zd(nxd,nyd)
REAL(dp), INTENT(IN)     :: pdd(3,nxd,nyd)
INTEGER, INTENT(IN)  :: nip
REAL(dp), INTENT(IN)     :: xi(nip)
REAL(dp), INTENT(IN)     :: yi(nip)
INTEGER, INTENT(IN)  :: inxi(nip)
INTEGER, INTENT(IN)  :: inyi(nip)
REAL(dp), INTENT(OUT)    :: zi(nip)

!     ..
!     .. Local Scalars ..
REAL(dp) :: a, b, c, d, dx, dxsq, dy, dysq, p00, p01, p02, p03, p10, p11,  &
        p12, p13, p20, p21, p22, p23, p30, p31, p32, p33, q0, q1, q2,  &
        q3, u, v, x0, xii, y0, yii, z00, z01, z0dx, z0dy, z10, z11,  &
        z1dx, z1dy, zdxdy, zii, zx00, zx01, zx0dy, zx10, zx11,  &
        zx1dy, zxy00, zxy01, zxy10, zxy11, zy00, zy01, zy0dx, zy10, zy11, zy1dx
INTEGER :: iip, ixd0, ixd1, ixdi, ixdipv, iyd0, iyd1, iydi, iydipv
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC        MAX
!     ..

! Calculation
! Outermost DO-loop with respect to the output point
DO  iip = 1,nip
  xii = xi(iip)
  yii = yi(iip)
  IF (iip == 1) THEN
    ixdipv = -1
    iydipv = -1
  ELSE
    ixdipv = ixdi
    iydipv = iydi
  END IF
  ixdi = inxi(iip)
  iydi = inyi(iip)

! Retrieves the z and partial derivative values at the origin of
! the coordinate for the rectangle.
  IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
    ixd0 = MAX(1,ixdi)
    iyd0 = MAX(1,iydi)
    x0 = xd(ixd0)
    y0 = yd(iyd0)
    z00 = zd(ixd0,iyd0)
    zx00 = pdd(1,ixd0,iyd0)
    zy00 = pdd(2,ixd0,iyd0)
    zxy00 = pdd(3,ixd0,iyd0)
  END IF

! Case 1.  When the rectangle is inside the data area in both the
! x and y directions.
  IF ((ixdi > 0 .AND. ixdi < nxd) .AND. (iydi > 0 .AND. iydi < nyd)) THEN
! Retrieves the z and partial derivative values at the other three
! vertices of the rectangle.
    IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
      ixd1 = ixd0 + 1
      dx = xd(ixd1) - x0
      dxsq = dx*dx
      iyd1 = iyd0 + 1
      dy = yd(iyd1) - y0
      dysq = dy*dy
      z10 = zd(ixd1,iyd0)
      z01 = zd(ixd0,iyd1)
      z11 = zd(ixd1,iyd1)
      zx10 = pdd(1,ixd1,iyd0)
      zx01 = pdd(1,ixd0,iyd1)
      zx11 = pdd(1,ixd1,iyd1)
      zy10 = pdd(2,ixd1,iyd0)
      zy01 = pdd(2,ixd0,iyd1)
      zy11 = pdd(2,ixd1,iyd1)
      zxy10 = pdd(3,ixd1,iyd0)
      zxy01 = pdd(3,ixd0,iyd1)
      zxy11 = pdd(3,ixd1,iyd1)
! Calculates the polynomial coefficients.
      z0dx = (z10-z00)/dx
      z1dx = (z11-z01)/dx
      z0dy = (z01-z00)/dy
      z1dy = (z11-z10)/dy
      zx0dy = (zx01-zx00)/dy
      zx1dy = (zx11-zx10)/dy
      zy0dx = (zy10-zy00)/dx
      zy1dx = (zy11-zy01)/dx
      zdxdy = (z1dy-z0dy)/dx
      a = zdxdy - zx0dy - zy0dx + zxy00
      b = zx1dy - zx0dy - zxy10 + zxy00
      c = zy1dx - zy0dx - zxy01 + zxy00
      d = zxy11 - zxy10 - zxy01 + zxy00
      p00 = z00
      p01 = zy00
      p02 = (2.0* (z0dy-zy00)+z0dy-zy01)/dy
      p03 = (-2.0*z0dy+zy01+zy00)/dysq
      p10 = zx00
      p11 = zxy00
      p12 = (2.0* (zx0dy-zxy00)+zx0dy-zxy01)/dy
      p13 = (-2.0*zx0dy+zxy01+zxy00)/dysq
      p20 = (2.0* (z0dx-zx00)+z0dx-zx10)/dx
      p21 = (2.0* (zy0dx-zxy00)+zy0dx-zxy10)/dx
      p22 = (3.0* (3.0*a-b-c)+d)/ (dx*dy)
      p23 = (-6.0*a+2.0*b+3.0*c-d)/ (dx*dysq)
      p30 = (-2.0*z0dx+zx10+zx00)/dxsq
      p31 = (-2.0*zy0dx+zxy10+zxy00)/dxsq
      p32 = (-6.0*a+3.0*b+2.0*c-d)/ (dxsq*dy)
      p33 = (2.0* (2.0*a-b-c)+d)/ (dxsq*dysq)
    END IF

! Evaluates the polynomial.
    u = xii - x0
    v = yii - y0
    q0 = p00 + v* (p01+v* (p02+v*p03))
    q1 = p10 + v* (p11+v* (p12+v*p13))
    q2 = p20 + v* (p21+v* (p22+v*p23))
    q3 = p30 + v* (p31+v* (p32+v*p33))
    zii = q0 + u* (q1+u* (q2+u*q3))
! End of Case 1

! Case 2.  When the rectangle is inside the data area in the x
! direction but outside in the y direction.
  ELSE IF ((ixdi > 0.AND.ixdi < nxd) .AND.  &
        (iydi <= 0.OR.iydi >= nyd)) THEN
! Retrieves the z and partial derivative values at the other
! vertex of the semi-infinite rectangle.
    IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
      ixd1 = ixd0 + 1
      dx = xd(ixd1) - x0
      dxsq = dx*dx
      z10 = zd(ixd1,iyd0)
      zx10 = pdd(1,ixd1,iyd0)
      zy10 = pdd(2,ixd1,iyd0)
      zxy10 = pdd(3,ixd1,iyd0)
! Calculates the polynomial coefficients.
      z0dx = (z10-z00)/dx
      zy0dx = (zy10-zy00)/dx
      p00 = z00
      p01 = zy00
      p10 = zx00
      p11 = zxy00
      p20 = (2.0* (z0dx-zx00)+z0dx-zx10)/dx
      p21 = (2.0* (zy0dx-zxy00)+zy0dx-zxy10)/dx
      p30 = (-2.0*z0dx+zx10+zx00)/dxsq
      p31 = (-2.0*zy0dx+zxy10+zxy00)/dxsq
    END IF
! Evaluates the polynomial.
    u = xii - x0
    v = yii - y0
    q0 = p00 + v*p01
    q1 = p10 + v*p11
    q2 = p20 + v*p21
    q3 = p30 + v*p31
    zii = q0 + u* (q1+u* (q2+u*q3))
! End of Case 2

! Case 3.  When the rectangle is outside the data area in the x
! direction but inside in the y direction.
  ELSE IF ((ixdi <= 0.OR.ixdi >= nxd) .AND.  &
        (iydi > 0 .AND. iydi < nyd)) THEN
! Retrieves the z and partial derivative values at the other
! vertex of the semi-infinite rectangle.
    IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
      iyd1 = iyd0 + 1
      dy = yd(iyd1) - y0
      dysq = dy*dy
      z01 = zd(ixd0,iyd1)
      zx01 = pdd(1,ixd0,iyd1)
      zy01 = pdd(2,ixd0,iyd1)
      zxy01 = pdd(3,ixd0,iyd1)
! Calculates the polynomial coefficients.
      z0dy = (z01-z00)/dy
      zx0dy = (zx01-zx00)/dy
      p00 = z00
      p01 = zy00
      p02 = (2.0*(z0dy-zy00)+z0dy-zy01)/dy
      p03 = (-2.0*z0dy+zy01+zy00)/dysq
      p10 = zx00
      p11 = zxy00
      p12 = (2.0*(zx0dy-zxy00) + zx0dy - zxy01)/dy
      p13 = (-2.0*zx0dy + zxy01 + zxy00)/dysq
    END IF

! Evaluates the polynomial.
    u = xii - x0
    v = yii - y0
    q0 = p00 + v* (p01 + v*(p02+v*p03))
    q1 = p10 + v* (p11 + v*(p12+v*p13))
    zii = q0 + u*q1
! End of Case 3

! Case 4.  When the rectangle is outside the data area in both the
! x and y direction.
  ELSE IF ((ixdi <= 0 .OR. ixdi >= nxd) .AND.  &
           (iydi <= 0 .OR. iydi >= nyd)) THEN
! Calculates the polynomial coefficients.
    IF (ixdi /= ixdipv .OR. iydi /= iydipv) THEN
      p00 = z00
      p01 = zy00
      p10 = zx00
      p11 = zxy00
    END IF
! Evaluates the polynomial.
    u = xii - x0
    v = yii - y0
    q0 = p00 + v*p01
    q1 = p10 + v*p11
    zii = q0 + u*q1
  END IF
! End of Case 4
  zi(iip) = zii
END DO

RETURN
END SUBROUTINE rgplnl

!!!!!!!!!!!!!!!!
! Routine for integration of a function by cubic interpolation
!!!!!!!!!!!!!!!!
subroutine cubint ( ntab, xtab, ftab, ia, ib, result, error )

!*****************************************************************************80
!
!! CUBINT approximates an integral using cubic interpolation of data.
!
!  Discussion:
!
!    The integral to be approximated is
! 
!      Integral ( XTAB(IB) <= X <= XTAB(IA) ) F(X) DX
!
!    The routine estimates the error in integration.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Philip Gill, GF Miller,
!    An algorithm for the integration of unequally spaced data,
!    The Computer Journal, 
!    Number 15, Number 1, 1972, pages 80-83.
!
!  Parameters:
!
!    Input, integer NTAB, the number of tabulated points.
!    NTAB must be at least 4.
!
!    Input, real ( dp ) XTAB(NTAB), contains the points at which the
!    function was tabulated.  XTAB should contain distinct
!    values, given in ascending order.
!
!    Input, real ( dp ) FTAB(NTAB), contains the tabulated function
!    values, FTAB(I) = F(XTAB(I)).
!
!    Input, integer IA, the entry of XTAB at which integration
!    is to begin.  IA must be no less than 1 and no greater
!    than NTAB.
!
!    Input, integer IB, the entry of XTAB at which integration
!    is to end.  IB must be no less than 1 and no greater than
!    NTAB.
!
!    Output, real ( dp ) RESULT, the approximate value of the
!    integral from XTAB(IA) to XTAB(IB) of the function.
!
!    Output, real ( dp ) ERROR, an estimate of the error in
!    integration.
!
  implicit none

  integer ntab

  real ( dp ) c
  real ( dp ) d1
  real ( dp ) d2
  real ( dp ) d3
  real ( dp ) error
  real ( dp ) ftab(ntab)
  real ( dp ) h1
  real ( dp ) h2
  real ( dp ) h3
  real ( dp ) h4
  integer i
  integer ia
  integer ib
  integer ind
  integer it
  integer j
  integer k
  real ( dp ) r1
  real ( dp ) r2
  real ( dp ) r3
  real ( dp ) r4
  real ( dp ) result
  real ( dp ) s
  real ( dp ) term
  real ( dp ) xtab(ntab)

  result = 0.0_dp
  error = 0.0_dp
 
  if ( ia == ib ) then
    return
  end if
 
  if ( ntab < 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  NTAB must be at least 4, but input NTAB = ', ntab
    stop
  end if
 
  if ( ia < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IA must be at least 1, but input IA = ', ia
    stop
  end if
 
  if ( ntab < ia ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IA must be <= NTAB, but input IA = ', ia
    stop
  end if
 
  if ( ib < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IB must be at least 1, but input IB = ', ib
    stop
  end if
 
  if ( ntab < ib ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IB must be <= NTAB, but input IB = ', ib
    stop
  end if
!
!  Temporarily switch IA and IB, and store minus sign in IND
!  so that, while integration is carried out from low X's
!  to high ones, the sense of the integral is preserved.
!
  if ( ib < ia ) then
    ind = -1
    it = ib
    ib = ia
    ia = it
  else
    ind = 1
  end if
 
  s = 0.0_dp
  c = 0.0_dp
  r4 = 0.0_dp
  j = ntab-2
  if ( ia < ntab-1 .or. ntab == 4 ) then
    j = max ( 3, ia )
  end if

  k = 4
  if ( 2 < ib .or. ntab == 4 ) then
    k = min ( ntab, ib + 2 ) - 1
  end if
 
  do i = j, k
 
    if ( i <= j ) then
 
      h2 = xtab(j-1) - xtab(j-2)
      d3 = ( ftab(j-1) - ftab(j-2) ) / h2
      h3 = xtab(j) - xtab(j-1)
      d1 = ( ftab(j) - ftab(j-1) ) / h3
      h1 = h2 + h3
      d2 = ( d1 - d3 ) / h1
      h4 = xtab(j+1) - xtab(j)
      r1 = ( ftab(j+1) - ftab(j) ) / h4
      r2 = ( r1 - d1 ) / ( h4 + h3 )
      h1 = h1 + h4
      r3 = (r2-d2) / h1
 
      if ( ia <= 1 ) then
        result = h2 * ( ftab(1) + h2 * ( 0.5_dp * d3 - h2 &
          * ( d2 / 6.0_dp -(h2+h3+h3)*r3/12.0_dp)))
        s = -h2**3 * (h2*(3.0_dp*h2+5.0_dp*h4)+10.0_dp*h3*h1) / 60.0_dp
      end if
 
    else
 
      h4 = xtab(i+1) - xtab(i)
      r1 = ( ftab(i+1) - ftab(i) ) / h4
      r4 = h4 + h3
      r2 = ( r1 - d1 ) / r4
      r4 = r4 + h2
      r3 = ( r2 - d2 ) / r4
      r4 = ( r3 - d3 ) / ( r4 + h1 )
 
    end if
 
    if ( ia < i .and. i <= ib ) then
 
      term = h3 * ( ( ftab(i) + ftab(i-1) ) * 0.5_dp &
        -h3 * h3 * ( d2 + r2 + ( h2 - h4 ) * r3 ) / 12.0_dp )
      result = result + term
      c = h3**3 * ( 2.0_dp * h3 * h3 &
        + 5.0_dp * ( h3 * ( h4 + h2 ) + 2.0_dp * h2 * h4 ) ) / 120.0_dp
      error = error + (c+s)*r4
 
      if ( i /= j ) then
        s = c
      else
        s = s + c + c
      end if
 
    else
 
      error = error + r4 * s
 
    end if
 
    if ( k <= i ) then
 
      if ( ntab <= ib ) then
        term = h4 * ( ftab(ntab) - h4 * ( 0.5 * r1 &
          + h4 * ( r2 / 6.0_dp + ( h3 + h3 + h4 ) * r3 / 12.0_dp )))
        result = result + term
        error = error - h4**3 * r4 * &
          ( h4 * ( 3.0_dp * h4 + 5.0_dp * h2 ) &
          + 10.0_dp * h3 * ( h2 + h3 + h4 ) ) / 60.0_dp
      end if
 
      if ( ntab-1 <= ib ) then
        error = error + s * r4
      end if

    else

      h1 = h2
      h2 = h3
      h3 = h4
      d1 = r1
      d2 = r2
      d3 = r3
    end if
 
  end do
!
!  Restore original values of IA and IB, reverse signs
!  of RESULT and ERROR, to account for integration
!  that proceeded from high X to low X.
!
  if ( ind /= 1 ) then
    it = ib
    ib = ia
    ia = it
    result = -result
    error = -error
  end if
 
  return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!
! Weddle's rule for numerical integration
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wedint ( ntab, h, ftab, result )

!*****************************************************************************80
!
!! WEDINT uses Weddle's rule to integrate data at equally spaced points.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NTAB, is the number of data points.  (NTAB-1) must be
!    divisible by 6.
!
!    Input, real ( dp ) H, is the spacing between the points at which
!    the data was evaluated.
!
!    Input, real ( dp ) FTAB(NTAB), contains the tabulated data values.
!
!    Output, real ( dp ) RESULT, is the approximation to the integral.
!
  implicit none

  integer ntab

  real ( dp ) ftab(ntab)
  real ( dp ) h
  integer i
  real ( dp ) result

  result = 0.0_dp
 
  if ( ntab <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDINT - Fatal error!'
    write ( *, '(a)' ) '  NTAB < 2'
    write ( *, '(a,i8)' ) '  NTAB = ', ntab
    stop
  end if
 
  if ( mod ( ntab, 6 ) /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDINT - Fatal error!'
    write ( *, '(a)' ) '  NTAB must equal 6*N+1 for some N!'
    stop
  end if
 
  do i = 1, ntab-6, 6
    result = result & 
      +           ftab(i)   &
      + 5.0_dp * ftab(i+1) &
      +           ftab(i+2) &
      + 6.0_dp * ftab(i+3) &
      +           ftab(i+4) &
      + 5.0_dp * ftab(i+5) &
      +           ftab(i+6)
  end do
 
  result = 3.0_dp * h * result / 10.0_dp
 
  return
end subroutine

!!!!!!!!!!!!!!!
! Adaptive simpsons rule
!!!!!!!!!!!!!!!

subroutine simp ( func, a, b, eps, result )

!*****************************************************************************80
!
!! SIMP approximates the integral of a function by an adaptive Simpson's rule.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    James Lyness,
!    Algorithm 379:
!    SQUANK - Simpson quadrature used adaptively, noise killed,
!    Communications of the ACM,
!    Volume 13, Number 4, April 1970, pages 260-263.
!
!    William McKeeman, Lawrence Tesler,
!    Algorithm 182:
!    Nonrecursive adaptive integration,
!    Communications of the ACM,
!    Volume 6, 1963, page 315.
!
!  Parameters:
!
!    Input, real ( dp ), external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      FUNCTION FUNC ( X )
!    which evaluates the function at the point X.
!
!    Input, real ( dp ) A, the lower limit of integration.
!
!    Input, real ( dp ) B, the upper limit of integration.
!
!    Input, real ( dp ) EPS, the requested error tolerance.
!
!    Output, real ( dp ) RESULT, the approximation to the integral.
!
  implicit none

  integer, parameter :: maxlev = 30

  real ( dp ) a
  real ( dp ) a1
  real ( dp ) absar
  real ( dp ) b
  real ( dp ) da
  real ( dp ) dx(maxlev)
  real ( dp ) ep
  real ( dp ) eps
  real ( dp ) epsp(maxlev)
  real ( dp ) est
  real ( dp ) est1
  real ( dp ) est2(maxlev)
  real ( dp ) est3(maxlev)
  real ( dp ) f1
  real ( dp ) f2(maxlev)
  real ( dp ) f3(maxlev)
  real ( dp ) f4(maxlev)
  real ( dp ) fa
  real ( dp ) fb
  real ( dp ) fbp(maxlev)
  real ( dp ) fm
  real ( dp ) fmp(maxlev)
  real ( dp ), external :: func
  integer l
  integer lvl
  integer nrtr(maxlev)
  real ( dp ) pval(maxlev,3)
  real ( dp ) result
  real ( dp ) sum1
  real ( dp ) sx
  real ( dp ) x2(maxlev)
  real ( dp ) x3(maxlev)

  result = 0.0_dp
 
  if ( a == b ) then
    return
  end if
 
  ep = eps
  a1 = a
  nrtr(1:maxlev) = 0
  pval(1:maxlev,1:3) = 0.0_dp
 
  lvl = 0
  absar = 0.0_dp
  est = 0.0_dp
  da = b - a1

  fa = func ( a1 )
  fm = 4.0_dp * func ( ( a1 + b ) * 0.5_dp )
  fb = func ( b )
!
!  1 = RECUR
!
   30 continue
 
  lvl = lvl + 1
  dx(lvl) = da / 3.0_dp
  sx = dx(lvl) / 6.0_dp
  f1 = 4.0_dp * func ( 0.5_dp * dx(lvl) + a1 )
  x2(lvl) = a1 + dx(lvl)
  f2(lvl) = func ( x2(lvl) )
  x3(lvl) = x2(lvl) + dx(lvl)
  f3(lvl) = func ( x3(lvl) )
  epsp(lvl) = ep
  f4(lvl) = 4.0_dp * func ( dx(lvl) * 0.5_dp + x3(lvl) )
  fmp(lvl) = fm
  est1 = sx * (fa+f1+f2(lvl))
  fbp(lvl) = fb
  est2(lvl) = sx * ( f2(lvl) + f3(lvl) + fm )
  est3(lvl) = sx * ( f3(lvl) + f4(lvl) + fb )
  sum1 = est1 + est2(lvl) + est3(lvl)
  absar = absar - abs ( est ) + abs ( est1 ) + abs ( est2(lvl) ) &
    + abs ( est3(lvl) )

  if ( abs ( est - sum1 ) <= epsp(lvl) * absar ) then
    go to 40
  end if

  if ( maxlev <= lvl ) then
    go to 50
  end if
!
!  2 = UP
!
40 continue
 
  if ( 1 < lvl ) then
    lvl = lvl-1
  end if

  l = nrtr(lvl)

  if ( l == 0 ) then
    go to 50
  end if

  pval(lvl,l) = sum1

  if ( l == 1 ) then
    go to 60
  end if

  if ( l == 2 ) then
    go to 70
  end if

  if ( l == 3 ) then
    go to 80
  end if
 
50 continue
 
  nrtr(lvl) = 1
  est = est1
  fm = f1
  fb = f2(lvl)
  ep = epsp(lvl) / 1.7_dp
  da = dx(lvl)
  go to 30
 
60 continue
 
  nrtr(lvl) = 2
  fa = f2(lvl)
  fm = fmp(lvl)
  fb = f3(lvl)
  est = est2(lvl)
  a1 = x2(lvl)
  ep = epsp(lvl) / 1.7_dp
  da = dx(lvl)
  go to 30
 
70 continue
 
  nrtr(lvl) = 3
  fa = f3(lvl)
  fm = f4(lvl)
  fb = fbp(lvl)
  est = est3(lvl)
  a1 = x3(lvl)
  ep = epsp(lvl) / 1.7_dp
  da = dx(lvl)
  go to 30
 
80 continue
 
  sum1 = pval(lvl,1) + pval(lvl,2) + pval(lvl,3)

  if ( 1 < lvl ) then
    go to 40
  end if
 
90 continue
 
  result = sum1
 
  return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INVERSE NORMAL DISTRIBUTION FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(dp) function r8_normal_01_cdf_inverse ( p )

!*****************************************************************************80
!
!! R8_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10**16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    The Percentage Points of the Normal Distribution,
!    Algorithm AS 241,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( dp ) P, the value of the cumulative probability 
!    densitity function.  0 < P < 1.  If P is outside this range,
!    an "infinite" value will be returned.
!
!    Output, real ( dp ) D_NORMAL_01_CDF_INVERSE, the normal deviate 
!    value with the property that the probability of a standard normal 
!    deviate being less than or equal to the value is P.
!
  implicit none
!  real( dp ) :: r8_normal_01_cdf_inverse 
  real( dp ), intent(in):: p

  real    ( dp ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080_dp, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real    ( dp ), parameter, dimension ( 8 ) :: b = (/ &
    1.0_dp, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real   ( dp ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734_dp, &
    4.63033784615654529590_dp, &
    5.76949722146069140550_dp, &
    3.64784832476320460504_dp, &
    1.27045825245236838258_dp, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real    ( dp ), parameter :: const1 = 0.180625_dp
  real    ( dp ), parameter :: const2 = 1.6_dp
  real    ( dp ), parameter, dimension ( 8 ) :: d = (/ &
    1.0_dp, &
    2.05319162663775882187_dp, &
    1.67638483018380384940_dp, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real    ( dp ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720_dp, &
    5.46378491116411436990_dp, &
    1.78482653991729133580_dp, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real    ( dp ), parameter, dimension ( 8 ) :: f = (/ &
    1.0_dp, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real    ( dp ) q
  real    ( dp ) r
!  real    ( dp ) r8poly_value
  real    ( dp ), parameter :: split1 = 0.425_dp
  real    ( dp ), parameter :: split2 = 5.0_dp

  if ( p <= 0.0_dp ) then
    r8_normal_01_cdf_inverse = - huge ( p )
    return
  end if

  if ( 1.0_dp <= p ) then
    r8_normal_01_cdf_inverse = huge ( p )
    return
  end if

  q = p - 0.5_dp

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    r8_normal_01_cdf_inverse = q * r8poly_value ( 8, a, r ) &
                                 / r8poly_value ( 8, b, r )

  else

    if ( q < 0.0_dp ) then
      r = p
    else
      r = 1.0_dp - p
    end if

    if ( r <= 0.0_dp ) then
      r8_normal_01_cdf_inverse = - 1.0_dp
      stop
    end if

    r = sqrt ( -log ( r ) )

    if ( r <= split2 ) then

      r = r - const2
      r8_normal_01_cdf_inverse = r8poly_value ( 8, c, r ) &
                               / r8poly_value ( 8, d, r )

    else

      r = r - split2
      r8_normal_01_cdf_inverse = r8poly_value ( 8, e, r ) &
                               / r8poly_value ( 8, f, r )
   
    end if

    if ( q < 0.0_dp ) then
      r8_normal_01_cdf_inverse = - r8_normal_01_cdf_inverse
    end if

  end if

  return
end function r8_normal_01_cdf_inverse

REAL(dp) function r8poly_value ( n, a, x )

!*****************************************************************************80
!
!! R8POLY_VALUE evaluates an R8POLY
!
!  Discussion:
!
!    For sanity's sake, the value of N indicates the NUMBER of 
!    coefficients, or more precisely, the ORDER of the polynomial,
!    rather than the DEGREE of the polynomial.  The two quantities
!    differ by 1, but cause a great deal of confusion.
!
!    Given N and A, the form of the polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( dp ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( dp ) X, the point at which the polynomial is 
!    to be evaluated.
!
!    Output, real ( dp ) R8POLY_VALUE, the value of the polynomial at X.
!
  implicit none

  integer( kind = 4 ),intent(in):: n
  real( dp ),intent(in):: a(n)
  real( dp ), intent(in):: x
!  real( dp ):: r8poly_value
  
  integer( kind = 4 ) i

  r8poly_value = 0.0_dp
  do i = n, 1, -1
    r8poly_value = r8poly_value * x + a(i)
  end do

  return
end function r8poly_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Normal cumulative density function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp) function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    David Hill
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( dp ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( dp ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none
  logical, intent(in):: upper
  real ( dp ),intent(in):: x

  real ( dp ), parameter :: a1 = 5.75885480458_dp
  real ( dp ), parameter :: a2 = 2.62433121679_dp
  real ( dp ), parameter :: a3 = 5.92885724438_dp
  real ( dp ), parameter :: b1 = -29.8213557807_dp
  real ( dp ), parameter :: b2 = 48.6959930692_dp
  real ( dp ), parameter :: c1 = -0.000000038052_dp
  real ( dp ), parameter :: c2 = 0.000398064794_dp
  real ( dp ), parameter :: c3 = -0.151679116635_dp
  real ( dp ), parameter :: c4 = 4.8385912808_dp
  real ( dp ), parameter :: c5 = 0.742380924027_dp
  real ( dp ), parameter :: c6 = 3.99019417011_dp
  real ( dp ), parameter :: con = 1.28_dp
  real ( dp ), parameter :: d1 = 1.00000615302_dp
  real ( dp ), parameter :: d2 = 1.98615381364_dp
  real ( dp ), parameter :: d3 = 5.29330324926_dp
  real ( dp ), parameter :: d4 = -15.1508972451_dp
  real ( dp ), parameter :: d5 = 30.789933034_dp
  real ( dp ), parameter :: ltone = 7.0_dp
  real ( dp ), parameter :: p = 0.398942280444_dp
  real ( dp ), parameter :: q = 0.39990348504_dp
  real ( dp ), parameter :: r = 0.398942280385_dp
  logical up
  real ( dp ), parameter :: utzero = 18.66_dp
  real ( dp ) y
  real ( dp ) z

  up = upper
  z = x

  if ( z < 0.0_dp ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0_dp
    else
      alnorm = 1.0_dp
    end if

    return

  end if

  y = 0.5_dp * z * z

  if ( z <= con ) then

    alnorm = 0.5_dp - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 & 
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0_dp - alnorm
  end if

  return
end function alnorm

end module SDP_routines

