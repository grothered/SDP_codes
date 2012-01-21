This repo contains routines for a version of the greedy forester type model that I was working on with BW. It can be run either as straight fortran code (in which case there is an openMP parallel version), or via R (which requires compiling SDP_take2.f90 as a shared library). From memory the R version has a similar speed to the single processor fortran version. 



Why I found this numerically difficult:
The total future return function has a discontinuity in it at x=mu_crit -- because if logging is permitted in the present time step, then the population will fall for the next time step, discontinuously affecting the future total return.

Mathematically, this discontinuity stops numerical integrators from rapidly converging to the correct value of the integral -- even supposedly high order methods will only be first order here. However, we know that the discontinuity occurs at mu_crit (I have checked with simulation too) -- and so if we integrate first up to mu_crit, and then above, we should get more accurate results.

It seems that we do -- very nice results with 1001 integration points. We still need a large number, because in places the TR function is steep -- perhaps some adaptive methods would work better here -- or even just placing more points near the steep sections? 

One question is whether we are better off having more integration points, or a bicubic interpolation of TR_j+1. The answer seems to be: More points for the integration.
