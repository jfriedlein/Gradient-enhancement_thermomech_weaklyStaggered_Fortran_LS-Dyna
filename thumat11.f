     subroutine thumat11(c1,c2,c3,cvl,dcvdtl,hsrcl,dhsrcdtl,
     1 hsv,hsvm,nmecon,r_matp,crv,nnpcrv,
     2 nel,nep,iep,eltype,dt,atime,ihsrcl,temp,ttconv,ttnext,
     3 tmconv,tmnext,inithis)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c For this thermal user material to properly work, the following two cards need to be used (in addition to others):
c *CONTROL_THERMAL_SOLVER
c $ eqheat&fwork: The generation of heat (temperature) from plastic deformation as eqheat and fwork are deactived by setting them to "1e-60" (numerically 0). Note that entering "0" would render the default value of "1", which does not deactive this.
c $#   atype     ptype    solver     cgtol       gpt    eqheat     fwork       sbc
c          1         1         11.00000E-4         0     1e-60     1e-60       0.0
c $#  msglvl    maxitr    abstol    reltol     omega    unused    unused       tsf
c          0       5001.0000E-101.00000E-4       1.0                           1.0
c $#   mxdmp      dtvf    varden      
c          0       0.0         0
c *CONTROL_THERMAL_TIMESTEP
c $ tip: use "1.0" to avoid intermediate steps, which would require a reduced source term
c $ its: the initial thermal time step size (its) is here set to "0.1". Note that you need to set this according to your simulation model (ideally equal to mechanical time step size).
c $ dtemp: use "1e6" to make sure that the thermal step size is never reduced, which would require a reduced source term; a reduced thermal step is not helpful to increase the accuracy, but would only be necessary in case of convergence failure
c $#      ts       tip       its      tmin      tmax     dtemp      tscp      lcts
c          0       1.0       0.1       0.0         1       1e6       0.5         0
c
c******************************************************************
c
      include 'nlqparm'
      include 'iounits.inc'
c
      character*(*) eltype
      logical inithis
      dimension hsv(*),hsvm(*),r_matp(*),crv(lq1,2,*)
      integer nnpcrv(*)    
c       
       real*8 L, d_loc
c
c Internal length as parameter 1 of thermal material card
       L = r_matp(8+ 1 )
c       
c Retrieve truly local damage (entry 1, see umat41) from mechanical subroutine history, herein named "hsvm". To be able to do this, we require IHVE=1, so we check this first:
       if ( r_matp(4) .NE. 1 ) then
            write(*,*) "thumat11<< ERROR. The gradient-enhancement
     &requires IHVE=1 on *MAT_THERMAL_USER_DEFINED"
            call cstop('E R R O R  T E R M I N A T I O N')
       endif
       d_loc = hsvm(7+ 1 )
c
c Thermal conductivity (parameter in front of Laplacian of temperature)
       c1 = L**2
c
c Other parameters of thermal conductivity are zero
       c2=0.
       c3=0.
       cvl=0.
c       
c Activate the source term
       ihsrcl = 1
c       
c Enter the source term as difference between temperature and truly local damage
       hsrcl = -temp + d_loc
c      
c Derivative of source term hsrcl with respect to temperature temp
        dhsrcdtl = -1.
c
       return
       end
