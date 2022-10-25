# Gradient-enhancement_thermomech_weaklyStaggered_Fortran_LS-Dyna
Repository containing Fortran source code for the implementation of gradient-enhancement by weakly staggered thermomechanical coupling especially suitable for commerical FEM-Software, exemplified for LS-Dyna

## Overview of implementation
The files `umat41.f` and `thumat11.f` (both from ls-dyna_smp_d_R12_1_0_x64_centos610_ifort160) show the key parts of the implementation and are detailed in the following.

<img src="https://github.com/jfriedlein/Gradient-enhancement_thermomech_weaklyStaggered_Fortran_LS-Dyna/blob/main/images/Gradient-enhancement%20-%20Numerical%20implementation%20-%20overview.png" width="1000">

By this way you can directly apply the gradient-enhancement in LS-Dyna


### Thermal user-defined material model
Here, we exemplarily use the thermal user-defined material model 11 with the default arguments.
```fortran
     subroutine thumat11(c1,c2,c3,cvl,dcvdtl,hsrcl,dhsrcdtl,
     1 hsv,hsvm,nmecon,r_matp,crv,nnpcrv,
     2 nel,nep,iep,eltype,dt,atime,ihsrcl,temp,ttconv,ttnext,
     3 tmconv,tmnext,inithis)
c
      include 'nlqparm'
      include 'iounits.inc'
c
      character*(*) eltype
      logical inithis
      dimension hsv(*),hsvm(*),r_matp(*),crv(lq1,2,*)
      integer nnpcrv(*)    
c
```
Next, we declare two additional variables, namely the internal length parameter `L` and the truly local damage variable `d_loc`.
```fortran


test



## How to make it better

<img src="https://github.com/jfriedlein/Gradient-enhancement_thermomech_weaklyStaggered_Fortran_LS-Dyna/blob/main/images/Gradient-enhancement%20-%20Numerical%20implementation%20II%20-%20overview.png" width="1000">

+ Thermal solver is again available for computing the temperature, etc.

- Not yet available in LS-Dyna (please motivate them to add this)

## Notes
R12.0 is goofy regarding the thermal solver


