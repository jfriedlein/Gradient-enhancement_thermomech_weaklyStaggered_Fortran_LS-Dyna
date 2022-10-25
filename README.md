# Gradient-enhancement_thermomech_weaklyStaggered_Fortran_LS-Dyna
Repository containing Fortran source code for the implementation of gradient-enhancement by weakly staggered thermomechanical coupling especially suitable for commerical FEM-Software, exemplified for LS-Dyna

Here, exemplified for a gradient-damage model, so extending a scalar damage variable by a gradient. Also applicable for other scalar non-local variables.

We focus on implicit time integration, but briefly refer to explicit time integration (see [Gradient enhanced damage: modelling, implementation and applications](https://www.researchgate.net/publication/357376191_Gradient_enhanced_damage_modelling_implementation_and_applications) for more details).

## Overview of implementation for implicit time integration
The files `umat41.f` and `thumat11.f` (both from ls-dyna_smp_d_R12_1_0_x64_centos610_ifort160) show the key parts of the implementation and are detailed in the following.

<img src="https://github.com/jfriedlein/Gradient-enhancement_thermomech_weaklyStaggered_Fortran_LS-Dyna/blob/main/images/Gradient-enhancement%20-%20Numerical%20implementation%20-%20overview.png" width="1000">

By this way you can directly apply the gradient-enhancement in LS-Dyna


### Thermal user-defined material model
A huge thanks to Marvin Nahrmann and Henning Schmidt for initially providing this file ([Details on the implementation](https://www.researchgate.net/publication/357376191_Gradient_enhanced_damage_modelling_implementation_and_applications), [Application and comparison](https://doi.org/10.1016/j.ijsolstr.2021.111166)).

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
       real*8 L, d_loc
```
We retrieve the internal length `L` from the thermal material parameters (entered as p1... on *MAT_THERMAL_USER_DEFINED). Here, the length parameter is the first parameter p1, so we extract the first entry "1". For details on the leading "8+", see docu in "dyn21tumat.f - thusrmat".
```fortran
       L = r_matp(8+ 1 )
```
The truly local damage variable `d_loc` is, in this example (see `umat41.f`), stored as the first mechanical history variable "1". In the thermal material model the mechanical history is named `hsvm`. To be able to access this, we require IHVE=1 (see LS-Dyna manual Vol. I Appendix H), so we check this first:
```fortran
       if ( r_matp(4) .NE. 1 ) then
            write(*,*) "thumat11<< ERROR. The gradient-enhancement
     &requires IHVE=1 on *MAT_THERMAL_USER_DEFINED"
            call cstop('E R R O R  T E R M I N A T I O N')
       endif
       d_loc = hsvm(7+ 1 )
```
The thermal conductivity named `c1` is set to our length parameter squared (see partial differential equation for regularisation in the figures).
```fortran
       c1 = L**2
```
All the other conductivity parameters are set to zero to deactivate their contributions.
```fortran
       c2=0.
       c3=0.
       cvl=0.
```
For the regularisation, we require the source term, thus we activate it by setting ihsrcl to 1.
```fortran
       ihsrcl = 1
```
Next, we enter the source term as the difference between the temperature `temp` (our non-local variable) and the truly local damage `d_loc`.
```fortran
       hsrcl = -temp + d_loc
```
The thermal time steps are solved implicitly, so we require the linearisation of the source term. Here this is trivial.
```fortran
       dhsrcdtl = -1.
```
And this is it.
```fortran
       return
       end
```
You can leave this subroutine unchanged for all your user material models. The only thing that might change is the location of your truly local variable `d_loc` in hsvm. (Recommendation: use a [history variable manager](https://github.com/jfriedlein/history_hsv-manager_LS-Dyna))


### User-defined material model
Concering the user-defined material model `umat`, we focus on the changes that are needed to incorporate the gradient-enhancement. By adding the following lines of code, you can easily transform a local damage material model to a gradient-enhanced model (necessary keyword cards are listed below).

We start in the umat41 subroutine with the default arguments.
```fortran
      subroutine umat41 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
```
We declare different damage variables, namely the truly local damage variable `d_loc` (what you compute with your damage law or evolution equation), the local damage variable `d` (that actually affects the material behaviour), and the non-local damage variable `d_NL` (which is equal to the pseudo temperature and the non-locally smoothed/distributed damage). Additionally, we use affix `_n` to indicate values from the last converged load step n, and `_n1` to describe the new/current values from this step.
```fortran
        real*8 d_loc_n1, d_loc_n
        real*8 d_n1, d_n
        real*8 d_NL_n1, d_NL_n
        real*8 K_NL
```
After computing everything as usual in your material model and now computing your damage by whatever damage law you like, you now assign this value to the current truly local damage `d_loc_n1`
```fortran
        d_loc_n1 = ...
```
The input argument `temper` contains our old non-local damage variable `d_NL_n`. A note on efficiency: Clearly it would be more efficient to directly use some of the variables or simplify parts of the code (tell me about improvements in the "Issues" sections of the repository). However, here the focus is on explaining the framework and being overly explicit.
```fortran
        d_NL_n = temper
```       
The old truly local damage variable needs to be stored in the history variables, here at position "1".
```fortran
        d_loc_n = hsv(1)
```  
Next, we compute the non-local factor `K_NL` as the ratio of the old non-local damage and the old truly local damage.
```fortran
        K_NL = d_NL_n / d_loc_n
```  
To avoid numerical difficulties, like division by zero or exploding values and fluctuations, we correct/limit the non-local factor. Again, here we try to explain the thought-process and not optimise the expressions.
```fortran
        ! if old non-local variable is zero, we have to use the truly local
        !  else e.g. the damage would always be delayed
         if ( d_NL_n < 1e-14 ) then
           K_NL = 1.
        ! if old truly local variable is zero, we have to use the truly local
         elseif ( d_loc_n < 1e-14 ) then
           K_NL = 1.
        ! If damage shall be spread (K_NL>1), better do not scale the truly local variable to avoid extremely large values. We rely on the below "max(max(..,phi_n1),...)" for spreading.
         elseif ( K_NL > 1. ) then
           K_NL = 1.
         endif
```  
After the non-local factor is possibly corrected, we mimic the change of the non-local variable to obtain its (approximate) current value `d_NL_n1` as
```fortran
        d_NL_n1 = K_NL * d_loc_n1
```  
We retrieve the old local damage `d_n` from the history variable to ensure that the new damage (to be computed next) is larger or equal to `d_n`. This ensures the irreversibility of damage ("damage can only increase").
```fortran
        d_n = hsv(2)
```  
We compute the new local damage variable `d_n1` as
```fortran
        d_n1 = max( max( d_n, d_NL_n ), d_NL_n1 )
```  
This ensures that `d_n1>=d_n` (no decrease of damage), that damage grows as corrected by the regularisation `d_NL_n1`, and damage also spreads by `d_NL_n` into areas without truly local damage (`d_loc_n1=1` therefore also `d_NL_n1=0`).
Now you can use `d_n1` as the damage variable that affects the material behaviour, e. g. decreases the stiffness.

Lastly, we store the current truly local damage `d_loc_n1` and the current local damage `d_n1` into their positions of the history variable. (Recommendation: use a [history variable manager](https://github.com/jfriedlein/history_hsv-manager_LS-Dyna))

### Further notes
Do not forget to alter the tangent in `utan` and if needed the linearisation of your local iterations inside the material model, here the factor K_NL can be handy.

## Overview of implementation for explicit time integration
The thermal user-defined material model is the same for explicit time integration. The changes in the user-defined material model are even simpler, because the explicit time steps are so small that hardly anything happens between them. Therefore, mimicking the change of the non-local variable by the non-local factor appears to be not necessary. Thus, you can directly use the non-local variable to affect the material behaviour. Details for explicit time integration see: [Gradient enhanced damage: modelling, implementation and applications](https://www.researchgate.net/publication/357376191_Gradient_enhanced_damage_modelling_implementation_and_applications).

## Application and keyword cards
Firstly, we need to activate a thermo-mechanical solution, that solves the mechanical (umat) and thermal (thumat) problem together, by setting soln=2.
```
*CONTROL_SOLUTION
$#    soln       nlq     isnan     lcint     lcacc     ncdcf     
         2         0         0       100         0         1
```
For the thermal user material thumat11, we need the card *MAT_THERMAL_USER_DEFINED with "ihve=1" (to be able to access the mechanical history `d_loc` in the thermal routine. Further, we set the parameter "p1" to 2, which is equal to our internal length L (depending on the unit system the internal length L is here for instance 2 mm).
```
*MAT_THERMAL_USER_DEFINED
$#    tmid        ro        mt       lmc       nhv      aopt    iortho      ihve
        11       0.0        11         8         0       0.0         0         1
$#      p1        p2        p3        p4        p5        p6        p7        p8
         2       0.0       0.0       0.0       0.0       0.0       0.0       0.0
```
To properly work for the regularisation, the following two cards need to be set up as follows:
The generation of heat (temperature) from plastic deformation as eqheat and fwork are deactived by setting them to "1e-60" (numerically 0). Note that entering "0" would render the default value of "1", which does not deactive this.
```
*CONTROL_THERMAL_SOLVER
$#   atype     ptype    solver     cgtol       gpt    eqheat     fwork       sbc
         1         1         11.00000E-4         0     1e-60     1e-60       0.0
$#  msglvl    maxitr    abstol    reltol     omega    unused    unused       tsf
         0       5001.0000E-101.00000E-4       1.0                           1.0
$#   mxdmp      dtvf    varden      
         0       0.0         0
```
Use "tip=1.0" to avoid intermediate steps, which would require a reduced source term.
The initial thermal time step size (its) is here set to "0.1". Note that you need to set this according to your simulation model (ideally equal to mechanical time step size).
Use "dtemp=1e6" to make sure that the thermal step size is never reduced, which would require a reduced source term; a reduced thermal step is not helpful to increase the accuracy, but would only be necessary in case of convergence failure.
```
*CONTROL_THERMAL_TIMESTEP
$#      ts       tip       its      tmin      tmax     dtemp      tscp      lcts
         0       1.0       0.1       0.0         1       1e6       0.5         0
```
Additionally, you need your user-defined material by the keyword card *MAT_USER_DEFINED_MATERIAL_MODELS_TITLE. Here, "itherm=1" is needed to make the temperature from the thermal solver available as "temper" in umat.
```
*MAT_USER_DEFINED_MATERIAL_MODELS_TITLE
User-defined material model
$#     mid        ro        mt       lmc       nhv    iortho     ibulk        ig
        107.83000E-6        41         8         2         0         3         4
$#   ivect     ifail    itherm    ihyper      ieos      lmca    unused    unused
         0         0         1         0         0         0                    
$#      p1        p2        p3        p4        p5        p6        p7        p8
         ...
```
Lastly, you need to link the mechanical umat as "mid=10" and the thermal thumat as "tmid=11" to your part in *PART.
```
*PART
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         1        10         0         0         0         0        11
```

## How to make it better

<img src="https://github.com/jfriedlein/Gradient-enhancement_thermomech_weaklyStaggered_Fortran_LS-Dyna/blob/main/images/Gradient-enhancement%20-%20Numerical%20implementation%20II%20-%20overview.png" width="1000">

+ Thermal solver is again available for computing the temperature, etc.

- Not yet available in LS-Dyna (please motivate them to add this)


## Notes
R12.0 is goofy regarding the thermal solver
In newer/older releases some variables might be named differently or located in different files, be aware of this!


