*KEYWORD
*CONTROL_SOLUTION
$#    soln       nlq     isnan     lcint     lcacc     ncdcf     
         2         0         0       100         0         1
*CONTROL_THERMAL_SOLVER
$ eqheat&fwork: The generation of heat (temperature) from plastic deformation as eqheat and fwork are deactived by setting them to "1e-60" (numerically 0). Note that entering "0" would render the default value of "1", which does not deactive this.
$#   atype     ptype    solver     cgtol       gpt    eqheat     fwork       sbc
         1         1         11.00000E-4         0     1e-60     1e-60       0.0
$#  msglvl    maxitr    abstol    reltol     omega    unused    unused       tsf
         0       5001.0000E-101.00000E-4       1.0                           1.0
$#   mxdmp      dtvf    varden      
         0       0.0         0
*CONTROL_THERMAL_TIMESTEP
$ tip: use "1.0" to avoid intermediate steps, which would require a reduced source term
$ its: the initial thermal time step size (its) is here set to "0.1". Note that you need to set this according to your simulation model (ideally equal to mechanical time step size).
$ dtemp: use "1e6" to make sure that the thermal step size is never reduced, which would require a reduced source term; a reduced thermal step is not helpful to increase the accuracy, but would only be necessary in case of convergence failurecase of convergence failure
$#      ts       tip       its      tmin      tmax     dtemp      tscp      lcts
         0       1.0         1       0.0         1       1e6       0.5         0
*DATABASE_TPRINT
$#      dt    binary      lcur     ioopt     
       1.0         3         0         1
*END
