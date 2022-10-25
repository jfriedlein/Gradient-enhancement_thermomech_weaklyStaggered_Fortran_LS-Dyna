      subroutine umat41 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c ...
        real*8 d_loc_n1, d_loc_n
        real*8 d_n1, d_n
        real*8 d_NL_n1, d_NL_n
        real*8 K_NL
c ...
c Compute the truly local damage "d_loc_n1" with whatever damage law you like, as before
        d_loc_n1 = ...
c
c The temperature "temper" contains the old non-local damage variable (this requires ITHERM=1)
        d_NL_n = temper
c Retrieve the old truly local damage "d_loc_n" from the history "hsv"
        d_loc_n = hsv(1)
c Compute the non-local factor "K_NL"
        K_NL = d_NL_n / d_loc_n
c To avoid issues with division by zero and large undesired fluctuations
c  the non-local factor K_NL is possibly corrected as follows:
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
c Update the non-local damage "d_NL_n1"
        d_NL_n1 = K_NL * d_loc_n1
c Retrieve the old local damage variable "d_n" from the history "hsv"
        d_n = hsv(2)
c Compute the new local damage "d_n1"
        d_n1 = max( max( d_n, d_NL_n ), d_NL_n1 )
c
c Use the new local damage "d_n1" to damage the material behaviour, as before
c e.g. stress = (1-d_n1) * stress
c ...
c
c Store as history variable "hsv":
c  the truly local damage "d_loc_n1"
        hsv(1) = d_loc_n1
c  the new local damage variable "d_n1"
        hsv(2) = d_n1
c
      return
      end

