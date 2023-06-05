$# LS-DYNA Keyword file created by LS-PrePost(R) V4.7.20 -24Aug2020
$# Created on Jan-11-2022 (13:25:20)
*KEYWORD
*PARAMETER
$
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
Rtmax     100                                                                 
*BOUNDARY_PRESCRIBED_MOTION_SET_ID
$#      id                                                               heading
         1BC_uy
$#    nsid       dof       vad      lcid        sf       vid     death     birth
         1         2         2         3       1.0         01.00000E28       0.0
*CONSTRAINED_GLOBAL
$#      tc        rc       dir         x         y         z       tol  
         1         0         2       0.0       0.0       0.0       0.0
         2         0         2       0.0       0.0       0.0       0.0 
         3         0         3       0.0       0.0       0.0       0.0
         3         0         3       0.0       0.0       1.0       0.0
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol     
&tmax              0       0.0       0.01.000000E8         0
*DATABASE_CROSS_SECTION_PLANE_ID
$#    csid                                                                 title
         2                                                                      
$#    psid       xct       yct       zct       xch       ych       zch    radius
         0       0.0      14.0       1.0       0.0      15.0       1.0       0.0
$#    xhev      yhev      zhev      lenl      lenm        id     itype     
       1.0      14.0       1.0       5.0      1.25         0         0
*DEFINE_CURVE_TITLE
Curve_load
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         3         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                 0.0
&tmax                                3.0
*PART
$#                                                                         title
auto created
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         1       244         0         0         0         0         1
*SECTION_SOLID_TITLE
tensileSpeci
$#   secid    elform       aet 
         1        -1         0
*SET_PART_LIST
$#     sid       da1       da2       da3       da4    solver      
         1       0.0       0.0       0.0       0.0MECH
$#    pid1      pid2      pid3      pid4      pid5      pid6      pid7      pid8
         1         0         0         0         0         0         0         0
*END
