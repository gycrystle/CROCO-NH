












! $Id: cppdefs.h 1628 2015-01-10 13:53:00Z marchesiello $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!




                      
                      
                      
                      

                      



! $Id: set_global_definitions.h 1618 2014-12-18 14:39:51Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!


 







 



 
  

















































 

 

 
!!#define NF_CLOBBER nf_clobber, nf_64bit_offset



! $Id: set_global_definitions.h 1618 2014-12-18 14:39:51Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!


























!-# define float dfloat
!-# define FLoaT dfloat
!-# define FLOAT dfloat
!-# define sqrt dsqrt
!-# define SQRT dsqrt
!-# define exp dexp
!-# define EXP dexp
!-# define dtanh dtanh
!-# define TANH dtanh



 




MODULE par_pisces
   !!======================================================================
   !!                        ***  par_pisces  ***
   !! TOP :   set the PISCES parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 2.0 , LOCEAN-IPSL (2007) 
   !! $Id: par_pisces.F90 1152 2008-06-26 14:11:13Z rblod $ 
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   IMPLICIT NONE
   PUBLIC

   INTEGER, PUBLIC, PARAMETER ::   jp_l      = 0      !: cumulative number of already defined TRC
   INTEGER, PUBLIC, PARAMETER ::   jp_l_2d   = 0   !:
   INTEGER, PUBLIC, PARAMETER ::   jp_l_3d   = 0   !:
   INTEGER, PUBLIC, PARAMETER ::   jp_l_trd  = 0  !:

   !!---------------------------------------------------------------------
   !!   Default                                   No CFC geochemical model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .FALSE.  !: CFC flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_kriest     = .FALSE.  !: Kriest flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     =  0       !: No CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  =  0       !: No CFC additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  =  0       !: No CFC additional 3d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_trd =  0       !: number of sms trends for PISCES

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jptra       = jp_pisces                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0     = jp_l + 1                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1     = jp_l + jp_pisces          !: Last  index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_2d  = jp_l_2d + 1               !: First index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_2d  = jp_l_2d + jp_pisces_2d    !: Last  index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_3d  = jp_l_3d + 1               !: First index of 3D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_3d  = jp_l_3d + jp_pisces_3d    !: Last  index of 3d diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_trd = jp_l_trd + 1              !: First index of bio diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_trd = jp_l_trd + jp_pisces_trd  !: Last  index of bio diag

   !!======================================================================
END MODULE par_pisces
