












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



 



 module module_parallel_nbq
  implicit none
  integer, parameter :: ouest=1,est=2,nord=3,sud=4,haut=5,bas=6
  integer, parameter :: sudouest=7,sudest=8,nordouest=9,nordest=10
  integer, parameter :: ouestest=1,nordsud=2

  integer,parameter :: MPI_PROC_NULL=-2

  type infopar_croco
         integer ::  comm2d                 !COMMUNICATEUR GLOBAL
         integer ::  rank
         integer,dimension(10)                      ::  tvoisin
  end type 

!  integer,dimension(8),parameter :: liste_voisin = &
!      (/ ouest, est, nord, sud, sudouest, sudest, nordouest, nordest /)
  integer :: ierr,mynode
  type(infopar_croco) :: par

 end module module_parallel_nbq
