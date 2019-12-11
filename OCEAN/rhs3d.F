! $Id: rhs3d.F 1458 2014-02-03 15:01:25Z gcambon $
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
#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine rhs3d (tile)
!
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call rhs3d_tile (Istr,Iend,Jstr,Jend,
     &                                    A3d(1,1,trd), A3d(1,2,trd),
# ifdef NHMG
     &                                    A3d(1,5,trd),
# endif
     &                      A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                      A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                      A2d(1,4,trd), A2d(1,5,trd), A2d(1,6,trd))
      return
      end

      subroutine rhs3d_tile (Istr,Iend,Jstr,Jend, ru,rv, 
# ifdef NHMG
     &                                                   rw,
# endif
     &                                                   CF,FC,DC,
     &                                 wrk1,wrk2, UFx,UFe, VFx,VFe)
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k
# ifdef MPI
     &                   ,imin,imax,jmin,jmax
# endif
# include "param.h"
# ifdef DIAGNOSTICS_UV
#  include "diagnostics.h"
# endif
      real ru(PRIVATE_2D_SCRATCH_ARRAY,N),
     &     rv(PRIVATE_2D_SCRATCH_ARRAY,N),
     &     CF(PRIVATE_1D_SCRATCH_ARRAY,0:N),  cff,cff1,cff2,
     &     FC(PRIVATE_1D_SCRATCH_ARRAY,0:N),  cffX,
     &     DC(PRIVATE_1D_SCRATCH_ARRAY,0:N),  cffE,      
     &     wrk1(PRIVATE_2D_SCRATCH_ARRAY),    curvX,
     &     wrk2(PRIVATE_2D_SCRATCH_ARRAY),    curvE,
     &      UFx(PRIVATE_2D_SCRATCH_ARRAY),
     &      UFe(PRIVATE_2D_SCRATCH_ARRAY),
     &      VFx(PRIVATE_2D_SCRATCH_ARRAY), 
     &      VFe(PRIVATE_2D_SCRATCH_ARRAY),    gamma 
# ifdef NHMG
     &    ,rw(PRIVATE_2D_SCRATCH_ARRAY,0:N)
# endif
      parameter (gamma = -0.25)
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# include "scalars.h"
# ifdef NBQ
#  include "nbq.h"
# endif
# ifdef MRL_WCI
      real vstu,dudx,dvdx,ustv,dude,dvde
# endif
!!!! 29/11/2019, @gyc : vadv
!--------------------------------------------------------------------
! Definition of flux operators: 1st, 2nd, 3rd, 4th, 5th or 6th order,
! used in UP5 and C6 advection schemes (and order degradation near 
! land masks). cdiff is part of laplacian diffusion in flux1 (used 
! near mask): 
!    0 --> flux1=flux2 (second order C2 advection scheme)
!    1 --> flux1 gives 1st order monotonic UP1 advection scheme
!--------------------------------------------------------------------
!
# if defined UV_HADV_WENO5 || defined UV_HADV_C6 || defined UV_HADV_UP5 \
  || defined UV_VADV_WENO5 || defined UV_VADV_C6 || defined UV_VADV_UP5

      REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2
      REAL    :: ua, vel, cdiff, cdif
      REAL    :: flux1, flux2, flux3, flux4, flux5, flux6
      REAL    :: flx2, flx3, flx4, flx5
      REAL    :: mask0, mask1, mask2, mask3

      flux2(q_im1, q_i, ua, cdiff) = 0.5*( q_i + q_im1 )

      flux1(q_im1, q_i, ua, cdiff) = flux2(q_im1, q_i, ua, cdiff) -
     &      0.5*cdiff*sign(1.,ua)*(q_i-q_im1)

      flux4(q_im2, q_im1, q_i, q_ip1, ua) = 
     &      ( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua) =
     &      flux4(q_im2, q_im1, q_i, q_ip1, ua) +
     &      sign(1.,ua)*((q_ip1 -
     &      q_im2)-3.*(q_i-q_im1))/12.0

      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      ( 37.*(q_i+q_im1) - 8.*(q_ip1+q_im2)
     &      +(q_ip2+q_im3) )/60.0

      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
     &      -sign(1.,ua)*( 
     &      (q_ip2-q_im3)-5.*(q_ip1-q_im2)+10.*(q_i-q_im1) )/60.0
# endif
# if  defined UV_HADV_WENO5 || defined UV_VADV_WENO5
      REAL    :: flux3_weno, flux5_weno
# endif
!
!!!! end of 29/11/2019, @gyc : vadv
!
# include "compute_auxiliary_bounds.h"
!
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
# if defined UV_HADV_TVD || defined UV_VADV_TVD
      if(nrhs.eq.3) then
        cdt =     dt                 !<-- Corrector
      elseif(FIRST_TIME_STEP) then
        cdt = 0.5*dt
      else
        cdt = (1.-1./6.)*dt        !<-- Predictor
      endif
# endif

# ifdef BODYFORCE
!====================================================================
!
!  Apply surface and bottom stresses as body forces if selected
!
!====================================================================
!
#  define Uwrk UFx
#  define Vwrk VFe
#  define wrk  UFe
!
!  Apply surface stress as a bodyforce: determine the thickness [m]
!  of the surface layer; then add in surface stress as a bodyfoce.
!
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          wrk(i,j)=0.
        enddo
      enddo
      do k=N,levsfrc,-1
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            wrk(i,j)=wrk(i,j)+Hz(i,j,k)
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          Uwrk(i,j)=sustr(i,j)*4./((wrk(i-1,j)+wrk(i,j))
     &                             *(pm (i-1,j)+pm (i,j))
     &                             *(pn (i-1,j)+pn (i,j)))
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          Vwrk(i,j)=svstr(i,j)*4./((wrk(i,j-1)+wrk(i,j))
     &                             *(pm (i,j-1)+pm (i,j))
     &                             *(pn (i,j-1)+pn (i,j)))
        enddo
      enddo             !--> discard wrk
      do k=levsfrc,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru(i,j,k)=ru(i,j,k)+Uwrk(i,j)*(Hz(i,j,k)+
     &                                     Hz(i-1,j,k))
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rv(i,j,k)=rv(i,j,k)+Vwrk(i,j)*(Hz(i,j,k)+
     &                                     Hz(i,j-1,k))
          enddo
        enddo
      enddo          !--> discard Uwrk,Vwrk
!
!  Apply bottom stress as a bodyforce: determine the thickness [m]
!  of the bottom layer; then add in bottom stress as a bodyfoce.
!
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          wrk(i,j)=0.
        enddo
      enddo
      do k=1,levbfrc
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            wrk(i,j)=wrk(i,j)+Hz(i,j,k)
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          Uwrk(i,j)=bustr(i,j)*4./((wrk(i-1,j)+wrk(i,j))
     &                             *(pm (i-1,j)+pm (i,j))
     &                             *(pn (i-1,j)+pn (i,j)))
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          Vwrk(i,j)=bvstr(i,j)*4./((wrk(i,j-1)+wrk(i,j))
     &                             *(pm (i,j-1)+pm (i,j))
     &                             *(pn (i,j-1)+pn (i,j)))
        enddo
      enddo             !--> discard wrk
      do k=1,levbfrc
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru(i,j,k)=ru(i,j,k)-Uwrk(i,j)*(Hz(i,j,k)+
     &                                               Hz(i-1,j,k))
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rv(i,j,k)=rv(i,j,k)-Vwrk(i,j)*(Hz(i,j,k)+
     &                                               Hz(i,j-1,k))
          enddo
        enddo
      enddo         !--> discard Uwrk,Vwrk
#  undef wrk
#  undef Vwrk
#  undef Uwrk
# endif /* BODYFORCE */

      do k=1,N
!
# if defined UV_COR || (defined CURVGRID && defined UV_ADV)
!
!====================================================================
!
!  Add in Coriolis, Stokes-Coriolis and curvilinear transformation 
!  terms, if any.
!
!====================================================================
!
        do j=JstrV-1,Jend+1
          do i=IstrU-1,Iend+1
            cff=0.5*Hz(i,j,k)*(
#  ifdef UV_COR
     &              fomn(i,j)
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     &             +0.5*( (v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*dndx(i,j)
     &                   -(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*dmde(i,j))
#  endif 
     &                                                             )
#  ifdef MRL_WCI
#    if defined CURVGRID && defined UV_ADV
            cff1=0.5*Hz(i,j,k)*(
     &           0.5*( dndx(i,j)*(vst(i,j,k)+vst(i,j+1,k))
     &                -dmde(i,j)*(ust(i,j,k)+ust(i+1,j,k)) ))
#    else
            cff1 = 0.0
#    endif
            UFx(i,j)=(cff+cff1)*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))  
            VFe(i,j)=(cff+cff1)*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
            UFe(i,j)=cff*(vst(i,j,k)+vst(i,j+1,k))
            VFx(i,j)=cff*(ust(i,j,k)+ust(i+1,j,k))
#  else
            UFx(i,j)=cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
#   if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_VRT || \
       defined DIAGNOSTICS_EK || defined DIAGNOSTICS_PV
            UFxCor(i,j)=0.5*Hz(i,j,k)*fomn(i,j)
     &                 *(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            VFeCor(i,j)=0.5*Hz(i,j,k)*fomn(i,j)
     &                 *(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
#   endif
#  endif
          enddo
        enddo

        do j=Jstr,Jend
          do i=IstrU,Iend+1
            ru(i,j,k)=ru(i,j,k)+0.5*(UFx(i,j)+UFx(i-1,j))
#  ifdef MRL_WCI
     &                         +0.5*(UFe(i,j)+UFe(i-1,j))
#  endif

# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
            if (nnew.ne.3) then
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL
              MCor(i,j,k,1) = 0.5*(UFxCor(i,j)+UFxCor(i-1,j))
     &                        SWITCH umask(i,j)
# endif
              MXadv(i,j,k,1) = (0.5*(UFx(i,j)+UFx(i-1,j))
     &                       - 0.5*(UFxCor(i,j)+UFxCor(i-1,j)))
     &                        SWITCH umask(i,j)
#   ifdef MRL_WCI
              MStCo(i,j,k,1) = 0.5*(UFe(i,j)+UFe(i-1,j))
     &                        SWITCH umask(i,j)
#   endif
            endif
# elif defined DIAGNOSTICS_EK
            if (nnew.ne.3) then
              if (k.eq.1) then
                ekwrkCor(i,j,1) = 0.5*(UFxCor(i,j)+UFxCor(i-1,j))
     &                    * u(i,j,k,nrhs) SWITCH umask(i,j)
                ekwrkHadv(i,j,1) = (0.5*(UFx(i,j)+UFx(i-1,j))
     &                       - 0.5*(UFxCor(i,j)+UFxCor(i-1,j)))
     &                    * u(i,j,k,nrhs) SWITCH umask(i,j)
              else
                ekwrkCor(i,j,1)=ekwrkCor(i,j,1)
     &                         +0.5*(UFxCor(i,j)+UFxCor(i-1,j))
     &                    * u(i,j,k,nrhs) SWITCH umask(i,j)
                ekwrkHadv(i,j,1)=ekwrkHadv(i,j,1)
     &                       +(0.5*(UFx(i,j)+UFx(i-1,j))
     &                       - 0.5*(UFxCor(i,j)+UFxCor(i-1,j)))
     &                    * u(i,j,k,nrhs) SWITCH umask(i,j)
              endif
# if defined DIAGNOSTICS_EK_MLD
              if (k.eq.kbl(i,j)) then
                ekwrkCor_mld(i,j,1)=0.5*(UFxCor(i,j)+UFxCor(i-1,j))
     &                    * u(i,j,k,nrhs) SWITCH umask(i,j)
                ekwrkHadv_mld(i,j,1)=(0.5*(UFx(i,j)+UFx(i-1,j))
     &                       - 0.5*(UFxCor(i,j)+UFxCor(i-1,j)))
     &                    * u(i,j,k,nrhs) SWITCH umask(i,j)
              elseif (k.gt.kbl(i,j)) then
                ekwrkCor_mld(i,j,1)=ekwrkCor_mld(i,j,1)
     &                         +0.5*(UFxCor(i,j)+UFxCor(i-1,j))
     &                    * u(i,j,k,nrhs) SWITCH umask(i,j)
                ekwrkHadv_mld(i,j,1)=ekwrkHadv_mld(i,j,1)
     &                       +(0.5*(UFx(i,j)+UFx(i-1,j))
     &                       - 0.5*(UFxCor(i,j)+UFxCor(i-1,j)))
     &                    * u(i,j,k,nrhs) SWITCH umask(i,j)
              endif
# endif
            endif
#  endif
#  if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
            if (nnew.ne.3) then
              if (k.eq.1) then
                wrkCor(i,j,1) = 0.5*(UFxCor(i,j)+UFxCor(i-1,j))
     &                        SWITCH umask(i,j)
                wrkXadv(i,j,1) = (0.5*(UFx(i,j)+UFx(i-1,j))
     &                       - 0.5*(UFxCor(i,j)+UFxCor(i-1,j)))
     &                        SWITCH umask(i,j)
              else
                wrkCor(i,j,1) = wrkCor(i,j,1) 
     &                        + 0.5*(UFxCor(i,j)+UFxCor(i-1,j))
     &                        SWITCH umask(i,j)
                wrkXadv(i,j,1) = wrkXadv(i,j,1)
     &                       +(0.5*(UFx(i,j)+UFx(i-1,j))
     &                       - 0.5*(UFxCor(i,j)+UFxCor(i-1,j)))
     &                        SWITCH umask(i,j)

              endif
            endif
#  endif
          enddo
        enddo

        do j=JstrV,Jend+1
          do i=Istr,Iend
            rv(i,j,k)=rv(i,j,k)-0.5*(VFe(i,j)+VFe(i,j-1))
#  ifdef MRL_WCI
     &                         -0.5*(VFx(i,j)+VFx(i,j-1))
#  endif

# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL
            if (nnew.ne.3) then
              MCor(i,j,k,2) = -0.5*(VFeCor(i,j)+VFeCor(i,j-1))
     &                        SWITCH vmask(i,j)
              MXadv(i,j,k,2) = (-0.5*(VFe(i,j)+VFe(i,j-1))
     &                       +0.5*(VFeCor(i,j)+VFeCor(i,j-1)))
     &                        SWITCH vmask(i,j)
#   ifdef MRL_WCI
              MStCo(i,j,k,2) = -0.5*(VFx(i,j)+VFx(i,j-1))
     &                        SWITCH vmask(i,j)
#   endif
            endif
# elif defined DIAGNOSTICS_EK
            if (nnew.ne.3) then
              if (k.eq.1) then
              ekwrkCor(i,j,2) = -0.5*(VFeCor(i,j)+VFeCor(i,j-1))
     &                     * v(i,j,k,nrhs) SWITCH vmask(i,j)
              ekwrkHadv(i,j,2) = (-0.5*(VFe(i,j)+VFe(i,j-1))
     &                       +0.5*(VFeCor(i,j)+VFeCor(i,j-1)))
     &                     * v(i,j,k,nrhs) SWITCH vmask(i,j)
              else
              ekwrkCor(i,j,2)=ekwrkCor(i,j,2)
     &                     -0.5*(VFeCor(i,j)+VFeCor(i,j-1))
     &                     * v(i,j,k,nrhs) SWITCH vmask(i,j)
              ekwrkHadv(i,j,2)=ekwrkHadv(i,j,2) + (
     &                        -0.5*(VFe(i,j)+VFe(i,j-1))
     &                       +0.5*(VFeCor(i,j)+VFeCor(i,j-1)))
     &                     * v(i,j,k,nrhs) SWITCH vmask(i,j)
              endif
# if defined DIAGNOSTICS_EK_MLD
              if (k.eq.kbl(i,j)) then
              ekwrkCor_mld(i,j,2)=
     &                     -0.5*(VFeCor(i,j)+VFeCor(i,j-1))
     &                     * v(i,j,k,nrhs) SWITCH vmask(i,j)
              ekwrkHadv_mld(i,j,2)=(
     &                        -0.5*(VFe(i,j)+VFe(i,j-1))
     &                       +0.5*(VFeCor(i,j)+VFeCor(i,j-1)))
     &                     * v(i,j,k,nrhs) SWITCH vmask(i,j)
              elseif (k.gt.kbl(i,j)) then
              ekwrkCor_mld(i,j,2)=ekwrkCor_mld(i,j,2)
     &                     -0.5*(VFeCor(i,j)+VFeCor(i,j-1))
     &                     * v(i,j,k,nrhs) SWITCH vmask(i,j)
              ekwrkHadv_mld(i,j,2)=ekwrkHadv_mld(i,j,2) + (
     &                        -0.5*(VFe(i,j)+VFe(i,j-1))
     &                       +0.5*(VFeCor(i,j)+VFeCor(i,j-1)))
     &                     * v(i,j,k,nrhs) SWITCH vmask(i,j)
              endif
# endif
            endif
#  endif
# if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
            if (nnew.ne.3) then
              if (k.eq.1) then
                wrkCor(i,j,2) = -0.5*(VFeCor(i,j)+VFeCor(i,j-1))
     &                        SWITCH vmask(i,j)
                wrkXadv(i,j,2) = (-0.5*(VFe(i,j)+VFe(i,j-1))
     &                       +0.5*(VFeCor(i,j)+VFeCor(i,j-1)))
     &                        SWITCH vmask(i,j)
              else
                wrkCor(i,j,2) = wrkCor(i,j,2) 
     &                        -0.5*(VFeCor(i,j)+VFeCor(i,j-1))
     &                        SWITCH vmask(i,j)
                wrkXadv(i,j,2) = wrkXadv(i,j,2) + (
     &                       -0.5*(VFe(i,j)+VFe(i,j-1))
     &                       +0.5*(VFeCor(i,j)+VFeCor(i,j-1)))
     &                        SWITCH vmask(i,j)
              endif
            endif
#  endif
          enddo
        enddo
!
#  ifdef NONTRAD_COR
         ! wrk1 is the vertical momentum at RHO-points
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
             wrk1(i,j) = 0.5*Hz(i,j,k)*(
     &            wz(i,j,k,nrhs)+wz(i,j,k-1,nrhs)
     &           + dzdxi(i,j,k)*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
     &           +dzdeta(i,j,k)*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs)) )
          enddo
        enddo        
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru(i,j,k)=ru(i,j,k)
     &          - 0.5*( wrk1(i,j)*fetaomn(i,j) + wrk1(i-1,j)*fetaomn(i-1,j) )
          enddo
        enddo        
        do j=JstrV,Jend
          do i=Istr,Iend
            rv(i,j,k)=rv(i,j,k)
     &          + 0.5*( wrk1(i,j)*fxiomn(i,j) + wrk1(i,j-1)*fxiomn(i,j-1) )
          enddo
        enddo        
#endif

!
# endif /* UV_COR || (CURVGRID && UV_ADV) */

# ifdef UV_ADV
!
!======================================================================
!
!  Add in horizontal advection of momentum.
!  Compute diagonal [UFx,VFe] and off-diagonal [UFe,VFx] components
!  of tensor of momentum flux due to horizontal advection; after that
!  add in horizontal advection terms.
!
!======================================================================
!
#  if !(defined UV_HADV_UP3 || defined UV_HADV_C4  || \
        defined UV_HADV_C2  || defined UV_HADV_UP5 || \
        defined UV_HADV_C6  || defined UV_HADV_WENO5 )
#   define UV_HADV_UP3
! why ?
#  endif
!!!! 01/12/2019 @gyc: adding 5th order hadv 
#  if defined UV_HADV_UP5 || defined UV_HADV_C6 || defined UV_HADV_WENO5
!
!----------------------------------------------------------------------
! 5th or 6th order or WENO5 advection schemes
!----------------------------------------------------------------------
!
#   ifdef UV_HADV_WENO5
#    define FLUX5 flux5_weno
#    define FLUX3 flux3_weno
#    define FLUX2 flux1
#    define UP5_MASKING
          cdif=1.
#    include "u_hadv_order5.h"
#    include "v_hadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
#    undef UP5_MASKING

#   else  /* C6 or UP5 */

        if (nnew.eq.3) then  ! predictor
#    define FLUX5 flux6
#    define FLUX3 flux4
#    define FLUX2 flux2
#    undef  UP5_MASKING
          cdif=0. 
#    include "u_hadv_order5.h"
#    include "v_hadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
        else                 ! corrector
#    if UV_HADV_C6
#     define FLUX5 flux6
#     define FLUX3 flux4
#     define FLUX2 flux2
#     undef  UP5_MASKING
          cdif=0.
#    else
#     define FLUX5 flux5
#     define FLUX3 flux3
#     define FLUX2 flux2
#     define UP5_MASKING
          cdif=0.5
#    endif

#    include "u_hadv_order5.h"
#    include "v_hadv_order5.h"

#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
#    undef UP5_MASKING

#    if defined UV_HADV_UP5
#      if defined DIAGNOSTICS_UV  || defined DIAGNOSTICS_EK  || \
        defined DIAGNOSTICS_VRT || defined DIAGNOSTICS_PV
#        define FLUX5 flux6
#        define FLUX3 flux4
#        define FLUX2 flux2
           cdif=0.
#          include "uhdiff_order5.h"
#          include "vhdiff_order5.h"
#        undef FLUX5
#        undef FLUX3
#        undef FLUX2
#      endif
#    endif



        endif
#   endif /* WENO5 */
!!!! end of 01/12/2019 @gyc: adding 5th order hadv 
#  elif defined UV_HADV_UP3 || defined UV_HADV_C4 
!
!----------------------------------------------------------
! Fourth or Third order advection scheme (default)
!----------------------------------------------------------
!
#   define uxx wrk1
#   define Huxx wrk2
#   ifdef EW_PERIODIC
#    define IU_EXT_RANGE IstrU-1,Iend+1
#   else
#    ifdef MPI
        if (WEST_INTER) then
          imin=IstrU-1
        else
          imin=max(IstrU-1,2)
        endif
        if (EAST_INTER) then
          imax=Iend+1
        else
          imax=min(Iend+1,Lmmpi)
        endif
#     define IU_EXT_RANGE imin,imax
#    else
#     define IU_EXT_RANGE max(IstrU-1,2),min(Iend+1,Lm)
#    endif
#   endif 
!/* EW_PERIODIC */
        do j=Jstr,Jend
          do i=IU_EXT_RANGE
            uxx(i,j)=(u(i-1,j,k,nrhs)-2.*u(i,j,k,nrhs)
     &               +u(i+1,j,k,nrhs))  SWITCH umask(i,j)
            Huxx(i,j)=(Huon(i-1,j,k)-2.*Huon(i,j,k)
     &                +Huon(i+1,j,k))   SWITCH umask(i,j)
          enddo
        enddo
#   undef IU_EXT_RANGE
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=Jstr,Jend
#    if defined OBC_MODIF_UP1
            uxx(1,j)=(u(1,j,k,nrhs)-u(2,j,k,nrhs))/gamma
            Huxx(1,j)=-Huxx(2,j)
#    elif defined OBC_MODIF_C2
            uxx(1,j)=0.
            Huxx(1,j)=-Huxx(2,j)
#    else
            uxx(1,j)=uxx(2,j)
            Huxx(1,j)=Huxx(2,j)
#    endif
          enddo
        endif
        if (EASTERN_EDGE) then
#    ifdef MPI        
          do j=Jstr,Jend
#    if defined OBC_MODIF_UP1
            uxx(Lmmpi+1,j)=(u(Lmmpi+1,j,k,nrhs)-u(Lmmpi,j,k,nrhs))/gamma
            Huxx(Lmmpi+1,j)=-Huxx(Lmmpi,j)
#    elif defined OBC_MODIF_C2
            uxx(Lmmpi+1,j)=0.
            Huxx(Lmmpi+1,j)=-Huxx(Lmmpi,j)
#    else
            uxx(Lmmpi+1,j)=uxx(Lmmpi,j)
            Huxx(Lmmpi+1,j)=Huxx(Lmmpi,j)
#    endif
          enddo
#    else
          do j=Jstr,Jend
#    if defined OBC_MODIF_UP1
            uxx(Lm+1,j)=(u(Lm+1,j,k,nrhs)-u(Lm,j,k,nrhs))/gamma
            Huxx(Lm+1,j)=-Huxx(Lm,j)
#    elif defined OBC_MODIF_C2
            uxx(Lm+1,j)=0.
            Huxx(Lm+1,j)=-Huxx(Lm,j)
#    else
            uxx(Lm+1,j)=uxx(Lm,j)
            Huxx(Lm+1,j)=Huxx(Lm,j)
#    endif
          enddo
#    endif          
        endif
#   endif
        do j=Jstr,Jend
          do i=IstrU-1,Iend
#   ifdef UV_HADV_C4
            UFx(i,j)=0.25*( u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
     &                         -0.125*(uxx(i,j)+uxx(i+1,j))
     &                  )*( Huon(i,j,k)+Huon(i+1,j,k)
     &                      -0.125*(Huxx(i,j)+Huxx(i+1,j)))
#   else
#     if defined DIAGNOSTICS_UV  || defined DIAGNOSTICS_EK || \
         defined DIAGNOSTICS_VRT || defined DIAGNOSTICS_PV
            if (nnew.ne.3) then
            UFx4(i,j)=0.25*( u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
     &                         -0.125*(uxx(i,j)+uxx(i+1,j))
     &                  )*( Huon(i,j,k)+Huon(i+1,j,k)
     &                      -0.125*(Huxx(i,j)+Huxx(i+1,j)))
            endif
#     endif
            cffX=u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
            if (cffX.gt.0.) then
              curvX=uxx(i,j)
            else
              curvX=uxx(i+1,j)
            endif
            UFx(i,j)=0.25*(cffX+gamma*curvX)*( Huon(i,j,k)+
     &            Huon(i+1,j,k)-0.125*(Huxx(i,j)+Huxx(i+1,j)))
#   endif
          enddo
        enddo
#   undef Huxx
#   undef uxx

#   define vee wrk1
#   define Hvee wrk2
      
#   ifdef NS_PERIODIC
#    define JV_EXT_RANGE JstrV-1,Jend+1
#   else
#    ifdef MPI
        if (SOUTH_INTER) then
          jmin=JstrV-1
        else
          jmin=max(JstrV-1,2)
        endif
        if (NORTH_INTER) then
          jmax=Jend+1
        else
          jmax=min(Jend+1,Mmmpi)
        endif
#     define JV_EXT_RANGE jmin,jmax
#    else
#     define JV_EXT_RANGE max(JstrV-1,2),min(Jend+1,Mm)
#    endif
#   endif
        do j=JV_EXT_RANGE
          do i=Istr,Iend
            vee(i,j)=(v(i,j-1,k,nrhs)-2.*v(i,j,k,nrhs)
     &               +v(i,j+1,k,nrhs))  SWITCH  vmask(i,j)
            Hvee(i,j)=(Hvom(i,j-1,k)-2.*Hvom(i,j,k)
     &                +Hvom(i,j+1,k))   SWITCH  vmask(i,j)
          enddo
        enddo
#   undef JV_EXT_RANGE
#   ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=Istr,Iend
#    if defined OBC_MODIF_UP1
            vee(i,1)=(v(i,1,k,nrhs)-v(i,2,k,nrhs))/gamma
            Hvee(i,1)=-Hvee(i,2)
#    elif defined OBC_MODIF_C2
            vee(i,1)=0.
            Hvee(i,1)=-Hvee(i,2)
#    else
            vee(i,1)=vee(i,2)
            Hvee(i,1)=Hvee(i,2)
#    endif
          enddo
        endif
        if (NORTHERN_EDGE) then
#    ifdef MPI
          do i=Istr,Iend
#    if defined OBC_MODIF_UP1
            vee(i,Mmmpi+1)=(v(i,Mmmpi+1,k,nrhs)-v(i,Mmmpi,k,nrhs))/gamma
            Hvee(i,Mmmpi+1)=-Hvee(i,Mmmpi)
#    elif defined OBC_MODIF_C2
            vee(i,Mmmpi+1)=0.
            Hvee(i,Mmmpi+1)=-Hvee(i,Mmmpi)
#    else
            vee(i,Mmmpi+1)=vee(i,Mmmpi)
            Hvee(i,Mmmpi+1)=Hvee(i,Mmmpi)
#    endif
          enddo
#    else        
          do i=Istr,Iend
#    if defined OBC_MODIF_UP1
            vee(i,Mm+1)=(v(i,Mm+1,k,nrhs)-v(i,Mm,k,nrhs))/gamma
            Hvee(i,Mm+1)=-Hvee(i,Mm)
#    elif defined OBC_MODIF_C2
            vee(i,Mm+1)=0.
            Hvee(i,Mm+1)=-Hvee(i,Mm)
#    else
            vee(i,Mm+1)=vee(i,Mm)
            Hvee(i,Mm+1)=Hvee(i,Mm)
#    endif
          enddo
#    endif          
        endif
#   endif
        do j=JstrV-1,Jend
          do i=Istr,Iend
#   ifdef UV_HADV_C4
            VFe(i,j)=0.25*( v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
     &                        -0.125*(vee(i,j)+vee(i,j+1))
     &                  )*( Hvom(i,j,k)+Hvom(i,j+1,k)
     &                      -0.125*(Hvee(i,j)+Hvee(i,j+1)))
#   else
#     if defined DIAGNOSTICS_UV  || defined DIAGNOSTICS_EK || \
         defined DIAGNOSTICS_VRT || defined DIAGNOSTICS_PV
            if (nnew.ne.3) then
            VFe4(i,j)=0.25*( v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
     &                        -0.125*(vee(i,j)+vee(i,j+1))
     &                  )*( Hvom(i,j,k)+Hvom(i,j+1,k)
     &                      -0.125*(Hvee(i,j)+Hvee(i,j+1)))
            endif
#     endif
            cffE=v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
            if (cffE.gt.0.) then
              curvE=vee(i,j)
            else
              curvE=vee(i,j+1)
            endif
            VFe(i,j)=0.25*(cffE+gamma*curvE)*( Hvom(i,j,k)+
     &             Hvom(i,j+1,k)-0.125*(Hvee(i,j)+Hvee(i,j+1)))
#   endif
          enddo
        enddo
#   undef Hvee
#   undef vee

#   define uee wrk1
#   define Hvxx wrk2
#   ifdef NS_PERIODIC
#    define JU_EXT_RANGE Jstr-1,Jend+1
#   else
#    ifdef MPI
        if (SOUTH_INTER) then
          jmin=Jstr-1
        else
          jmin=max(Jstr-1,1)
        endif
        if (NORTH_INTER) then
          jmax=Jend+1
        else
          jmax=min(Jend+1,Mmmpi)
        endif
#     define JU_EXT_RANGE jmin,jmax
#    else
#     define JU_EXT_RANGE max(Jstr-1,1),min(Jend+1,Mm)
#    endif
#   endif
        do j=JU_EXT_RANGE
          do i=IstrU,Iend
            uee(i,j)=(u(i,j-1,k,nrhs)-u(i,j,k,nrhs)) SWITCH pmask(i,j  )
     &              +(u(i,j+1,k,nrhs)-u(i,j,k,nrhs)) SWITCH pmask(i,j+1)
          enddo
        enddo
#   undef JU_EXT_RANGE
#   ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=IstrU,Iend
#    if defined OBC_MODIF_UP1
            uee(i,0)=(u(i,0,k,nrhs)-u(i,1,k,nrhs))/gamma
#    elif defined OBC_MODIF_C2
            uee(i,0)=0.
#    else
            uee(i,0)=uee(i,1)
#    endif
          enddo
        endif
        if (NORTHERN_EDGE) then
#    ifdef MPI
          do i=IstrU,Iend
#    if defined OBC_MODIF_UP1
            uee(i,Mmmpi+1)=(u(i,Mmmpi+1,k,nrhs)-u(i,Mmmpi,k,nrhs))/gamma
#    elif defined OBC_MODIF_C2
            uee(i,Mmmpi+1)=0.
#    else
            uee(i,Mmmpi+1)=uee(i,Mmmpi)
#    endif
          enddo
#    else        
          do i=IstrU,Iend
#    if defined OBC_MODIF_UP1
            uee(i,Mm+1)=(u(i,Mm+1,k,nrhs)-u(i,Mm,k,nrhs))/gamma
#    elif defined OBC_MODIF_C2
            uee(i,Mm+1)=0.
#    else
            uee(i,Mm+1)=uee(i,Mm)
#    endif
          enddo
#    endif          
        endif
#   endif
        do j=Jstr,Jend+1
          do i=IstrU-1,Iend
           Hvxx(i,j)=Hvom(i-1,j,k)-2.*Hvom(i,j,k)+Hvom(i+1,j,k)
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=IstrU,Iend
#   if defined UV_HADV_C4
            UFe(i,j)=0.25*( u(i,j,k,nrhs)+u(i,j-1,k,nrhs)
     &                        -0.125*(uee(i,j)+uee(i,j-1))
     &                  )*( Hvom(i,j,k)+Hvom(i-1,j,k)
     &                     -0.125*(Hvxx(i,j)+Hvxx(i-1,j)))
#   else
#     if defined DIAGNOSTICS_UV  || defined DIAGNOSTICS_EK  || \
         defined DIAGNOSTICS_VRT || defined DIAGNOSTICS_PV
            if (nnew.ne.3) then
            UFe4(i,j)=0.25*( u(i,j,k,nrhs)+u(i,j-1,k,nrhs)
     &                        -0.125*(uee(i,j)+uee(i,j-1))
     &                  )*( Hvom(i,j,k)+Hvom(i-1,j,k)
     &                     -0.125*(Hvxx(i,j)+Hvxx(i-1,j)))
            endif
#     endif
            cffX=u(i,j,k,nrhs)+u(i,j-1,k,nrhs)
            cffE=Hvom(i,j,k)+Hvom(i-1,j,k)
            if (cffE.gt.0.) then
              curvX=uee(i,j-1)
            else
              curvX=uee(i,j)
            endif
            UFe(i,j)=0.25*(cffX+gamma*curvX)*(cffE-0.125*(
     &                             Hvxx(i,j)+Hvxx(i-1,j) ))
#   endif
          enddo
        enddo
#   undef Hvxx
#   undef uee

#   define vxx wrk1
#   define Huee wrk2
#   ifdef EW_PERIODIC 
#    define IV_EXT_RANGE Istr-1,Iend+1
#   else
#    ifdef MPI
        if (WEST_INTER) then
          imin=Istr-1
        else
          imin=max(Istr-1,1)
        endif
        if (EAST_INTER) then
          imax=Iend+1
        else
          imax=min(Iend+1,Lmmpi)
        endif
#     define IV_EXT_RANGE imin,imax
#    else
#     define IV_EXT_RANGE max(Istr-1,1),min(Iend+1,Lm)
#    endif
#   endif
        do j=JstrV,Jend
          do i=IV_EXT_RANGE
            vxx(i,j)=(v(i-1,j,k,nrhs)-v(i,j,k,nrhs)) SWITCH pmask(i  ,j)
     &              +(v(i+1,j,k,nrhs)-v(i,j,k,nrhs)) SWITCH pmask(i+1,j)
          enddo
        enddo
#   undef IV_EXT_RANGE
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=JstrV,Jend
#    if defined OBC_MODIF_UP1
            vxx(0,j)=(v(0,j,k,nrhs)-v(1,j,k,nrhs))/gamma
#    elif defined OBC_MODIF_C2
            vxx(0,j)=0.
#    else
            vxx(0,j)=vxx(1,j)
#    endif
          enddo
        endif
        if (EASTERN_EDGE) then
#    ifdef MPI
          do j=JstrV,Jend
#    if defined OBC_MODIF_UP1
            vxx(Lmmpi+1,j)=(v(Lmmpi+1,j,k,nrhs)-v(Lmmpi,j,k,nrhs))/gamma
#    elif defined OBC_MODIF_C2
            vxx(Lmmpi+1,j)=0.
#    else
            vxx(Lmmpi+1,j)=vxx(Lmmpi,j)
#    endif
          enddo
#    else        
          do j=JstrV,Jend
#    if defined OBC_MODIF_UP1
            vxx(Lm+1,j)=(v(Lm+1,j,k,nrhs)-v(Lm,j,k,nrhs))/gamma
#    elif defined OBC_MODIF_C2
            vxx(Lm+1,j)=0.
#    else
            vxx(Lm+1,j)=vxx(Lm,j)
#    endif
          enddo
#    endif          
        endif
#   endif
        do j=JstrV-1,Jend
          do i=Istr,Iend+1
           Huee(i,j)=Huon(i,j-1,k)-2.*Huon(i,j,k)+Huon(i,j+1,k)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend+1
#   if defined  UV_HADV_C4
            VFx(i,j)=0.25*( v(i,j,k,nrhs)+v(i-1,j,k,nrhs)
     &                        -0.125*(vxx(i,j)+vxx(i-1,j))
     &                  )*( Huon(i,j,k)+Huon(i,j-1,k)
     &                     -0.125*(Huee(i,j)+Huee(i,j-1)))
#   else
#     if defined DIAGNOSTICS_UV  || defined DIAGNOSTICS_EK  || \
         defined DIAGNOSTICS_VRT || defined DIAGNOSTICS_PV
            if (nnew.ne.3) then
            VFx4(i,j)=0.25*( v(i,j,k,nrhs)+v(i-1,j,k,nrhs)
     &                        -0.125*(vxx(i,j)+vxx(i-1,j))
     &                  )*( Huon(i,j,k)+Huon(i,j-1,k)
     &                     -0.125*(Huee(i,j)+Huee(i,j-1)))
            endif
#     endif
            cffE=v(i,j,k,nrhs)+v(i-1,j,k,nrhs)
            cffX=Huon(i,j,k)+Huon(i,j-1,k)
            if (cffX.gt.0.) then
              curvE=vxx(i-1,j)
            else
              curvE=vxx(i,j)
            endif
            VFx(i,j)=0.25*(cffE+gamma*curvE)*(cffX-0.125*(
     &                             Huee(i,j)+Huee(i,j-1) ))
#   endif
          enddo
        enddo
#   undef Huee
#   undef vxx

#  elif defined UV_HADV_C2 || defined UV_HADV_TVD
!
!----------------------------------------------------------
! Second order advection scheme with TVD correction
!----------------------------------------------------------
!
        do j=Jstr,Jend
          do i=IstrU-1,Iend
            UFx(i,j)=0.25*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
     &                      *(Huon(i,j,k)+Huon(i+1,j,k))
          enddo
        enddo
        do j=JstrV-1,Jend
          do i=Istr,Iend
            VFe(i,j)=0.25*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
     &                       *(Hvom(i,j,k)+Hvom(i,j+1,k))
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=IstrU,Iend
            UFe(i,j)=0.25*(u(i,j,k,nrhs)+u(i,j-1,k,nrhs))
     &                      *( Hvom(i,j,k)+Hvom(i-1,j,k))
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend+1
            VFx(i,j)=0.25*(v(i,j,k,nrhs)+v(i-1,j,k,nrhs))
     &                       *(Huon(i,j,k)+Huon(i,j-1,k))
          enddo
        enddo
!
! Correction with Total Variance Diminishing method
!
#  ifdef UV_HADV_TVD
#   include "uvhadv_tvd.h"
#  endif

#  endif /* UV_HADV_UP3 */
!
!----------------------------------------------------------
! Finalize horizontal advection: compute flux divergences
!----------------------------------------------------------
!
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru(i,j,k)=ru(i,j,k)-UFx(i,j  )+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i  ,j)
#  if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
            if (nnew.ne.3) then
              MXadv(i,j,k,1) = MXadv(i,j,k,1)
     &                       + (-UFx(i,j)+UFx(i-1,j))
     &                         SWITCH umask(i,j)
              MYadv(i,j,k,1) = (-UFe(i,j+1)+UFe(i,j))
     &                         SWITCH umask(i,j)
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
              MHdiff(i,j,k,1) = (UFx4(i-1,j)-UFx4(i,j)
     &                        +  UFe4(i,j)-UFe4(i,j+1))
     &                         SWITCH umask(i,j)
#    endif
            endif
#  elif defined DIAGNOSTICS_EK
          if (nnew.ne.3) then
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
            if (k.eq.1) then
             ekwrkHdiff(i,j,1) = 0
            endif
#    endif
             ekwrkHadv(i,j,1)=ekwrkHadv(i,j,1)
     &                       +(UFx(i-1,j)-UFx(i,j)
     &                       +UFe(i,j)-UFe(i,j+1))
     &                       * u(i,j,k,nrhs) SWITCH umask(i,j)
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
             ekwrkHdiff(i,j,1) = ekwrkHdiff(i,j,1)
     &                         + (UFx4(i-1,j)-UFx4(i,j)
     &                         +  UFe4(i,j)-UFe4(i,j+1))
     &                         *  u(i,j,k,nrhs) SWITCH umask(i,j)
#    endif
#    if defined DIAGNOSTICS_EK_MLD
              if (k.eq.kbl(i,j)) then
             ekwrkHadv_mld(i,j,1)=
     &                       +(UFx(i-1,j)-UFx(i,j)
     &                       +UFe(i,j)-UFe(i,j+1))
     &                       * u(i,j,k,nrhs) SWITCH umask(i,j)
#      if defined UV_HADV_UP3 || defined UV_HADV_UP5
             ekwrkHdiff_mld(i,j,1) =
     &                         + (UFx4(i-1,j)-UFx4(i,j)
     &                         +  UFe4(i,j)-UFe4(i,j+1))
     &                         *  u(i,j,k,nrhs) SWITCH umask(i,j)
#      endif
              elseif (k.gt.kbl(i,j)) then
             ekwrkHadv_mld(i,j,1)=ekwrkHadv_mld(i,j,1)
     &                       +(UFx(i-1,j)-UFx(i,j)
     &                       +UFe(i,j)-UFe(i,j+1))
     &                       * u(i,j,k,nrhs) SWITCH umask(i,j)
#      if defined UV_HADV_UP3 || defined UV_HADV_UP5
             ekwrkHdiff_mld(i,j,1) = ekwrkHdiff_mld(i,j,1)
     &                         + (UFx4(i-1,j)-UFx4(i,j)
     &                         +  UFe4(i,j)-UFe4(i,j+1))
     &                         *  u(i,j,k,nrhs) SWITCH umask(i,j)
#      endif
              endif
#    endif
           endif
#  endif
#  if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
            if (nnew.ne.3) then
             if (k.eq.1) then
              wrkYadv(i,j,1) = (-UFe(i,j+1)+UFe(i,j))
     &                         SWITCH umask(i,j)
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
              wrkHdiff(i,j,1) = (UFx4(i-1,j)-UFx4(i,j)
     &                          +  UFe4(i,j)-UFe4(i,j+1))
     &                         SWITCH umask(i,j)
#    endif
             else
              wrkYadv(i,j,1) = wrkYadv(i,j,1)
     &                       + (-UFe(i,j+1)+UFe(i,j))
     &                         SWITCH umask(i,j)
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
              wrkHdiff(i,j,1) = wrkHdiff(i,j,1)
     &                        + (UFx4(i-1,j)-UFx4(i,j)
     &                        +  UFe4(i,j)-UFe4(i,j+1))
     &                         SWITCH umask(i,j)
#    endif
             endif
              wrkXadv(i,j,1) = wrkXadv(i,j,1)
     &                       + (-UFx(i,j)+UFx(i-1,j))
     &                         SWITCH umask(i,j)
            endif
#  endif
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rv(i,j,k)=rv(i,j,k)-VFx(i+1,j)+VFx(i,j  )
     &                         -VFe(i  ,j)+VFe(i,j-1)
#  if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
            if (nnew.ne.3) then
              MXadv(i,j,k,2) = MXadv(i,j,k,2)
     &                         + (-VFx(i+1,j)+VFx(i,j))
     &                        SWITCH vmask(i,j)
              MYadv(i,j,k,2) = (-VFe(i,j)+VFe(i,j-1))
     &                         SWITCH vmask(i,j)
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
              MHdiff(i,j,k,2) = (VFx4(i,j)-VFx4(i+1,j)
     &                        +  VFe4(i,j-1)-VFe4(i,j))
     &                         SWITCH vmask(i,j)
#    endif
            endif
#  elif defined DIAGNOSTICS_EK
            if (nnew.ne.3) then
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
              if (k.eq.1) then
               ekwrkHdiff(i,j,2) = 0
              endif
#    endif
                ekwrkHadv(i,j,2) = ekwrkHadv(i,j,2)
     &                           +(VFx(i,j)-VFx(i+1,j)
     &                           +  VFe(i,j-1)-VFe(i,j))
     &                           *  v(i,j,k,nrhs) SWITCH vmask(i,j)
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
                ekwrkHdiff(i,j,2) = ekwrkHdiff(i,j,2) +
     &                           (VFx4(i,j)-VFx4(i+1,j)
     &                           + VFe4(i,j-1)-VFe4(i,j))
     &                           * v(i,j,k,nrhs) SWITCH vmask(i,j)
#    endif
#    if defined DIAGNOSTICS_EK_MLD
              if (k.eq.kbl(i,j)) then
                ekwrkHadv_mld(i,j,2) =
     &                           +(VFx(i,j)-VFx(i+1,j)
     &                           +  VFe(i,j-1)-VFe(i,j))
     &                           *  v(i,j,k,nrhs) SWITCH vmask(i,j)
#      if defined UV_HADV_UP3 || defined UV_HADV_UP5
                ekwrkHdiff_mld(i,j,2) =
     &                           (VFx4(i,j)-VFx4(i+1,j)
     &                           + VFe4(i,j-1)-VFe4(i,j))
     &                           * v(i,j,k,nrhs) SWITCH vmask(i,j)
#      endif
              elseif (k.gt.kbl(i,j)) then
                ekwrkHadv_mld(i,j,2) = ekwrkHadv_mld(i,j,2)
     &                           +(VFx(i,j)-VFx(i+1,j)
     &                           +  VFe(i,j-1)-VFe(i,j))
     &                           *  v(i,j,k,nrhs) SWITCH vmask(i,j)
#      if defined UV_HADV_UP3 || defined UV_HADV_UP5
                ekwrkHdiff_mld(i,j,2) = ekwrkHdiff_mld(i,j,2) +
     &                           (VFx4(i,j)-VFx4(i+1,j)
     &                           + VFe4(i,j-1)-VFe4(i,j))
     &                           * v(i,j,k,nrhs) SWITCH vmask(i,j)
#      endif
              endif
#    endif
            endif
#  endif
#  if defined DIAGNOSTICS_VRT && ! defined DIAGNOSTICS_UV
            if (nnew.ne.3) then
              if (k.eq.1) then
                wrkYadv(i,j,2) = (-VFe(i,j)+VFe(i,j-1))
     &                         SWITCH vmask(i,j)
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
                wrkHdiff(i,j,2) = (VFx4(i,j)-VFx4(i+1,j)
     &                        +  VFe4(i,j-1)-VFe4(i,j))
     &                         SWITCH vmask(i,j)
#    endif
              else
                wrkYadv(i,j,2) = wrkYadv(i,j,2)
     &                         + (-VFe(i,j)+VFe(i,j-1))
     &                        SWITCH vmask(i,j)
#    if defined UV_HADV_UP3 || defined UV_HADV_UP5
                wrkHdiff(i,j,2) = wrkHdiff(i,j,2)
     &                        +  (VFx4(i,j)-VFx4(i+1,j)
     &                        +  VFe4(i,j-1)-VFe4(i,j))
     &                         SWITCH vmask(i,j)
#    endif
              endif
                wrkXadv(i,j,2) = wrkXadv(i,j,2)
     &                         + (-VFx(i+1,j)+VFx(i,j))
     &                        SWITCH vmask(i,j)
            endif
#  endif
          enddo
        enddo
# endif /* UV_ADV */
      enddo    ! <-- k

      do j=Jstr,Jend

# if defined UV_COR & defined NHMG
!     Vertical component of Coriolis on vertical velocity
!     We are re-using FC as a work array
         do k = 1,N
            do i=Istr,Iend
              FC(i,k) = 0.5*fomn(i,j)*Hz(i,j,k)*( 
     &                 dzdeta(i,j,k)*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
     &                - dzdxi(i,j,k)*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs)) )
            enddo
         enddo

         do k = 1,N-1
            do i=Istr,Iend
              rw(i,j,k)= rw(i,j,k) + 0.5*(FC(i,k+1)+FC(i,k))
            enddo
         enddo
         do i=Istr,Iend
            rw(i,j,N) = rw(i,j,N) + FC(i,N)
         enddo
#  ifdef NONTRAD_COR
         ! FC is the vertical momentum at W-points
         do k = 1,N
            do i=Istr,Iend
               FC(i,k) = 0.5*Hz(i,j,k)*(
     &  +(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
     &    * ( fetaomn(i,j) * (1+dzdxi(i,j,k)**2)
     &      - fxiomn(i,j)*dzdxi(i,j,k)*dzdeta(i,j,k) )
     &  -(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
     &    * ( fxiomn(i,j) * (1+dzdeta(i,j,k)**2)
     &       +fetaomn(i,j)*dzdxi(i,j,k)*dzdeta(i,j,k) ))
            enddo
         enddo
         do k = 1,N-1
            do i=Istr,Iend
               rw(i,j,k) = rw(i,j,k) + 0.5*(FC(i,k)+FC(i,k+1))
     & +wz(i,j,k,nrhs)*0.5*Hz(i,j,k)*(
     &   -fxiomn(i,j)*(dzdeta(i,j,k)+dzdeta(i,j,k+1))
     &   +fetaomn(i,j)*( dzdxi(i,j,k)+ dzdxi(i,j,k+1)) )
            enddo
         enddo
         do i=Istr,Iend
            rw(i,j,N) = rw(i,j,N) + FC(i,N)
     & +wz(i,j,N,nrhs)*Hz(i,j,N)*(
     & -fxiomn(i,j)*dzdeta(i,j,N)
     & +fetaomn(i,j)* dzdxi(i,j,N) )
         enddo

#  endif

# endif

# ifdef UV_ADV
!
!======================================================================
!
!  Compute and add in vertical advection terms
!
!======================================================================
!
!----------------------------------------------------------
!  U-component VERTICAL ADVECTION:
!
!  Construct conservative parabolic splines for the 
!  vertical derivatives "CF" of u-momentum.
!----------------------------------------------------------
!
#  ifdef UV_VADV_SPLINES
        do k=1,N
          do i=IstrU,Iend
            DC(i,k)=0.5625*(HZR(i  ,j,k)+HZR(i-1,j,k))
     &             -0.0625*(HZR(i+1,j,k)+HZR(i-2,j,k))
          enddo
        enddo
        do i=IstrU,Iend     ! Construct parabolic splines: here
          FC(i,0)=0.        ! CF is the set off vertical derivatives
          CF(i,0)=0.        ! of the velocity field u(:,:,:,nrhs),
        enddo               ! FC is an auxiliary scratch variable.
        do k=1,N-1,+1
          do i=IstrU,Iend
            cff=1./(2.*DC(i,k+1)+DC(i,k)*(2.-FC(i,k-1)))
            FC(i,k)=cff*DC(i,k+1)
            CF(i,k)=cff*(6.*(u(i,j,k+1,nrhs)-u(i,j,k,nrhs))
     &                                  -DC(i,k)*CF(i,k-1))
          enddo
        enddo
        do i=IstrU,Iend
          CF(i,N)=0.
        enddo
        do k=N-1,1,-1
          do i=IstrU,Iend
            CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
          enddo
        enddo               !--> discard FC
!
! Compute vertical advective fluxes 
! FC=W*[spline-interpolated velocity]
!
        do k=1,N-1
          do i=IstrU,Iend
            FC(i,k)=( 0.5625*(We(i  ,j,k)+We(i-1,j,k))-0.0625*(
     &                                 We(i+1,j,k)+We(i-2,j,k) ))
     &             *( u(i,j,k,nrhs)+DC(i,k)*(
     &                              0.33333333333333*CF(i,k  )
     &                             +0.16666666666667*CF(i,k-1)
     &                                                       ))
          enddo
        enddo               !--> discard CF,DC
        do i=IstrU,Iend
          FC(i,N)=0.
#   ifdef MOVING_BATHY
          FC(i,0)=0.5*u(i,j,1,nrhs)*
     &            (We(i  ,j,0)+
     &             We(i-1,j,0))
#   else
          FC(i,0)=0.
#   endif
        enddo
!!!! 29/11/2019 @gyc: vadv
#  elif defined UV_VADV_WENO5 || defined UV_VADV_UP5 \
                              || defined UV_VADV_C6
!
!----------------------------------------------------------
! Compute vertical advective fluxes 
! using 5th-order WENO scheme
!----------------------------------------------------------
!
#   ifdef UV_VADV_WENO5
#    define FLUX5 flux5_weno
#    define FLUX3 flux3_weno
#    define FLUX2 flux1
          cdif=1.
#    include "u_vadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
#   else
        if (nnew.eq.3) then  ! predictor
#    define FLUX5 flux6
#    define FLUX3 flux4
#    define FLUX2 flux2
          cdif=0.
#    include "u_vadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
        else                 ! corrector
#    ifdef UV_VADV_C6
#     define FLUX5 flux6
#     define FLUX3 flux4
#     define FLUX2 flux2
          cdif=0.
#    else
#     define FLUX5 flux5
#     define FLUX3 flux3
#     define FLUX2 flux1
          cdif=0.5
#    endif
#    include "u_vadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
         endif
#   endif /* WENO5 */
!!!! end of 29/11/2019 @gyc: vadv
#  elif defined UV_VADV_C2 || defined UV_VADV_TVD
!
!----------------------------------------------------------
!  Compute vertical advective fluxes using
!  second order advection scheme with TVD correction
!----------------------------------------------------------
!
        do k=1,N-1
          do i=IstrU,Iend
            FC(i,k)=0.25*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))
     &                            *(We(i,j,k)+We(i-1,j,k))
          enddo
        enddo
        do i=IstrU,Iend
#   ifdef MOVING_BATHY
          FC(i,0)=0.5*u(i,j,1,nrhs)*
     &            (We(i  ,j,0)+
     &             We(i-1,j,0))
#   else
          FC(i,0)=0.
#   endif
          FC(i,N)=0.
        enddo
!
!  Correction with TVD method
!
#   ifdef UV_VADV_TVD
        do k=1,N-1
          kp2=MIN(N,k+2)
          kp1=MIN(N,k+1)
          km1=MAX(1,k-1)
          do i=IstrU,Iend
            cff = 0.5 *(We(i,j,k)+We(i-1,j,k) )
            cu = cff*pm_u(i,j)*pn_u(i,j)
     &           * cdt*2/(-z_r(i  ,j,k)+z_r(i  ,j,k+1)
     &                    -z_r(i-1,j,k)+z_r(i-1,j,k+1))      
            rrp= (u(i,j,kp2,nstp) -u(i,j,kp1,nstp)) 
            rr = (u(i,j,kp1,nstp) -u(i,j,k  ,nstp)) 
            rrm= (u(i,j,k  ,nstp) -u(i,j,km1,nstp))                 
            cff1=0.5*( cff*(u(i,j,k,nstp)+u(i,j,kp1,nstp))
     &                                      -abs(cff)*rr )
            Cr=limiter(cu,FC(i,k),cff1,rrm,rr,rrp)
            FC(i,k) = cff1*(1-cr) + Cr*FC(i,k)
          enddo
        enddo
#   endif

#  else
!
!----------------------------------------------------------
!  Compute vertical advective fluxes using
!  Fourth order advection scheme
!----------------------------------------------------------
!
        do k=2,N-2
          do i=IstrU,Iend
            FC(i,k)=( 0.5625*(u(i,j,k  ,nrhs)+u(i,j,k+1,nrhs))
     &                 -0.0625*(u(i,j,k-1,nrhs)+u(i,j,k+2,nrhs)))
     &                         *( 0.5625*(We(i  ,j,k)+We(i-1,j,k))
     &                           -0.0625*(We(i+1,j,k)+We(i-2,j,k)))
          enddo
        enddo
        do i=IstrU,Iend
          FC(i,N)=0.
          FC(i,N-1)=( 0.5625*(u(i,j,N-1,nrhs)+u(i,j,N,nrhs))
     &                   -0.0625*(u(i,j,N-2,nrhs)+u(i,j,N,nrhs)))
     &                     *( 0.5625*(We(i  ,j,N-1)+We(i-1,j,N-1))
     &                       -0.0625*(We(i+1,j,N-1)+We(i-2,j,N-1)))

          FC(i,  1)=( 0.5625*(u(i,j,  1,nrhs)+u(i,j,2,nrhs))
     &                   -0.0625*(u(i,j,  1,nrhs)+u(i,j,3,nrhs)))
     &                         *( 0.5625*(We(i  ,j,1)+We(i-1,j,1))
     &                           -0.0625*(We(i+1,j,1)+We(i-2,j,1)))
#   ifdef MOVING_BATHY
          FC(i,0)=2.*
     &     (0.5625* u(i,j,1,nrhs)-
     &     -0.0625* u(i,j,2,nrhs)
     &     )*
     &     (0.5625*(We(i  ,j,0)+
     &              We(i-1,j,0)
     &             )
     &      -0.0625*(We(i+1,j,0)+
     &                   We(i-2,j,0)
     &              )
     &      )
#   else
          FC(i,0)=0.
#   endif
        enddo
#  endif
!
!----------------------------------------------------------
! Finalize U vertical advection: compute flux divergences
!----------------------------------------------------------
!
#  if defined DIAGNOSTICS_EK && ! defined DIAGNOSTICS_EK_FULL && ! defined DIAGNOSTICS_UV
            if (nnew.ne.3) then
             do i=istrU,iend
              ekwrkVadv(i,j,1) = 0
#    if defined DIAGNOSTICS_EK_MLD
              ekwrkVadv_mld(i,j,1) = 0
#    endif
             enddo         
            endif
#  endif

        do k=1,N
          do i=IstrU,Iend
            ru(i,j,k)=ru(i,j,k)-FC(i,k)+FC(i,k-1)
#  if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL
            if (nnew.ne.3) then
              MVadv(i,j,k,1) = (-FC(i,k)+FC(i,k-1)) SWITCH umask(i,j)
            endif
#  else
#    ifdef DIAGNOSTICS_EK
            if (nnew.ne.3) then
              ekwrkVadv(i,j,1) = ekwrkVadv(i,j,1)
     &                         + (-FC(i,k)+FC(i,k-1))
     &                      * u(i,j,k,nrhs)  SWITCH umask(i,j)
#      if defined DIAGNOSTICS_EK_MLD
              if (k.ge.kbl(i,j)) then
                ekwrkVadv_mld(i,j,1) = ekwrkVadv_mld(i,j,1)
     &                         + (-FC(i,k)+FC(i,k-1))
     &                      * u(i,j,k,nrhs)  SWITCH umask(i,j)
              endif
#      endif
            endif
#    endif
#  endif
          enddo
        enddo               !--> discard FC

!
!----------------------------------------------------------
!  V-component VERTICAL ADVECTION:
!
!  Construct conservative parabolic splines for the 
!  vertical derivatives "CF" of u-momentum.
!----------------------------------------------------------
!
        if (j.ge.JstrV) then
!
#  ifdef UV_VADV_SPLINES
          do k=1,N
            do i=Istr,Iend
              DC(i,k)=0.5625*(HZR(i  ,j,k)+HZR(i,j-1,k))
     &               -0.0625*(HZR(i,j+1,k)+HZR(i,j-2,k))
            enddo             
          enddo             
          do i=Istr,Iend    ! Construct parabolic splines: here
            FC(i,0)=0.      ! CF is the set off vertical derivatives
            CF(i,0)=0.      ! of the velocity field v(:,:,:,nrhs),
          enddo             ! FC is an auxiliary scratch variable.
          do k=1,N-1,+1
            do i=Istr,Iend
              cff=1./(2.*DC(i,k+1)+DC(i,k)*(2.-FC(i,k-1)))
              FC(i,k)=cff*DC(i,k+1)
              CF(i,k)=cff*(6.*(v(i,j,k+1,nrhs)-v(i,j,k,nrhs))
     &                                    -DC(i,k)*CF(i,k-1))
            enddo
          enddo
          do i=Istr,Iend
            CF(i,N)=0.
          enddo
          do k=N-1,1,-1
            do i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            enddo
          enddo             !--> discard FC
!
! Compute vertical advective fluxes
! FC=W*[spline-interpolated velocity]
!
          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=( 0.5625*(We(i,j  ,k)+We(i,j-1,k))-0.0625*(
     &                                   We(i,j+1,k)+We(i,j-2,k) ))
     &               *( v(i,j,k,nrhs)+DC(i,k)*(
     &                                0.33333333333333*CF(i,k  )
     &                               +0.16666666666667*CF(i,k-1)
     &                                                         ))
            enddo
          enddo             !--> discard CF,DC
          do i=Istr,Iend
            FC(i,N)=0.
#   ifdef MOVING_BATHY
            FC(i,0)=0.5*v(i,j,1,nrhs)*
     &              (We(i,j  ,0)+
     &               We(i,j-1,0))
#   else
            FC(i,0)=0.
#   endif
          enddo

#  elif defined UV_VADV_WENO5 || defined UV_VADV_C6 \
                              || defined UV_VADV_UP5
!
!----------------------------------------------------------
! Compute vertical advective fluxes 
! using 5th-order WENO scheme
!----------------------------------------------------------
!
#   ifdef UV_VADV_WENO5
#    define FLUX5 flux5_weno
#    define FLUX3 flux3_weno
#    define FLUX2 flux1
          cdif=1.
#    include "v_vadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
#   else
        if (nnew.eq.3) then  ! predictor
#    define FLUX5 flux6
#    define FLUX3 flux4
#    define FLUX2 flux2
          cdif=0.
#    include "v_vadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
        else                 ! corrector
#    ifdef UV_VADV_C6
#     define FLUX5 flux6
#     define FLUX3 flux4
#     define FLUX2 flux2
          cdif=0.
#    else
#     define FLUX5 flux5
#     define FLUX3 flux3
#     define FLUX2 flux1
          cdif=0.5
#    endif
#    include "v_vadv_order5.h"
#    undef FLUX5
#    undef FLUX3
#    undef FLUX2
        endif
#   endif /* WENO5 */

#  elif defined UV_VADV_C2 || defined UV_VADV_TVD
!
!----------------------------------------------------------
!  Second order advection scheme with TVD correction
!----------------------------------------------------------
!
          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=0.25*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))
     &                              *(We(i,j,k)+We(i,j-1,k))
            enddo
          enddo
          do i=Istr,Iend
#   ifdef MOVING_BATHY
            FC(i,0)=0.5*v(i,j,1,nrhs)*
     &              (We(i,j  ,0)+
     &               We(i,j-1,0))
#   else
            FC(i,0)=0.
#   endif
            FC(i,N)=0.
          enddo
!
!  Correction using TVD method
!
#   ifdef UV_VADV_TVD
          do k=1,N-1
            kp2=min(N,k+2)
            kp1=min(N,k+1)
            km1=max(1,k-1)
            do i=Istr,Iend
              cff = 0.5*(We(i,j,k)+We(i,j-1,k))
              cu = cff*pm_v(i,j)*pn_v(i,j)
     &             * cdt*2/(-z_r(i,j  ,k)+z_r(i,j  ,k+1)
     &                      -z_r(i,j-1,k)+z_r(i,j-1,k+1))      
              rrp= v(i,j,kp2,nstp)-v(i,j,kp1,nstp) 
              rr = v(i,j,kp1,nstp)-v(i,j,k  ,nstp) 
              rrm= v(i,j,k  ,nstp)-v(i,j,km1,nstp)                    
              cff1=0.5*( cff*(v(i,j,k,nstp)+v(i,j,kp1,nstp))
     &                                        -abs(cff)*rr )
              Cr=limiter(cu,FC(i,k),cff1,rrm,rr,rrp)
              FC(i,k) = cff1*(1-cr) + Cr*FC(i,k)
            enddo
          enddo         
#   endif

#  else
!
!----------------------------------------------------------
!  Fourth order advection scheme
!----------------------------------------------------------
!
          do k=2,N-2
            do i=Istr,Iend
              FC(i,k)=( 0.5625*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))
     &                   -0.0625*(v(i,j,k-1,nrhs)+v(i,j,k+2,nrhs)))
     &                           *( 0.5625*(We(i,j  ,k)+We(i,j-1,k))
     &                             -0.0625*(We(i,j+1,k)+We(i,j-2,k)))
            enddo
          enddo
          do i=Istr,Iend
            FC(i,N)=0.
            FC(i,N-1)=( 0.5625*(v(i,j,N-1,nrhs)+v(i,j,N,nrhs))
     &                   -0.0625*(v(i,j,N-2,nrhs)+v(i,j,N,nrhs)))
     &                     *( 0.5625*(We(i,j  ,N-1)+We(i,j-1,N-1))
     &                       -0.0625*(We(i,j+1,N-1)+We(i,j-2,N-1)))

            FC(i,  1)=(+0.5625*(v(i,j,  1,nrhs)+v(i,j,2,nrhs))
     &                   -0.0625*(v(i,j,  1,nrhs)+v(i,j,3,nrhs)))
     &                         *( 0.5625*(We(i,j  ,1)+We(i,j-1,1))
     &                           -0.0625*(We(i,j+1,1)+We(i,j-2,1)))
#   ifdef MOVING_BATHY
          FC(i,0)=2.*
     &     (0.5625* v(i,j,1,nrhs)-
     &     -0.0625* v(i,j,2,nrhs)
     &     )*
     &     (0.5625*(We(i,j  ,0)+
     &              We(i,j-1,0)
     &             )
     &      -0.0625*(We(i,j+1,0)+
     &               We(i,j-2,0)
     &              )
     &      )
#   else
          FC(i,0)=0.
#   endif
           enddo
#  endif
!
!----------------------------------------------------------
! Finalize V vertical advection: compute flux divergences
!----------------------------------------------------------
!
#  if defined DIAGNOSTICS_EK && ! defined DIAGNOSTICS_EK_FULL  && ! defined DIAGNOSTICS_UV
            if (nnew.ne.3) then
             do i=Istr,Iend
              ekwrkVadv(i,j,2) = 0
#    if defined DIAGNOSTICS_EK_MLD
              ekwrkVadv_mld(i,j,2) = 0
#    endif
             enddo         
            endif
#  endif
          do k=1,N
            do i=Istr,Iend
              rv(i,j,k)=rv(i,j,k)-FC(i,k)+FC(i,k-1)
#  if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL
            if (nnew.ne.3) then
            MVadv(i,j,k,2) = (-FC(i,k)+FC(i,k-1)) SWITCH vmask(i,j)
            endif
#  else
#  ifdef DIAGNOSTICS_EK
            if (nnew.ne.3) then
             ekwrkVadv(i,j,2) = ekwrkVadv(i,j,2)
     &                         + (-FC(i,k)+FC(i,k-1))
     &                         * v(i,j,k,nrhs)
     &                           SWITCH vmask(i,j)
# if defined DIAGNOSTICS_EK_MLD
              if (k.ge.kbl(i,j)) then
                ekwrkVadv_mld(i,j,2) = ekwrkVadv_mld(i,j,2)
     &                         + (-FC(i,k)+FC(i,k-1))
     &                         * v(i,j,k,nrhs)
     &                           SWITCH vmask(i,j)
              endif
# endif
            endif
#  endif
#  endif
            enddo
          enddo
        endif
# endif /* UV_ADV */

!
!======================================================================
!
! Add combined vortex-force and advection terms & breaking terms
!
!======================================================================
!
# ifdef MRL_WCI
        do k=1,N
          do i=IstrU,iend
            cff  =0.5*on_u(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))
            vstu =0.25*( vst(i  ,j,k) +vst(i  ,j+1,k)
     &                  +vst(i-1,j,k) +vst(i-1,j+1,k) )
            dudx =0.5*( u(i+1,j,k,nrhs) -u(i-1,j,k,nrhs) )
            dvdx =0.5*( v(i  ,j,k,nrhs) -v(i-1,j,k,nrhs)
     &                 +v(i,j+1,k,nrhs) -v(i-1,j+1,k,nrhs) )
            ru(i,j,k) =ru(i,j,k) + cff
     &                           *(ust(i,j,k)*dudx+vstu*dvdx)
#  ifndef WAVE_SFC_BREAK
     &                           +cff*om_u(i,j)*brk3dx(i,j,k)
#  endif
#  if defined WAVE_STREAMING && defined WAVE_BODY_STREAMING
     &                           +cff*om_u(i,j)*frc3dx(i,j,k)
#  endif
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              Mvf(i,j,k,1)  = cff*(ust(i,j,k)*dudx +vstu*dvdx)
#   ifndef WAVE_SFC_BREAK
              Mbrk(i,j,k,1) = cff*om_u(i,j)*brk3dx(i,j,k)
#   endif
#   if defined WAVE_STREAMING && defined WAVE_BODY_STREAMING
              Mfrc(i,j,k,1) = cff*om_u(i,j)*frc3dx(i,j,k)
#   endif
            endif
#  endif
          enddo


          do i=istr,iend
            cff  =0.5*om_v(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
            ustv =0.25*( ust(i,j  ,k) +ust(i+1,j  ,k)
     &                  +ust(i,j-1,k) +ust(i+1,j-1,k) )
            dude =0.5*( u(i  ,j,k,nrhs) -u(i  ,j-1,k,nrhs)
     &                 +u(i+1,j,k,nrhs) -u(i+1,j-1,k,nrhs) )
ccc            dvde =0.5*( v(i,j-1,k,nrhs) -v(i,j+1,k,nrhs) )
            dvde =0.5*( v(i,j+1,k,nrhs) -v(i,j-1,k,nrhs) )
            rv(i,j,k) =rv(i,j,k) + cff
     &                            *(ustv*dude+vst(i,j,k)*dvde)
#  ifndef WAVE_SFC_BREAK
     &                            +cff*on_v(i,j)*brk3de(i,j,k)
#  endif
#  if defined WAVE_STREAMING && defined WAVE_BODY_STREAMING
     &                            +cff*on_v(i,j)*frc3de(i,j,k)
#  endif
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              Mvf(i,j,k,2)  = cff*(ustv*dude + vst(i,j,k)*dvde)
#   ifndef WAVE_SFC_BREAK
              Mbrk(i,j,k,2) = cff*on_v(i,j)*brk3de(i,j,k)
#   endif
#   if defined WAVE_STREAMING && defined WAVE_BODY_STREAMING
              Mfrc(i,j,k,2) = cff*on_v(i,j)*frc3de(i,j,k)
#   endif
            endif
#  endif
         enddo
      enddo

#  if defined WAVE_STREAMING && !defined WAVE_BODY_STREAMING
        do i=IstrU,iend
          ru(i,j,1) =ru(i,j,1) + om_u(i,j)*on_u(i,j)*frc2dx(i,j)
        enddo
        do i=istr,iend
          rv(i,j,1) =rv(i,j,1) + om_v(i,j)*on_v(i,j)*frc2de(i,j)
        enddo
#  endif

# endif   /* MRL_WCI */

!
!======================================================================
!
!  Start computation of the forcing terms for the 2D (barotropic
!  mode) momentum equations: vertically integrate the just computed
!  r.h.s "ru" and "rv". Also, if so prescribed, add in the
!  difference between surface (wind) and bottom (drag) stresses.
!
!  The computation of the 2D forcing terms will be finalized in
!  "rhs2d" during the first barotropic time step, when the
!  barotropically computed r.h.ss "rubar", "rvbar" will be subtracted
!  from the vertically integrated (here) "rufrc", "rvfrc".
!
!======================================================================
!

        do i=IstrU,Iend
          rufrc(i,j)=ru(i,j,1)
# ifndef BODYFORCE
     &        +(sustr(i,j)-bustr(i,j))*om_u(i,j)*on_u(i,j)
# endif
# if defined MRL_WCI && defined WAVE_SFC_BREAK
     &                    +brk2dx(i,j)*om_u(i,j)*on_u(i,j)
# endif
        enddo
        do k=2,N
          do i=IstrU,Iend
            rufrc(i,j)=rufrc(i,j)+ru(i,j,k)
          enddo
        enddo

        if (j.ge.JstrV) then
          do i=Istr,Iend
            rvfrc(i,j)=rv(i,j,1)
# ifndef BODYFORCE
     &       +(svstr(i,j)-bvstr(i,j))*om_v(i,j)*on_v(i,j)
# endif
# if defined MRL_WCI && defined WAVE_SFC_BREAK
     &                   +brk2de(i,j)*om_v(i,j)*on_v(i,j)
# endif
          enddo
          do k=2,N
            do i=Istr,Iend
              rvfrc(i,j)=rvfrc(i,j)+rv(i,j,k)
            enddo
          enddo
        endif

      enddo      !<-- j loop

      return
      end
!
!=======================================================================
!
# if defined UV_VADV_TVD || defined UV_HADV_TVD
 
      function limiter(nu,hf,lf,rrm,rr,rrp)      
      implicit none
      real    :: cff,limiter,Rj,nu,hf,lf,rrm,rr,rrp
!
!  Upwind        Limiter(Rj)=0.
!  Lax-Wendroff  Limiter(Rj)=1.
!  Min-Mod       Limiter(Rj)=max(0.,min(1.,Rj))
!  Suberbee      Limiter(Rj)=max(0.,max(min(1.,2*Rj),min(2.,Rj)))
!  Van Leer      Limiter(Rj)=(Rj+abs(Rj))/(1+abs(Rj))  => default
!
      if (abs(rr).gt.1.e-20) THEN
        if (nu .gt. 1.e-20) THEN
          Rj=rrm/rr
        else
          Rj=rrp/rr
        endif
      else
        if (nu.gt.1.e-20) then
          Rj=Rrm*1.E20
        else
          Rj=Rrp*1.E20
        endif
      endif

#  ifdef SUPERBEE
       limiter=max(0.,max(min(1.,2.*Rj),min(2.,Rj)))
#  elif defined MINMOD
       limiter=max(0.,min(1.,Rj))
#  else
       limiter=(Rj+abs(Rj))/(1+abs(Rj))  ! Van Leer
#  endif

      end function

# endif /* UV_VADV_TVD || UV_VADV_TVD */
!
!=======================================================================
!
#else
      subroutine rhs3d_empty
      end
#endif /* SOLVE3D */

