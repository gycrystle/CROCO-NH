! $Id: gls_mixing.F 1524 2014-04-14 17:00:06Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
#include "cppdefs.h"

#if defined SOLVE3D && defined GLS_MIX2017 
   !!======================================================================
   !!                 ***  SUBROUTINE  gls_mixing  ***
   !!
   !!   Vertical mixing coefficients computed from the gls turbulent  
   !!  closure parameterization (3 schemes: k-epsilon, k-omega, k-gen)
   !!
   !! *** References :
   !!  
   !! Umlauf, L., and H. Burchard, A generic length-scale equation 
   !! for geophysical turbulence models, J. Mar. Res., 61, 235-265, 2003
   !!
   !! Umlauf, L., and H. Burchard, Second-order turbulence closure models 
   !!                        Cont. Shelf. Res., 25, 795-827, 2005
   !!
   !! Duhaut, T., Notes sur les schemas de turbulence a deux equations, 
   !!                     Technical Report. IFREMER report , 2009   
   !!
   !! *** History :
   !!      !  2016-11 (F. Lemarié) : original code inspired by the Mars3D 
   !!                                implementation by S. Petton, A.-C. Bennis,  
   !!                                T. Duhaut, V. Garnier, F. Dumas 
   !!
   !!======================================================================

# define PSI_SMOOTH
# define TKE_SMOOTH
c++ # define GLS_HADV_UP1
c++ # undef GLS_KOMEGA
c++ # define GLS_KEPSILON
c++ # define CANUTO_A

   
#if   defined GLS_KOMEGA  
#elif defined GLS_KEPSILON
#else
#undef GLS_KOMEGA  
#undef GLS_KEPSILON
#endif

# if   defined CANUTO_A  
# elif defined GibLau_78
# elif defined MelYam_82
# elif defined KanCla_94
# elif defined Luyten_96
# elif defined CANUTO_B 
# elif defined Cheng_02
# else
# define CANUTO_A
#endif


!=======================================================================
      SUBROUTINE gls_mix2017 (tile)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER         :: tile, trd
      INTEGER         :: omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "ocean3d.h"
# include "compute_tile_bounds.h"
      
      trd=omp_get_thread_num()
      
      call gls_mix2017_tile ( Istr, Iend, Jstr, Jend,
     &                      A3d(1, 1,trd), A2d(1,  2,trd), 
     &                      A2d(1, 3,trd), A2d(1, 4,trd), A2d(1,  5,trd)
     &                    , A2d(1, 6,trd), A2d(1, 7,trd), A2d(1,  8,trd) )
      
      RETURN
!-----------------------------------------------------------------------
      END
!=======================================================================






!=======================================================================
       SUBROUTINE gls_mix2017_tile ( Istr, Iend, Jstr, Jend, 
     &                            shear2    ,   diss, ustar_sfc_sq, 
     &                      ustar_bot_sq, DC,     FC,      CF,  RH  ) 
!-----------------------------------------------------------------------                                                                
      IMPLICIT NONE
# include "param.h"
      ! Local integers
      INTEGER         ::   Istr, Iend, Jstr, Jend
      INTEGER         ::   i,       j,    k, tind
      INTEGER         ::   imin, imax, jmin, jmax
      INTEGER         ::   ig,    ig1,  ig2
      ! Local arrays
      REAL            ::  shear2      (PRIVATE_2D_SCRATCH_ARRAY,0:N-1)
      REAL            ::  diss        (PRIVATE_1D_SCRATCH_ARRAY,1:N-1)
      REAL            ::  ustar_sfc_sq(PRIVATE_2D_SCRATCH_ARRAY      )
      REAL            ::  ustar_bot_sq(PRIVATE_2D_SCRATCH_ARRAY      )
      REAL            ::  DC          (PRIVATE_1D_SCRATCH_ARRAY,0:N  )     
      REAL            ::  FC          (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL            ::  CF          (PRIVATE_1D_SCRATCH_ARRAY,1:N-1)
      REAL            ::  RH          (PRIVATE_1D_SCRATCH_ARRAY,1:N-1)     
      ! Local scalars
      REAL            ::  cff        , cff1       , cff2  , cff3m    , cff3p
      REAL            ::  invk       , invG       , Bprod , Sprod    , epsilon
      REAL            ::  alpha_n    , alpha_m    , c_mu  , c_mu_prim, Denom
      REAL            ::  alpha_n_min, alpha_m_max, cm0   , cm0inv2  , gls
      REAL            ::  flux_top   , flux_bot   , lgthsc, L_lim, dU,  dV    
      REAL            ::  trb_sfc    , trb_bot    ,   z0_s, z0_b     , gls_min  
      REAL            ::  HUon_w     , HVom_w     ,  trb_min(2)         
      ! Parameter values      
      REAL, PARAMETER ::  eps_min =  1.0E-12             ! min value for dissipation rate
      REAL, PARAMETER ::  tke_min =  1.0E-06             ! min value for TKE    
      REAL, PARAMETER ::  eps     =  1.0E-14             ! min stable stratification        
      REAL, PARAMETER ::  galp    =  0.53                ! parameter for the Galperin mixing length limitation
      REAL, PARAMETER ::  nuws    =  0.1E-04             ! background value for diffusivities 
      REAL, PARAMETER ::  nuwm    =  1.0E-04             ! background value for viscosities
# include "scalars.h"
      REAL, PARAMETER ::  chk     =  1400./g             ! charnock coefficient   
      !REAL, PARAMETER ::  chk     =  2.E+05/g           
      ! Choice of GLS model (!!! warning, it requires rn<0 !!!)
      REAL            :: rp,    rm,    rn                !<-- n,m and p exponents
      REAL            :: beta1, beta2, beta3m, beta3p    !<-- beta terms for the psi equation
      REAL            :: OneOverSig(2)                !<-- inverse of Schmidt number for tke and psi
# if defined GLS_KOMEGA                            /* K-omega model */
      PARAMETER( rp    = -1.0 , rm    = 0.5  , rn     = -1.0               )
      PARAMETER( beta1 = 0.555, beta2 = 0.833, beta3m = -0.6, beta3p = 1.0 )
      PARAMETER( OneOverSig = (/ 0.5, 0.5 /) )       
# elif defined GLS_KEPSILON                      /* K-epsilon model */
      PARAMETER( rp    = 3.0 , rm    = 1.5 , rn     = -1.0               )
      PARAMETER( beta1 = 1.44, beta2 = 1.92, beta3m = -0.4, beta3p = 1.0 )
      PARAMETER( OneOverSig = (/ 1.0, 0.7692 /) )  
# else                                                 /* GEN model */
      PARAMETER( rp    = 0.0, rm    = 1.0 , rn     = -0.67               )     
      PARAMETER( beta1 = 1.0, beta2 = 1.22, beta3m =  0.05, beta3p = 1.0 )
      PARAMETER( OneOverSig = (/ 1.25, 0.9345 /) ) 
# endif
# ifdef MRL_WCI
      REAL :: cmu_fac, cb_wallE, gls_sigp_cb
      REAL, PARAMETER :: zalp = 0.05     ! fraction of wave energy
# endif

      REAL, PARAMETER :: e1 =  3.0 + 1.*rp / rn
      REAL, PARAMETER :: e2 =  1.5 + 1.*rm / rn
      REAL, PARAMETER :: e3 = -1.0 / rn
# if defined PSI_SMOOTH || defined TKE_SMOOTH
      ! Hanning filter      
      !REAL, PARAMETER :: smth_a = 1.D0/8.D0
      !REAL, PARAMETER :: smth_b = 1.D0/4.D0     
      ! 9-point isotropic laplacian filter
      REAL, PARAMETER :: smth_a = 1.D0/12.D0
      REAL, PARAMETER :: smth_b = 3.D0/16.D0 
# endif      
      ! Choice of Stability function (Default : Canuto A)
      !--------------------------------------------------
# include "gls_stab_func.h"     
      !--------------------------------------------------
# include "ocean3d.h"
# include "forces.h"
# include "grid.h"
# include "mixing.h"
      !--------------------------------------------------
      ! Compute extended bounds (necessary because of 
      !                          the TKE/GLS smoothing)
# ifdef EW_PERIODIC
#  define I_EXT_RANGE Istr-1,Iend+1
      imin = Istr-1
      imax = Iend+1
# else
      if (WESTERN_EDGE) then
        imin=Istr
      else
        imin=Istr-1
      endif
      if (EASTERN_EDGE) then
        imax=Iend
      else
        imax=Iend+1
      endif
#  define I_EXT_RANGE imin,imax
# endif
# ifdef NS_PERIODIC
#  define J_EXT_RANGE Jstr-1,Jend+1
      jmin = Jstr-1
      jmax = Jend+1
# else
      if (SOUTHERN_EDGE) then
        jmin=Jstr
      else
        jmin=Jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend
      else
        jmax=Jend+1
      endif
#  define J_EXT_RANGE jmin,jmax
# endif    
      !--------------------------------------------------
      ! Initialization of various constants
      cm0     =  ( (a2*a2 - 3.0*a3*a3 + 3.0*a1*nn)/(3.0*nn*nn) )**0.25  ! Compute cmu0 
      cm0inv2 = 1./cm0**2                                               ! inverse of cmu0 squared   
      ! minimum value of alpha_n to ensure that alpha_m is positive
      alpha_n_min = 0.5*( - ( sf_d1 + sf_nb0 )  + sqrt(  ( sf_d1 + sf_nb0 )**2    
     &            - 4. * sf_d0 *( sf_d4 + sf_nb1 ) ) ) / ( sf_d4 + sf_nb1 )  
      cff     = (cm0**3 )*(tke_min**1.5) / eps_min                      ! Compute gls_min consistently
      gls_min = (cm0**rp)*(tke_min**rm ) * ( cff**rn )                  !  with eps_min/tke_min         
      
      trb_min(itke) = tke_min
      trb_min(igls) = gls_min

# ifdef MRL_WCI
      cmu_fac=(1.5/OneOverSig(1))**(1./3.)
     &             /(cm0**(4./3.))
      cb_wallE=1.
      gls_sigp_cb=vonKar**2/(cm0**2*beta2*cb_wallE)
     &            *(rn**2-cff*rn/3.*(4.*rm+1.)
     &                      +cff**2*rm/9.*(2.+4.*rm))
# endif

      !--------------------------------------------------
      IF (FIRST_TIME_STEP) THEN
         DO k=0,N
            DO j=jmin,jmax
               DO i=imin,imax    
                  trb( i, j, k, nstp, itke ) = trb_min(itke)
                  trb( i, j, k, nstp, igls ) = trb_min(igls)             
               ENDDO
            ENDDO  
         ENDDO
      ENDIF 


      !--------------------------------------------------
      ! Horizontal advection (first-order upwind, not yet finalized)
#   define FX ustar_sfc_sq      
#   define FE ustar_bot_sq      
      DO ig = 1,ngls      
         DO k=1,N-1
!# ifdef GLS_HADV_UP1
!#include "gls_compute_horiz_fluxes.h"
!            DO j=jmin,jmax
!               DO i=imin,imax          
!                  cff    = 2. / ( Hz(i,j,k)+ Hz(i,j,k+1) )
!                  trb(i,j,k,nnew,ig)= trb(i,j,k,nstp,ig) 
!     &                           -dt*cff*pm(i,j)*pn(i,j)
!     &                             *(  FX(i+1,j)-FX(i,j)
!     &                               + FE(i,j+1)-FE(i,j) ) 
!               ENDDO
!            ENDDO
!# else
            DO j=jmin,jmax
               DO i=imin,imax  
                  trb(i,j,k,nnew,ig)=trb(i,j,k,nstp,ig)
               ENDDO
            ENDDO   
!# endif            
         ENDDO
      ENDDO         
# undef FX
# undef FE
      !--------------------------------------------------
      ! Compute the vertical shear
      tind  = nrhs
      DO k=1,N-1
         DO j=jmin,jmax
            DO i=imin,imax 
               cff = 1. / ( Hz( i, j, k ) + Hz( i, j, k+1 ) )
               du  = cff*( u(i, j, k+1, tind)+u(i+1, j, k+1, tind)
     &                    -u(i, j, k  , tind)-u(i+1, j, k  , tind))
               dv  = cff*( v(i, j, k+1, tind)+v(i, j+1, k+1, tind)
     &                    -v(i, j, k  , tind)-v(i, j+1, k  , tind)) 
               shear2(i,j,k) = du*du + dv*dv
            ENDDO
         ENDDO
      ENDDO

      !--------------------------------------------------   
      ! Compute ustar squared at the surface and at the bottom
      DO j=jmin,jmax
         DO i=imin,imax 
            ustar_sfc_sq( i, j ) = sqrt( (0.5*(sustr(i,j)+sustr(i+1,j)))**2
     &                                  +(0.5*(svstr(i,j)+svstr(i,j+1)))**2 )
            ustar_bot_sq( i, j ) = sqrt( (0.5*(bustr(i,j)+bustr(i+1,j)))**2
     &                                  +(0.5*(bvstr(i,j)+bvstr(i,j+1)))**2 )     
         ENDDO
      ENDDO 

      !--------------------------------------------------
      DO j=jmin,jmax         !<-- j-outer loop
      !--------------------------------------------------
      
         ! Compute the dissipation rate   
         DO i=imin,imax        
            DO k=1,N-1
               cff       = (cm0**e1) * ( trb( i,j,k,nstp,itke )**e2 ) 
     &                               * ( trb( i,j,k,nstp,igls )**e3 )              
               diss(i,k) = MAX( cff , eps_min )  
            ENDDO  
         ENDDO !<-- terminate i-loop               
         
         !--------------------------------------------------
         DO ig = 1,ngls     ! ig = 2 for gls and = 1 for tke
         !--------------------------------------------------
         
            ! Off-diagonal terms for the tridiagonal problem 
            cff=-0.5*dt
            DO k=2,N-1
               DO i=imin,imax
                  FC(i,k)= cff* OneOverSig(ig)* 
     &                   ( Akv(i,j,k)+Akv(i,j,k-1) ) / Hz(i,j,k)
               ENDDO
            ENDDO

            DO i=imin,imax
               FC(i,1)=0.
               FC(i,N)=0.
            ENDDO   
           
            ! Production/Dissipation terms and diagonal term
            DO k=1,N-1
               DO i=imin,imax            
                  ig1   = (igls-ig); ig2 = (ig-itke)           
                  invk  =     1. / trb( i,j,k,nstp,itke )
                  gls   =          trb( i,j,k,nstp,igls )
                  ! invG = 1 for tke invg=1/psi for gls                     
                  invG  =  ig1+ig2*(1./gls)             
                  cff1  =  ig1+ig2*beta1   * invk*gls
                  cff2  = (ig1+ig2*beta2 ) * invk          
                  cff3m =  ig1+ig2*beta3m  * invk*gls
                  cff3p =  ig1+ig2*beta3p  * invk*gls
                  ! Shear and buoyancy production
                  Sprod =  cff1*Akv(i,j,k) * shear2(i,j,k)    
                  Bprod = -Akt(i,j,k,itemp)*( cff3m*MAX(bvf(i,j,k),0.)       
     &                                     +  cff3p*MIN(bvf(i,j,k),0.) )
                  ! Patankar trick to ensure non-negative solutions
                  cff   =       0.5*(Hz(i,j,k)+Hz(i,j,k+1))
                  IF( (Bprod + Sprod) .gt. 0.) THEN
                     RH(i,k) = cff*( trb(i,j,k,nnew,ig) + dt*(Bprod+Sprod) )       
                     DC(i,k) = cff*(1.+dt*cff2*diss(i,k))-FC(i,k)-FC(i,k+1)
                  ELSE              
                     RH(i,k) = cff*( trb(i,j,k,nnew,ig) + dt*       Sprod  )
                     DC(i,k) = cff*(1.+dt*(cff2*diss(i,k)  
     &                              -invG*Bprod)) - FC(i,k) - FC(i,k+1) 
                  ENDIF
               ENDDO
            ENDDO  
            
            ! Boundary conditions
            IF( ig == itke ) THEN
               DO i=imin,imax            
                  ! surface
# if defined MRL_WCI
                  trb_sfc      = MAX( tke_min, cmu_fac*
     &                       (zalp*wepb(i,j))**(2./3.) )
                  flux_top=zalp*wepb(i,j)
# else
                  trb_sfc      = MAX( tke_min, cm0inv2*ustar_sfc_sq(i,j) ) 
                  flux_top     = 0.
# endif
                  ! bottom                  
                  trb_bot      = MAX( tke_min, cm0inv2*ustar_bot_sq(i,j) )
                  flux_bot     = 0.              
                  ! finalize
                  RH(  i,1   ) = RH(i,  1) + dt*flux_bot
                  RH(  i,N-1 ) = RH(i,N-1) + dt*flux_top                  
                  trb( i,j,N,nnew,ig ) = trb_sfc
                  trb( i,j,0,nnew,ig ) = trb_bot 
               
               ENDDO
            ELSE
               DO i=imin,imax
                  ! surface
# if !defined MRL_WCI
                  z0_s    = MAX( 1.e-2 , chk*ustar_sfc_sq(i,j) )   !<-- Charnock
# else
                  z0_s=MAX(whrm(i,j)/2.,0.01)
# endif
!                  cff     = 30.*tanh( 0.6 / (28.*sqrt( ustar_sfc_sq(i,j) )) )
!                  z0_s    = MAX( 1.e-2 , 
!     &                       1.3*( 782.353/g )*ustar_sfc_sq(i,j)*(cff**1.5) )                 
                  cff         = 0.5*( trb( i,j,N-1,nnew,itke )
     &                            +   trb( i,j,N  ,nnew,itke ) )
                  lgthsc      = vonKar*(0.5*Hz(i,j,N)+z0_s)

                  trb_sfc     = MAX(gls_min,(cm0**rp)*(lgthsc**rn)
     &                                                 *(cff**rm))
                  flux_top    = -rn*cm0**(rp+1.)
     &                             *vonKar*OneOverSig(igls)
     &                             *(cff**(rm+0.5))*(lgthsc**rn)

# ifdef MRL_WCI
                  flux_top=flux_top-
     &             OneOverSig(itke)*OneOverSig(igls)
     &                 *   cm0**(rp+1.)
     &                 *   rm
     &                 *   cff**(rm-0.5)
     &                 *   lgthsc** (rn+1.)
     &                 *   zalp*   wepb(i,j)
# endif 
                  ! bottom
                  z0_b        = MAX( Zob , 1.E-04 )
                  cff         = 0.5*( trb( i,j,  1,nnew,itke )
     &                            +   trb( i,j,  0,nnew,itke ) )
                  lgthsc      = vonKar*(0.5*Hz(i,j,1)+z0_b)
                  trb_bot     = MAX(gls_min,(cm0**rp)*(lgthsc**rn)
     &                                               *(cff**rm)) 
                  flux_bot    =-rn*cm0**(rp+1.)
     &                          *vonKar*OneOverSig(igls)
     &                          *(cff**(rm+0.5))*(lgthsc**rn)
                  ! finalize
                  RH( i,  1 ) = RH(i,  1) + dt*flux_bot
                  RH( i,N-1 ) = RH(i,N-1) + dt*flux_top
                  trb( i,j,N,nnew,ig ) = trb_sfc
                  trb( i,j,0,nnew,ig ) = trb_bot 
               
               ENDDO
            ENDIF               
            
            ! tridiagonal resolution 
            DO i=imin,imax             
               cff       =  1./DC(i,N-1)                
               CF(i,N-1) = cff*FC(i,N-1)           
               RH(i,N-1) = cff*RH(i,N-1)
            ENDDO   
                  
            DO k=N-2,1,-1
               DO i=imin,imax      
                  cff     =   1./(DC(i,k)-CF(i,k+1)*FC(i,k+1))      
                  CF(i,k) = cff*FC(i,k)
                  RH(i,k) = cff*( RH(i,k)-FC(i,k+1)*RH(i,k+1))
               ENDDO
            ENDDO   

            DO i=imin,imax           
               trb( i,j,1,nnew,ig ) = MAX( RH(i,1), trb_min(ig) )
            ENDDO  
              
            DO k=2,N-1
               DO i=imin,imax
                  RH(i,k) = RH(i,k)-CF(i,k)*RH(i,k-1)
                  trb( i,j,k,nnew,ig ) = MAX( RH(i,k), trb_min(ig) )
               ENDDO     
            ENDDO  

         !--------------------------------------------------
         ENDDO     ! ig loop
         !--------------------------------------------------             
    
      !--------------------------------------------------
      ENDDO     ! j loop
      !-------------------------------------------------- 
      !
      !    GLS Smoothing
      !
# ifdef PSI_SMOOTH
      ig = igls
#  define FX ustar_sfc_sq      
#  define FE ustar_bot_sq 
#  define FE1 shear2       
      DO k=0,N
#        include "gls2017_smooth.h"
      ENDDO
#  undef FX
#  undef FE      
#  undef FE1
# endif

# ifdef TKE_SMOOTH
        ig = itke
#  define FX ustar_sfc_sq      
#  define FE ustar_bot_sq 
#  define FE1 shear2       
      DO k=0,N
#        include "gls2017_smooth.h"
      ENDDO
#  undef FX
#  undef FE      
#  undef FE1
# endif

# undef I_EXT_RANGE
# undef J_EXT_RANGE

      !--------------------------------------------------        
      DO j=jstr,jend  !<-- j-outer loop
      !--------------------------------------------------         
         DO k=1,N-1         
            DO i=istr,iend  
               !
               ! Galperin limitation : l <= l_lim
               L_lim = galp * sqrt( 2.* trb(i,j,k,nnew,itke)) / 
     &                            ( sqrt(max(eps, bvf(i,j,k)))  )
               !
               ! Limitation on psi (use MAX because rn is negative)
               cff = (cm0**rp) * (L_lim**rn) * (trb(i,j,k,nnew,itke)**rm)
               trb( i,j,k,nnew,igls ) = MAX( trb( i,j,k,nnew,igls ),cff )
               !
               ! Dissipation rate
               epsilon = (cm0**e1) * ( trb( i,j,k,nnew,itke )**e2 ) 
     &                             * ( trb( i,j,k,nnew,igls )**e3 )   
               epsilon = MAX(epsilon,eps_min)     
               ! 
               ! Compute alpha_n and alpha_m       
               cff     = ( trb(i,j,k,nnew,itke)/epsilon )**2                        
               alpha_m     = cff*  shear2(i,j,k)
               alpha_n     = cff*     bvf(i,j,k) 
               !
               ! Limitation of alpha_n and alpha_m
               alpha_n     = MIN(  MAX( 0.73*alpha_n_min , alpha_n ) , 1.0e10 )
               alpha_m_max = ( lim_am0 + lim_am1 * alpha_n 
     &                      +  lim_am2 * alpha_n**2 + lim_am3 * alpha_n**3) / 
     &                ( lim_am4 + lim_am5 * alpha_n + lim_am6 * alpha_n**2)                                  
               alpha_m = MIN(alpha_m , alpha_m_max) 
               !
               ! Compute stability functions
               Denom = sf_d0  +  sf_d1*alpha_n +  sf_d2*alpha_m 
     &           + sf_d3*alpha_n*alpha_m + sf_d4*alpha_n**2 + sf_d5*alpha_m**2               
               c_mu      = (sf_n0  +  sf_n1*alpha_n +  sf_n2*alpha_m)/Denom
               c_mu_prim = (sf_nb0 + sf_nb1*alpha_n + sf_nb2*alpha_m)/Denom
               !
               ! Finalize the computation of Akv and Akt
               cff = trb( i,j,k,nnew,itke )**2 / epsilon                 
               Akv(i,j,k      )= MAX( cff*c_mu     ,nuwm )
#ifdef MASKING               
     &               *rmask(i,j)
#endif               
               Akt(i,j,k,itemp)= MAX( cff*c_mu_prim,nuws )
#ifdef MASKING    
     &               *rmask(i,j)        
#endif 
#ifdef SALINITY
               Akt(i,j,k,isalt)= Akt(i,j,k,itemp) 
#endif  
               Lscale( i, j , k ) =  cm0 * cm0 * cm0 * cff 
     &                            / sqrt( trb( i,j,k,nnew,itke ) ) 
#ifdef MASKING     
               Lscale( i, j , k ) = Lscale( i, j , k ) * rmask(i,j)
#endif     
            ENDDO
         ENDDO   
         
         DO i=istr,iend 
             Akv(i,j,N)  = MAX( 1.5*Akv(i,j,N-1)
     &                         -0.5*Akv(i,j,N-2), nuwm)          
             Akv(i,j,0)  = MAX( 1.5*Akv(i,j,  1)
     &                         -0.5*Akv(i,j,  2), nuwm)
             Akt(i,j,N,itemp)=  MAX(  1.5*Akt(i,j,N-1,itemp)
     &                          -0.5*Akt(i,j,N-2,itemp), nuws )
             Akt(i,j,0,itemp)=  MAX(  1.5*Akt(i,j,  1,itemp)
     &                          -0.5*Akt(i,j,  2,itemp), nuws )      
# ifdef SALINITY
             Akt(i,j,N,isalt)=Akt(i,j,N,itemp)
             Akt(i,j,0,isalt)=Akt(i,j,0,itemp)
# endif
         ENDDO    
      !--------------------------------------------------        
      ENDDO  !<-- end j-outer loop
      !--------------------------------------------------  

!
!------------------------------------------------
!  Apply boundary conditions
!------------------------------------------------
!
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          do k=0,N
            trb(istr-1,j,k,nnew,itke)=
     &                    trb(istr,j,k,nnew,itke)            
            trb(istr-1,j,k,nnew,igls)=
     &                    trb(istr,j,k,nnew,igls)                         
            Akv(istr-1,j,k      )=Akv(istr,j,k      )
            Akt(istr-1,j,k,itemp)=Akt(istr,j,k,itemp)
#  ifdef SALINITY
            Akt(istr-1,j,k,isalt)=Akt(istr,j,k,isalt)
#  endif
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          do k=0,N
            trb(iend+1,j,k,nnew,itke)=
     &                    trb(iend,j,k,nnew,itke)            
            trb(iend+1,j,k,nnew,igls)=
     &                    trb(iend,j,k,nnew,igls)  
            Akv(iend+1,j,k      )=Akv(iend,j,k      )
            Akt(iend+1,j,k,itemp)=Akt(iend,j,k,itemp)
#  ifdef SALINITY
            Akt(iend+1,j,k,isalt)=Akt(iend,j,k,isalt)
#  endif
          enddo
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          do k=0,N
            trb(i,jstr-1,k,nnew,itke)=
     &                    trb(i,jstr,k,nnew,itke)            
            trb(i,jstr-1,k,nnew,igls)=
     &                    trb(i,jstr,k,nnew,igls)  
            Akv(i,jstr-1,k      )=Akv(i,jstr,k      )
            Akt(i,jstr-1,k,itemp)=Akt(i,jstr,k,itemp)
#  ifdef SALINITY
            Akt(i,jstr-1,k,isalt)=Akt(i,jstr,k,isalt)
#  endif
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          do k=0,N
            trb(i,jend+1,k,nnew,itke)=
     &                    trb(i,jend,k,nnew,itke)            
            trb(i,jend+1,k,nnew,igls)=
     &                    trb(i,jend,k,nnew,igls) 
            Akv(i,jend+1,k)=Akv(i,jend,k)
            Akt(i,jend+1,k,itemp)=Akt(i,jend,k,itemp)
#  ifdef SALINITY
            Akt(i,jend+1,k,isalt)=Akt(i,jend,k,isalt)
#  endif
          enddo
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=0,N
          trb(istr-1,jstr-1,k,nnew,itke)=
     &                    trb(istr,jstr,k,nnew,itke)            
          trb(istr-1,jstr-1,k,nnew,igls)=
     &                    trb(istr,jstr,k,nnew,igls)  
          Akv(istr-1,jstr-1,k      )=Akv(istr,jstr,k      )
          Akt(istr-1,jstr-1,k,itemp)=Akt(istr,jstr,k,itemp)
#  ifdef SALINITY
          Akt(istr-1,jstr-1,k,isalt)=Akt(istr,jstr,k,isalt)
#  endif
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=0,N
          trb(istr-1,jend+1,k,nnew,itke)=
     &                    trb(istr,jend,k,nnew,itke)            
          trb(istr-1,jend+1,k,nnew,igls)=
     &                    trb(istr,jend,k,nnew,igls) 
          Akv(istr-1,jend+1,k      )=Akv(istr,jend,k      )
          Akt(istr-1,jend+1,k,itemp)=Akt(istr,jend,k,itemp)
#  ifdef SALINITY
          Akt(istr-1,jend+1,k,isalt)=Akt(istr,jend,k,isalt)
#  endif
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=0,N
          trb(iend+1,jstr-1,k,nnew,itke)=
     &                    trb(iend,jstr,k,nnew,itke)            
          trb(iend+1,jstr-1,k,nnew,igls)=
     &                    trb(iend,jstr,k,nnew,igls)
          Akv(iend+1,jstr-1,k      )=Akv(iend,jstr,k      )
          Akt(iend+1,jstr-1,k,itemp)=Akt(iend,jstr,k,itemp)
#  ifdef SALINITY
          Akt(iend+1,jstr-1,k,isalt)=Akt(iend,jstr,k,isalt)
#  endif
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=0,N
          trb(iend+1,jend+1,k,nnew,itke)=
     &                    trb(iend,jend,k,nnew,itke)            
          trb(iend+1,jend+1,k,nnew,igls)=
     &                    trb(iend,jend,k,nnew,igls)
          Akv(iend+1,jend+1,k)=Akv(iend,jend,k)
          Akt(iend+1,jend+1,k,itemp)=Akt(iend,jend,k,itemp)
#  ifdef SALINITY
          Akt(iend+1,jend+1,k,isalt)=Akt(iend,jend,k,isalt)
#  endif
        enddo
      endif
#  endif
# endif
!
! Exchange computational margines and/or periodic boundaries:
!--------- ------------- -------- ------ -------- ----------- 
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        trb(START_2D_ARRAY,0,nnew,itke))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        trb(START_2D_ARRAY,0,nnew,igls))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend, Akv)
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(START_2D_ARRAY,0,itemp))
#   ifdef SALINITY
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(START_2D_ARRAY,0,isalt))
#   endif
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Lscale(START_2D_ARRAY,0))
#  endif
!====================================================================== 
      return
      end
!======================================================================       
#endif  /* SOLVE3D & GLS_MIX2017 */
