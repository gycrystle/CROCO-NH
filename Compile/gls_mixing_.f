      SUBROUTINE gls_mixing (tile)
      IMPLICIT NONE
      INTEGER*4         :: tile, trd
      INTEGER*4         :: omp_get_thread_num
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=192,  MMm0=192,  N=512)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=4,  NP_ETA=4,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6, Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, itemp
      integer*4   ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      parameter (itemp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_sed=0)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      integer*4 NGLS
      parameter(NGLS=2)
      integer*4 itke
      parameter(itke=1)
      integer*4 igls
      parameter(igls=2)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_surf=0)
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,5,0:NPP-1)
      common /private_scratch/ A2d,A3d
      real u(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real v(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real t(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hz_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_w(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Huon(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hvom(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Wi(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      common /grid_Wi/Wi
      real wz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,3)
      real nhdu(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N,2)
      real nhdv(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N,2)
      real nhdw(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2)
      real dzdxi(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N)
      real dzdeta(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N)
      real Hz_half(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Pnh(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_wz/wz
      common /ocean_nhdu/nhdu
      common /ocean_nhdv/nhdv
      common /ocean_nhdw/nhdw
      common /ocean_dzdxi/dzdxi
      common /ocean_dzdeta/dzdeta
      common /grid_Hz_half/Hz_half
      common /ocean_pnh/Pnh
      real rho1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real rho(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lmmpi+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lmmpi)/2
      chunk_size_E=(Mmmpi+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mmmpi)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lmmpi)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mmmpi)
      trd=omp_get_thread_num()
      call gls_mixing_tile ( Istr, Iend, Jstr, Jend,
     &                    A3d(1, 1,trd), A2d(1, 2,trd),
     &                    A2d(1, 3,trd), A2d(1, 4,trd), A2d(1, 5,trd),
     &                    A2d(1, 6,trd), A2d(1, 7,trd), A2d(1, 8,trd) )
      RETURN
      END
      SUBROUTINE gls_mixing_tile ( Istr, Iend, Jstr, Jend,
     &                             shear2, diss, ustar_sfc_sq,
     &                             ustar_bot_sq, DC, FC, CF, RH )
      IMPLICIT NONE
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=192,  MMm0=192,  N=512)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=4,  NP_ETA=4,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6, Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, itemp
      integer*4   ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      parameter (itemp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_sed=0)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      integer*4 NGLS
      parameter(NGLS=2)
      integer*4 itke
      parameter(itke=1)
      integer*4 igls
      parameter(igls=2)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_surf=0)
      INTEGER*4         ::   Istr, Iend, Jstr, Jend
      INTEGER*4         ::   i,       j,    k, tind, kref
      INTEGER*4         ::   imin, imax, jmin, jmax
      INTEGER*4         ::   ig,    ig1,  ig2
      REAL            ::  shear2      
     &                               (Istr-2:Iend+2,Jstr-2:Jend+2,0:N-1)
      REAL            ::  diss        (Istr-2:Iend+2,1:N-1)
      REAL            ::  ustar_sfc_sq(Istr-2:Iend+2,Jstr-2:Jend+2      
     &                                                                 )
      REAL            ::  ustar_bot_sq(Istr-2:Iend+2,Jstr-2:Jend+2      
     &                                                                 )
      REAL            ::  DC          (Istr-2:Iend+2,0:N  )
      REAL            ::  FC          (Istr-2:Iend+2,0:N  )
      REAL            ::  CF          (Istr-2:Iend+2,1:N-1)
      REAL            ::  RH          (Istr-2:Iend+2,1:N-1)
      REAL            ::  cff , cff1 , cff2, cff3m, cff3p
      REAL            ::  invk, invG, Bprod, Sprod, epsilon
      REAL            ::  alpha_n, alpha_m, c_mu, c_mu_prim
      REAL            ::  alpha_n_min, alpha_m_max, cm0, cm0inv2, gls
      REAL            ::  flux_top , flux_bot , lgthsc, L_lim, du,dv,dw
      REAL            ::  trb_sfc  , trb_bot  , z0_s, z0_b , gls_min
      REAL            ::  HUon_w   , HVom_w   , trb_min(2), Denom
      REAL, PARAMETER ::  eps_min =  1.0D-12
      REAL, PARAMETER ::  tke_min =  1.0D-06
      REAL, PARAMETER ::  eps     =  1.0D-04
      REAL, PARAMETER ::  nuws    =  1.0D-06
      REAL, PARAMETER ::  nuwm    =  1.0D-06
      REAL, PARAMETER ::  galp    =  0.53D0
      real dt, dtfast, time, time2, time_start, tdays
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
     &      , iprec1, iprec2
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &                       ndtfast, iic, kstp, krhs, knew, next_kstp,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       iprec1, iprec2,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zob
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real R0,T0,S0, Tcoef, Scoef
      real weight(6,0:NWEIGHT)
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
      logical ldefhis
      logical got_tini(NT)
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zob,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1,       tnu2,    tnu4
     &                      , R0,T0,S0,  Tcoef,   Scoef
     &                      , weight
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                      , got_tini
     &                      , ldefhis
      real Akv_bak
      real Akt_bak(NT)
      common /scalars_akt/ Akv_bak, Akt_bak
      real charnok_alpha, zos_hsig_alpha, sz_alpha, crgban_cw,
     &     Zos, Akk_bak, Akp_bak, gls_diff2
      real gls_p, gls_m, gls_n, gls_Kmin, gls_Pmin, gls_c1, gls_c2,
     &     gls_c3m, gls_c3p, gls_sigk, gls_sigp, gls_cmu0,
     &     gls_Gh0, gls_Ghcri, gls_Ghmin, gls_E2
      real my_A1, my_A2, my_B1, my_B2, my_C1, my_C2, my_C3,
     &     my_E1, my_E2, my_Gh0, my_Sq, my_dtfac, my_lmax, my_qmin,
     &     my_B1p2o3, my_B1pm1o3, my_E1o2, my_Sh1, my_Sh2, my_Sm1,
     &     my_Sm2, my_Sm3, my_Sm4
      common /gls_par1/ charnok_alpha, zos_hsig_alpha, sz_alpha,
     &     crgban_cw, Zos, Akk_bak, Akp_bak, gls_diff2
      common /gls_par2/ gls_p, gls_m, gls_n, gls_Kmin, gls_Pmin, gls_c1,
     &     gls_c2,  gls_c3m, gls_c3p, gls_sigk, gls_sigp, gls_cmu0,
     &     gls_Gh0, gls_Ghcri, gls_Ghmin, gls_E2
      common /gls_par3/ my_A1, my_A2, my_B1, my_B2, my_C1, my_C2, my_C3,
     &     my_E1, my_E2, my_Gh0, my_Sq, my_dtfac, my_lmax, my_qmin,
     &     my_B1p2o3, my_B1pm1o3, my_E1o2, my_Sh1, my_Sh2, my_Sm1,
     &     my_Sm2, my_Sm3, my_Sm4
      real gls_s0, gls_s1, gls_s2, gls_s3, gls_s4, gls_s5, gls_s6,
     &     gls_b0, gls_b1, gls_b2, gls_b3, gls_b4, gls_b5, gls_L1,
     &     gls_L2, gls_L3, gls_L4, gls_L5, gls_L6, gls_L7, gls_L8
      common /gls_par4/ gls_s0, gls_s1, gls_s2, gls_s3, gls_s4, gls_s5,
     &     gls_s6,
     &     gls_b0, gls_b1, gls_b2, gls_b3, gls_b4, gls_b5, gls_L1,
     &     gls_L2, gls_L3, gls_L4, gls_L5, gls_L6, gls_L7, gls_L8
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
     &        , avg_vol, avg_rho
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
     &        , avg_vol, avg_rho
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      integer*4 mynode, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE, p_NW,p_NE
      common /comm_setup/ mynode, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     &  p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  day2sec=86400.D0,
     &           sec2day=1.D0/86400.D0, jul_off=2440000.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-9999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      REAL, PARAMETER ::  chk     =  1400.D0/g
      REAL            :: rp,    rm,    rn
      REAL            :: beta1, beta2, beta3m, beta3p
      REAL            :: OneOverSig(2)
      PARAMETER( rp    = 3.0D0 , rm    = 1.5D0 , rn     = -1.0D0        
     &                                                                 )
      PARAMETER( beta1 = 1.44D0, beta2 = 1.92D0, beta3m = -0.4D0, 
     &                                                   beta3p = 1.0D0)
      PARAMETER( OneOverSig = (/ 1.0D0, 0.7692D0 /) )
      REAL, PARAMETER :: e1 =  3.0D0 + 1.D0*rp / rn
      REAL, PARAMETER :: e2 =  1.5D0 + 1.D0*rm / rn
      REAL, PARAMETER :: e3 = -1.0D0 / rn
      REAL, PARAMETER :: smth_a = 1.D0/12.D0
      REAL, PARAMETER :: smth_b = 3.D0/16.D0
       REAL ::  c1   ,c2    ,c3    ,c4    ,c5    , c6
       REAL :: cb1   ,cb2   ,cb3   ,cb4   ,cb5   ,cbb
       REAL :: a1    ,a2    ,a3    ,a5    ,nn
       REAL :: ab1   ,ab2   ,ab3   ,ab5   ,nb
       REAL :: sf_d0 ,sf_d1 ,sf_d2 ,sf_d3 ,sf_d4 , sf_d5
       REAL :: sf_n0 ,sf_n1 ,sf_n2
       REAL :: sf_nb0,sf_nb1,sf_nb2
       REAL :: lim_am0,lim_am1,lim_am2,lim_am3,lim_am4,lim_am5,lim_am6
       PARAMETER(c1=5.D0   ,
     &           c2=0.8D0  ,
     &           c3=1.968D0,
     &           c4=1.136D0,
     &           c5=0.D0   ,
     &           c6=0.4D0   )
       PARAMETER(cb1=5.95D0  ,
     &           cb2=0.6D0   ,
     &           cb3=1.D0    ,
     &           cb4=0.D0    ,
     &           cb5=0.3333D0,
     &           cbb=0.72D0   )
       PARAMETER(  a1 = 0.66666666667D0 - 0.5D0*c2 )
       PARAMETER(  a2 = 1.D0            - 0.5D0*c3 )
       PARAMETER(  a3 = 1.D0            - 0.5D0*c4 )
       PARAMETER(  a5 = 0.5D0           - 0.5D0*c6 )
       PARAMETER( ab1 = 1.D0 - cb2               )
       PARAMETER( ab2 = 1.D0 - cb3               )
       PARAMETER( ab3 = 2.D0*(1.D0-cb4)            )
       PARAMETER( ab5 = 2.D0*cbb*(1.D0-cb5)        )
       PARAMETER( nn  = 0.5D0*c1                 )
       PARAMETER( nb  = cb1                    )
       PARAMETER( sf_d0 = 36.0D0*nn*nn*nn*nb*nb                         
     &                                                                 )
       PARAMETER( sf_d1 = 84.0D0*a5*ab3*nn*nn*nb+36.0D0*ab5*nn*nn*nn*nb 
     &                                                                 )
       PARAMETER( sf_d2 = 9.0D0*(ab2*ab2-ab1*ab1)*nn*nn*nn
     &                  - 12.0D0*(a2*a2-3.D0*a3*a3)*nn*nb*nb)
       PARAMETER( sf_d3 = 12.0D0*a5*ab3*(a2*ab1-3.0D0*a3*ab2)* nn
     &                    + 12.0D0*a5*ab3*(    a3*a3-a2*a2)* nb
     &                    + 12.0D0*   ab5*(3.0D0*a3*a3-a2*a2)*nn*nb     
     &                                                                 )
       PARAMETER( sf_d4 = 48.0D0*a5*a5*ab3*ab3*nn + 
     &                                         36.0D0*a5*ab3*ab5*nn*nn )
       PARAMETER( sf_d5 = 3.0D0*(a2*a2-3.0D0*a3*a3)
     &                       *(ab1*ab1-ab2*ab2)*nn    )
       PARAMETER( sf_n0  = 36.0D0*a1*nn*nn*nb*nb )
       PARAMETER( sf_n1  = - 12.0D0*a5*ab3*(ab1+ab2)*nn*nn
     &                    + 8.0D0*a5*ab3*(6.0D0*a1-a2-3.0D0*a3)*nn*nb
     &                    + 36.0D0*a1*ab5*nn*nn*Nb )
       PARAMETER( sf_n2  = 9.0D0*a1*(ab2*ab2-ab1*ab1)*nn*nn )
       PARAMETER( sf_nb0 = 12.0D0*ab3*nn*nn*nn*nb  )
       PARAMETER( sf_nb1 = 12.0D0*a5*ab3*ab3*nn*nn )
       PARAMETER( sf_nb2 = 9.0D0*a1*ab3*(ab1-ab2)*nn*nn + ( 
     &                                            6.0D0*a1*(a2-3.0D0*a3)
     &                               - 4.0D0*(a2*a2-3.0D0*a3*a3) )*ab3 
     &                                                        * nn * nb)
       PARAMETER( lim_am0 = sf_d0*sf_n0               )
       PARAMETER( lim_am1 = sf_d0*sf_n1 + sf_d1*sf_n0 )
       PARAMETER( lim_am2 = sf_d1*sf_n1 + sf_d4*sf_n0 )
       PARAMETER( lim_am3 = sf_d4*sf_n1               )
       PARAMETER( lim_am4 = sf_d2*sf_n0               )
       PARAMETER( lim_am5 = sf_d2*sf_n1+sf_d3*sf_n0   )
       PARAMETER( lim_am6 = sf_d3*sf_n1               )
      real u(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real v(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real t(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hz_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_w(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Huon(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hvom(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Wi(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      common /grid_Wi/Wi
      real wz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,3)
      real nhdu(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N,2)
      real nhdv(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N,2)
      real nhdw(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2)
      real dzdxi(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N)
      real dzdeta(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N)
      real Hz_half(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Pnh(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_wz/wz
      common /ocean_nhdu/nhdu
      common /ocean_nhdv/nhdv
      common /ocean_nhdw/nhdw
      common /ocean_dzdxi/dzdxi
      common /ocean_dzdeta/dzdeta
      common /grid_Hz_half/Hz_half
      common /ocean_pnh/Pnh
      real rho1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real rho(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real sustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real svstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real bustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real bvstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real bvstrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen,    bms_tstart,   bms_tend,  tsbms,   
     &                                                            sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real stflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_stflx/stflx
      real btflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_btflx/btflx
      real srflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_srflx/srflx
      real h(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real hinv(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real f(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real fomn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real xp(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real xr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real yp(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real yr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_xr/xr /grid_xp/xp /grid_yp/yp /grid_yr/yr
      real pm(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real pmon_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real grdscl(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real diff2_sponge(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /mixing_diff2_sponge/diff2_sponge
      real diff4_sponge(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real diff4(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /mixing_diff4_sponge/diff4_sponge
      common /mixing_diff4/diff4
      real diff3d_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real diff3d_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /mixing_diff3d_u/diff3d_u
      common /mixing_diff3d_v/diff3d_v
      real dRdx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real dRde(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real idRz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /mixing_dRdx/dRdx
      common /mixing_dRde/dRde
      common /mixing_idRz/idRz
      real Akv(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Akt(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real bvf(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /mixing_bvf/ bvf
      real trb(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2,NGLS)
      common /gls_trb/trb
      real Lscale(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /gls_lsc/Lscale
      integer*4 kbl(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /gls_kbl/ kbl
      real hbl(-2:Lm+3+padd_X,-2:Mm+3+padd_E  )
      common /gls_hbl/ hbl
      imin = Istr-1
      imax = Iend+1
      jmin = Jstr-1
      jmax = Jend+1
      cm0     =  ( (a2*a2 - 3.0D0*a3*a3 + 3.0D0*a1*nn)/(3.0D0*nn*nn) 
     &                                                         )**0.25D0
      cm0inv2 = 1.D0/cm0**2
      alpha_n_min = 0.5D0*( - ( sf_d1 + sf_nb0 )
     &             + sqrt(  ( sf_d1 + sf_nb0 )**2
     &            - 4.D0 * sf_d0 *( sf_d4 + sf_nb1 ) ) )
     &                            / ( sf_d4 + sf_nb1 )
      cff     = (cm0**3 )*(tke_min**1.5D0) / eps_min
      gls_min = (cm0**rp)*(tke_min**rm ) * ( cff**rn )
      trb_min(itke) = tke_min
      trb_min(igls) = gls_min
      IF (iic.eq.ntstart) THEN
        DO k=0,N
          DO j=jmin,jmax
            DO i=imin,imax
              trb( i, j, k, nstp, itke ) = trb_min(itke)
              trb( i, j, k, nstp, igls ) = trb_min(igls)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      DO ig = 1,ngls
         DO k=1,N-1
            DO j=jmin,jmax
               DO i=imin,imax
                  trb(i,j,k,nnew,ig)=trb(i,j,k,nstp,ig)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      tind  = nrhs
      DO k=1,N-1
         DO j=jmin,jmax
            DO i=imin,imax
               cff = 1.D0 / ( Hz( i, j, k ) + Hz( i, j, k+1 ) )
               du  = cff*( u(i, j, k+1,tind)+u(i+1, j, k+1,tind)
     &                    -u(i, j, k  ,tind)-u(i+1, j, k  ,tind))
               dv  = cff*( v(i, j, k+1,tind)+v(i, j+1, k+1,tind)
     &                    -v(i, j, k  ,tind)-v(i, j+1, k  ,tind))
               shear2(i,j,k) = du*du + dv*dv
               du  = 0.5D0*pm(i,j)
     &                  *( u(i+1, j, k+1,tind)-u(i, j, k+1,tind)
     &                    +u(i+1, j, k  ,tind)-u(i, j, k  ,tind))
               dv  = 0.5D0*pn(i,j)
     &                  *( v(i, j+1, k+1,tind)-v(i, j, k+1,tind)
     &                    +v(i, j+1, k  ,tind)-v(i, j, k  ,tind))
               dw = 0.D0
               shear2(i,j,k) = 2.D0*(du*du + dv*dv + dw*dw)
     &                                              + shear2(i,j,k)
               du  = 0.125D0*pn(i,j)*
     &                        (u(i,j+1,k  ,tind)+u(i+1,j+1,k  ,tind)
     &                        -u(i,j-1,k  ,tind)-u(i+1,j-1,k  ,tind)
     &                        +u(i,j+1,k+1,tind)+u(i+1,j+1,k+1,tind)
     &                        -u(i,j-1,k+1,tind)-u(i+1,j-1,k+1,tind))
               dv  = 0.125D0*pm(i,j)*
     &                        (v(i+1,j,k  ,tind)+v(i+1,j+1,k  ,tind)
     &                        -v(i-1,j,k  ,tind)-v(i-1,j+1,k  ,tind)
     &                        +v(i+1,j,k+1,tind)+v(i+1,j+1,k+1,tind)
     &                        -v(i-1,j,k+1,tind)-v(i-1,j+1,k+1,tind))
               shear2(i,j,k) = (du+dv)*(du+dv) + shear2(i,j,k)
           ENDDO
         ENDDO
      ENDDO
      DO j=jmin,jmax
         DO i=imin,imax
            ustar_sfc_sq( i, j ) =
     &                        sqrt( (0.5D0*(sustr(i,j)+sustr(i+1,j)))**2
     &                             
     &                           +(0.5D0*(svstr(i,j)+svstr(i,j+1)))**2 )
            ustar_bot_sq( i, j ) =
     &                        sqrt( (0.5D0*(bustr(i,j)+bustr(i+1,j)))**2
     &                             
     &                           +(0.5D0*(bvstr(i,j)+bvstr(i,j+1)))**2 )
         ENDDO
      ENDDO
      DO j=jmin,jmax
         DO i=imin,imax
            DO k=1,N-1
               cff       = (cm0**e1) * ( trb( i,j,k,nstp,itke )**e2 )
     &                               * ( trb( i,j,k,nstp,igls )**e3 )
               diss(i,k) = MAX( cff , eps_min )
            ENDDO
         ENDDO
         DO ig = 1,ngls
            cff=-0.5D0*dt
            DO k=2,N-1
               DO i=imin,imax
                  FC(i,k)= cff* OneOverSig(ig)*
     &                   ( Akv(i,j,k)+Akv(i,j,k-1) ) / Hz(i,j,k)
               ENDDO
            ENDDO
            DO i=imin,imax
               FC(i,1)=0.D0
               FC(i,N)=0.D0
            ENDDO
            DO k=1,N-1
               DO i=imin,imax
                  ig1   = (igls-ig); ig2 = (ig-itke)
                  invk  =     1.D0 / trb( i,j,k,nstp,itke )
                  gls   =          trb( i,j,k,nstp,igls )
                  invG  =  ig1+ig2*(1.D0/gls)
                  cff1  =  ig1+ig2*beta1   * invk*gls
                  cff2  = (ig1+ig2*beta2 ) * invk
                  cff3m =  ig1+ig2*beta3m  * invk*gls
                  cff3p =  ig1+ig2*beta3p  * invk*gls
                  Sprod =  cff1*Akv(i,j,k) * shear2(i,j,k)
                  Bprod = -Akt(i,j,k,itemp)*( cff3m*MAX(bvf(i,j,k),0.D0)
     &                                     +  
     &                                      cff3p*MIN(bvf(i,j,k),0.D0) )
                  cff   =       0.5D0*(Hz(i,j,k)+Hz(i,j,k+1))
                  IF( (Bprod + Sprod) .gt. 0.D0) THEN
                     RH(i,k) = cff*( trb(i,j,k,nnew,ig)
     &                                              + dt*(Bprod+Sprod) )
                     DC(i,k) = cff*(1.D0+dt*cff2*diss(i,k))
     &                                                -FC(i,k)-FC(i,k+1)
                  ELSE
                     RH(i,k) = cff*( trb(i,j,k,nnew,ig) + dt*Sprod  )
                     DC(i,k) = cff*(1.D0+dt*(cff2*diss(i,k)
     &                              -invG*Bprod)) - FC(i,k) - FC(i,k+1)
                  ENDIF
               ENDDO
            ENDDO
            IF( ig == itke ) THEN
               DO i=imin,imax
                  trb_sfc = MAX( tke_min, cm0inv2*ustar_sfc_sq(i,j) )
                  flux_top = 0.D0
                  trb_bot = MAX( tke_min, cm0inv2*ustar_bot_sq(i,j) )
                  flux_bot = 0.D0
                  RH(  i,1   ) = RH(i,  1) + dt*flux_bot
                  RH(  i,N-1 ) = RH(i,N-1) + dt*flux_top
                  trb( i,j,N,nnew,ig ) = trb_sfc
                  trb( i,j,0,nnew,ig ) = trb_bot
               ENDDO
            ELSE
               DO i=imin,imax
                  z0_s=MAX( 1.D-2 , chk*ustar_sfc_sq(i,j) )
                  cff         = 0.5D0*( trb( i,j,N-1,nnew,itke )
     &                            +   trb( i,j,N  ,nnew,itke ) )
                  lgthsc      = vonKar*(0.5D0*Hz(i,j,N)+z0_s)
                  trb_sfc     = MAX(gls_min,(cm0**rp)*(lgthsc**rn)
     &                                                 *(cff**rm))
                  flux_top    = -rn*cm0**(rp+1.D0)
     &                             *vonKar*OneOverSig(igls)
     &                             *(cff**(rm+0.5D0))*(lgthsc**rn)
                  z0_b        = MAX( Zob, 1.D-04 )
                  cff         = 0.5D0*( trb( i,j,  1,nnew,itke )
     &                            +   trb( i,j,  0,nnew,itke ) )
                  lgthsc      = vonKar*(0.5D0*Hz(i,j,1)+z0_b)
                  trb_bot     = MAX(gls_min,(cm0**rp)*(lgthsc**rn)
     &                                               *(cff**rm))
                  flux_bot    =-rn*cm0**(rp+1.D0)
     &                          *vonKar*OneOverSig(igls)
     &                          *(cff**(rm+0.5D0))*(lgthsc**rn)
                  RH( i,  1 ) = RH(i,  1) + dt*flux_bot
                  RH( i,N-1 ) = RH(i,N-1) + dt*flux_top
                  trb( i,j,N,nnew,ig ) = trb_sfc
                  trb( i,j,0,nnew,ig ) = trb_bot
               ENDDO
            ENDIF
            DO i=imin,imax
               cff       =  1.D0/DC(i,N-1)
               CF(i,N-1) = cff*FC(i,N-1)
               RH(i,N-1) = cff*RH(i,N-1)
            ENDDO
            DO k=N-2,1,-1
               DO i=imin,imax
                  cff     =   1.D0/(DC(i,k)-CF(i,k+1)*FC(i,k+1))
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
         ENDDO
      ENDDO
      ig = igls
      DO k=0,N
         DO j=jstr-1,jend+1
            DO i=istr,iend+1
               ustar_sfc_sq (i,j  )=( trb(i  ,j,k,nnew,ig)
     &                   -  trb(i-1,j,k,nnew,ig) )
            ENDDO
         ENDDO
         DO j=jstr,jend+1
            DO i=istr-1,iend+1
               shear2(i,j,0)=( trb(i,j  ,k,nnew,ig)
     &                    - trb(i,j-1,k,nnew,ig) )
            ENDDO
            DO i=istr,iend
              ustar_bot_sq(i,j)=shear2(i,j,0)
     &         +  smth_a*( ustar_sfc_sq(i+1,j)+ustar_sfc_sq(i  ,j-1)
     &                 -ustar_sfc_sq(i  ,j)-ustar_sfc_sq(i+1,j-1))
            ENDDO
         ENDDO
         DO j=jstr,jend
            DO i=istr,iend+1
              ustar_sfc_sq(i,j)=ustar_sfc_sq(i,j  )
     &         + smth_a*( shear2(i,j+1,0)+shear2(i-1,j  ,0)
     &                -shear2(i,j  ,0)-shear2(i-1,j+1,0))
            ENDDO
            DO i=istr,iend
               trb(i,j,k,nnew,ig)=trb(i,j,k,nnew,ig)
     &         + smth_b*( ustar_sfc_sq(i+1,j)-ustar_sfc_sq(i,j)
     &                 +ustar_bot_sq(i,j+1)-ustar_bot_sq(i,j) )
            ENDDO
         ENDDO
      ENDDO
        ig = itke
      DO k=0,N
         DO j=jstr-1,jend+1
            DO i=istr,iend+1
               ustar_sfc_sq (i,j  )=( trb(i  ,j,k,nnew,ig)
     &                   -  trb(i-1,j,k,nnew,ig) )
            ENDDO
         ENDDO
         DO j=jstr,jend+1
            DO i=istr-1,iend+1
               shear2(i,j,0)=( trb(i,j  ,k,nnew,ig)
     &                    - trb(i,j-1,k,nnew,ig) )
            ENDDO
            DO i=istr,iend
              ustar_bot_sq(i,j)=shear2(i,j,0)
     &         +  smth_a*( ustar_sfc_sq(i+1,j)+ustar_sfc_sq(i  ,j-1)
     &                 -ustar_sfc_sq(i  ,j)-ustar_sfc_sq(i+1,j-1))
            ENDDO
         ENDDO
         DO j=jstr,jend
            DO i=istr,iend+1
              ustar_sfc_sq(i,j)=ustar_sfc_sq(i,j  )
     &         + smth_a*( shear2(i,j+1,0)+shear2(i-1,j  ,0)
     &                -shear2(i,j  ,0)-shear2(i-1,j+1,0))
            ENDDO
            DO i=istr,iend
               trb(i,j,k,nnew,ig)=trb(i,j,k,nnew,ig)
     &         + smth_b*( ustar_sfc_sq(i+1,j)-ustar_sfc_sq(i,j)
     &                 +ustar_bot_sq(i,j+1)-ustar_bot_sq(i,j) )
            ENDDO
         ENDDO
      ENDDO
      DO j=jstr,jend
         DO k=1,N-1
            DO i=istr,iend
               L_lim = galp * sqrt( 2.D0* trb(i,j,k,nnew,itke)) /
     &                            ( sqrt(max(eps, bvf(i,j,k)))  )
               cff = (cm0**rp)*(L_lim**rn)*(trb(i,j,k,nnew,itke)**rm)
               trb( i,j,k,nnew,igls ) = MAX( trb( i,j,k,nnew,igls ),cff)
               epsilon = (cm0**e1) * ( trb( i,j,k,nnew,itke )**e2 )
     &                             * ( trb( i,j,k,nnew,igls )**e3 )
               epsilon = MAX(epsilon,eps_min)
               cff     = ( trb(i,j,k,nnew,itke)/epsilon )**2
               alpha_m     = cff*  shear2(i,j,k)
               alpha_n     = cff*     bvf(i,j,k)
               alpha_n     = MIN(  MAX( 0.73D0*alpha_n_min , alpha_n ) 
     &                                                        , 1.0D10 )
               alpha_m_max = ( lim_am0 + lim_am1 * alpha_n
     &                                 + lim_am2 * alpha_n**2
     &                                 + lim_am3 * alpha_n**3) /
     &                         ( lim_am4 + lim_am5 * alpha_n
     &                                   + lim_am6 * alpha_n**2 )
               alpha_m = MIN(alpha_m , alpha_m_max)
               Denom = sf_d0  + sf_d1*alpha_n +  sf_d2*alpha_m
     &                        + sf_d3*alpha_n*alpha_m
     &                        + sf_d4*alpha_n**2 + sf_d5*alpha_m**2
               c_mu      = (sf_n0  + sf_n1*alpha_n  + sf_n2*alpha_m)
     &                                                        /Denom
               c_mu_prim = (sf_nb0 + sf_nb1*alpha_n + sf_nb2*alpha_m)
     &                                                        /Denom
               cff = trb( i,j,k,nnew,itke )**2 / epsilon
               Akv(i,j,k      )= MAX( cff*c_mu , nuwm )
               Akt(i,j,k,itemp)= MAX( cff*c_mu_prim,nuws )
               Akt(i,j,k,isalt)= Akt(i,j,k,itemp)
               Lscale( i, j , k ) =  cm0 * cm0 * cm0 * cff
     &                            / sqrt( trb( i,j,k,nnew,itke ) )
            ENDDO
         ENDDO
         DO i=istr,iend
             Akv(i,j,N)  = MAX( 1.5D0*Akv(i,j,N-1)
     &                         -0.5D0*Akv(i,j,N-2), nuwm)
             Akv(i,j,0)  = MAX( 1.5D0*Akv(i,j,  1)
     &                         -0.5D0*Akv(i,j,  2), nuwm)
             Akt(i,j,N,itemp)=  MAX(  1.5D0*Akt(i,j,N-1,itemp)
     &                          -0.5D0*Akt(i,j,N-2,itemp), nuws )
             Akt(i,j,0,itemp)=  MAX(  1.5D0*Akt(i,j,  1,itemp)
     &                          -0.5D0*Akt(i,j,  2,itemp), nuws )
             Akt(i,j,N,isalt)=Akt(i,j,N,itemp)
             Akt(i,j,0,isalt)=Akt(i,j,0,itemp)
         ENDDO
      ENDDO
      do j=Jstr,Jend
        do i=Istr,Iend
          kbl(i,j)=1
          hbl(i,j)=z_w(i,j,N)-z_r(i,j,1)
        enddo
      enddo
      kref=max(1,N-3)
      do k=1,kref
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=rho1(i,j,k)-rho1(i,j,kref)
            if (cff.gt.0.01D0) then
              kbl(i,j)=k
              hbl(i,j)=z_w(i,j,N)-z_r(i,j,k)
            endif
          enddo
        enddo
      enddo
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        trb(-2,-2,0,nnew,itke))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        trb(-2,-2,0,nnew,igls))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend, Akv)
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(-2,-2,0,itemp))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(-2,-2,0,isalt))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Lscale(-2,-2,0))
      return
      end
