      subroutine check_switches1 (ierr)
      implicit none
      integer*4 ierr, is,ie, iexample
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
      integer*4 max_opt_size
      parameter (max_opt_size=3400)
      character*3400 Coptions,srcs
      common /strings/ Coptions,srcs
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
      if (mynode.eq.0) write(stdout,'(/1x,A/)')
     &      'Activated C-preprocessing Options:'
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      iexample=0
      is=1
      iexample=iexample+1
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TAYLOR_INST'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TAYLOR_INST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'MPI'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NEW_S_COORD'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NEW_S_COORD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SOLVE3D'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SOLVE3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_COR'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_COR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_ADV'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_ADV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SALINITY'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SALINITY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TS_DIF4'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_DIF4'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'VADV_ADAPT_IMP'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_ADAPT_IMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLS_MIXING'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_MIXING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLS_KEPSILON'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_KEPSILON'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'CANUTO_A'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CANUTO_A'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_BSFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BSFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_BTFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BTFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_GRID'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_GRID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_INITIAL'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_INITIAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_SMFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SMFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_STFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_STFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'EW_PERIODIC'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='EW_PERIODIC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_SRFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SRFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_SSFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SSFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_MIX_S'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_MIX_S'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NHMG'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NHMG'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NS_PERIODIC'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NS_PERIODIC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TS_HADV_RSUP5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_RSUP5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_HADV_UP5'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_UP5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLS_MIXING_3D'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_MIXING_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TS_MIX_GEO'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_GEO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'M2FILTER_POWER'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_POWER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'VAR_RHO_2D'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VAR_RHO_2D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'RESET_RHO0'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RESET_RHO0'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'HZR'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HZR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_HADV_UP3'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'DIF_COEF_3D'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIF_COEF_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_VADV_SPLINES'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VADV_SPLINES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TS_HADV_C6'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_C6'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TS_MIX_IMP'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_IMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NTRA_T3DMIX'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NTRA_T3DMIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TS_VADV_SPLINES'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_VADV_SPLINES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'LIMIT_BSTRESS'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LIMIT_BSTRESS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NF_CLOBBER'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NF_CLOBBER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(/)')
      if (iexample.eq.0) then
        if (mynode.eq.0) write(stdout,'(1x,A)')
     & 'ERROR in "cppdefs.h": no configuration is specified.'
        ierr=ierr+1
      elseif (iexample.gt.1) then
        if (mynode.eq.0) write(stdout,'(1x,A)')
     & 'ERROR: more than one configuration in "cppdefs.h".'
        ierr=ierr+1
      endif
      return
  99  if (mynode.eq.0) write(stdout,'(/1x,A,A/14x,A)')
     &  'CHECKDEFS -- ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
