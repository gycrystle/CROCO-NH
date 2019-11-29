      subroutine t3dmix (tile)
      implicit none
      integer*4 tile, itrc, trd, omp_get_thread_num
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
      do itrc=1,NT
        call t3dmix_tile (istr,iend,jstr,jend, itrc,   A3d(1,1,trd),
     &                                                 A3d(1,2,trd),
     &                                    A3d(1,3,trd),A3d(1,4,trd),
     &                      A2d(1,1,trd), A2d(1,2,trd),A2d(1,3,trd),
     &                      A2d(1,5,trd), A2d(1,7,trd),A2d(1,9,trd)
     &                   ,A2d(1,10,trd),A2d(1,11,trd), A2d(1,12,trd),
     &                                   A2d(1,13,trd),A2d(1,14,trd)
     &                    )
      enddo
      return
      end
      subroutine t3dmix_tile (istr,iend,jstr,jend, itrc, LapT,
     &                                                          Akz,
     &                                                diff3u,diff3v,
     &                                      FX,FE,FC,dTdr, dTdx,dTde
     &                                              ,FFC,CF,BC,CD,DC
     &                        )
      implicit none
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
      integer*4 istr,iend,jstr,jend, itrc, i,j,k,k1,k2, kmld,
     &        imin,imax,jmin,jmax, indx, idx,ide
      real LapT (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &   diff3u (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &   diff3v (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &      Akz (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &       CF (Istr-2:Iend+2,0:N),
     &       DC (Istr-2:Iend+2,0:N),
     &       CD (Istr-2:Iend+2,0:N),
     &       BC (Istr-2:Iend+2,0:N),
     &       FFC(Istr-2:Iend+2,0:N),
     &       FX (Istr-2:Iend+2,Jstr-2:Jend+2),
     &       FE (Istr-2:Iend+2,Jstr-2:Jend+2),
     &       FC (Istr-2:Iend+2,Jstr-2:Jend+2,2),
     &     dTdr (Istr-2:Iend+2,Jstr-2:Jend+2,2),
     &     dTdx (Istr-2:Iend+2,Jstr-2:Jend+2,2),
     &     dTde (Istr-2:Iend+2,Jstr-2:Jend+2,2)
       real cff,cff1,
     &      TRIADS1,TRIADS2,TRIADS3,TRIADS4,sumX,sumE,sig,
     &      SLOPEXQ1,SLOPEXQ2,SLOPEXQ3,SLOPEXQ4,
     &      SLOPEYQ1,SLOPEYQ2,SLOPEYQ3,SLOPEYQ4
       real wgt(0:4)
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
      integer*4 IstrR,IendR,JstrR,JendR
      if (.not.WEST_INTER) then
        IstrR=Istr-2
      else
        IstrR=Istr
      endif
      if (.not.EAST_INTER) then
        IendR=Iend+2
      else
        IendR=Iend
      endif
      if (.not.SOUTH_INTER) then
        JstrR=Jstr-2
      else
        JstrR=Jstr
      endif
      if (.not.NORTH_INTER) then
        JendR=Jend+2
      else
        JendR=Jend
      endif
      wgt=(/0.D0,1.D0,0.5D0,0.33333333333D0,0.25D0/)
      imin=istr-1
      imax=iend+1
      jmin=jstr-1
      jmax=jend+1
      do k=1,N
        do j=jmin,jmax
          do i=imin,imax+1
            diff3u(i,j,k)=
     &                     +sqrt(
     &                      0.5D0*(diff4(i,j,itrc)+diff4(i-1,j,itrc))
     &                     +diff3d_u(i,j,k)
     &                           )
          enddo
        enddo
        do j=jmin,jmax+1
          do i=imin,imax
            diff3v(i,j,k)=
     &                     +sqrt(
     &                      0.5D0*(diff4(i,j,itrc)+diff4(i,j-1,itrc))
     &                     +diff3d_v(i,j,k)
     &                           )
          enddo
        enddo
      enddo
      k2=1
      do k=0,N,+1
       k1=k2
       k2=3-k1
        if (k.lt.N) then
          do j=jmin,jmax
            do i=imin,imax+1
              cff=0.5D0*(pm(i,j)+pm(i-1,j))
              dTdx(i,j,k2)=cff*( t(i  ,j,k+1,nstp,itrc)
     &                          -t(i-1,j,k+1,nstp,itrc)
     &                         )
            enddo
          enddo
          do j=jmin,jmax+1
            do i=imin,imax
              cff=0.5D0*(pn(i,j)+pn(i,j-1))
              dTde(i,j,k2)=cff*( t(i,j  ,k+1,nstp,itrc)
     &                          -t(i,j-1,k+1,nstp,itrc)
     &                         )
            enddo
          enddo
        endif
        if (k.eq.0 .or. k.eq.N) then
          do j=jmin-1,jmax+1
            do i=imin-1,imax+1
               FC  (i,j,k2) = 0.0D0
               Akz (i,j,k )= 0.0D0
             enddo
          enddo
          if (k.eq.0) then
            do j=jmin-1,jmax+1
              do i=imin-1,imax+1
                dTdr(i,j,k2)= idRz(i,j,1)*( t(i,j,2,nstp,itrc)
     &                                    - t(i,j,1,nstp,itrc)
     &                                     )
              enddo
            enddo
          endif
        else
          do j=jmin-1,jmax+1
            do i=imin-1,imax+1
              FC(i,j,k2)  = idRz(i,j,k)*( z_r (i,j,k+1)-z_r (i,j,k) )
              dTdr(i,j,k2)= idRz(i,j,k)*( t(i,j,k+1,nstp,itrc)
     &                                  - t(i,j,k  ,nstp,itrc)
     &                                  )
            enddo
          enddo
        endif
        if (k.gt.0) then
          cff=0.5D0
          do j=jmin,jmax
            do i=imin,imax+1
              FX(i,j)=cff*diff3u(i,j,k)*(Hz(i,j,k)+Hz(i-1,j,k))
     &               *on_u(i,j)*(   dTdx(i,j,k1) -
     &       0.5D0*( min(dRdx(i,j,k),0.D0)*(dTdr(i-1,j,k1)+dTdr(i,j,k2))
     &            +max(dRdx(i,j,k),0.D0)*(dTdr(i-1,j,k2)+dTdr(i,j,k1)))
     &                                                              )
            enddo
          enddo
          do j=jmin,jmax+1
            do i=imin,imax
              FE(i,j)=cff*diff3v(i,j,k)*(Hz(i,j,k)+Hz(i,j-1,k))
     &               *om_v(i,j)*(   dTde(i,j,k1) -
     &       0.5D0*( min(dRde(i,j,k),0.D0)*(dTdr(i,j-1,k1)+dTdr(i,j,k2))
     &            +max(dRde(i,j,k),0.D0)*(dTdr(i,j-1,k2)+dTdr(i,j,k1)))
     &                                                              )
            enddo
          enddo
          if (k.lt.N) then
            do j=jmin,jmax
              do i=imin,imax
                TRIADS1=dRdx(i  ,j,k  )*dTdr(i,j,k2)-dTdx(i  ,j,k1)
                TRIADS2=dRdx(i  ,j,k+1)*dTdr(i,j,k2)-dTdx(i  ,j,k2)
                TRIADS3=dRdx(i+1,j,k+1)*dTdr(i,j,k2)-dTdx(i+1,j,k2)
                TRIADS4=dRdx(i+1,j,k  )*dTdr(i,j,k2)-dTdx(i+1,j,k1)
                sumX=0.D0
                idx=0
                if (dRdx(i  ,j,k  ) .lt. 0.D0) then
                 sumX=     diff3u(i  ,j,k  )*dRdx(i  ,j,k  )*TRIADS1
                 idx=idx+1
                endif
                if (dRdx(i  ,j,k+1) .gt. 0.D0) then
                 sumX=sumX+diff3u(i  ,j,k+1)*dRdx(i  ,j,k+1)*TRIADS2
                 idx=idx+1
                endif
                if (dRdx(i+1,j,k+1) .lt. 0.D0) then
                 sumX=sumX+diff3u(i+1,j,k+1)*dRdx(i+1,j,k+1)*TRIADS3
                 idx=idx+1
                endif
                if (dRdx(i+1,j,k  ) .gt. 0.D0) then
                 sumX=sumX+diff3u(i+1,j,k  )*dRdx(i+1,j,k  )*TRIADS4
                 idx=idx+1
                endif
                TRIADS1=dRde(i,j  ,k  )*dTdr(i,j,k2)-dTde(i,j  ,k1)
                TRIADS2=dRde(i,j  ,k+1)*dTdr(i,j,k2)-dTde(i,j  ,k2)
                TRIADS3=dRde(i,j+1,k+1)*dTdr(i,j,k2)-dTde(i,j+1,k2)
                TRIADS4=dRde(i,j+1,k  )*dTdr(i,j,k2)-dTde(i,j+1,k1)
                sumE=0.D0
                ide=0
                if (dRde(i,j  ,k  ) .lt. 0.D0) then
                 sumE=     diff3v(i,j  ,k  )*dRde(i,j  ,k  )*TRIADS1
                 ide=ide+1
                endif
                if (dRde(i,j  ,k+1) .gt. 0.D0) then
                 sumE=sumE+diff3v(i,j  ,k+1)*dRde(i,j  ,k+1)*TRIADS2
                 ide=ide+1
                endif
                if (dRde(i,j+1,k+1) .lt. 0.D0) then
                 sumE=sumE+diff3v(i,j+1,k+1)*dRde(i,j+1,k+1)*TRIADS3
                 ide=ide+1
                endif
                if (dRde(i,j+1,k  ) .gt. 0.D0) then
                 sumE=sumE+diff3v(i,j+1,k  )*dRde(i,j+1,k  )*TRIADS4
                 ide=ide+1
                endif
                FC(i,j,k2)=(sumX*wgt(idx)+sumE*wgt(ide))*FC(i,j,k2)
              enddo
            enddo
          endif
          do j=jmin,jmax
            do i=imin,imax
              LapT(i,j,k)=( pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                       +FE(i,j+1)-FE(i,j))
     &                     +FC(i,j,k2)-FC(i,j,k1)    )/Hz(i,j,k)
            enddo
          enddo
        endif
      enddo
      k2=1
      do k=0,N,+1
       k1=k2
       k2=3-k1
        if (k.lt.N) then
          do j=jstr,jend
            do i=istr,iend+1
              cff=0.5D0*(pm(i,j)+pm(i-1,j))
              dTdx(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i-1,j,k+1))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              cff=0.5D0*(pn(i,j)+pn(i,j-1))
              dTde(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i,j-1,k+1))
            enddo
          enddo
        endif
        if (k.eq.0 .or. k.eq.N) then
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              FC(i,j,k2)=0.0D0
              dTdr(i,j,k2)=0.0D0
              Akz (i,j,k )= 0.0D0
            enddo
          enddo
          if (k.eq.0) then
            do j=jstr-1,jend+1
              do i=istr-1,iend+1
                dTdr(i,j,k2)= idRz(i,j,1)
     &                    *( LapT(i,j,2)-LapT(i,j,1) )
              enddo
            enddo
          endif
        else
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              FC(i,j,k2)  = idRz(i,j,k)*( z_r (i,j,k+1)-z_r (i,j,k) )
              dTdr(i,j,k2)= idRz(i,j,k)*( LapT(i,j,k+1)-LapT(i,j,k) )
            enddo
          enddo
        endif
        if (k.gt.0) then
          cff=0.5D0
          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=-cff*diff3u(i,j,k)*(Hz(i,j,k)+Hz(i-1,j,k))
     &         *on_u(i,j)*(   dTdx(i  ,j,k1) -
     &       0.5D0*( min(dRdx(i,j,k),0.D0)*(dTdr(i-1,j,k1)+dTdr(i,j,k2))
     &            +max(dRdx(i,j,k),0.D0)*(dTdr(i-1,j,k2)+dTdr(i,j,k1)))
     &                                                              )
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=-cff*diff3v(i,j,k)*(Hz(i,j,k)+Hz(i,j-1,k))
     &        *om_v(i,j)*(  dTde(i,j  ,k1) -
     &       0.5D0*( min(dRde(i,j,k),0.D0)*(dTdr(i,j-1,k1)+dTdr(i,j,k2))
     &            +max(dRde(i,j,k),0.D0)*(dTdr(i,j-1,k2)+dTdr(i,j,k1)))
     &                                                              )
            enddo
          enddo
          if (k.lt.N) then
            do j=jstr,jend
              do i=istr,iend
                TRIADS1=dRdx(i  ,j,k  )*dTdr(i,j,k2)-dTdx(i  ,j,k1)
                TRIADS2=dRdx(i  ,j,k+1)*dTdr(i,j,k2)-dTdx(i  ,j,k2)
                TRIADS3=dRdx(i+1,j,k+1)*dTdr(i,j,k2)-dTdx(i+1,j,k2)
                TRIADS4=dRdx(i+1,j,k  )*dTdr(i,j,k2)-dTdx(i+1,j,k1)
                sumX=0.D0
                idx=0
                if (dRdx(i  ,j,k  ) .lt. 0.D0) then
                 sumX=     diff3u(i  ,j,k  )*dRdx(i  ,j,k  )*TRIADS1
                 idx=idx+1
                endif
                if (dRdx(i  ,j,k+1) .gt. 0.D0) then
                 sumX=sumX+diff3u(i  ,j,k+1)*dRdx(i  ,j,k+1)*TRIADS2
                 idx=idx+1
                endif
                if (dRdx(i+1,j,k+1) .lt. 0.D0) then
                 sumX=sumX+diff3u(i+1,j,k+1)*dRdx(i+1,j,k+1)*TRIADS3
                 idx=idx+1
                endif
                if (dRdx(i+1,j,k  ) .gt. 0.D0) then
                 sumX=sumX+diff3u(i+1,j,k  )*dRdx(i+1,j,k  )*TRIADS4
                 idx=idx+1
                endif
                TRIADS1=dRde(i,j  ,k  )*dTdr(i,j,k2)-dTde(i,j  ,k1)
                TRIADS2=dRde(i,j  ,k+1)*dTdr(i,j,k2)-dTde(i,j  ,k2)
                TRIADS3=dRde(i,j+1,k+1)*dTdr(i,j,k2)-dTde(i,j+1,k2)
                TRIADS4=dRde(i,j+1,k  )*dTdr(i,j,k2)-dTde(i,j+1,k1)
                sumE=0.D0
                ide=0
                if (dRde(i,j  ,k  ) .lt. 0.D0) then
                 sumE=     diff3v(i,j  ,k  )*dRde(i,j  ,k  )*TRIADS1
                 ide=ide+1
                endif
                if (dRde(i,j  ,k+1) .gt. 0.D0) then
                 sumE=sumE+diff3v(i,j  ,k+1)*dRde(i,j  ,k+1)*TRIADS2
                 ide=ide+1
                endif
                if (dRde(i,j+1,k+1) .lt. 0.D0) then
                 sumE=sumE+diff3v(i,j+1,k+1)*dRde(i,j+1,k+1)*TRIADS3
                 ide=ide+1
                endif
                if (dRde(i,j+1,k  ) .gt. 0.D0) then
                 sumE=sumE+diff3v(i,j+1,k  )*dRde(i,j+1,k  )*TRIADS4
                 ide=ide+1
                endif
                SLOPEXQ1=(FC(i,j,k2)*dRdx(i  ,j,k  ))**2
                SLOPEXQ2=(FC(i,j,k2)*dRdx(i  ,j,k+1))**2
                SLOPEXQ3=(FC(i,j,k2)*dRdx(i+1,j,k+1))**2
                SLOPEXQ4=(FC(i,j,k2)*dRdx(i+1,j,k  ))**2
                SLOPEYQ1=(FC(i,j,k2)*dRde(i,j  ,k  ))**2
                SLOPEYQ2=(FC(i,j,k2)*dRde(i,j  ,k+1))**2
                SLOPEYQ3=(FC(i,j,k2)*dRde(i,j+1,k+1))**2
                SLOPEYQ4=(FC(i,j,k2)*dRde(i,j+1,k  ))**2
                cff = 1.D0/(z_r(i,j,k+1)-z_r(i,j,k))
                Akz(i,j,k) = 8.D0*( max(
     &                       diff3u(i  ,j,k  )*SLOPEXQ1,
     &                       diff3u(i  ,j,k+1)*SLOPEXQ2,
     &                       diff3u(i+1,j,k+1)*SLOPEXQ3,
     &                       diff3u(i+1,j,k  )*SLOPEXQ4)
     &                           +max(
     &                       diff3v(i,j  ,k  )*SLOPEYQ1,
     &                       diff3v(i,j  ,k+1)*SLOPEYQ2,
     &                       diff3v(i,j+1,k+1)*SLOPEYQ3,
     &                       diff3v(i,j+1,k  )*SLOPEYQ4)
     &                          )*( max(
     &                diff3u(i  ,j,k  )*(pm(i  ,j)**2+SLOPEXQ1*cff**2),
     &                diff3u(i  ,j,k+1)*(pm(i  ,j)**2+SLOPEXQ2*cff**2),
     &                diff3u(i+1,j,k+1)*(pm(i+1,j)**2+SLOPEXQ3*cff**2),
     &                diff3u(i+1,j,k  )*(pm(i+1,j)**2+SLOPEXQ4*cff**2))
     &                             +max(
     &                diff3v(i,j  ,k  )*(pn(i,j  )**2+SLOPEYQ1*cff**2),
     &                diff3v(i,j  ,k+1)*(pn(i,j  )**2+SLOPEYQ2*cff**2),
     &                diff3v(i,j+1,k+1)*(pn(i,j+1)**2+SLOPEYQ3*cff**2),
     &                diff3v(i,j+1,k  )*(pn(i,j+1)**2+SLOPEYQ4*cff**2))
     &                            )
                FC(i,j,k2)=-(sumX*wgt(idx)+sumE*wgt(ide))*FC(i,j,k2)
             enddo
            enddo
          endif
          do j=jstr,jend
            do i=istr,iend
              t(i,j,k,nnew,itrc)=Hz(i,j,k)*t(i,j,k,nnew,itrc)
     &         + dt*(
     &                   pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                    +FE(i,j+1)-FE(i,j))
     &                  +FC(i,j,k2)-FC(i,j,k1)    )
            enddo
          enddo
        endif
      enddo
       do i=Istr,Iend
            DC(i,0)=dt*pn(i,j)*pm(i,j)
       enddo
      do j=jstr,jend
          indx=min(itrc,isalt)
          do i=istr,iend
            do k=1,N-1
              CD(i,k) = Akz(i,j,k)*(
     &           t(i,j,k+1,nstp,itrc)-t(i,j,k,nstp,itrc)
     &           )  / ( z_r(i,j,k+1)-z_r(i,j,k) )
            enddo
            CD(i,0) = 0.D0
            CD(i,N) = 0.D0
          enddo
          do i=istr,iend
            FFC(i,1)=dt*(Akt(i,j,1,indx)+Akz(i,j,1))
     &                               /( z_r(i,j,2)-z_r(i,j,1) )
            BC(i,1)=DC(i,0)*Wi(i,j,1)
            cff=1.D0/(Hz(i,j,1)      +FFC(i,1)+max(BC(i,1),0.D0))
            CF(i,1)=cff*(           FFC(i,1)-min(BC(i,1),0.D0))
            DC(i,1)= cff*(t(i,j,1,nnew,itrc)-dt*(CD(i,1)-CD(i,0)))
          enddo
          do k=2,N-1,+1
            do i=istr,iend
              FFC(i,k)=dt*(Akt(i,j,k,indx)+Akz(i,j,k))
     &                              /( z_r(i,j,k+1)-z_r(i,j,k) )
              BC(i,k)=DC(i,0)*Wi(i,j,k)
              cff=1.D0/(      Hz(i,j,k) +FFC(i,k)+max(BC(i,k),0.D0)
     &                              +FFC(i,k-1)-min(BC(i,k-1),0.D0)
     &                   -CF(i,k-1)*(FFC(i,k-1)+max(BC(i,k-1),0.D0))
     &                                                          )
              CF(i,k)=cff*(FFC(i,k)-min(BC(i,k),0.D0))
              DC(i,k)=cff*( t(i,j,k,nnew,itrc) +DC(i,k-1)*(
     &                          FFC(i,k-1)+max(BC(i,k-1),0.D0) )
     &                                    -dt*(CD(i,k)-CD(i,k-1))
     &                                                          )
            enddo
          enddo
          do i=istr,iend
            t(i,j,N,nnew,itrc)=( t(i,j,N,nnew,itrc)
     &                           -dt*(CD(i,N)-CD(i,N-1))
     &                                           +DC(i,N-1)*(
     &                                FFC(i,N-1)+max(BC(i,N-1),0.D0) )
     &               )/( Hz(i,j,N) +FFC(i,N-1)-min(BC(i,N-1),0.D0)
     &                      -CF(i,N-1)*(FFC(i,N-1)+max(BC(i,N-1),0.D0))
     &                                                            )
        enddo
          do k=N-1,1,-1
            do i=istr,iend
              t(i,j,k,nnew,itrc)=DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc)
            enddo
          enddo
      enddo
        call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                               t(-2,-2,1,nnew,itrc))
      return
      end
