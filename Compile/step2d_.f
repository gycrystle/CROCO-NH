      subroutine step2d (tile)
      implicit none
      integer*4 tile, trd
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
C$    integer*4 omp_get_thread_num
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
      trd=0
C$    trd=omp_get_thread_num()
      call step2D_FB_tile ( Istr,Iend,Jstr,Jend, A2d(1,1,trd)
     &                    , A2d(1, 2,trd), A2d(1, 3,trd), A2d(1, 4,trd)
     &                    , A2d(1, 5,trd), A2d(1, 6,trd), A2d(1, 7,trd)
     &                    , A2d(1, 8,trd), A2d(1, 9,trd)
     &                    , A2d(1,10,trd), A2d(1,11,trd)
     &                    , A2d(1,12,trd), A2d(1,13,trd)
     &                    )
      return
      end
      subroutine step2D_FB_tile (Istr,Iend,Jstr,Jend, zeta_new,
     &                           Dnew,rubar,rvbar,
     &                           Drhs, UFx,UFe,
     &                           VFx,VFe
     &                          ,urhs,vrhs
     &                          ,DUon,DVom
     &                          )
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
      integer*4 Istr,Iend,Jstr,Jend, i,j, kbak, kold,
     &         err,
     &        imin,imax,jmin,jmax
      real sum_c
      real    VMAX,VMAXL
      real zeta_new(Istr-2:Iend+2,Jstr-2:Jend+2),  cff,
     &         Dnew(Istr-2:Iend+2,Jstr-2:Jend+2),  cff0,
     &        rubar(Istr-2:Iend+2,Jstr-2:Jend+2),  cff1,
     &        rvbar(Istr-2:Iend+2,Jstr-2:Jend+2),  cff2,
     &         Drhs(Istr-2:Iend+2,Jstr-2:Jend+2),  cff3,
     &          UFx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &          UFe(Istr-2:Iend+2,Jstr-2:Jend+2),  DUnew,
     &          VFx(Istr-2:Iend+2,Jstr-2:Jend+2),  DVnew,
     &          VFe(Istr-2:Iend+2,Jstr-2:Jend+2)
      real     urhs(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         vrhs(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         DUon(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         DVom(Istr-2:Iend+2,Jstr-2:Jend+2)
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
      real zeta(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      real ubar(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      real vbar(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar
      real nh_ubar(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real nh_vbar(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real nh_wcor(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nh_wcor/nh_ubar,nh_vbar,nh_wcor
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
      real rhoA(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rhoS(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rvfrc(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rufrc_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real rvfrc_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real DU_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,5)
      real DV_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,5)
      real DU_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real DV_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      include 'mpif.h'
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
      if (iif.eq.1) then
        kbak=kstp
        kold=kstp
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      elseif (iif.eq.1+1) then
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      else
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak-1
        if (kold.lt.1) kold=4
        cff1= 1.781105D0
        cff2=-1.06221D0
        cff3= 0.281105D0
      endif
      imin=Istr-2
      imax=Iend+1
      jmin=Jstr-2
      jmax=Jend+1
      do j=jmin,jmax
        do i=imin,imax
          Drhs(i,j)=cff1*zeta(i,j,kstp)+cff2*zeta(i,j,kbak)
     &                                 +cff3*zeta(i,j,kold)
     &                                             + h(i,j)
        enddo
      enddo
      do j=Jstr-1,Jend+1
        do i=imin+1,imax
          urhs(i,j)=cff1*ubar(i,j,kstp) +cff2*ubar(i,j,kbak)
     &                                  +cff3*ubar(i,j,kold)
          DUon(i,j)=0.5D0*(Drhs(i,j)+Drhs(i-1,j))*on_u(i,j)*( urhs(i,j)
     &                                                              )
        enddo
      enddo
      do j=jmin+1,jmax
        do i=Istr-1,Iend+1
          vrhs(i,j)=cff1*vbar(i,j,kstp) +cff2*vbar(i,j,kbak)
     &                                  +cff3*vbar(i,j,kold)
          DVom(i,j)=0.5D0*(Drhs(i,j)+Drhs(i,j-1))*om_v(i,j)*( vrhs(i,j)
     &                                                              )
        enddo
      enddo
      if (iif.eq.1) then
        cff0=0.D0
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      elseif (iif.eq.1+1) then
        cff0= 1.0833333333333D0
        cff1=-0.1666666666666D0
        cff2= 0.0833333333333D0
        cff3= 0.D0
      else
        cff0=0.614D0
        cff1=0.285D0
        cff2=0.088D0
        cff3=0.013D0
      endif
      do j=Jstr-1,Jend
        do i=Istr-1,Iend
          zeta_new(i,j)=zeta(i,j,kstp) + dtfast*pm(i,j)*pn(i,j)
     &                                   *(DUon(i,j)-DUon(i+1,j  )
     &                                    +DVom(i,j)-DVom(i  ,j+1))
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr-1,Iend
          zeta_new(i,j)=zeta_new(i,j)
          Dnew(i,j)=zeta_new(i,j)+h(i,j)
          UFx(i,j)=cff0*zeta_new(i,j) +cff1*zeta(i,j,kstp)
     &             +cff2*zeta(i,j,kbak)+cff3*zeta(i,j,kold)
          UFe(i,j)=(1.D0+rhoS(i,j))*UFx(i,j)
          VFe(i,j)=UFe(i,j)*UFx(i,j)
          VFx(i,j)=UFx(i,j)*(rhoS(i,j)-rhoA(i,j))
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr-1,Iend
          zeta(i,j,knew)=zeta_new(i,j)
        enddo
      enddo
      call zetabc_tile (Istr,Iend,Jstr,Jend)
      cff1=weight(1,iif)
      cff2=weight(2,iif)
      if (iif.eq.1) then
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=cff1*zeta(i,j,knew)
            DU_avg1(i,j,nnew)=0.D0
            DV_avg1(i,j,nnew)=0.D0
            DU_avg2(i,j)=cff2*DUon(i,j)
            DV_avg2(i,j)=cff2*DVom(i,j)
          enddo
        enddo
      else
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
            DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
            DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
          enddo
        enddo
      endif
      cff=0.5D0*g
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=cff*on_u(i,j)*( (h(i-1,j)+h(i,j))*(UFe(i-1,j)
     &                        -UFe(i,j)) +VFe(i-1,j)-VFe(i,j)
     &              +(h(i-1,j)-h(i,j))*( VFx(i-1,j)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i-1,j)-rhoA(i,j))
     &                                      *(UFx(i-1,j)-UFx(i,j)))
     &                                                              )
          rvbar(i,j)=cff*om_v(i,j)*( (h(i,j-1)+h(i,j))*(UFe(i,j-1)
     &                        -UFe(i,j)) +VFe(i,j-1)-VFe(i,j)
     &              +(h(i,j-1)-h(i,j))*( VFx(i,j-1)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i,j-1)-rhoA(i,j))
     &                                      *(UFx(i,j-1)-UFx(i,j)))
     &                                                              )
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr-1,Iend
          UFx(i,j)=0.25D0*(DUon(i,j)+DUon(i+1,j))
     &                 *(urhs(i,j)+urhs(i+1,j))
          VFx(i+1,j)=0.25D0*(DUon(i+1,j)+DUon(i+1,j-1))
     &                   *(vrhs(i+1,j)+vrhs(i,j))
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr,Iend
          VFe(i,j)=0.25D0*(DVom(i,j)+DVom(i,j+1))
     &                 *(vrhs(i,j)+vrhs(i,j+1))
          UFe(i,j+1)=0.25D0*(DVom(i,j+1)+DVom(i-1,j+1))
     &                   *(urhs(i,j+1)+urhs(i,j))
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)-UFx(i,j)+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i,j)
          rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j)+VFe(i,j-1)
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr-1,Iend
          cff=Drhs(i,j)*(
     &                   fomn(i,j)
     &                   )
          UFx(i,j)=cff*(vrhs(i,j)+vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i,j)+urhs(i+1,j))
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)+0.25D0*(UFx(i,j)+UFx(i-1,j))
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rvbar(i,j)=rvbar(i,j)-0.25D0*(VFe(i,j)+VFe(i,j-1))
        enddo
      enddo
      if (rdrg2.gt.0.D0) then
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=0.25D0*( vbar(i  ,j,kstp)+vbar(i  ,j+1,kstp)
     &                +vbar(i-1,j,kstp)+vbar(i-1,j+1,kstp))
            rubar(i,j)=rubar(i,j) - ubar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,kstp)*ubar(i,j,kstp)+cff*cff)
     &                               )*om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=0.25D0*( ubar(i,j  ,kstp)+ubar(i+1,j  ,kstp)
     &                +ubar(i,j-1,kstp)+ubar(i+1,j-1,kstp))
            rvbar(i,j)=rvbar(i,j) - vbar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,kstp)*vbar(i,j,kstp))
     &                               )*om_v(i,j)*on_v(i,j)
          enddo
        enddo
      else if (rdrg.gt.0.0D0) then
        do j=Jstr,Jend
          do i=Istr,Iend
            rubar(i,j)=rubar(i,j) - rdrg*ubar(i,j,kstp)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - rdrg*vbar(i,j,kstp)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
      endif
      if (iif.eq.1) then
        if (iic.eq.ntstart) then
          cff3=0.D0
          cff2=0.D0
          cff1=1.D0
        elseif (iic.eq.ntstart+1) then
          cff3=0.D0
          cff2=-0.5D0
          cff1=1.5D0
        else
          cff3=0.281105D0
          cff2=-0.5D0-2.D0*cff3
          cff1=1.5D0+cff3
        endif
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=rufrc(i,j)-rubar(i,j)
            rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                          + cff3*rufrc_bak(i,j,nstp)
            rufrc_bak(i,j,nstp)=cff
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=rvfrc(i,j)-rvbar(i,j)
            rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                          + cff3*rvfrc_bak(i,j,nstp)
            rvfrc_bak(i,j,nstp)=cff
          enddo
        enddo
        do j=Jstr-1,Jend
          do i=Istr-1,Iend
            UFx(i,j)=zeta_new(i,j)-zeta(i,j,kstp)
            UFe(i,j)=(1.D0+rhoS(i,j))*UFx(i,j)
            VFe(i,j)=UFe(i,j)*(zeta_new(i,j)+zeta(i,j,kstp))
            VFx(i,j)=UFx(i,j)*(rhoS(i,j)-rhoA(i,j))
          enddo
        enddo
        cff=0.5D0*g
        do j=Jstr,Jend
          do i=Istr,Iend
            rubar(i,j)=rubar(i,j) +cff*on_u(i,j)*( (h(i-1,j)+h(i,j))
     &          *(UFe(i-1,j)-UFe(i,j)) +VFe(i-1,j)-VFe(i,j)
     &              +(h(i-1,j)-h(i,j))*( VFx(i-1,j)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(UFx(i-1,j)-UFx(i,j)) )
     &                                                              )
            rvbar(i,j)=rvbar(i,j) +cff*om_v(i,j)*( (h(i,j-1)+h(i,j))
     &          *(UFe(i,j-1)-UFe(i,j)) +VFe(i,j-1)-VFe(i,j)
     &              +(h(i,j-1)-h(i,j))*( VFx(i,j-1)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(UFx(i,j-1)-UFx(i,j)) )
     &                                                              )
          enddo
        enddo
      endif
      do j=Jstr-1,Jend
        do i=Istr-1,Iend
          DUon(i,j)=zeta(i,j,kstp)+h(i,j)
        enddo
      enddo
      cff=0.5D0*dtfast
      cff1=0.5D0*weight(1,iif)
      do j=Jstr,Jend
        do i=Istr,Iend
          DUnew=( (DUon(i,j)+DUon(i-1,j))*ubar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
     &                             *(rubar(i,j)+rufrc(i,j))
     &                                                    )
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*on_u(i,j)*( DUnew
     &                                                   )
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          DVnew=( (DUon(i,j)+DUon(i,j-1))*vbar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
     &                             *(rvbar(i,j)+rvfrc(i,j))
     &                                                    )
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*om_v(i,j)*(DVnew
     &                                                   )
        enddo
      enddo
      call u2dbc_tile (Istr,Iend,Jstr,Jend, UFx)
      call v2dbc_tile (Istr,Iend,Jstr,Jend, UFx)
      cff1=0.5D0*weight(1,iif)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                   zeta(-2,-2,knew))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   ubar(-2,-2,knew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vbar(-2,-2,knew))
      VMAXL=100.D0
      VMAX=0.D0
      do j=Jstr,Jend
        do i=Istr,Iend
          cff1=ubar(i,j,knew)
          cff2=vbar(i,j,knew)
          cff=max(abs(cff1),abs(cff2))
          IF (cff.GE.VMAX .or. cff1.ne.cff1 .or. cff2.ne.cff2) THEN
            IF (cff.GE.VMAX .and. cff1.eq.cff1 .and. cff2.eq.cff2) THEN
              VMAX=cff
            ELSE
              VMAX=666.D0
            ENDIF
            imax=i+iminmpi-1
            jmax=j+jminmpi-1
          ENDIF
        enddo
      enddo
      IF (VMAX.GT.VMAXL) THEN
        write(stdout,'(9(A/))')
     &     '                                         ',
     &     '                                         ',
     &     ' ======================================= ',
     &     ' =                                     = ',
     &     ' =   STEP2D:   ABNORMAL JOB END        = ',
     &     ' =                 BLOW UP             = ',
     &     ' =                                     = ',
     &     ' ======================================= ',
     &     '                                         '
        if (VMAX.eq.666.D0) then
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =   NaN'
        else
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =',VMAX
        endif
        write(stdout,'(A,2I6)')
     &                                       '  IMAX JMAX  =',imax,jmax
        write(stdout,'(A,2I6/)')
     &                                         '  IINT IEXT  =',iic,iif
        may_day_flag=1
        call mpi_abort (MPI_COMM_WORLD, err)
      ENDIF
      return
      end
