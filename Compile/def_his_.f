      subroutine def_his (ncid, total_rec, ierr)
      implicit none
      logical create_new_file
      integer*4 ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  u2dgrd(3), v2dgrd(3), auxil(2), checkdims
     &      , b3dgrd(4)
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4), w3dgrd(4), itrc
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
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diabio, filetype_diabio_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diabio=7,filetype_diabio_avg=8)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxU, indxV, indxT
      parameter (indxU=6, indxV=7, indxT=8)
      integer*4 indxS
      parameter (indxS=indxT+1)
      integer*4 indxBSD, indxBSS
      parameter (indxBSD=indxT+ntrc_salt+ntrc_pas+ntrc_bio+1,
     &           indxBSS=101)
      integer*4 indxO, indxW, indxPnh, indxR, indxVisc, indxDiff,
     &        indxAkv, indxAkt
      parameter (indxO=indxT+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed
     &                      +ntrc_diats+ntrc_diauv+ntrc_diabio+1,
     &           indxW=indxO+1,
     &           indxPnh=indxO+2, indxR=indxO+3, indxVisc=indxO+4,
     &           indxDiff=indxO+5,indxAkv=indxO+6, indxAkt=indxO+7)
      integer*4 indxAks
      parameter (indxAks=indxAkt+4)
      integer*4 indxHbl
      parameter (indxHbl=indxAkt+5)
      integer*4 indxAkk
      parameter (indxAkk=indxAkt+7)
      integer*4 indxAkp
      parameter (indxAkp=indxAkt+8)
      integer*4 indxTke
      parameter (indxTke=indxAkt+9)
      integer*4 indxGls
      parameter (indxGls=indxAkt+10)
      integer*4 indxLsc
      parameter (indxLsc=indxAkt+11)
      integer*4 indxSSH
      parameter (indxSSH=indxAkt+12)
      integer*4 indxbvf
      parameter (indxbvf=indxSSH+1)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+2, indxSVSTR=indxSSH+3)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+4)
      integer*4 indxShflx, indxShflx_rsw
      parameter (indxShflx=indxSSH+5)
      integer*4 indxSwflx
      parameter (indxSwflx=indxShflx+1, indxShflx_rsw=indxShflx+2)
      integer*4 indxSST, indxdQdSST
      parameter (indxSST=indxShflx_rsw+1, indxdQdSST=indxShflx_rsw+2)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+21)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+22)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+23)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+24)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB
      parameter (indxWWA=indxSUSTR+32, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      parameter (xi_rho=LLm+2,  xi_u=xi_rho-1,
     &           eta_rho=MMm+2, eta_v=eta_rho-1)
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms ,
     &        ntsrf,  ntssh,  ntsst, ntsss, ntuclm,
     &        ntbulk, ncidqbar, ntqbar, ntww
      integer*4 nttclm(NT), ntstf(NT), nttsrc(NT)
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
      integer*4 rstT(NT)
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisU,   hisV,   hisR,    hisHbl, hisHbbl
     &      , hisO,   hisW,   hisVisc, hisDiff
     &      , hisPnh
     &      , hisAkv, hisAkt, hisAks
     &      , hisbvf
     &      , hisAkk, hisAkp, hisTke, hisGls, hisLsc
      integer*4 hisT(NT)
      logical wrthis(500+NT)
      common/incscrum/
     &        ncidfrc, ncidbulk,ncidclm, ntsms, ntsrf, ntssh, ntsst
     &      , ntuclm, ntsss, ntbulk, ncidqbar, ntqbar, ntww
     &                        ,  nttclm, ntstf, nttsrc
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV,   rstT
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisU,    hisV,     hisT,    hisR
     &      , hisO,    hisW,     hisVisc, hisDiff
     &      , hisPnh
     &      , hisAkv,  hisAkt,   hisAks
     &      , hisHbl,  hisHbbl
     &      , hisbvf
     &      , hisAkk, hisAkp, hisTke, hisGls, hisLsc
     &      , wrthis
      character*80 date_str, title, start_date
      character*80 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
      character*75  vname(20, 500)
      common /cncscrum/       date_str,   title,  start_date
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &                      ,  vname
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      integer*4 nf_diskless
      integer*4 nf_mmap
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      integer*4 nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4 nf_inq_path
      external nf_inq_path
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_inq_format
      external        nf_inq_format
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nf_ubyte
      integer*4 nf_ushort
      integer*4 nf_uint
      integer*4 nf_int64
      integer*4 nf_uint64
      integer*4 nf_string
      integer*4 nf_vlen
      integer*4 nf_opaque
      integer*4 nf_enum
      integer*4 nf_compound
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)
      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)
      integer*4           nf_fill_ubyte
      integer*4           nf_fill_ushort
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)
      integer*4 nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)
      integer*4 nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)
      integer*4 nf_netcdf4
      parameter (nf_netcdf4 = 4096)
      integer*4 nf_classic_model
      parameter (nf_classic_model = 256)
      integer*4 nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer*4 nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer*4 nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)
      integer*4 nf_endian_native
      parameter (nf_endian_native = 0)
      integer*4 nf_endian_little
      parameter (nf_endian_little = 1)
      integer*4 nf_endian_big
      parameter (nf_endian_big = 2)
      integer*4 nf_chunked
      parameter (nf_chunked = 0)
      integer*4 nf_contiguous
      parameter (nf_contiguous = 1)
      integer*4 nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer*4 nf_fletcher32
      parameter (nf_fletcher32 = 1)
      integer*4 nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer*4 nf_shuffle
      parameter (nf_shuffle = 1)
      integer*4 nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer*4 nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)
      integer*4 nf_mpiio
      parameter (nf_mpiio = 8192)
      integer*4 nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer*4 nf_pnetcdf
      parameter (nf_pnetcdf = 32768)
      integer*4 nf_independent
      parameter (nf_independent = 0)
      integer*4 nf_collective
      parameter (nf_collective = 1)
      integer*4 nf_ehdferr
      parameter (nf_ehdferr = -101)
      integer*4 nf_ecantread
      parameter (nf_ecantread = -102)
      integer*4 nf_ecantwrite
      parameter (nf_ecantwrite = -103)
      integer*4 nf_ecantcreate
      parameter (nf_ecantcreate = -104)
      integer*4 nf_efilemeta
      parameter (nf_efilemeta = -105)
      integer*4 nf_edimmeta
      parameter (nf_edimmeta = -106)
      integer*4 nf_eattmeta
      parameter (nf_eattmeta = -107)
      integer*4 nf_evarmeta
      parameter (nf_evarmeta = -108)
      integer*4 nf_enocompound
      parameter (nf_enocompound = -109)
      integer*4 nf_eattexists
      parameter (nf_eattexists = -110)
      integer*4 nf_enotnc4
      parameter (nf_enotnc4 = -111)
      integer*4 nf_estrictnc3
      parameter (nf_estrictnc3 = -112)
      integer*4 nf_enotnc3
      parameter (nf_enotnc3 = -113)
      integer*4 nf_enopar
      parameter (nf_enopar = -114)
      integer*4 nf_eparinit
      parameter (nf_eparinit = -115)
      integer*4 nf_ebadgrpid
      parameter (nf_ebadgrpid = -116)
      integer*4 nf_ebadtypid
      parameter (nf_ebadtypid = -117)
      integer*4 nf_etypdefined
      parameter (nf_etypdefined = -118)
      integer*4 nf_ebadfield
      parameter (nf_ebadfield = -119)
      integer*4 nf_ebadclass
      parameter (nf_ebadclass = -120)
      integer*4 nf_emaptype
      parameter (nf_emaptype = -121)
      integer*4 nf_elatefill
      parameter (nf_elatefill = -122)
      integer*4 nf_elatedef
      parameter (nf_elatedef = -123)
      integer*4 nf_edimscale
      parameter (nf_edimscale = -124)
      integer*4 nf_enogrp
      parameter (nf_enogrp = -125)
      integer*4 nf_create_par
      external nf_create_par
      integer*4 nf_open_par
      external nf_open_par
      integer*4 nf_var_par_access
      external nf_var_par_access
      integer*4 nf_inq_ncid
      external nf_inq_ncid
      integer*4 nf_inq_grps
      external nf_inq_grps
      integer*4 nf_inq_grpname
      external nf_inq_grpname
      integer*4 nf_inq_grpname_full
      external nf_inq_grpname_full
      integer*4 nf_inq_grpname_len
      external nf_inq_grpname_len
      integer*4 nf_inq_grp_parent
      external nf_inq_grp_parent
      integer*4 nf_inq_grp_ncid
      external nf_inq_grp_ncid
      integer*4 nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid
      integer*4 nf_inq_varids
      external nf_inq_varids
      integer*4 nf_inq_dimids
      external nf_inq_dimids
      integer*4 nf_def_grp
      external nf_def_grp
      integer*4 nf_rename_grp
      external nf_rename_grp
      integer*4 nf_def_var_deflate
      external nf_def_var_deflate
      integer*4 nf_inq_var_deflate
      external nf_inq_var_deflate
      integer*4 nf_def_var_fletcher32
      external nf_def_var_fletcher32
      integer*4 nf_inq_var_fletcher32
      external nf_inq_var_fletcher32
      integer*4 nf_def_var_chunking
      external nf_def_var_chunking
      integer*4 nf_inq_var_chunking
      external nf_inq_var_chunking
      integer*4 nf_def_var_fill
      external nf_def_var_fill
      integer*4 nf_inq_var_fill
      external nf_inq_var_fill
      integer*4 nf_def_var_endian
      external nf_def_var_endian
      integer*4 nf_inq_var_endian
      external nf_inq_var_endian
      integer*4 nf_inq_typeids
      external nf_inq_typeids
      integer*4 nf_inq_typeid
      external nf_inq_typeid
      integer*4 nf_inq_type
      external nf_inq_type
      integer*4 nf_inq_user_type
      external nf_inq_user_type
      integer*4 nf_def_compound
      external nf_def_compound
      integer*4 nf_insert_compound
      external nf_insert_compound
      integer*4 nf_insert_array_compound
      external nf_insert_array_compound
      integer*4 nf_inq_compound
      external nf_inq_compound
      integer*4 nf_inq_compound_name
      external nf_inq_compound_name
      integer*4 nf_inq_compound_size
      external nf_inq_compound_size
      integer*4 nf_inq_compound_nfields
      external nf_inq_compound_nfields
      integer*4 nf_inq_compound_field
      external nf_inq_compound_field
      integer*4 nf_inq_compound_fieldname
      external nf_inq_compound_fieldname
      integer*4 nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex
      integer*4 nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset
      integer*4 nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype
      integer*4 nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims
      integer*4 nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes
      integer*4 nf_def_vlen
      external nf_def_vlen
      integer*4 nf_inq_vlen
      external nf_inq_vlen
      integer*4 nf_free_vlen
      external nf_free_vlen
      integer*4 nf_def_enum
      external nf_def_enum
      integer*4 nf_insert_enum
      external nf_insert_enum
      integer*4 nf_inq_enum
      external nf_inq_enum
      integer*4 nf_inq_enum_member
      external nf_inq_enum_member
      integer*4 nf_inq_enum_ident
      external nf_inq_enum_ident
      integer*4 nf_def_opaque
      external nf_def_opaque
      integer*4 nf_inq_opaque
      external nf_inq_opaque
      integer*4 nf_put_att
      external nf_put_att
      integer*4 nf_get_att
      external nf_get_att
      integer*4 nf_put_var
      external nf_put_var
      integer*4 nf_put_var1
      external nf_put_var1
      integer*4 nf_put_vara
      external nf_put_vara
      integer*4 nf_put_vars
      external nf_put_vars
      integer*4 nf_get_var
      external nf_get_var
      integer*4 nf_get_var1
      external nf_get_var1
      integer*4 nf_get_vara
      external nf_get_vara
      integer*4 nf_get_vars
      external nf_get_vars
      integer*4 nf_put_var1_int64
      external nf_put_var1_int64
      integer*4 nf_put_vara_int64
      external nf_put_vara_int64
      integer*4 nf_put_vars_int64
      external nf_put_vars_int64
      integer*4 nf_put_varm_int64
      external nf_put_varm_int64
      integer*4 nf_put_var_int64
      external nf_put_var_int64
      integer*4 nf_get_var1_int64
      external nf_get_var1_int64
      integer*4 nf_get_vara_int64
      external nf_get_vara_int64
      integer*4 nf_get_vars_int64
      external nf_get_vars_int64
      integer*4 nf_get_varm_int64
      external nf_get_varm_int64
      integer*4 nf_get_var_int64
      external nf_get_var_int64
      integer*4 nf_get_vlen_element
      external nf_get_vlen_element
      integer*4 nf_put_vlen_element
      external nf_put_vlen_element
      integer*4 nf_set_chunk_cache
      external nf_set_chunk_cache
      integer*4 nf_get_chunk_cache
      external nf_get_chunk_cache
      integer*4 nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer*4 nf_get_var_chunk_cache
      external nf_get_var_chunk_cache
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      character*70 text
      ierr=0
      lstr=lenstr(hisname)
      if (nrpfhis.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfhis))
        call insert_time_index (hisname, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefhis
      if (ncid.ne.-1) create_new_file=.false.
      if (mynode.gt.0) create_new_file=.false.
  10  if (create_new_file) then
        ierr=nf_create(hisname(1:lstr),nf_64bit_offset, ncid)
        if (ierr .ne. nf_noerr) then
          write(stdout,'(/3(1x,A)/)') 'ERROR in def_his/avg:',
     &           'Cannot create netCDF file:', hisname(1:lstr)
          goto 99
        endif
        if (nrpfhis.eq.0) total_rec=0
        call put_global_atts (ncid, ierr)
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        b3dgrd(1)=r2dgrd(1)
        b3dgrd(2)=r2dgrd(2)
        b3dgrd(4)=timedim
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2)
        r3dgrd(4)=timedim
        u3dgrd(1)=u2dgrd(1)
        u3dgrd(2)=r2dgrd(2)
        u3dgrd(3)=r3dgrd(3)
        u3dgrd(4)=timedim
        v3dgrd(1)=r2dgrd(1)
        v3dgrd(2)=v2dgrd(2)
        v3dgrd(3)=r3dgrd(3)
        v3dgrd(4)=timedim
        w3dgrd(1)=r2dgrd(1)
        w3dgrd(2)=r2dgrd(2)
        w3dgrd(4)=timedim
      if (total_rec.le.1) then
         call def_grid_3d(ncid, r2dgrd, u2dgrd, v2dgrd
     &                   ,r3dgrd, w3dgrd)
      endif
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 hisTstep)
        ierr=nf_put_att_text (ncid, hisTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_REAL, 1, timedim, hisTime)
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, hisTime, 'long_name', lvar,
     &                                vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, hisTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, hisTime, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))
        call nf_add_attribute(ncid, hisTime, indxTime, 5,
     &       NF_REAL, ierr)
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                            NF_REAL, 1, timedim, hisTime2)
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, hisTime2, 'long_name', lvar,
     &                                vname(2,indxTime2)(1:lvar))
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, hisTime2, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, hisTime2, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))
        call nf_add_attribute(ncid, hisTime2, indxTime2, 5,
     &       NF_REAL, ierr)
        if (wrthis(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_def_var (ncid, vname(1,indxZ)(1:lvar),
     &                              NF_REAL, 3, r2dgrd, hisZ)
          lvar=lenstr(vname(2,indxZ))
          ierr=nf_put_att_text (ncid, hisZ, 'long_name', lvar,
     &                                  vname(2,indxZ)(1:lvar))
          lvar=lenstr(vname(3,indxZ))
          ierr=nf_put_att_text (ncid, hisZ, 'units',     lvar,
     &                                  vname(3,indxZ)(1:lvar))
          lvar=lenstr(vname(4,indxZ))
          ierr=nf_put_att_text (ncid, hisZ, 'field',     lvar,
     &                                  vname(4,indxZ)(1:lvar))
          call nf_add_attribute(ncid, hisZ, indxZ, 5, NF_REAL, ierr)
       endif
        if (wrthis(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_def_var (ncid, vname(1,indxUb)(1:lvar),
     &                             NF_REAL, 3, u2dgrd, hisUb)
          lvar=lenstr(vname(2,indxUb))
          ierr=nf_put_att_text (ncid, hisUb, 'long_name', lvar,
     &                                  vname(2,indxUb)(1:lvar))
          lvar=lenstr(vname(3,indxUb))
          ierr=nf_put_att_text (ncid, hisUb, 'units',     lvar,
     &                                  vname(3,indxUb)(1:lvar))
          lvar=lenstr(vname(4,indxUb))
          ierr=nf_put_att_text (ncid, hisUb, 'field',    lvar,
     &                                  vname(4,indxUb)(1:lvar))
          call nf_add_attribute(ncid, hisUb, indxUb, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_def_var (ncid, vname(1,indxVb)(1:lvar),
     &                              NF_REAL, 3, v2dgrd, hisVb)
          lvar=lenstr(vname(2,indxVb))
          ierr=nf_put_att_text (ncid, hisVb, 'long_name', lvar,
     &                                  vname(2,indxVb)(1:lvar))
          lvar=lenstr(vname(3,indxVb))
          ierr=nf_put_att_text (ncid, hisVb, 'units',     lvar,
     &                                  vname(3,indxVb)(1:lvar))
          lvar=lenstr(vname(4,indxVb))
          ierr=nf_put_att_text (ncid, hisVb, 'field',     lvar,
     &                                  vname(4,indxVb)(1:lvar))
          call nf_add_attribute(ncid, hisVb, indxVb, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_def_var (ncid, vname(1,indxU)(1:lvar),
     &                             NF_REAL, 4, u3dgrd, hisU)
          lvar=lenstr(vname(2,indxU))
          ierr=nf_put_att_text (ncid, hisU, 'long_name', lvar,
     &                                  vname(2,indxU)(1:lvar))
          lvar=lenstr(vname(3,indxU))
          ierr=nf_put_att_text (ncid, hisU, 'units',     lvar,
     &                                  vname(3,indxU)(1:lvar))
          lvar=lenstr(vname(4,indxU))
          ierr=nf_put_att_text (ncid, hisU, 'field',     lvar,
     &                                  vname(4,indxU)(1:lvar))
          call nf_add_attribute(ncid, hisU, indxU, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_def_var (ncid, vname(1,indxV)(1:lvar),
     &                             NF_REAL, 4, v3dgrd, hisV)
          lvar=lenstr(vname(2,indxV))
          ierr=nf_put_att_text (ncid, hisV, 'long_name', lvar,
     &                                  vname(2,indxV)(1:lvar))
          lvar=lenstr(vname(3,indxV))
          ierr=nf_put_att_text (ncid, hisV, 'units',     lvar,
     &                                  vname(3,indxV)(1:lvar))
          lvar=lenstr(vname(4,indxV))
          ierr=nf_put_att_text (ncid, hisV, 'field',     lvar,
     &                                  vname(4,indxV)(1:lvar))
          call nf_add_attribute(ncid, hisV, indxV, 5, NF_REAL, ierr)
        endif
        do itrc=1,NT
          if (wrthis(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_def_var (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, hisT(itrc))
            lvar=lenstr(vname(2,indxT+itrc-1))
            ierr=nf_put_att_text (ncid, hisT(itrc), 'long_name',
     &                         lvar, vname(2,indxT+itrc-1)(1:lvar))
            lvar=lenstr(vname(3,indxT+itrc-1))
            ierr=nf_put_att_text (ncid, hisT(itrc), 'units', lvar,
     &                               vname(3,indxT+itrc-1)(1:lvar))
            lvar=lenstr(vname(4,indxT+itrc-1))
            ierr=nf_put_att_text (ncid, hisT(itrc), 'field', lvar,
     &                               vname(4,indxT+itrc-1)(1:lvar))
            call nf_add_attribute(ncid,hisT(itrc),indxT+itrc-1,5,
     &           NF_REAL, ierr)
            endif
      enddo
        if (wrthis(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_def_var (ncid, vname(1,indxR)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, hisR)
          lvar=lenstr(vname(2,indxR))
          ierr=nf_put_att_text (ncid, hisR, 'long_name', lvar,
     &                                  vname(2,indxR)(1:lvar))
          lvar=lenstr(vname(3,indxR))
          ierr=nf_put_att_text (ncid, hisR, 'units',     lvar,
     &                                  vname(3,indxR)(1:lvar))
          lvar=lenstr(vname(4,indxR))
          ierr=nf_put_att_text (ncid, hisR, 'field',     lvar,
     &                                  vname(4,indxR)(1:lvar))
          call nf_add_attribute(ncid, hisR, indxR, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxbvf)) then
          lvar=lenstr(vname(1,indxbvf))
          ierr=nf_def_var (ncid, vname(1,indxbvf)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, hisbvf)
          lvar=lenstr(vname(2,indxbvf))
          ierr=nf_put_att_text (ncid, hisbvf, 'long_name', lvar,
     &                                  vname(2,indxbvf)(1:lvar))
          lvar=lenstr(vname(3,indxbvf))
          ierr=nf_put_att_text (ncid, hisbvf, 'units',     lvar,
     &                                  vname(3,indxbvf)(1:lvar))
          lvar=lenstr(vname(4,indxbvf))
          ierr=nf_put_att_text (ncid, hisR, 'field',     lvar,
     &                                  vname(4,indxbvf)(1:lvar))
          call nf_add_attribute(ncid, hisbvf, indxbvf, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_def_var (ncid, vname(1,indxO)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, hisO)
          lvar=lenstr(vname(2,indxO))
          ierr=nf_put_att_text (ncid, hisO, 'long_name', lvar,
     &                                  vname(2,indxO)(1:lvar))
          lvar=lenstr(vname(3,indxO))
          ierr=nf_put_att_text (ncid, hisO, 'units',     lvar,
     &                                  vname(3,indxO)(1:lvar))
          lvar=lenstr(vname(4,indxO))
          ierr=nf_put_att_text (ncid, hisO, 'field',     lvar,
     &                                  vname(4,indxO)(1:lvar))
          call nf_add_attribute(ncid, hisO, indxO, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_def_var (ncid, vname(1,indxW)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, hisW)
          lvar=lenstr(vname(2,indxW))
          ierr=nf_put_att_text (ncid, hisW, 'long_name', lvar,
     &                                  vname(2,indxW)(1:lvar))
          lvar=lenstr(vname(3,indxW))
          ierr=nf_put_att_text (ncid, hisW, 'units',     lvar,
     &                                  vname(3,indxW)(1:lvar))
          lvar=lenstr(vname(4,indxW))
          ierr=nf_put_att_text (ncid, hisW, 'field',     lvar,
     &                                  vname(4,indxW)(1:lvar))
          call nf_add_attribute(ncid, hisW, indxW, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxPnh)) then
          lvar=lenstr(vname(1,indxPnh))
          ierr=nf_def_var (ncid, vname(1,indxPnh)(1:lvar),
     &                          NF_REAL, 4, r3dgrd, hisPnh)
          lvar=lenstr(vname(2,indxPnh))
          ierr=nf_put_att_text (ncid, hisPnh, 'long_name', lvar,
     &                                  vname(2,indxPnh)(1:lvar))
          lvar=lenstr(vname(3,indxPnh))
          ierr=nf_put_att_text (ncid, hisPnh, 'units',     lvar,
     &                                  vname(3,indxPnh)(1:lvar))
          lvar=lenstr(vname(4,indxPnh))
          ierr=nf_put_att_text (ncid, hisPnh, 'field',     lvar,
     &                                  vname(4,indxPnh)(1:lvar))
          call nf_add_attribute(ncid, hisPnh, indxPnh, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxBostr)) then
          lvar=lenstr(vname(1,indxBostr))
          ierr=nf_def_var (ncid, vname(1,indxBostr)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, hisBostr)
          lvar=lenstr(vname(2,indxBostr))
          ierr=nf_put_att_text (ncid, hisBostr, 'long_name', lvar,
     &                                 vname(2,indxBostr)(1:lvar))
          lvar=lenstr(vname(3,indxBostr))
          ierr=nf_put_att_text (ncid, hisBostr, 'units',     lvar,
     &                                 vname(3,indxBostr)(1:lvar))
          call nf_add_attribute(ncid, hisBostr,indxBostr,5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxWstr)) then
          lvar=lenstr(vname(1,indxWstr))
          ierr=nf_def_var (ncid, vname(1,indxWstr)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, hisWstr)
          lvar=lenstr(vname(2,indxWstr))
          ierr=nf_put_att_text (ncid, hisWstr, 'long_name', lvar,
     &                                 vname(2,indxWstr)(1:lvar))
          lvar=lenstr(vname(3,indxWstr))
          ierr=nf_put_att_text (ncid, hisWstr, 'units',     lvar,
     &                                 vname(3,indxWstr)(1:lvar))
          call nf_add_attribute(ncid, hisWstr, indxWstr, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxUWstr)) then
          lvar=lenstr(vname(1,indxUWstr))
          ierr=nf_def_var (ncid, vname(1,indxUWstr)(1:lvar),
     &                             NF_REAL, 3, u2dgrd, hisUWstr)
          lvar=lenstr(vname(2,indxUWstr))
          ierr=nf_put_att_text (ncid, hisUWstr, 'long_name', lvar,
     &                                 vname(2,indxUWstr)(1:lvar))
          lvar=lenstr(vname(3,indxUWstr))
          ierr=nf_put_att_text (ncid, hisUWstr, 'units',     lvar,
     &                                 vname(3,indxUWstr)(1:lvar))
          call nf_add_attribute(ncid, hisUWstr, indxUWstr, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxVWstr)) then
          lvar=lenstr(vname(1,indxVWstr))
          ierr=nf_def_var (ncid, vname(1,indxVWstr)(1:lvar),
     &                             NF_REAL, 3, v2dgrd, hisVWstr)
          lvar=lenstr(vname(2,indxVWstr))
          ierr=nf_put_att_text (ncid, hisVWstr, 'long_name', lvar,
     &                                 vname(2,indxVWstr)(1:lvar))
          lvar=lenstr(vname(3,indxVWstr))
          ierr=nf_put_att_text (ncid, hisVWstr, 'units',     lvar,
     &                                 vname(3,indxVWstr)(1:lvar))
          call nf_add_attribute(ncid, hisVWstr, indxVWstr, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxDiff)) then
          lvar=lenstr(vname(1,indxDiff))
          ierr=nf_def_var (ncid, vname(1,indxDiff)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, hisDiff)
          lvar=lenstr(vname(2,indxDiff))
          ierr=nf_put_att_text (ncid, hisDiff, 'long_name', lvar,
     &                                  vname(2,indxDiff)(1:lvar))
          lvar=lenstr(vname(4,indxDiff))
          ierr=nf_put_att_text (ncid, hisDiff, 'field',     lvar,
     &                                  vname(4,indxDiff)(1:lvar))
          call nf_add_attribute(ncid, hisDiff, indxDiff, 5,
     &                          NF_REAL, ierr)
       endif
        if (wrthis(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_def_var (ncid, vname(1,indxAkv)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, hisAkv)
          lvar=lenstr(vname(2,indxAkv))
          ierr=nf_put_att_text (ncid, hisAkv, 'long_name', lvar,
     &                                  vname(2,indxAkv)(1:lvar))
          lvar=lenstr(vname(3,indxAkv))
          ierr=nf_put_att_text (ncid, hisAkv, 'units',     lvar,
     &                                  vname(3,indxAkv)(1:lvar))
          lvar=lenstr(vname(4,indxAkv))
          ierr=nf_put_att_text (ncid, hisAkv, 'field',     lvar,
     &                                  vname(4,indxAkv)(1:lvar))
          call nf_add_attribute(ncid, hisAkv, indxAkv, 5, NF_REAL,
     &                                                          ierr)
        endif
        if (wrthis(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_def_var (ncid, vname(1,indxAkt)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisAkt)
          lvar=lenstr(vname(2,indxAkt))
          ierr=nf_put_att_text (ncid, hisAkt, 'long_name', lvar,
     &                                  vname(2,indxAkt)(1:lvar))
          lvar=lenstr(vname(3,indxAkt))
          ierr=nf_put_att_text (ncid, hisAkt, 'units',     lvar,
     &                                  vname(3,indxAkt)(1:lvar))
          lvar=lenstr(vname(4,indxAkt))
          ierr=nf_put_att_text (ncid, hisAkt, 'field',     lvar,
     &                                  vname(4,indxAkt)(1:lvar))
          call nf_add_attribute(ncid, hisAkt, indxAkt, 5, NF_REAL,
     &                                                           ierr)
        endif
        if (wrthis(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_def_var (ncid, vname(1,indxAks)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisAks)
          lvar=lenstr(vname(2,indxAks))
          ierr=nf_put_att_text (ncid, hisAks, 'long_name', lvar,
     &                                  vname(2,indxAks)(1:lvar))
          lvar=lenstr(vname(3,indxAks))
          ierr=nf_put_att_text (ncid, hisAks, 'units',     lvar,
     &                                  vname(3,indxAks)(1:lvar))
          lvar=lenstr(vname(4,indxAks))
          ierr=nf_put_att_text (ncid, hisAks, 'field',     lvar,
     &                                  vname(4,indxAks)(1:lvar))
          call nf_add_attribute(ncid, hisAks, indxAks, 5, NF_REAL,
     &                                                           ierr)
        endif
        if (wrthis(indxAkk)) then
          lvar=lenstr(vname(1,indxAkk))
          ierr=nf_def_var (ncid, vname(1,indxAkk)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisAkk)
          lvar=lenstr(vname(2,indxAkk))
          ierr=nf_put_att_text (ncid, hisAkk, 'long_name', lvar,
     &                                  vname(2,indxAkk)(1:lvar))
          lvar=lenstr(vname(3,indxAkk))
          ierr=nf_put_att_text (ncid, hisAkk, 'units',     lvar,
     &                                  vname(3,indxAkk)(1:lvar))
          lvar=lenstr(vname(4,indxAkk))
          ierr=nf_put_att_text (ncid, hisAkk, 'field',     lvar,
     &                                  vname(4,indxAkk)(1:lvar))
          call nf_add_attribute(ncid, hisAkk, indxAkk, 5, NF_REAL,
     &                                                             ierr)
        endif
        if (wrthis(indxAkp)) then
          lvar=lenstr(vname(1,indxAkp))
          ierr=nf_def_var (ncid, vname(1,indxAkp)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisAkp)
          lvar=lenstr(vname(2,indxAkp))
          ierr=nf_put_att_text (ncid, hisAkp, 'long_name', lvar,
     &                                  vname(2,indxAkp)(1:lvar))
          lvar=lenstr(vname(3,indxAkp))
          ierr=nf_put_att_text (ncid, hisAkp, 'units',     lvar,
     &                                  vname(3,indxAkp)(1:lvar))
          lvar=lenstr(vname(4,indxAkp))
          ierr=nf_put_att_text (ncid, hisAkp, 'field',     lvar,
     &                                  vname(4,indxAkp)(1:lvar))
          call nf_add_attribute(ncid, hisAkp, indxAkp, 5, NF_REAL,
     &                                                             ierr)
        endif
        if (wrthis(indxTke)) then
          lvar=lenstr(vname(1,indxTke))
          ierr=nf_def_var (ncid, vname(1,indxTke)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisTke)
          lvar=lenstr(vname(2,indxTke))
          ierr=nf_put_att_text (ncid, hisTke, 'long_name', lvar,
     &                                  vname(2,indxTke)(1:lvar))
          lvar=lenstr(vname(3,indxTke))
          ierr=nf_put_att_text (ncid, hisTke, 'units',     lvar,
     &                                  vname(3,indxTke)(1:lvar))
          lvar=lenstr(vname(4,indxTke))
          ierr=nf_put_att_text (ncid, hisTke, 'field',     lvar,
     &                                  vname(4,indxTke)(1:lvar))
          call nf_add_attribute(ncid, hisTke, indxTke, 5, NF_REAL,
     &                                                            ierr)
        endif
        if (wrthis(indxGls)) then
          lvar=lenstr(vname(1,indxGls))
          ierr=nf_def_var (ncid, vname(1,indxGls)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisGls)
          lvar=lenstr(vname(2,indxGls))
          ierr=nf_put_att_text (ncid, hisGls, 'long_name', lvar,
     &                                  vname(2,indxGls)(1:lvar))
          lvar=lenstr(vname(3,indxGls))
          ierr=nf_put_att_text (ncid, hisGls, 'units',     lvar,
     &                                  vname(3,indxGls)(1:lvar))
          lvar=lenstr(vname(4,indxGls))
          ierr=nf_put_att_text (ncid, hisGls, 'field',     lvar,
     &                                  vname(4,indxGls)(1:lvar))
          call nf_add_attribute(ncid, hisGls, indxGls, 5, NF_REAL,
     &                                                            ierr)
        endif
        if (wrthis(indxLsc)) then
          lvar=lenstr(vname(1,indxLsc))
          ierr=nf_def_var (ncid, vname(1,indxLsc)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisLsc)
          lvar=lenstr(vname(2,indxLsc))
          ierr=nf_put_att_text (ncid, hisLsc, 'long_name', lvar,
     &                                  vname(2,indxLsc)(1:lvar))
          lvar=lenstr(vname(3,indxLsc))
          ierr=nf_put_att_text (ncid, hisLsc, 'units',     lvar,
     &                                  vname(3,indxLsc)(1:lvar))
          lvar=lenstr(vname(4,indxLsc))
          ierr=nf_put_att_text (ncid, hisLsc, 'field',     lvar,
     &                                  vname(4,indxLsc)(1:lvar))
          call nf_add_attribute(ncid, hisLsc, indxLsc, 5, NF_REAL,
     &                                                             ierr)
        endif
        if (wrthis(indxShflx)) then
          lvar=lenstr(vname(1,indxShflx))
          ierr=nf_def_var (ncid, vname(1,indxShflx)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, hisShflx)
          lvar=lenstr(vname(2,indxShflx))
          ierr=nf_put_att_text (ncid, hisShflx, 'long_name', lvar,
     &                                 vname(2,indxShflx)(1:lvar))
          lvar=lenstr(vname(3,indxShflx))
          ierr=nf_put_att_text (ncid, hisShflx, 'units',     lvar,
     &                                 vname(3,indxShflx)(1:lvar))
          call nf_add_attribute(ncid, hisShflx, indxShflx, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxSwflx)) then
          lvar=lenstr(vname(1,indxSwflx))
          ierr=nf_def_var (ncid, vname(1,indxSwflx)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, hisSwflx)
          lvar=lenstr(vname(2,indxSwflx))
          ierr=nf_put_att_text (ncid, hisSwflx, 'long_name', lvar,
     &                                 vname(2,indxSwflx)(1:lvar))
          lvar=lenstr(vname(3,indxSwflx))
          ierr=nf_put_att_text (ncid, hisSwflx, 'units',     lvar,
     &                                 vname(3,indxSwflx)(1:lvar))
          call nf_add_attribute(ncid, hisSwflx, indxSwflx, 5,
     &                          NF_REAL, ierr)
        endif
      if (wrthis(indxShflx_rsw)) then
        lvar=lenstr(vname(1,indxShflx_rsw))
        ierr=nf_def_var (ncid, vname(1,indxShflx_rsw)(1:lvar),
     &                           NF_REAL, 3, r2dgrd, hisShflx_rsw)
        lvar=lenstr(vname(2,indxShflx_rsw))
        ierr=nf_put_att_text (ncid, hisShflx_rsw, 'long_name', lvar,
     &                               vname(2,indxShflx_rsw)(1:lvar))
        lvar=lenstr(vname(3,indxShflx_rsw))
        ierr=nf_put_att_text (ncid, hisShflx_rsw, 'units',     lvar,
     &                               vname(3,indxShflx_rsw)(1:lvar))
        call nf_add_attribute(ncid, hisShflx_rsw, indxShflx_rsw,5,
     &                                                   NF_REAL, ierr)
      endif
        ierr=nf_enddef(ncid)
        write(stdout,'(6x,4A,1x,A,i4)') 'DEF_HIS/AVG - Created ',
     &                'new netCDF file ''', hisname(1:lstr), '''.'
     &                 ,' mynode =', mynode
      elseif (ncid.eq.-1) then
        ierr=nf_open (hisname(1:lstr), nf_write, ncid)
        if (ierr. eq. nf_noerr) then
          ierr=checkdims (ncid, hisname, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfhis.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfhis))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &            ) 'DEF_HIS/AVG WARNING: Actual number of records',
     &               rec,  'in netCDF file',  '''',  hisname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpfhis.eq.0) then
              total_rec=rec+1
              if (mynode.gt.0) total_rec=total_rec-1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr. ne. nf_noerr) then
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
            write(stdout,'(/1x,4A,2x,A,I4/)') 'DEF_HIS/AVG ERROR: ',
     &                  'Cannot open file ''', hisname(1:lstr), '''.'
     &                   ,' mynode =', mynode
            goto 99
          endif
        endif
        ierr=nf_inq_varid (ncid, 'time_step', hisTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', hisname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),hisTime)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), hisname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),hisTime2)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), hisname(1:lstr)
          goto 99
        endif
        if (wrthis(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), hisZ)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxZ)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), hisUb)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxUb)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), hisVb)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxVb)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxBostr)) then
          lvar=lenstr(vname(1,indxBostr))
          ierr=nf_inq_varid (ncid,vname(1,indxBostr)(1:lvar),
     &                                                   hisBostr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxBostr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxWstr)) then
          lvar=lenstr(vname(1,indxWstr))
          ierr=nf_inq_varid (ncid,vname(1,indxWstr)(1:lvar),
     &                                                   hisWstr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxWstr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxUWstr)) then
          lvar=lenstr(vname(1,indxUWstr))
          ierr=nf_inq_varid (ncid,vname(1,indxUWstr)(1:lvar),
     &                                                   hisUWstr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxUWstr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxVWstr)) then
          lvar=lenstr(vname(1,indxVWstr))
          ierr=nf_inq_varid (ncid,vname(1,indxVWstr)(1:lvar),
     &                                                   hisVWstr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxVWstr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_inq_varid (ncid, vname(1,indxU)(1:lvar), hisU)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxU)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_inq_varid (ncid, vname(1,indxV)(1:lvar), hisV)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxV)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        do itrc=1,NT
          if (wrthis(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_inq_varid (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                                                 hisT(itrc))
            if (ierr .ne. nf_noerr) then
              write(stdout,1) vname(1,indxT+itrc-1)(1:lvar),
     &                                       hisname(1:lstr)
              goto 99
            endif
          endif
        enddo
        if (wrthis(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_inq_varid (ncid, vname(1,indxR)(1:lvar), hisR)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxR)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxbvf)) then
          lvar=lenstr(vname(1,indxbvf))
          ierr=nf_inq_varid (ncid, vname(1,indxbvf)(1:lvar), hisbvf)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxbvf)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_inq_varid (ncid, vname(1,indxO)(1:lvar), hisO)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxO)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_inq_varid (ncid, vname(1,indxW)(1:lvar), hisW)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxW)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxPnh)) then
          lvar=lenstr(vname(1,indxPnh))
          ierr=nf_inq_varid (ncid, vname(1,indxPnh)(1:lvar), hisPnh)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxPnh)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxDiff)) then
          lvar=lenstr(vname(1,indxDiff))
          ierr=nf_inq_varid (ncid, vname(1,indxDiff)(1:lvar), hisDiff)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxDiff)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_inq_varid (ncid, vname(1,indxAkv)(1:lvar), hisAkv)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAkv)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_inq_varid (ncid,vname(1,indxAkt)(1:lvar), hisAkt)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAkt)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_inq_varid (ncid,vname(1,indxAks)(1:lvar), hisAks)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAks)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxAkk)) then
          lvar=lenstr(vname(1,indxAkk))
          ierr=nf_inq_varid (ncid,vname(1,indxAkk)(1:lvar), hisAkk)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAkk)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxAkp)) then
          lvar=lenstr(vname(1,indxAkp))
          ierr=nf_inq_varid (ncid,vname(1,indxAkp)(1:lvar), hisAkp)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAkp)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxTke)) then
          lvar=lenstr(vname(1,indxTke))
          ierr=nf_inq_varid (ncid,vname(1,indxTke)(1:lvar), hisTke)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTke)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxGls)) then
          lvar=lenstr(vname(1,indxGls))
          ierr=nf_inq_varid (ncid,vname(1,indxGls)(1:lvar), hisGls)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxGls)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxLsc)) then
          lvar=lenstr(vname(1,indxLsc))
          ierr=nf_inq_varid (ncid,vname(1,indxLsc)(1:lvar), hisLsc)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxLsc)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxShflx)) then
          lvar=lenstr(vname(1,indxShflx))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx)(1:lvar),
     &                                                   hisShflx)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxShflx)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxSwflx)) then
          lvar=lenstr(vname(1,indxSwflx))
          ierr=nf_inq_varid (ncid,vname(1,indxSwflx)(1:lvar),
     &                                                   hisSwflx)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxSwflx)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxShflx_rsw)) then
         lvar=lenstr(vname(1,indxShflx_rsw))
         ierr=nf_inq_varid (ncid,vname(1,indxShflx_rsw)(1:lvar),
     &                                               hisShflx_rsw)
         if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxShflx_rsw)(1:lvar), 
     &                                                   hisname(1:lstr)
          goto 99
         endif
        endif
      if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'DEF_HIS/AVG -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
      else
        ierr=nf_open (hisname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          if (mynode.eq.0) write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'DEF_HIS/AVG ERROR: ',
     &                'Cannot open file ''', hisname(1:lstr), '''.'
     &                 ,' mynode =', mynode
          goto 99
        endif
      endif
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_HIS/AVG ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr
      endif
   1  format(/1x,'DEF_HIS/AVG ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)
        if (total_rec.le.1) call wrt_grid (ncid, hisname, lstr)
  99  return
      end
