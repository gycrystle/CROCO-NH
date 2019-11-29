      subroutine MessPass3D_tile (Istr,Iend,Jstr,Jend, A, nmax)
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      include 'mpif.h'
      integer*4 Npts,ipts,jpts
      parameter (Npts=2)
      integer*4 nmax
      real A(-2:Lm+3+padd_X,-2:Mm+3+padd_E,nmax)
      integer*4 Istr,Iend,Jstr,Jend, i,j,k, isize,jsize,ksize,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr
      integer*4 iter, mdii, mdjj
      integer*4 sub_X,size_X, sub_E,size_E, size_Z
      parameter (size_Z=Npts*Npts*(N+1),
     &     sub_X=(Lm+NSUB_X-1)/NSUB_X,  size_X=Npts*(N+1)*(4+sub_X),
     &     sub_E=(Mm+NSUB_E-1)/NSUB_E,  size_E=Npts*(N+1)*(4+sub_E))
      real buf_snd4(size_Z),  ibuf_sndN(size_X),  buf_snd2(size_Z),
     &     buf_rev4(size_Z),  ibuf_revN(size_X),  buf_rev2(size_Z),
     &    jbuf_sndW(size_E),                      jbuf_sndE(size_E),
     &    jbuf_revW(size_E),                      jbuf_revE(size_E),
     &     buf_snd1(size_Z),  ibuf_sndS(size_X),  buf_snd3(size_Z),
     &     buf_rev1(size_Z),  ibuf_revS(size_X),  buf_rev3(size_Z)
      integer*4 imin,imax,ishft, jmin,jmax,jshft
      if (NP_XI.eq.1) then
        imin=Istr-2
        imax=Iend+2
      else
        imin=Istr
        imax=Iend
      endif
      ishft=imax-imin+1
      if (NP_ETA.eq.1) then
        jmin=Jstr-2
        jmax=Jend+2
      else
        jmin=Jstr
        jmax=Jend
      endif
      jshft=jmax-jmin+1
      ksize=Npts*Npts*nmax
      isize=Npts*ishft*nmax
      jsize=Npts*jshft*nmax
      do iter=0,1
        mdii=mod(ii+iter,2)
        mdjj=mod(jj+iter,2)
        if (mdii.eq.0) then
          if ((WEST_INTER)) then
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndW(k+nmax*(j-jmin+(ipts-1)*jshft))=A(ipts,j,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (jbuf_revW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 1, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (EAST_INTER) then
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndE(k+nmax*(j-jmin+(ipts-1)*jshft))=
     &                                       A(Lmmpi-Npts+ipts,j,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (jbuf_revE, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 2, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdjj.eq.0) then
          if (SOUTH_INTER) then
            ibuf_snds = 0.D0
            do k=1,nmax
              do i=imin,imax
                do jpts=1,Npts
                  ibuf_sndS(k+nmax*(i-imin+(jpts-1)*ishft))=A(i,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (ibuf_revS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER) then
            ibuf_sndn = 0.D0
            do k=1,nmax
              do i=imin,imax
                do jpts=1,Npts
                  ibuf_sndN(k+nmax*(i-imin+(jpts-1)*ishft))=
     &                                         A(i,Mmmpi-Npts+jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (ibuf_revN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdii.eq.0) then
          if (SOUTH_INTER .and. WEST_INTER) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd1(k+nmax*(ipts-1+Npts*(jpts-1)))=A(ipts,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev1, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                5, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (NORTH_INTER .and. EAST_INTER) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd2(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                            A(Lmmpi+ipts-Npts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev2, ksize, MPI_DOUBLE_PRECISION, p_NE,
     &                                  5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2, ksize, MPI_DOUBLE_PRECISION, p_NE,
     &                                  6, MPI_COMM_WORLD,        ierr)
          endif
        endif
        if (mdii.eq.1) then
          if (SOUTH_INTER .and. EAST_INTER) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd3(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(Lmmpi+ipts-Npts,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev3, ksize, MPI_DOUBLE_PRECISION, p_SE,
     &                                8, MPI_COMM_WORLD, req(7),ierr)
            call MPI_Send  (buf_snd3, ksize, MPI_DOUBLE_PRECISION, p_SE,
     &                                7, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (NORTH_INTER .and. WEST_INTER) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd4(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(ipts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                                7, MPI_COMM_WORLD, req(8),ierr)
            call MPI_Send  (buf_snd4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                                8, MPI_COMM_WORLD,        ierr)
          endif
        endif
      enddo
      if (WEST_INTER) then
        call MPI_Wait (req(1),status(1,1),ierr)
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
             A(ipts-Npts,j,k)=jbuf_revW(k+nmax*(j-jmin+(ipts-1)*jshft))
            enddo
          enddo
        enddo
      endif
      if (EAST_INTER) then
        call MPI_Wait (req(2),status(1,2),ierr)
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(Lmmpi+ipts,j,k)=
     &                         jbuf_revE(k+nmax*(j-jmin+(ipts-1)*jshft))
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER) then
        call MPI_Wait (req(3),status(1,3),ierr)
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,jpts-Npts,k)=ibuf_revS(k+nmax*(i-imin+(jpts-1)*ishft))
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER) then
        call MPI_Wait (req(4),status(1,4),ierr)
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,Mmmpi+jpts,k)=
     &                         ibuf_revN(k+nmax*(i-imin+(jpts-1)*ishft))
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(5),status(1,5),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,jpts-Npts,k)=
     &                           buf_rev1(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER.and.EAST_INTER) then
        call MPI_Wait (req(6),status(1,6),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,Mmmpi+jpts,k)=
     &                           buf_rev2(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(7),status(1,7),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,jpts-Npts,k)=
     &                           buf_rev3(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(8),status(1,8),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,Mmmpi+jpts,k)=
     &                           buf_rev4(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      return
      end
      subroutine MessPass3D_3pts_tile (Istr,Iend,Jstr,Jend, A, nmax)
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      include 'mpif.h'
      integer*4 Npts,ipts,jpts
      parameter (Npts=3)
      integer*4 nmax
      real A(-2:Lm+3+padd_X,-2:Mm+3+padd_E,nmax)
      integer*4 Istr,Iend,Jstr,Jend, i,j,k, isize,jsize,ksize,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr
      integer*4 iter, mdii, mdjj
      integer*4 sub_X,size_X, sub_E,size_E, size_Z
      parameter (size_Z=Npts*Npts*(N+1),
     &     sub_X=(Lm+NSUB_X-1)/NSUB_X,  size_X=Npts*(N+1)*(4+sub_X),
     &     sub_E=(Mm+NSUB_E-1)/NSUB_E,  size_E=Npts*(N+1)*(4+sub_E))
      real buf_snd4(size_Z),  ibuf_sndN(size_X),  buf_snd2(size_Z),
     &     buf_rev4(size_Z),  ibuf_revN(size_X),  buf_rev2(size_Z),
     &    jbuf_sndW(size_E),                      jbuf_sndE(size_E),
     &    jbuf_revW(size_E),                      jbuf_revE(size_E),
     &     buf_snd1(size_Z),  ibuf_sndS(size_X),  buf_snd3(size_Z),
     &     buf_rev1(size_Z),  ibuf_revS(size_X),  buf_rev3(size_Z)
      integer*4 imin,imax,ishft, jmin,jmax,jshft
      if (NP_XI.eq.1) then
        imin=Istr-2
        imax=Iend+2
      else
        imin=Istr
        imax=Iend
      endif
      ishft=imax-imin+1
      if (NP_ETA.eq.1) then
        jmin=Jstr-2
        jmax=Jend+2
      else
        jmin=Jstr
        jmax=Jend
      endif
      jshft=jmax-jmin+1
      ksize=Npts*Npts*nmax
      isize=Npts*ishft*nmax
      jsize=Npts*jshft*nmax
      do iter=0,1
        mdii=mod(ii+iter,2)
        mdjj=mod(jj+iter,2)
        if (mdii.eq.0) then
          if ((WEST_INTER)) then
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndW(k+nmax*(j-jmin+(ipts-1)*jshft))=A(ipts,j,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (jbuf_revW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 1, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (EAST_INTER) then
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndE(k+nmax*(j-jmin+(ipts-1)*jshft))=
     &                                       A(Lmmpi-Npts+ipts,j,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (jbuf_revE, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 2, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdjj.eq.0) then
          if (SOUTH_INTER) then
            ibuf_snds = 0.D0
            do k=1,nmax
              do i=imin,imax
                do jpts=1,Npts
                  ibuf_sndS(k+nmax*(i-imin+(jpts-1)*ishft))=A(i,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (ibuf_revS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER) then
            ibuf_sndn = 0.D0
            do k=1,nmax
              do i=imin,imax
                do jpts=1,Npts
                  ibuf_sndN(k+nmax*(i-imin+(jpts-1)*ishft))=
     &                                         A(i,Mmmpi-Npts+jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (ibuf_revN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdii.eq.0) then
          if (SOUTH_INTER .and. WEST_INTER) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd1(k+nmax*(ipts-1+Npts*(jpts-1)))=A(ipts,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev1, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                5, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (NORTH_INTER .and. EAST_INTER) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd2(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                            A(Lmmpi+ipts-Npts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev2, ksize, MPI_DOUBLE_PRECISION, p_NE,
     &                                  5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2, ksize, MPI_DOUBLE_PRECISION, p_NE,
     &                                  6, MPI_COMM_WORLD,        ierr)
          endif
        endif
        if (mdii.eq.1) then
          if (SOUTH_INTER .and. EAST_INTER) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd3(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(Lmmpi+ipts-Npts,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev3, ksize, MPI_DOUBLE_PRECISION, p_SE,
     &                                8, MPI_COMM_WORLD, req(7),ierr)
            call MPI_Send  (buf_snd3, ksize, MPI_DOUBLE_PRECISION, p_SE,
     &                                7, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (NORTH_INTER .and. WEST_INTER) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd4(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(ipts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                                7, MPI_COMM_WORLD, req(8),ierr)
            call MPI_Send  (buf_snd4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                                8, MPI_COMM_WORLD,        ierr)
          endif
        endif
      enddo
      if (WEST_INTER) then
        call MPI_Wait (req(1),status(1,1),ierr)
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
             A(ipts-Npts,j,k)=jbuf_revW(k+nmax*(j-jmin+(ipts-1)*jshft))
            enddo
          enddo
        enddo
      endif
      if (EAST_INTER) then
        call MPI_Wait (req(2),status(1,2),ierr)
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(Lmmpi+ipts,j,k)=
     &                         jbuf_revE(k+nmax*(j-jmin+(ipts-1)*jshft))
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER) then
        call MPI_Wait (req(3),status(1,3),ierr)
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,jpts-Npts,k)=ibuf_revS(k+nmax*(i-imin+(jpts-1)*ishft))
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER) then
        call MPI_Wait (req(4),status(1,4),ierr)
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,Mmmpi+jpts,k)=
     &                         ibuf_revN(k+nmax*(i-imin+(jpts-1)*ishft))
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(5),status(1,5),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,jpts-Npts,k)=
     &                           buf_rev1(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER.and.EAST_INTER) then
        call MPI_Wait (req(6),status(1,6),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,Mmmpi+jpts,k)=
     &                           buf_rev2(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(7),status(1,7),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,jpts-Npts,k)=
     &                           buf_rev3(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(8),status(1,8),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,Mmmpi+jpts,k)=
     &                           buf_rev4(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      return
      end
