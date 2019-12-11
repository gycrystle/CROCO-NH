!
!----------------------------------------------------------
! Compute vertical advective fluxes 
! using 5th-order WENO scheme
!----------------------------------------------------------
!
        do k=3,N-3
          do i=IstrU,Iend
            vel=0.5*(We(i,j,k)+We(i-1,j,k))
            FC(i,k)=vel*
     &            FLUX5(
     &           u(i,j,k-2,nrhs), u(i,j,k-1,nrhs), 
     &           u(i,j,k  ,nrhs), u(i,j,k+1,nrhs),
     &           u(i,j,k+2,nrhs), u(i,j,k+3,nrhs), vel)
          enddo
        enddo

        do i=IstrU,Iend
          vel=0.5*(We(i,j,2)+We(i-1,j,2))
          FC(i,2)=vel*
     &            FLUX3(
     &         u(i,j,1,nrhs), u(i,j,2,nrhs), 
     &         u(i,j,3,nrhs), u(i,j,4,nrhs), vel)

          vel=0.5*(We(i,j,N-2)+We(i-1,j,N-2))
          FC(i,N-2)=vel*
     &            FLUX3(
     &         u(i,j,N-3,nrhs), u(i,j,N-2,nrhs), 
     &         u(i,j,N-1,nrhs), u(i,j,N  ,nrhs), vel)

          vel=0.5*(We(i,j,1)+We(i-1,j,1))
          FC(i,1)=vel*
     &            FLUX2(
     &            u(i,j,1,nrhs), u(i,j,2,nrhs), vel, cdif)

          vel=0.5*(We(i,j,N-1)+We(i-1,j,N-1))
          FC(i,N-1)=vel*
     &              FLUX2(
     &              u(i,j,N-1,nrhs), u(i,j,N,nrhs), vel, cdif)

#   ifdef MOVING_BATHY
          FC(i,0)=0.5*u(i,j,1,nrhs)*
     &            (We(i  ,j,0)+
     &             We(i-1,j,0))
#   else
          FC(i,0)=0.
#   endif
          FC(i,N)=0.
        enddo
