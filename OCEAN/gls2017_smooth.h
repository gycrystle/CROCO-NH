#   ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=J_EXT_RANGE
          trb(Istr-1,j,k,nnew,ig)=trb(Istr,j,k,nnew,ig)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=J_EXT_RANGE
          trb(Iend+1,j,k,nnew,ig)=trb(Iend,j,k,nnew,ig)
        enddo
      endif
#   endif
#   ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=I_EXT_RANGE
          trb(i,Jstr-1,k,nnew,ig)=trb(i,Jstr,k,nnew,ig)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=I_EXT_RANGE
          trb(i,Jend+1,k,nnew,ig)=trb(i,Jend,k,nnew,ig)
        enddo
      endif
#    ifndef EW_PERIODIC
      if (WESTERN_EDGE.and.SOUTHERN_EDGE) then
        trb(Istr-1,Jstr-1,k,nnew,ig)=trb(Istr,Jstr,k,nnew,ig)
      endif
      if (WESTERN_EDGE.and.NORTHERN_EDGE) then
        trb(Istr-1,Jend+1,k,nnew,ig)=trb(Istr,Jend,k,nnew,ig)
      endif
      if (EASTERN_EDGE.and.SOUTHERN_EDGE) then
        trb(Iend+1,Jstr-1,k,nnew,ig)=trb(Iend,Jstr,k,nnew,ig)
      endif
      if (EASTERN_EDGE.and.NORTHERN_EDGE) then
        trb(Iend+1,Jend+1,k,nnew,ig)=trb(Iend,Jend,k,nnew,ig)
      endif
#    endif
#   endif

         DO j=jstr-1,jend+1              
            DO i=istr,iend+1          
               FX (i,j  )=( trb(i  ,j,k,nnew,ig) 
     &                   -  trb(i-1,j,k,nnew,ig) ) 
#  ifdef MASKING
     &                                 *umask(i,j)  
#  endif
            ENDDO                    
         ENDDO
         DO j=jstr,jend+1                 
            DO i=istr-1,iend+1
               FE1(i,j,0)=( trb(i,j  ,k,nnew,ig) 
     &                    - trb(i,j-1,k,nnew,ig) ) 
#  ifdef MASKING
     &                                 *vmask(i,j)
#  endif
            ENDDO
            DO i=istr,iend
              FE(i,j)=FE1(i,j,0) 
     &         +  smth_a*( FX(i+1,j)+FX(i  ,j-1)
     &                 -FX(i  ,j)-FX(i+1,j-1))
            ENDDO
         ENDDO
         DO j=jstr,jend
            DO i=istr,iend+1
              FX(i,j)=FX(i,j  ) 
     &         + smth_a*( FE1(i,j+1,0)+FE1(i-1,j  ,0)
     &                -FE1(i,j  ,0)-FE1(i-1,j+1,0))
            ENDDO
            DO i=istr,iend
               trb(i,j,k,nnew,ig)=trb(i,j,k,nnew,ig) 
     &         + smth_b*( FX(i+1,j)-FX(i,j)
     &                 +FE(i,j+1)-FE(i,j) )
#  ifdef MASKING
     &                            *rmask(i,j)
#  endif
            ENDDO
         ENDDO              !--> discard FX,FE,FE1
