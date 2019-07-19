           ! first-order upwind 
           DO j=jmin,jmax
               DO i=imin,imax+1
                  HUon_w = 0.5*(Huon(i,j,k)+Huon(i,j,k+1))
# ifdef MASKING
     &                        *umask(i,j)   
# endif
                  FX(i,j)=
     &              trb(i-1,j,k,nstp,ig)*max(HUon_w,0.)
     &             +trb(i  ,j,k,nstp,ig)*min(HUon_w,0.)   
               ENDDO
            ENDDO
      
            DO j=jmin,jmax+1
               DO i=imin,imax  
                  HVom_w = 0.5*(Hvom(i,j,k)+Hvom(i,j,k+1))
# ifdef MASKING
     &                        *vmask(i,j)           
# endif
                  FE(i,j)=
     &              trb(i,j-1,k,nstp,ig)*max(HVom_w,0.)
     &             +trb(i,j  ,k,nstp,ig)*min(HVom_w,0.)
               ENDDO
            ENDDO
