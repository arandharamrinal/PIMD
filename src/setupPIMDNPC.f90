!****************************************************************************************************#  
subroutine  setupPIMDNPC()
      use allvars;use constants
      implicit none
      integer :: i,j,k,l,ii,dum
      real(8)::rn
      PreToNonpreNM = sqrtNbeadsInv  
      PreToNonpreCart = sqrtNbeads
      NonpreToPreNM   = 1.d0
      NonpreToPreCart = 1.d0 
      if (writebinary==1) then
        open(unit=iounitTraj,file=trim(adjustl(trajFile)),status="unknown",form='unformatted')
      else
        open(unit=iounitTraj,file=trim(adjustl(trajFile)),status="unknown")
      endif
      open(unit=iounitThermo,file=trim(adjustl(thermoFile)),status="unknown")
      open(unit=iounitCentTraj,file=trim(adjustl(centroidtrajFile)),status="unknown")
      open(unit=iounitXYZ,file=trim(adjustl(XYZFile)),status="unknown")
      open(unit=iounitCentXYZ,file=trim(adjustl(centroidXYZFile)),status="unknown")
      do i =1,1000
        call random_number(rn)
      enddo
      do j =0,Ncarts-1
           do k =0,Nbeads-1
              call random_number(rn)
              phi=pi*2.d0*rn
              if (k==0) then
                    u(j,k) = Xinit(j)*sqrtNbeads
                    up(j,k) = 0.0d0
              else
                    u(j,k)  = dsqrt(1.d0*Ethermal*massKprimeInv(j,k))*(dcos(phi)/OmegaK(k))
                    up(j,k) = dsqrt(1.d0*massKprime(j,k)*Ethermal)*dsin(phi)  
              endif 
            enddo
      enddo  
      !compute kinetic energies
      call getRPkinE()
      instTemp = NbeadsInv * 2.d0 * Ek *totDofInv * AuToKelvin
      call NCtrans()
      call NCtransp()
      if (Cayley==1) then
          call CayleyPolyMat()
      else
          call polyMat()
      endif 
      if (thermostattype==PILE) then	
          call pimdPIGLET()
      endif 
end subroutine 
!****************************************************************************************************#  
