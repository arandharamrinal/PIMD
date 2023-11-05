subroutine  setupPIMDPC()
      use allvars!, only:
      !Nbeads,PreToNonpreNM,PreToNonpreCart,NonpreToPreNM,NonpreToPreCart,sqrtNbeadsInv,sqrtNbeads,Ndim,Natoms,restart,q,p 
      implicit none
      integer :: i,j,k,l,ii,dum
      integer :: atomCtr,beadCtr,atomNum
      real(8) :: ke,usigma,vsigma,rn
      write(nbeadsChar,*)Nbeads
      888 format (A6,3E24.12)
      PreToNonpreNM   = 1.d0
      PreToNonpreCart   = 1.d0
      NonpreToPreNM   = sqrtNbeadsInv
      NonpreToPreCart = sqrtNbeads
      initStep = 0
      if (writebinary==1) then
        open(unit=iounitTraj,file=trim(adjustl(trajFile)),status="unknown",form='unformatted')
      else
        open(unit=iounitTraj,file=trim(adjustl(trajFile)),status="unknown")
      endif
      open(unit=iounitThermo,file=trim(adjustl(thermoFile)),status="unknown")
      open(unit=iounitCentTraj,file=trim(adjustl(centroidtrajFile)),status="unknown")
      open(unit=iounitXYZ,file=trim(adjustl(XYZFile)),status="unknown")
      open(unit=iounitCentXYZ,file=trim(adjustl(centroidXYZFile)),status="unknown")
      if (restart==1) then
          open( unit=iounitRestInTraj,file=inprestartTrajFile,status ="old" )
          open( unit=iounitRestInMNHC,file=inprestartmnhcTrajFile,status ="old" )
          do k = 0 , Nbeads-1
               do j = 0 , Natoms-1
                   read(iounitRestInTraj,*)initStep,(q(3*j+i,k),i=0,Ndim-1),(p(3*j+ii,k),ii=0,Ndim-1)         
               enddo
          enddo
          call CNTrans()
          call CNTransP()
          call NCTransP()
          qNose=0.d0
          VelNose=0.d0
          do i = 0,Nchains-1
             do k = 0 , Nbeads-1
                  do j = 0 , Natoms-1
                     read(iounitRestInMNHC,*)initStep,qNose(3*j:3*j+2,k,i),VelNose(3*j:3*j+2,k,i)     
                  enddo
               enddo
          enddo      
      else 
          do j =0,Ncarts-1
               do k =0,Nbeads-1
                  if (k==0) then
                        u(j,0) = Xinit(j)
                  else
                    !rms displacement at beta   
                    usigma = dsqrt(1.d0/beta/omegaP2/massK(j,k))
                    !rms displacement at beta_beadspread    
                    usigma = usigma*dsqrt(TbeadsInv/beta)
                    call getNormalRand(rn)
                    u(j,k) = usigma*rn
                  endif 
             enddo
          enddo 
          do j =0,Ncarts-1
               do k =0,Nbeads-1
                  if (k==0) then
                        up(j,k) = 0.0d0
                  else
                    call getNormalRand(rn)
                    vsigma = dsqrt(1.d0/beta/massKprime(j,k))
                    up(j,k) = vsigma*rn*massKprime(j,k)
                  endif 
             enddo
          enddo 
          call NCTrans()
          call NCTransP()
      endif
      !call CNTransP()
      call getCoM()
      call removeCoMCart()
      call CNTrans()
      call removeCoMVel()
      call NCTransP()
	  call IntpimdPC()
end subroutine 

