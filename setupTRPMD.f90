!****************************************************************************************************#  
subroutine setupTRPMD()
   use allvars;use constants
   implicit none 
   integer :: i,j,jj,k,l,ii,dum,frameNumber
   integer :: atomctr,beadctr,atomNum
   real(8) :: dumr,vp(0:Ndim-1,0:Natoms-1,0:Nbeads-1)
   real(8) :: cmtmp(0:Ndim-1)
   character(len=20)::framenamechar
   vp = 0.d0
   PreToNonpreNM   = sqrtNbeadsInv
   PreToNonpreCart   = sqrtNbeads
   NonpreToPreNM   = 1.d0
   NonpreToPreCart = 1.d0 
   !open PIMD trajectory file. The coordinates for each bead is written one after
   !the another
   !Read PIMD traj file frame by frame and start rpmd trajectories from each frame.
   do frameNumber= 1,NtotFrames
      open(unit=iounitTrajReadfile,file=trim(adjustl(inptrajFile)),status="unknown")
      !NlinesPerFrame is the total number of lines in each frame(Nbeads*Natoms)
      NlinesPerFrame=Natoms*Nbeads
      if (frameNumber>=JinitFrame) then
         write(framenamechar,*)frameNumber
         if (writebinary==1) then
                open(unit=iounitTraj,file=trim(adjustl(trajFile))//trim(adjustl(framenamechar)),status="unknown",form='unformatted')
         else
                open(unit=iounitTraj,file=trim(adjustl(trajFile))//trim(adjustl(framenamechar)),status="unknown")
         endif
         open(unit=iounitTraj,file=trim(adjustl(trajFile))//trim(adjustl(framenamechar)),status="unknown")
         open(unit=iounitThermo,file=trim(adjustl(thermoFile))//trim(adjustl(framenamechar)),status="unknown")
         open(unit=iounitCentTraj,file=trim(adjustl(centroidtrajFile))//trim(adjustl(framenamechar)),status="unknown")
         open(unit=iounitXYZ,file=trim(adjustl(XYZFile))//trim(adjustl(framenamechar)),status="unknown")
         open(unit=iounitCentXYZ,file=trim(adjustl(centroidXYZFile))//trim(adjustl(framenamechar)),status="unknown")
         atomNum = 0
         atomctr = 0
         beadctr = 0
         do j=0,NlinesPerFrame-1
             read(iounitTrajReadfile,*),dum,(q(3*atomNum+i,beadctr),i=0,Ndim-1),(p(3*atomNum+ii,beadctr),ii=0,Ndim-1)
             if (atomNum == Natoms-1) then
                  atomNum = 0
                  beadctr = beadctr + 1
             else 
                  atomNum = atomNum + 1
             endif
         enddo
         !print*,'p1=',p(0,:)
         call CNtransP()
         up = up * sqrtNbeadsInv
         do k = 0,Nbeads-1
            do j = 0,Ncarts-1
               up(j,k) = up(j,k) *  pimdTorpmdMass(j,k) 
            enddo
         enddo
         call NCtransP()
         p = p * sqrtNbeads
         !transform cartesian to normal mode coordinates
         !call get_center_of_mass_cart()
         !call removeCoM_cart()
         call CNtrans()
         call CNtransP()
         !do k = 0,Nbeads-1
         !   write(*,88)(up(j,k),j=0,Ncarts-1)
         !read(*,*)
         !enddo
         !call get_center_of_mass()
         call removeCoM()
                
         call NCtrans()
         call removeCoMvel()
         call NCtransP()
         !call get_center_of_mass_cart()
         call removeRotation()
         call NCtransP()
         !***********************Compute Eergies***************************! 
         !Harmonic energies of the ring polymer
         !call get_ringPolymer_energy()
         !compute kinetic energies
         !call get_kinetic_energy()
         call getTotalEnergy()
         !compute the external potential 
         !do k = 0,Nbeads-1
         !   call  get_2fePot(q(:,k),ext_V(k))
         !enddo
         !call get_OHPot()
         !call get_H2OPot()
         !call ext_V_nh3()
         !print*,E_tot*NbeadsInv,Ek,Ering,sum(ext_V)/Nbeads,instTemp
         !call the matrix for free ring polymer evolution 
         if (Cayley==1) then
             call CayleyPolyMat()
         else
             call polyMat()
         endif 
         !Run rpmd 
         if (method==RPMD) then
             call rpmdstep()
         elseif (method==TRPMD) then
             call trpmdstep()
         endif
         do jj = 0,skipFrames-1
            atomNum = 0
            atomctr = 0
            beadctr = 0
            do l=0,NlinesPerFrame-1
                read(iounitTrajReadfile,*),dum,(q(3*atomNum+i,beadctr),i=0,Ndim-1),(p(3*atomNum+ii,beadctr),ii=0,Ndim-1)
                if (atomNum == Natoms-1) then
                     atomNum = 0
                     beadctr = beadctr + 1
                else 
                     atomNum = atomNum + 1
                endif
            enddo
         enddo
      else 
         atomNum = 0
         atomctr = 0
         beadctr = 0
         do j=0,NlinesPerFrame-1
             read(iounitTrajReadfile,*),dum,(q(3*atomNum+i,beadctr),i=0,Ndim-1),(p(3*atomNum+ii,beadctr),ii=0,Ndim-1)
             if (atomNum == Natoms-1) then
                  atomNum = 0
                  beadctr = beadctr + 1
             else 
                  atomNum = atomNum + 1
             endif
         enddo
      endif
   enddo
end subroutine 
!****************************************************************************************************#  

