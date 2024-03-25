subroutine  setupClassicalMD()
	use sgdmlvars!use Gaussparam!use eg_pot_force
    use allvars!, only:
    use brngvars 
    implicit none
    integer :: i,j,k,l,ii,dum,jdim
    integer :: atomCtr,atomNum
    real(8) :: rn,meanv(0:2),sigmav(0:2,0:2),randomv(0:2)
	  real(8),dimension(0:Nbeads-1)::Vt
    initStep = 0
    if (writebinary==1) then
      open(unit=iounitTraj,file=trim(adjustl(trajFile)),status="unknown",form='unformatted')
    else
      open(unit=iounitTraj,file=trim(adjustl(trajFile)),status="unknown")
    endif
    open(unit=iounitThermo,file=trim(adjustl(thermoFile)),status="unknown")
    open(unit=iounitXYZ,file=trim(adjustl(XYZFile)),status="unknown")
    if (restart==1) then
        open( unit=iounitRestInTraj,file=inprestartTrajFile,status ="old" )

        do j = 0 , Natoms-1
            read(iounitRestInTraj,*)initStep,(q(3*j+i,k),i=0,Ndim-1),(p(3*j+ii,k),ii=0,Ndim-1)         
        enddo
        !If Nose-Hoover thermostattype is chosen
        if (thermostattype==1) then
            open( unit=iounitRestInMNHC,file=inprestartmnhcTrajFile,status ="old" )
          
            qNose=0.d0
            VelNose=0.d0
            do i = 0,Nchains-1
               do k = 0 , Nbeads-1
                    do j = 0 , Natoms-1
                       read(iounitRestInMNHC,*)initStep,qNose(3*j:3*j+2,k,i),VelNose(3*j:3*j+2,k,i)     
                    enddo
                 enddo
            enddo      
        endif
    else 
        do i =1,initRandSeed
          call random_number(rn)
        enddo
        do j =0,Natoms-1
              q(3*j:3*j+2,0) = Xinit(3*j:3*j+2)
              meanv = 0.5d0 * betaInv; 
              sigmav = 0.d0
              do jdim=0,2
                  sigmav(jdim,jdim) = dsqrt(1.d0 * betaInv * massInv(j))     
              enddo
              errcode = vdRngGaussianMV(brngmethod,stream,3,randomv,1,VSL_MATRIX_STORAGE_FULL,meanv,sigmav)
              do jdim = 0,2
                  p(3*j+jdim,0) = mass(j) * randomv(jdim)
              enddo
        enddo 
    endif
    call getCoM()
    call removeCoMCart()
    call removeCoMVel()
    call classicalMD() 
end subroutine 
