!************************************************************************************************!
!Integration Steps  
!************************************************************************************************!
subroutine trpmdOBCBO()
    use sgdmlvars!use Gaussparam!use eg_pot_force
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    use allvars, only: c1PIGLET,c2PIGLET,sqrtMbyBetan,q,p,u,up,dVdq,dt,psiMat,&
	cartmassInv,extV,Nbeads,Ncarts,sqrtMbyBetan
    use omp_lib
    !Updates velocity and momentum of rpmd simulations using verlet algorithm
    implicit none
    integer :: j,k,flag
    psiMat = 0.d0;
    !APPLY PILE THERMOSTAT
    !Generate gaussian random number array (psiMat) with zero mean and unit variance 
    !Cartesian to normal mode momentum transformation 
    call CNtransP
    !Generate gaussian random number array (psiMat) with zero mean and unit
    !variance
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_FULL,mean_arr,covmat) 
    !Apply Langevin thermostat 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLET(:)*up(j,:) + c2PIGLET(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo

    !Normal mode to Cartesian momentum transformation 
    call NCtransP
    
    !MOMENTUM UPDATE
    !update momentum(Half time step) due to the external force.
    p = p - 0.5d0 * dt *  dVdq

    !POSITION UPDATE 
    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do

    else
        call freeRPCayleyFullStep()
    end if

    !Compute External Forces

    !Get the external gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
        !    call getForce(q(:,k),dVdq(:,k),extV(k))
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !REMOVE NET EXTERNAL TORQUE

    !Transform dVdq to normal mode space 
    call CNpotTrans()
    !Remove external torque on the system
    call removeForce()
    !Convert Normal mode to cartesian gradient 
    call NCpotTrans()

    !MOMENTUM UPDATE
    !update momentum(Half time step) due to the external force.
    p = p - 0.5d0 * dt *  dVdq

    !APPLY PILE THERMOSTAT
    !Cartesian to normal mode momentum transformation 
    call CNtransP
    !Generate gaussian random number array (psiMat) with zero mean and unit
    !variance
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_FULL,mean_arr,covmat) 
    !Apply Langevin thermostat 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLET(:)*up(j,:) + c2PIGLET(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    call NCtransP
end subroutine trpmdOBCBO 
!************************************************************************************************!
!************************************************************************************************!
subroutine trpmdBCOCB()
    use sgdmlvars!use Gaussparam!use eg_pot_force
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    use allvars, only: c1PIGLETFull,c2PIGLETFull,sqrtMbyBetan,q,p,u,up,dVdq,dt,psiMat,cartmassInv, &
	extV,Nbeads,Ncarts,sqrtMbyBetan
    use omp_lib
    !Updates velocity and momentum of rpmd simulations using verlet algorithm
    implicit none
    integer :: j,k,flag
    psiMat = 0.d0;

    !**********Compute External Forces***************
    !Get the external potential and the cartesian gradient(dVdq) 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !**********REMOVE NET EXTERNAL TORQUE***********
    !Transform cartesian gradient (dVdq) to normal mode space 
    call CNpotTrans()
    !Remove external torque on the system
    call removeForce()
    !Convert Normal mode to cartesian gradient 
    call NCpotTrans()

    !*************MOMENTUM UPDATE ******************
    !update momentum(Half time step) due to the external force.
    p = p - 0.5d0 * dt *  dVdq

    !*************POSITION UPDATE ******************
    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do

    else
        call freeRPCayleyHalfStep()
    end if

    !**********APPLY PILE THERMOSTAT***********
    !Generate gaussian random number array (psiMat) with zero mean and unit variance 
    !Cartesian to normal mode momentum transformation 
    call CNtransP
    !Generate gaussian random number array (psiMat) with zero mean and unit
    !variance
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_FULL,mean_arr,covmat) 
    !Apply Langevin thermostat 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLETFull(:)*up(j,:) + c2PIGLETFull(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    !Normal mode to Cartesian momentum transformation 
    call NCtransP
    

    !*************POSITION UPDATE half timestep ******************
    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do

    else
        call freeRPCayleyHalfStep()
    end if
    !**********Compute External Forces***************

    !Get the external gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !**********REMOVE NET EXTERNAL TORQUE***********

    !Transform dVdq to normal mode space 
    call CNpotTrans()
    !Remove external torque on the system
    call removeForce()
    !Convert Normal mode to cartesian gradient 
    call NCpotTrans()

    !*************MOMENTUM UPDATE ******************
    !update momentum(Half time step) due to the external force.
    p = p - 0.5d0 * dt *  dVdq
end subroutine trpmdBCOCB 
!************************************************************************************************!
subroutine trpmdOBABO()
    use sgdmlvars!use Gaussparam!use eg_pot_force
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    use allvars, only: c1PIGLET,c2PIGLET,sqrtMbyBetan,q,p,u,up,dVdq,dt,psiMat,cartmassInv, &
	extV,Nbeads,Ncarts,sqrtMbyBetan
    use omp_lib
    !Updates velocity and momentum of rpmd simulations using verlet algorithm
    implicit none
    integer :: j,k,flag
    psiMat = 0.d0;
    !**********APPLY PILE THERMOSTAT***********
    !Generate gaussian random number array (psiMat) with zero mean and unit variance 
    !Cartesian to normal mode momentum transformation 
    call CNtransP
    !Generate gaussian random number array (psiMat) with zero mean and unit
    !variance
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_FULL,mean_arr,covmat) 
    !Apply Langevin thermostat 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLET(:)*up(j,:) + c2PIGLET(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo


    !Normal mode to Cartesian momentum transformation 
    call NCtransP
    
    !*************MOMENTUM UPDATE ******************
    !update momentum(Half time step) due to the external force.
    p = p - 0.5d0 * dt *  dVdq

    !*************POSITION UPDATE ******************
    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do

    else
        call freeRPFullStep()
    end if

    !**********Compute External Forces***************

    !Get the external gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
        !    call getForce(q(:,k),dVdq(:,k),extV(k),0)
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !**********REMOVE NET EXTERNAL TORQUE***********

    !Transform dVdq to normal mode space 
    call CNpotTrans()
    !Remove external torque on the system
    call removeForce()
    !Convert Normal mode to cartesian gradient 
    call NCpotTrans()

    !*************MOMENTUM UPDATE ******************
    !update momentum(Half time step) due to the external force.
    p = p - 0.5d0 * dt *  dVdq

    !APPLY PILE THERMOSTAT
    !Cartesian to normal mode momentum transformation 
    call CNtransP
    !Generate gaussian random number array (psiMat) with zero mean and unit
    !variance
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_FULL,mean_arr,covmat) 
    !Apply Langevin thermostat 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLET(:)*up(j,:) + c2PIGLET(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    call NCtransP
end subroutine trpmdOBABO 
!************************************************************************************************!
!************************************************************************************************!
subroutine trpmdBAOAB()
    use sgdmlvars!use Gaussparam!use eg_pot_force
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    use allvars, only: c1PIGLETFull,c2PIGLETFull,sqrtMbyBetan,q,p,u,up,dVdq,dt,psiMat,cartmassInv, &
	extV,Nbeads,Ncarts,sqrtMbyBetan
    use omp_lib
    !Updates velocity and momentum of rpmd simulations using verlet algorithm
    implicit none
    integer :: j,k,flag
    psiMat = 0.d0;

    !**********Compute External Forces***************
    !Get the external potential and the cartesian gradient(dVdq) 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !**********REMOVE NET EXTERNAL TORQUE***********
    !Transform cartesian gradient (dVdq) to normal mode space 
    call CNpotTrans()
    !Remove external torque on the system
    call removeForce()
    !Convert Normal mode to cartesian gradient 
    call NCpotTrans()

    !*************MOMENTUM UPDATE ******************
    !update momentum(Half time step) due to the external force.
    p = p - 0.5d0 * dt *  dVdq

    !*************POSITION UPDATE ******************
    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do

    else
        call freeRPHalfStep()
    end if

    !**********APPLY PILE THERMOSTAT***********
    !Generate gaussian random number array (psiMat) with zero mean and unit variance 
    !Cartesian to normal mode momentum transformation 
    call CNtransP
    !Generate gaussian random number array (psiMat) with zero mean and unit
    !variance
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_FULL,mean_arr,covmat) 
    !Apply Langevin thermostat 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLETFull(:)*up(j,:) + c2PIGLETFull(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    !Normal mode to Cartesian momentum transformation 
    call NCtransP
    

    !*************POSITION UPDATE half timestep ******************
    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do

    else
        call freeRPHalfStep()
    end if
    !**********Compute External Forces***************

    !Get the external gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
        !    call getForce(q(:,k),dVdq(:,k),extV(k),0)
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !**********REMOVE NET EXTERNAL TORQUE***********

    !Transform dVdq to normal mode space 
    call CNpotTrans()
    !Remove external torque on the system
    call removeForce()
    !Convert Normal mode to cartesian gradient 
    call NCpotTrans()

    !*************MOMENTUM UPDATE ******************
    !update momentum(Half time step) due to the external force.
    p = p - 0.5d0 * dt *  dVdq
end subroutine trpmdBAOAB 
!************************************************************************************************!
subroutine trpmdstep()
   use allvars, only : t,q,dVdq,extV,currStep,totSteps,nStepSave,dt,Ndim,Natoms,Nbeads,Cayley,trpmdIntscheme, &
	BCOCB,OBCBO,BAOAB,OBABO
   use sgdmlvars!use Gaussparam!use eg_pot_force
   implicit none
   integer :: i,j,k,l,chain   
   integer :: atomCtr,beadCtr,atomNum
   t=0.d0
   currStep = 0
   call removeCoMvel()
   call NCtrans()
   call NCtransP()
   call removeRotation()
   call NCtransP()
   if  (BCOCB==1) then 
            do  while (currStep<totSteps)
                 call trpmdBCOCB()
                 call CNtrans() 
                 call CNtransP() 
                 currStep = currStep + 1
                 if (mod(currStep,nStepSave)==0) then 
                    call writeToFiles() 
                 endif
            enddo    
   elseif (OBCBO==1) then
    		!**********Compute External Forces***************
    		!Get the external potential and the cartesian gradient(dVdq) 
    		!$OMP PARALLEL DO PRIVATE(k)
    		do k = 0,Nbeads-1
    		    call getForce(q(:,k),dVdq(:,k),extV(k))
				!call getForce(q(:,k),dVdq(:,k),extV(k))
    		enddo
    		!$OMP END PARALLEL DO

    		!**********REMOVE NET EXTERNAL TORQUE***********
    		!Transform cartesian gradient (dVdq) to normal mode space 
    		call CNpotTrans()
    		!Remove external torque on the system
    		call removeForce()
    		!Convert Normal mode to cartesian gradient 
    		call NCpotTrans()

            do  while (currStep<totSteps)
                 call trpmdOBCBO()
                 call CNtrans() 
                 call CNtransP() 
                 currStep = currStep + 1
                 if (mod(currStep,nStepSave)==0) then 
                    call writeToFiles() 
                 endif
            enddo    
    elseif (BAOAB==1) then
            do  while (currStep<totSteps)
                 call trpmdBAOAB()
                 call CNtrans() 
                 call CNtransP() 
                 currStep = currStep + 1
                 if (mod(currStep,nStepSave)==0) then 
                    call writeToFiles() 
                 endif
            enddo    
    elseif (OBABO==1) then
    		!**********Compute External Forces***************
    		!Get the external potential and the cartesian gradient(dVdq) 
    		!$OMP PARALLEL DO PRIVATE(k)
    		do k = 0,Nbeads-1
				call getForce(q(:,k),dVdq(:,k),extV(k))
    		enddo
    		!$OMP END PARALLEL DO

    		!**********REMOVE NET EXTERNAL TORQUE***********
    		!Transform cartesian gradient (dVdq) to normal mode space 
    		call CNpotTrans()
    		!Remove external torque on the system
    		call removeForce()
    		!Convert Normal mode to cartesian gradient 
    		call NCpotTrans()
            do  while (currStep<totSteps)
                 call trpmdOBABO()
                 call CNtrans() 
                 call CNtransP() 
                 currStep = currStep + 1    
                 if (mod(currStep,nStepSave)==0) then 
                    call writeToFiles() 
                 endif
            enddo    
   	endif
   	if (currStep==totSteps) then
   	  call writeRestartfiles()
   	  call writeToFiles()
   endif    
   Close(6666)
   Close(7777)
   Close(8888)
   Close(9999)
end subroutine
!****************************************************************************************************#  
