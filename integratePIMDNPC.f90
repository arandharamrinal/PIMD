!************************************************************************************************!
!Integration Steps  
!************************************************************************************************!
subroutine verletStep()
    use sgdmlvars!use Gaussparam!use egPot_force
    use allvars, only: q,p,dVdq,dt,cartmassInv,extV,NbeadsInv,Nbeads,Ncarts
    implicit none
    integer :: k, j
    !Get the external potential and the cartesian gradient(dVdq) 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !update momentum(Half time step) due to external force.
    p = p - 0.5d0 * dt * dVdq

    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do
    else
        call freeRPFullStep()
    end if

    !Get the external potential and the gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !Update momentum (half time step)
    p = p - 0.5d0 * dt * dVdq

end subroutine verletStep

!************************************************************************************************!
subroutine  pigletOBABO()
    use sgdmlvars!use Gaussparam!use egPot_force
    use allvars, only: up,sqrtMbyBetan,psiMat,c1PIGLET,c2PIGLET,q,p,dVdq,dt,cartmassInv,extV,NbeadsInv,Nbeads,Ncarts
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    !Updates the velocity and momentums of PIMD simulation using verlet
    !algorithm. 
    !Everything is done in the cartesian space and does not remove the external
    !torque  
    !since thermostat is attached. 
    implicit none
    integer :: k, j

    call CNTransP()
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_DIAGONAL,mean_arr,covmat)
    !Apply Langevin thermostate 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLET(:)*up(j,:) + c2PIGLET(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    call NCTransP()


    !Get the external potential and the cartesian gradient(dVdq) 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !update momentum(Half time step) due to external force.
    p = p - 0.5d0 * dt * dVdq

    ! Update position (full time step)
    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD

    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do
    else
        call freeRPFullStep()
    end if

    !Get the external potential and the gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !Update momentum (half time step)
    p = p - 0.5d0 * dt * dVdq

    call CNTransP()
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_DIAGONAL,mean_arr,covmat)
    !Apply Langevin thermostate 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLET(:)*up(j,:) + c2PIGLET(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    call NCTransP()

end subroutine pigletOBABO
!************************************************************************************************!
subroutine  pigletOBCBO()
    use sgdmlvars!use Gaussparam!use egPot_force
    use allvars, only: up,sqrtMbyBetan,psiMat,c1PIGLET,c2PIGLET,q,p,dVdq,dt,cartmassInv,extV,NbeadsInv,Nbeads,Ncarts
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    !Updates the velocity and momentums of PIMD simulation using verlet
    !algorithm. 
    !Everything is done in the cartesian space and does not remove the external
    !torque  
    !since thermostat is attached. 
    implicit none
    integer :: k, j

    call CNTransP()
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_DIAGONAL,mean_arr,covmat)
    !Apply Langevin thermostate 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLET(:)*up(j,:) + c2PIGLET(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    call NCTransP()


    !Get the external potential and the cartesian gradient(dVdq) 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !update momentum(Half time step) due to external force.
    p = p - 0.5d0 * dt * dVdq

    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do
    else
        call freeRPCayleyFullStep()
    end if

    !Get the external potential and the gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !Update momentum (half time step)
    p = p - 0.5d0 * dt * dVdq

    call CNTransP()
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_DIAGONAL,mean_arr,covmat)
    !Apply Langevin thermostate 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLET(:)*up(j,:) + c2PIGLET(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    call NCTransP()

end subroutine pigletOBCBO
!************************************************************************************************!
!************************************************************************************************!
subroutine  pigletBAOAB()
    use sgdmlvars!use Gaussparam!use egPot_force
    use allvars, only: up,sqrtMbyBetan,psiMat,c1PIGLETFull,c2PIGLETFull,q,p,dVdq,dt,cartmassInv,extV,NbeadsInv,Nbeads,Ncarts
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    !Updates the velocity and momentums of PIMD simulation using verlet
    !algorithm. 
    !Everything is done in the cartesian space and does not remove the external
    !torque  
    !since thermostat is attached. 
    implicit none
    integer :: k, j
    
    !B evolution half time step
    !Get the external potential and the cartesian gradient(dVdq) 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
        !       call get_force(q(:,k),dVdq(:,k),extV(k),0)
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !update momentum(Half time step) due to external force.
    p = p - 0.5d0 * dt * dVdq

    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do
    else
        call freeRPHalfStep()
    end if

    call CNTransP()
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_DIAGONAL,mean_arr,covmat)
    !Apply Langevin thermostate 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLETFull(:)*up(j,:) + c2PIGLETFull(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    call NCTransP()

    ! Update position (full time step)
    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD

    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do
    else
        call freeRPHalfStep()
    end if

    !Get the external potential and the gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
        !       call get_force(q(:,k),dVdq(:,k),extV(k),0)
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !Update momentum (half time step)
    p = p - 0.5d0 * dt * dVdq

end subroutine pigletBAOAB
!************************************************************************************************!
!************************************************************************************************!
subroutine  pigletBCOCB()
    use sgdmlvars!use Gaussparam!use egPot_force
    use allvars, only: up,sqrtMbyBetan,psiMat,c1PIGLETFull,c2PIGLETFull,q,p,dVdq,dt,cartmassInv,extV,NbeadsInv,Nbeads,Ncarts
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    !Updates the velocity and momentums of PIMD simulation using verlet
    !algorithm. 
    !Everything is done in the cartesian space and does not remove the external
    !torque  
    !since thermostat is attached. 
    implicit none
    integer :: k, j
    !Get the external potential and the cartesian gradient(dVdq) 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
        !       call get_force(q(:,k),dVdq(:,k),extV(k),0)
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !update momentum(Half time step) due to external force.
    p = p - 0.5d0 * dt * dVdq

    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do
    else
        call freeRPCayleyHalfStep()
    end if

    call CNTransP()
    errcode = vdRngGaussianMV(brngmethod,stream,Ncarts,psiMat,Nbeads,VSL_MATRIX_STORAGE_DIAGONAL,mean_arr,covmat)
    !Apply Langevin thermostate 
    do j = 0,Ncarts-1
        up(j,:) = c1PIGLETFull(:)*up(j,:) + c2PIGLETFull(:) * psiMat(j,:) * sqrtMbyBetan(j)
    enddo
    call NCTransP()

    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD
    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do
    else
        call freeRPCayleyHalfStep()
    end if

    !Get the external potential and the gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !Update momentum (half time step)
    p = p - 0.5d0 * dt * dVdq

end subroutine pigletBCOCB
!************************************************************************************************!
subroutine pimd_mnhc()
    use allvars, only : t,currStep,totSteps,Nbeads,Nchains,NstepSave,dt
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    implicit none
    integer :: i,j,k,chain   
    integer :: atom_ctr,bead_ctr,atom_num
    t=0.d0
    currStep = 0
    do  while (currStep<totSteps)
         call CNTransP()
         if (Nchains>0) then
            do k = 0,Nbeads-1
               if (k==0) then
                        call NHthermostatCent()
           else
                        call NHthermostat(k)
           endif
        enddo
         endif
         call NCTransP()
         call verletStep() 
         call CNTransP()
         if (Nchains>0) then
            do k = 0,Nbeads-1
           if (k==0) then
                           call NHthermostatCent()
           else
                           call NHthermostat(k)
           endif
            enddo
         endif   
         call NCTransP()
         if (mod(currStep,NstepSave)==0) then
            call writeToFiles()
         endif 
         currStep = currStep + 1    
         enddo    
    Close(6666)
    Close(7777)
    Close(8888)
    Close(9999)
end subroutine
!****************************************************************************************************#  
subroutine pimdPiglet()
    use allvars, only : Cayley,trpmdIntScheme,t,currStep,totSteps,Nbeads,NstepSave,dt,psiMat, &
	BCOCB,OBCBO,BAOAB,OBABO
    use brngvars
    use mkl_vsl_type
    use mkl_vsl
    implicit none
    integer :: i,j,k,l,chain   
    integer :: atom_ctr,bead_ctr,atom_num
    t=0.d0
    currStep = 0
    if  (BCOCB==1) then 
        do  while (currStep<totSteps)
             call pigletBCOCB()
             call CNTrans() 
             call CNTransP() 
             currStep = currStep + 1    
             if (mod(currStep,NstepSave)==0) then 
                call writeToFiles() 
             endif
        enddo    
	elseif (OBCBO==1) then
        do  while (currStep<totSteps)
             call pigletOBCBO()
             call CNTrans() 
             call CNTransP() 
             currStep = currStep + 1    
             if (mod(currStep,NstepSave)==0) then 
                call writeToFiles() 
             endif
        enddo    
	elseif (BAOAB==0) then
        do  while (currStep<totSteps)
             call pigletBAOAB()
             call CNTrans() 
             call CNTransP() 
             currStep = currStep + 1    
             if (mod(currStep,NstepSave)==0) then 
                call writeToFiles() 
             endif
        enddo    
    elseif (OBABO==1) then
        do  while (currStep<totSteps)
             call pigletOBABO()
             call CNTrans() 
             call CNTransP() 
             currStep = currStep + 1    
             if (mod(currStep,NstepSave)==0) then 
                call writeToFiles() 
             endif
        enddo    
    endif
	if (currStep==totSteps) then
		call writeRestartFiles()
		call writeToFiles()
	endif    
    Close(6666)
    Close(7777)
    Close(8888)
    Close(9999)
end subroutine
!****************************************************************************************************#  
