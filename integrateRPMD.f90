!************************************************************************************************!
!Integration Steps  
!************************************************************************************************!
subroutine updaterpmd()
    use sgdmlvars!use Gaussparam!use egPot_force
    use allvars, only: q,p,u,dVdq,dt,massKprime,cartmassInv,extV,Nbeads,Ncarts
    use omp_lib
    !Updates velocity and momentum of rpmd simulations using verlet algorithm
    implicit none
    integer :: j,k,flag
    !Get the external cartesian gradient(dVdq) 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
        !    call get_force(q(:,k),dVdq(:,k),extV(k),0)
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO
    !Transform cartesian gradient (dVdq) to normal mode space 
    call CNPottrans()

    !Remove external torque on the system
    call removeForce()

    !Convert Normal mode to cartesian gradient 
    call NCPottrans()

    !update momentum(Half time step) due to the external force.
    p = p - 0.5d0 * dt * dVdq

    !If Nbeads is 1 run classical MD. Else run RPMD/PIMD

    if (Nbeads .eq. 1) then
        do j = 0, Ncarts-1
                q(j,0) = q(j,0) + p(j,0) * dt * cartmassInv(j)
        end do

    else
        call freeRPFullStep()
    end if
        
    !Get the external gradient(cartesian) dVdq 
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO
    !Transform dVdq to normal mode space 
    call CNPottrans()

    !Remove external torque on the system
    call removeForce()
    !Convert Normal mode to cartesian gradient 
    call NCPottrans()

    ! Update momentum (half time step)
    p = p - 0.5d0 * dt * dVdq
end subroutine 
!************************************************************************************************!
subroutine rpmdstep()
   use allvars, only : t,currStep,totSteps,NstepSave,dt,Ndim,Natoms,Nbeads
   implicit none
   integer :: i,j,k,l,chain   
   integer :: atomctr,beadctr,atom_num
   real(8) :: vp(0:Ndim-1,0:Natoms-1,0:Nbeads-1)
   vp = 0.d0
   t=0.d0
   currStep = 0
   call removeCoMvel()
   call NCtrans()
   call NCtransP()
   call removeRotation()
   call NCtransP()
   do  while (currStep<totSteps)
        call updaterpmd()
        call CNtrans() 
        call CNtransP() 
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
