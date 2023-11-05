!****************************************************************************************************#  
subroutine classicalNVE()
use allvars, only:t,dt,currStep,plumed,Nbeads,Ncarts,mass,q,p,u,up,dVdu,dVdq,extV,NbeadsInv,totSteps, &
        Nref,dtRef,massKprimeInv,NrestartOut,NstepSave,initStep,thermostat,thermostattype
use sgdmlvars!use Gaussparam!use eg_pot_force
use omp_lib
implicit none
integer :: i,j,k,chain,numth 
real(8)::t1,t2
currStep = initStep
totSteps = totSteps + initStep

!Compute cartesian external force
call getForce(q(:,0),dVdq(:,0),extV(0))
call removeRotationClassical()
call removeForceClassical()

do  while (currStep<totSteps)
    !The following code updates the position and the momentum of the system by a time step of dt
    !using the RESPA algorithm

    !Update momentum due to external force(inter-atomic forces) [half-timestep]
    p = p -  0.5d0 * dt * dVdq

    !Update position [full-timestep]
    q = q + dt * p * massKprimeInv
    
    !getForce 
	call getForce(q(:,0),dVdq(:,0),extV(0))
    call removeForceClassical()

    !Update momentum due to external force(inter-atomic forces) [half-timestep]
    p = p -  0.5d0 * dt * dVdq

    !Count current step 
    currStep = currStep + 1    

    if (mod(currStep,NrestartOut)==0) then
       call writeRestartFiles()
    endif 

    if (mod(currStep,NstepSave)==0) then
       call writeToFiles()
    endif 
enddo

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
