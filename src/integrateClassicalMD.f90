!****************************************************************************************************#  
subroutine classicalMD()
use allvars, only:t,dt,currStep,plumed,Nbeads,Ncarts,mass,q,p,dVdq,extV,NbeadsInv,totSteps, &
        Nref,dtRef,massKprimeInv,NrestartOut,NstepSave,initStep,thermostat,thermostattype,&
		MNHC,WNLangevin,us,usPotCent,usGradCent	
use sgdmlvars
use plumedvars
use omp_lib
implicit none
integer :: i,j,k,chain,numth 
real(8)::t1,t2
currStep = initStep
totSteps = totSteps + initStep
!if (plumed==1) then
!    call InitPlumed()
!endif
 
!Compute cartesian external force
call getForce(q(:,0),dVdq(:,0),extV(0))
!Compute biased force 
!if (plumed==1) then
!     call plumedGetForce(q(:,0))
!     dVdq(:,0) = dVdq(:,0) - plmdForce 
!endif
if (us==1) then 
	usPotCent = 0.d0; usGradCent = 0.d0
	call computeUSgrad(q(:,0),usGradCent,usPotCent)
	dVdq(:,0) = dVdq(:,0) + usGradCent
endif 
!

do  while (currStep<totSteps)
    !The following code updates the position and the momentum of the system by a time step of dt
    !using the RESPA algorithm

    !Centroid Nose-Hoover thermostattype update(Half-timestep) 
    if (thermostat==1) then
        if (thermostattype==MNHC) then 
                call NHthermostatCent()
        else if (thermostattype==WNLangevin) then
            call Langevinthermostat()
        endif
    endif

    !Update momentum due to external force(inter-atomic forces) [half-timestep]
    p = p -  0.5d0 * dt * dVdq

    !Update position [full-timestep]
    q = q + dt * p * massKprimeInv
    
    !getForce 
	call getForce(q(:,0),dVdq(:,0),extV(0))

    !!Compute biased force 
    !if (plumed==1) then
    !     call plumedGetForce(q(:,0))
    !     dVdq(:,0) = dVdq(:,0) - plmdForce 
    !endif
	if (us==1) then 
		usPotCent = 0.d0; usGradCent = 0.d0
		call computeUSgrad(q(:,0),usGradCent,usPotCent)
		dVdq(:,0) = dVdq(:,0) + usGradCent
	endif 

    !Update momentum due to external force(inter-atomic forces) [half-timestep]
    p = p -  0.5d0 * dt * dVdq

    if (thermostattype==MNHC) then 
        !Centroid Nose-Hoover thermostattype update(Half-timestep) 
            call NHthermostatCent()
    else if (thermostattype==WNLangevin) then
        call LangevinThermostat()
    endif

    if (mod(currStep,NrestartOut)==0) then
       call writeRestartFiles()
    endif 

    if (mod(currStep,NstepSave)==0) then
       call writeToFiles()
    endif 
    !Count current step 
    currStep = currStep + 1    
enddo

if (mod(totSteps,NstepSave).ne.0) then
    call writeRestartFiles()
    call writeToFiles()
endif    

!if (plumed==1) then 
!    call plumed_f_gfinalize()
!endif

Close(6666)
Close(7777)
Close(8888)
Close(9999)
end subroutine
!****************************************************************************************************#  
