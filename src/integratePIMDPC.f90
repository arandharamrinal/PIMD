!****************************************************************************************************#  
subroutine IntpimdPC()
use allvars, only:t,dt,currStep,plumed,Nbeads,Ncarts,mass,q,p,u,up,dVdq,dVDu,extV,NbeadsInv,totSteps, &
        Nref,dtRef,massK,massKprimeInv,omegaP2,NrestartOut,NstepSave,initStep,us,usPotCent,usGradCent
use sgdmlvars!use Gaussparam!use egpot_force
use plumedvars
use omp_lib
implicit none
integer :: i,j,k,chain,numth 
real(8)::t1,t2
currStep = initStep
totSteps = totSteps + initStep

!if (plumed==1) then
!    call initPlumed()
!endif

!Convert Normal-mode to cartesian position
call NCtrans()
!Compute cartesian external force
!$OMP PARALLEL DO PRIVATE(k)
do k = 0,Nbeads-1
    !call get_force(q(:,k),dVdq(:,k),extV(k),0)
	call getForce(q(:,k),dVdq(:,k),extV(k))
enddo
!$OMP END PARALLEL DO

!Convert cartesian force to Normal mode force 
dVdq = NbeadsInv * dVdq
call CNpottrans()
!!Compute biased force 
!if (plumed==1) then
!	plmdpot = 0.d0; plmdForce = 0.d0 
!    call plumedGetForce(u(:,0))
!    dVdu(:,0) = dVdu(:,0) - plmdForce 
!endif
if (us==1) then 
	usPotCent = 0.d0; usGradCent = 0.d0
	call computeUSgrad(u(:,0),usGradCent,usPotCent)
	dVdu(:,0) = dVdu(:,0) + usGradCent
endif 

do  while (currStep<totSteps)
    !The following code updates the position and the momentum of the system by a time step of dt
    !using the RESPA algorithm

    !Centroid Nose-Hoover thermostat update(Half-timestep) 
    call NHthermostatCent()
    !openMP
    !Update momentum due to external force(inter-atomic forces) [half-timestep]
    up = up -  0.5d0 * dt * dVdu
    !MTS Step
    !OpenMP loop
    !$OMP PARALLEL DO PRIVATE(k,i)
    do k = 0,Nbeads-1
        !numth =  OMP_get_num_threads() 
        !loop over Nref 
        do i = 1,Nref

            !Non-Centroid Nose-Hoover thermostat update(Half-timestep) 
            if (k/=0) then 
                    call NHthermostat(k)
            endif

            !Update momentum due to harmonic forces of the ring-polymer
            up(:,k) = up(:,k) - 0.5d0  * dtRef  * massK(:,k) * omegaP2 * u(:,k)
            !Update positions of all the modes.
            u(:,k) =   u(:,k) + dtRef * up(:,k) * massKprimeInv(:,k)

            !Update momentum due to harmonic forces of the ring-polymer
            up(:,k) = up(:,k) - 0.5d0  * dtRef  * massK(:,k) * omegaP2 * u(:,k)

            !Non-Centroid Nose-Hoover thermostat update(Half-timestep ) 
            if (k/=0) then 
               call NHthermostat(k)
            endif

        enddo
    enddo
    !$OMP END PARALLEL DO
    !Convert Normal-mode to cartesian position
    call NCtrans()

    !Compute cartesian external force
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 0,Nbeads-1
		call getForce(q(:,k),dVdq(:,k),extV(k))
    enddo
    !$OMP END PARALLEL DO

    !Convert cartesian force to Normal mode force 
    dVdq = NbeadsInv * dVdq
    call CNpottrans()

    !!Compute biased force 
    !if (plumed==1) then
	!    plmdpot = 0.d0; plmdForce = 0.d0 
    !    call plumedGetForce(u(:,0))
    !    dVdu(:,0) = dVdu(:,0) - plmdForce 
    !endif
	if (us==1) then 
		usPotCent = 0.d0; usGradCent = 0.d0
		call computeUSgrad(u(:,0),usGradCent,usPotCent)
		!print*,usGradCent;read(*,*)
		dVdu(:,0) = dVdu(:,0) + usGradCent
	endif 
    !Update momentum due to external force(inter-atomic forces) [half-timestep]
    up = up - 0.5d0 * dt *  dVdu
    !Centroid Nose-Hoover thermostat update(Half-timestep ) 
    call NHThermostatCent()

    call NCtransP()

    if (mod(currStep,NrestartOut)==0) then
       call writeRestartFiles()
    endif 
    if (mod(currStep,NstepSave)==0) then
       call writeToFiles()
    endif 
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
