program PI_calculation
!************************************************************************************************!
use allvars, only : Ndim,Nbeads,plumed,us,method,Cayley,thermostat,thermostattype,&
				 CLASSICAL,PIMDNPC,PIMDPC,RPMD,TRPMD,MNHC,WNLangevin,initRandSeed
use constants
!use eg_pot_force
use omp_lib
!use potvars
use sgdmlvars!use Gaussparam
use mkl_vsl_type
use mkl_vsl
use brngvars
use plumedvars
implicit none
!************************************************************************************************!
!Initialize all the variables.
call initializeAllvars()

!Write all input variables to a file.
call writelogfile()

!Initialize Random seed
call initRandomSeed(initRandSeed)
!call initializePotvars()
!call initializeGaussParam()
call loadSgdmlModel()
if (us==1) then
	call initializeUSparam()
endif
if (thermostattype == MNHC) then
	call initializeMNHC()
endif

!Get the matrix for converting cartesian to normal mode
if (method.ne.CLASSICAL) then
    call getCmat()
    call eigvals()
    call getOmegaK()
endif

call getMass()

if (plumed==1) then
    call initializePlumedVars()
endif

if (method==CLASSICAL) then
    call initBrngVars(Ncarts)
elseif ((method==PIMDNPC).or.(method==TRPMD)) then
    call initBrngVars(Nbeads)
endif

if (method==TRPMD) then
    call trpmdVars()
endif
!************************************************************************************************!
!pimd step  
!(momentum p) mass m of each bead in PIMD is not the physical mass(momentum).This fact allows us to choose
! beads' masses so that numerical calculations are less expensive . Based on these ideas,
!there are many variants for the numerical implementation of PIMD which differ by the choice of masses. In 
!this code two such variants are implemented, one is given by Manolopoulos(setup_pimd), where masses for 
!the beads are chosen to be the same as the physical mass of that particular atom and, another is the conventional
!algorithm given by Parinello(setup_pimd_tm),Where masses are so chosen such that  each normal mode has the same frequency 
!(This makes the thermal sampling more efficient). setup_rpmd implements the ring polymer molcular dynamics, 
!where noo thermostattype is attached and masses are the physical masses.

if (method==CLASSICAL) then
   print*,"*******************************************"
   print*,"**       Starting Classical MD simulation        **"
   print*,"*******************************************"
   if (thermostat==1) then
        call setupClassicalMD()
   else
        call setupClassicalNVE()
   endif
elseif (method==PIMDNPC) then
   print*,"*******************************************"
   print*,"**       Starting PIMD simulation        **"
   print*,"**       non-preconditioned Algorithm    **"
   print*,"*******************************************"
   call setupPIMDNPC()
else if (method==PIMDPC) then
   print*,"*******************************************"
   print*,"**       Starting PIMD simulation        **"
   print*,"**       pre-conditioned  Algorithm      **"
   print*,"*******************************************"
   call setupPIMDPC()
!************************************************************************************************!
!RPMD step
else if ((method==RPMD).or.(method==TRPMD)) then
   print*,"*******************************************"
   print*,"**       Starting PIMD simulation        **"
   print*,"**       non-preconditioned Algorithm    **"
   print*,"*******************************************"
   call setupTRPMD()
endif  
if (((method==CLASSICAL).or.(method==PIMDNPC)).or.(method==TRPMD)) then
    !**** Deinitialize *****
    errcode=vsldeletestream( stream )
endif
end program
