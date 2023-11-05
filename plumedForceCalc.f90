subroutine plumedGetForce(x)
use plumedVars
use allvars, only: Ncarts,CurrStep,mass
implicit none
real(8),intent(in),dimension(0:Ncarts-1)::x
	!plmdPot = 0.d0
    !plmdForce = 0.d0
    !Calls to pass data to plumed
    call  plumed_f_gcmd("setStep"//char(0),CurrStep);                    !Pass a pointer to the current timestep to plumed
    !*** The way that you pass positions will depend on how they are stored in your code.  If the x, y and z position are all stored in a single array you may use:
    call  plumed_f_gcmd("setPositions"//char(0),x);               !Pass a pointer to the first element in the atomic positions array to plumed assuming they are stored in a x1,y1,z1,x2,y2,z2 ... kind of ordering
    !! *** Othersize if you pass the three separate vectors of x, y and z positions using:
    call  plumed_f_gcmd("setMasses"//char(0),mass);                     ! Pass a pointer to the first element in the masses array to plumed
    call  plumed_f_gcmd("setCharges"//char(0),charge);                      ! Pass a pointer to the first element in the charges array to plumed
    call  plumed_f_gcmd("setBox"//char(0),plmdBox);                          ! Pass a pointer to the first element in the box share array to plumed
    call  plumed_f_gcmd("setEnergy"//char(0),plmdPot);                   ! Pass a pointer to the current value of the potential energy to plumed?
    !call  plumed_f_gcmd("setEnergy"//char(0),extV(0));                   ! Pass a pointer to the current value of the potential energy to plumed?
	plmdForce = 0.d0	
    !*** The way that you pass forces will depend on how they are stored in your code.  If the x, y and z force are all stored in a single array you may use:
    call  plumed_f_gcmd("setForces"//char(0),plmdForce);                 ! Pass a pointer to the first element in the foces array to plumed
    call  plumed_f_gcmd("setVirial"//char(0),plmdVirForce);                 ! Pass a pointer to the first element in the virial array to plumed
    
    ! Calls to do actual calculations
    call  plumed_f_gcmd("calc"//char(0),0);                              !Calculate and apply forces from the biases defined in the plumed input
end subroutine     
