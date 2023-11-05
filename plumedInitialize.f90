subroutine InitPlumed()
use plumedvars
use allvars, only:dt
implicit none 
	call  plumed_f_gcreate();                 		! Create the plumed object
	call  plumed_f_gcmd("setRealPrecision"//char(0),real_precision);  ! Pass a pointer to an integer containing the size of a real number (4 or 8)
	call  plumed_f_gcmd("setMDEnergyUnits"//char(0),energyUnits);     ! Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
	call  plumed_f_gcmd("setMDLengthUnits"//char(0),lengthUnits);     ! Pass a pointer to the conversion factor between the length unit used in your code and nm 
	call  plumed_f_gcmd("setMDTimeUnits"//char(0),timeUnits);         ! Pass a pointer to the conversion factor between the time unit used in your code and ps
	
	call plumed_f_gcmd("getApiVersion"//char(0),plumed_version)
	
	if (plumed_version>3) then
		!This is valid only if API VERSION > 3
		call  plumed_f_gcmd("setMDChargeUnits"//char(0),chargeUnits);     ! Pass a pointer to the conversion factor between the charge unit used in your code and e
		!This is valid only if API VERSION > 3
		call  plumed_f_gcmd("setMDMassUnits"//char(0),massUnits);         ! Pass a pointer to the conversion factor between the mass unit used in your code and amu
	endif
	call  plumed_f_gcmd("setPlumedDat"//char(0),&
									trim(adjustl(plumedInput))//char(0));         ! Pass the name of the plumed input file from the md code to plumed
#ifdef MPI
	call  plumed_f_gcmd("setMPIFComm"//char(0),0);        ! Pass a pointer to the MPI communicator to plumed
# endif
	!notice that from fortran the command "setMPIFComm" should be used instead
	call  plumed_f_gcmd("setNatoms"//char(0),Natoms);                 ! Pass a pointer to the number of atoms in the system to plumed
	call  plumed_f_gcmd("setMDEngine"//char(0),"pimd_mrinal");             !/ Pass the name of your md engine to plumed (now it is just a label) 
	call  plumed_f_gcmd("setTimestep"//char(0),dt);              ! Pass a pointer to the molecular dynamics timestep to plumed
	if (plumed_version>1) then
		! This is valid only if API VERSION > 1
		call  plumed_f_gcmd("setKbT"//char(0),plmdkbT);                       ! Pointer to a real containing the value of kbT
	endif
	!Calls to do the actual initialization (all the above commands must appear before this call)
	call plumed_f_gcmd("init"//char(0),0);
endsubroutine 
