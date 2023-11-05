module plumedvars
use allvars, only:tempKelvin,temp,Natoms,Ncarts
implicit none
real(8)::energyUnits,lengthUnits,timeUnits,chargeUnits,massUnits,plmdKbT,plmdPot
real(8),allocatable::charge(:),plmdVirForce(:,:),plmdForce(:),plmdbox(:,:)
integer:: real_precision,res,plumed_version
CHARACTER(LEN=32):: plumedInput,fplog

contains 
subroutine initializePlumedVars()
allocate(charge(0:Natoms-1),plmdVirForce(0:2,0:2),plmdbox(0:2,0:2),plmdForce(0:Ncarts-1))
charge = 0.d0
plmdVirForce = 0.d0
plmdbox = 5.d0
plmdForce = 0.d0
plmdPot = 0.d0
real_precision = 8! integer containing the size of a real number (4 or 8)
energyUnits    = 4.3597447222071d-21*6.02214076d+23!conversion factor between the energy unit used in your code and kJ mol-1
lengthUnits    = 0.0529177210903 !conversion factor between the length unit used in your code and nm
timeUnits      = 2.4188843265857/1.0d+5    !conversion factor between the time unit used in your code and ps
chargeUnits    = 1.d0 !conversion factor between the charge unit used in your code and e
massUnits      = 1.d0/1822.8884572362088 !conversion factor between the mass unit used in your code and amu
plumedInput    =  "plumed.dat"!name of the plumed input file from the md code to plumed
fplog          = 'plumed.log'
!MPI_COMM_WORLD =  !MPI communicator to plumed
fplog          =  "plumed.log"!file on which to write out the plumed log
plmdkbT        =  1.d0 * temp!real containing the value of kbT
!plmdkbT        =   1.38064852d-23/4.359744650d-18 * temp!real containing the value of kbT
res            = 0!an integer saying if we are restarting (zero means no, one means yes)
end subroutine 
end module 
