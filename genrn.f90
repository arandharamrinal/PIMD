!******************************************************************************************************!
subroutine initRandomSeed(seed)
implicit none
integer,intent(in)   :: seed 
integer,dimension(1) :: seedarr
seedarr = seed
call random_seed(put=seedarr)
end subroutine 
!******************************************************************************************************!
subroutine getNormalRand(rn)
use constants, only:pi
  implicit none
  real(8), intent(out) :: rn 
  real(8) :: u1, u2, z2
  
  ! Generate two uniform random variables in the range (0, 1)
  call random_number(u1)
  call random_number(u2)
  
  ! Use Box-Muller transform to generate normal random variables
  rn = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * pi * u2)
  !z2 = sqrt(-2.0d0 * log(u1)) * sin(2.0d0 * pi * u2)
  
endsubroutine 

!******************************************************************************************************!
