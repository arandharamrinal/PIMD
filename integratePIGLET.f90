!****************************************************************************************************#  
subroutine PIGLET(ibead)
  use allvars, only: Ncarts, zetaPIGLET, up, c1PIGLET, sqrtMbyBetan, c2PIGLET,psiMat
  implicit none 
  integer,intent(in) :: ibead
  integer  :: j
  do j = 0,Ncarts-1 
     up(j,ibead) = c1PIGLET(ibead) * up(j,ibead) + sqrtMbyBetan(j) * c2PIGLET(ibead) *   psiMat(j,ibead)
  enddo

endsubroutine 
!****************************************************************************************************#  
