!****************************************************************************************************#  
subroutine LangevinThermostat()
  use allvars, only: Ncarts, sqrtMbyBetan,p,CLT2,CLT1
  use brngvars
  implicit none 
  integer  :: j
  real(8)  :: randomx(0:Ncarts-1)
  randomx = 0.d0
  errcode = vdRngGaussianMV(brngmethod,stream,1,randomx,Ncarts,VSL_MATRIX_STORAGE_FULL,mean_arr,covmat)
  p(:,0) = CLT1 * p(:,0) + sqrtMbyBetan * CLT2 * randomx  

endsubroutine 
!****************************************************************************************************#  
