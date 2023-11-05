subroutine gaussianrndist()

use mkl_vsl_type
use mkl_vsl
implicit none 
integer,intent(in)::Ncarts,Nbeads
real(8),allocatable,intent(in)::mean_arr(:),covmat(:,:)
real(8),allocatable,intent(out)::gamma_mat
allocate(mean_arr(0:Nbeads-1),covmat(0:Nbeads-1,0:Nbeads-1),gamma_mat(0:Ncarts-1,0:Nbeads-1))

mean_arr = 0.d0;covmat = 0.d0;gamma_mat=0.d0
do i=0,Nbeads-1 
        covmat(i,i) =0.d0
enddo 
vdRngGaussianMV(VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2,stream,Ncarts,gamma_mat, Nbeads,VSL_MATRIX_STORAGE_DIAGONAL,mean_arr,  )

