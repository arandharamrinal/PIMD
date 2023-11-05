module brngvars
use mkl_vsl_type
use mkl_vsl
use allvars, only: Ncarts,Nbeads
implicit none
integer(kind=4):: errcode
integer::brng,brngmethod
real(8),allocatable :: mean_arr(:),covmat(:,:)
TYPE (VSL_STREAM_STATE) :: stream
contains
subroutine  initBrngVars(rvdim)
        use mkl_vsl_type
        use mkl_vsl
        integer::i,seed
		integer,intent(in)::rvdim
        brng   = VSL_BRNG_MT19937
        brngmethod = VSL_RNG_METHOD_GAUSSIAN_ICDF
		call system_clock(count=seed)
        allocate(mean_arr(0:rvdim-1),covmat(0:rvdim-1,0:rvdim-1))
        
        mean_arr = 0.d0;covmat = 0.d0;
        errcode=vslnewstream( stream, brng,  seed )
        do i=0,rvdim-1
                covmat(i,i) =1.d0
        enddo
endsubroutine 
endmodule 
