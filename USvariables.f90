module usvars
implicit none 
real(8),allocatable :: usK(:),ucEq(:)
real(8) :: dx = 0.002,dx_inv = 500.d0
real(8) :: five_pt_1st_deriv(0:4)
integer :: numUc,iounitUCout,iounitUSparam,iounitUSucdef,nStepSaveUC
integer,allocatable::ucIdx(:,:),uctype(:)  
integer,allocatable::UCuniqIdxarr(:)
end module
