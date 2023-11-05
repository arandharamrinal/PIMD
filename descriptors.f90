subroutine getdesc(r,desc,rpdiff)
use allvars, only: Ncarts,Natoms
use sgdml_vars,only:dim_d
implicit none 
real(8),intent(in)::r(0:Ncarts-1)
real(8),intent(out)::rpdiff(0:Natoms-1,0:Natoms-1)
real(8),intent(out)::desc(0:dim_d-1)
integer::ctr
ctr = 0
rpdiff = 0.d0
do i = 0,Natoms-2
	do j = i+1,Natoms-1
		rpdiff(i,j) = r(3*i:3*i+2)  - r(3*j:3*j+2)
		rpdiff(j,i) = -rpdiff(i,j) 
		desc(ctr) = norm2(rdiff(i,j))
		ctr = ctr + 1
	enddo
enddo
desc = 1.d0/desc
end subroutine


subroutine r_to_d_desc(r,desc,d_desc)
use allvars, only: Ncarts,Natoms
use sgdml_vars,only:dim_d
implicit none
real(8),intent(in)::r(0:Ncarts-1)
real(8),intent(in)::rpdiff(0:Natoms-1,0:Natoms-1)
real(8),intent(out)::desc(0:dim_d-1)
real(8),intent(out)::d_desc(0:dim_d-1,0:Ncarts-1)
integer::ctr
pdist = 1.d0
ctr = 0
call getdesc(r,desc,rpdiff) 
do i = 0,Natoms-2
    do j = i+1,Natoms-1
		d_desc(ctr,3*i:3*i+2) = rpdiff(3*i:3*i+2) / desc(ctr)**3
		d_desc(ctr,3*j:3*j+2) = -d_desc(ctr,3*i:3*i+2) 
		ctr = ctr + 1
    enddo
enddo
end subroutine 
