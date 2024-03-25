module sgdmlvars
use allvars,only :Natoms
implicit none 
integer::use_E
real(8),allocatable::E_sgdml,F_sgdml(:),r_desc_test(:),r_d_desc_test(:,:),r_desc_perms_train(:,:),r_d_desc_alpha_perms_train(:,:),alphas_E_lin(:)
integer::io_sgdml,ioparam_sgdml,iomodel_sgdml,iostatparam_sgdml,iostatmodel_sgdml,Nperms,Ntrains,Natoms_sgdml,dim_i,dim_d,dim_c
real(8)::c_sgdml,std_sgdml,sqrt5,sigma_sgdml,sigma_sgdml_inv,mat52_base_fact,diag_scale_fact

contains 

!****************************************************************************************************# 
subroutine loadSgdmlModel()
implicit none 
integer::i,j
character(len=50)::modelfile
ioparam_sgdml =5000 
iomodel_sgdml = 5001
open(unit=io_sgdml,file='sgdml.param',status='old',iostat=iostatparam_sgdml)
if (iostatparam_sgdml.eq.0) then
	read(io_sgdml,*)modelfile
	open(unit=iomodel_sgdml,file=trim(adjustl(modelfile)),status='old',iostat=iostatmodel_sgdml)
	if (iostatmodel_sgdml.eq.0) then
		!Read Number of atoms(Natoms), Number of permutations(N_perm), And number of training set (N_train)
		read(iomodel_sgdml,*)
		read(iomodel_sgdml,*)Natoms_sgdml,Nperms,Ntrains
		if (Natoms.ne.Natoms_sgdml) then
			write(*,*)"Number of atoms does not match between MD input and SGDML model" 
			stop
		endif
		dim_i = 3*Natoms_sgdml
		!Read lattice vectors
 
		!dim_d : Number of descriptors, dim_i: Number of cartesian coordinates 
		dim_d = aint(Natoms * (Natoms-1) * 0.5d0)
		dim_c = Ntrains*Nperms
		read(iomodel_sgdml,*)
		read(iomodel_sgdml,*)sigma_sgdml,c_sgdml,std_sgdml
		!Read permutation indexes
		read(iomodel_sgdml,*)
		read(iomodel_sgdml,*)
		
		!Allocate arrays
		allocate(F_sgdml(0:dim_i-1)) 
		allocate(r_desc_test(0:dim_d-1),r_d_desc_test(0:dim_d-1,0:dim_i-1))
		allocate(r_desc_perms_train(0:dim_c-1,0:dim_d-1),r_d_desc_alpha_perms_train(0:dim_c-1,0:dim_d-1))

		!Initialize all arrays to zero
		!Read Training set descriptors 
		do i = 0,dim_c -1
			read(iomodel_sgdml,*)(r_desc_perms_train(i,j),j=0,dim_d-1)
		enddo
		read(iomodel_sgdml,*)
		read(iomodel_sgdml,*)
		!Read fitted coefficients of Jacobians 
		do i = 0,dim_c -1
			read(iomodel_sgdml,*)(r_d_desc_alpha_perms_train(i,j),j=0,dim_d-1)
		enddo
		!Read fitted Energy coefficents if 
		read(iomodel_sgdml,*)
		read(iomodel_sgdml,*)
		read(iomodel_sgdml,*)use_E
		if (use_E==1) then
			read(iomodel_sgdml,*)alphas_E_lin
		endif
	else
		write(*,*)"Error Reading SGDML model file!! "
	endif
	
else
	write(*,*)"Error Reading SGDML param file!! "
	stop
endif
!Define a few parameters 
sigma_sgdml_inv =  1.d0 / sigma_sgdml
mat52_base_fact = 5.d0 / (3.d0 * sigma_sgdml**3)
diag_scale_fact = 5.d0 / sigma_sgdml
sqrt5 = dsqrt(5.d0)
end subroutine 

!****************************************************************************************************# 
!Prediction function 
subroutine getForce(Xinp,F_sgdml,E_sgdml)
use allvars, only : Ncarts
implicit none 
integer::i
real(8),intent(in)::Xinp(0:Ncarts-1)
real(8),intent(out)::E_sgdml,F_sgdml(0:Ncarts-1)
real(8),allocatable :: r_desc_test(:),r_d_desc_test(:,:)
real(8),allocatable::diff_ab_perms(:,:),norm_ab_perms(:),mat52_base(:),a_x2(:),ax2_mat52_base(:)
real(8),allocatable::F_d_sgdml(:)
allocate(diff_ab_perms(0:dim_c-1,0:dim_d-1),norm_ab_perms(0:dim_c-1),mat52_base(0:dim_c-1),a_x2(0:dim_c-1))
allocate(ax2_mat52_base(0:dim_c-1),F_d_sgdml(0:dim_d-1))
allocate(r_desc_test(0:dim_d-1),r_d_desc_test(0:dim_d-1,0:Ncarts-1))
r_desc_test = 0.d0;r_d_desc_test = 0.d0
call r_to_desc_n_d_desc(Xinp,r_desc_test,r_d_desc_test)
!print*,'inside=',r_desc_test,r_d_desc_test
!Initialize 
diff_ab_perms = 0.d0;mat52_base=0.d0;mat52_base=0.d0

!Initialize force and energies to zeros
E_sgdml=0.d0;F_sgdml=0.d0

!subtract r_desc_test from R_d_desc_alpha_perms_train and store the result to diff_ab_perms
do i = 0,dim_c-1
	diff_ab_perms(i,:) =  r_desc_test - r_desc_perms_train(i,:)
enddo
!Take the norm of descriptors for each permuted points and multiply it by sqrt5
norm_ab_perms = sqrt5 * norm2(diff_ab_perms,dim = 2)
mat52_base =dexp(-norm_ab_perms*sigma_sgdml_inv)
mat52_base =  mat52_base *  mat52_base_fact

!Perform column wise dot product between diff_ab_perms and r_d_desc_alpha_perms_train
do i = 0,dim_c-1
	a_x2(i) = dot_product(diff_ab_perms(i,:),r_d_desc_alpha_perms_train(i,:)) 
enddo

!Get force in descriptors and transform to cartesians
do i = 0,dim_d-1
	F_d_sgdml(i) = dot_product((a_x2 * mat52_base),diff_ab_perms(:,i)) * diag_scale_fact
enddo
mat52_base = mat52_base*(norm_ab_perms + sigma_sgdml)

do i = 0,dim_d-1
	F_d_sgdml(i) = F_d_sgdml(i) - dot_product(mat52_base,r_d_desc_alpha_perms_train(:,i))
enddo
E_sgdml = dot_product(a_x2,mat52_base)*std_sgdml+c_sgdml  
F_sgdml = -std_sgdml*matmul(transpose(r_d_desc_test),F_d_sgdml) !'r_d_desc.T.dot(F)' for our special representation of 'r_d_desc'
end subroutine 


!****************************************************************************************************# 
!DESCRPTOR 
subroutine getdesc(r,desc,pddist,rpdiff)
use allvars, only: Ncarts,Natoms
!use sgdml_vars,only:dim_d
implicit none
real(8),intent(in)::r(0:Ncarts-1)
real(8),intent(out)::rpdiff(0:Natoms-1,0:Natoms-1,0:2)
real(8),intent(out)::pddist(0:dim_d-1)
real(8),intent(out)::desc(0:dim_d-1)
integer::i,j,ctr
rpdiff = 0.d0;desc=0.d0; pddist=0.d0;rpdiff=0.d0
ctr = 0
do i = 1,Natoms-1
    do j = 0,i-1
        rpdiff(i,j,:) = r(3*i:3*i+2)  - r(3*j:3*j+2)
        rpdiff(j,i,:) = -rpdiff(i,j,:)
        pddist(ctr) = norm2(rpdiff(i,j,:))
        ctr = ctr + 1
    enddo
enddo
desc = 1.d0/pddist
end subroutine

!****************************************************************************************************# 
subroutine r_to_desc_n_d_desc(r,desc,d_desc)
use allvars, only: Ncarts,Natoms
!use sgdml_vars,only:dim_d
implicit none
real(8),intent(in)::r(0:Ncarts-1)
real(8)::rpdiff(0:Natoms-1,0:Natoms-1,0:2),pddist(0:dim_d-1)
real(8),intent(out)::desc(0:dim_d-1)
real(8),intent(out)::d_desc(0:dim_d-1,0:Ncarts-1)
integer::i,j,ctr
d_desc = 0.d0 ; desc =0.d0; rpdiff =0.d0
ctr = 0
call getdesc(r,desc,pddist,rpdiff)
do i = 1,Natoms-1
    do j = 0,i-1 
        d_desc(ctr,3*i:3*i+2) = -rpdiff(i,j,:) / pddist(ctr)**3
        d_desc(ctr,3*j:3*j+2) = -d_desc(ctr,3*i:3*i+2)
        ctr = ctr + 1
    enddo
enddo
end subroutine

!****************************************************************************************************# 
subroutine matmul_elements(A,B,m,n,C)
implicit none 
integer,intent(in)::m,n
real(8),intent(in) :: A(0:m-1,0:n-1), B(0:m-1,0:n-1)
real(8),intent(out) :: C(0:m-1,0:n-1)
integer::i,j
C = 0.d0
! Using FORALL construct
forall (i = 0:m-1, j = 0:n-1)
    C(i,j) = A(i,j) * B(i,j)
end forall
end subroutine  
!****************************************************************************************************# 
endmodule
