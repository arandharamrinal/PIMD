module POTVARS
!use allvars, only : 
use allvars
implicit none 
integer::Mmax,Nmax,Lmax,n_eb,n_ea,n_ed, n_ob,n_oa,n_od,n_fix_tors
real(8),parameter,dimension(0:2)::tors_idx =[17,22,23]
real(8),parameter :: sqrt_2 = dsqrt(2.d0),sqrt_2_inv = 1.d0/dsqrt(2.d0)
real(8) :: dx = 0.002,dx_inv = 500.d0
real(8),allocatable::sym_mat(:,:),sym_mat_transpose(:,:)
character(len=200)::potParamDir
!Numerical derivative array
real(8) :: five_pt_1st_deriv(0:4)
integer,allocatable::atom_ind(:,:)
real(8),parameter:: degree_to_radian  = 4.d0*datan(1.d0)/180.d0,radian_to_degree = 180.d0 /4.d0/datan(1.d0)
integer,parameter,dimension(0:7,0:3)::sym_pair = reshape((/ 0,-1,-1,-1,1,2,-1,-1,3,4,5,6,7,8,-1,-1,9,10,-1,-1,11,12,13,14,15,16,-1,-1,18,19,20,21/),shape(sym_pair),order=(/2,1/))
!integer,dimension(0:7,0:3)::sym_pair != [[0,-1,-1,-1],[1,2,-1,-1],[3,4,5,6],[7,8,-1,-1],[9,10,-1,-1],[11,12,13,14],[15,16,-1,-1],[18,19,20,21]] 
integer,parameter,dimension(0:7)::es_fn_idx = [0,1,3,7,9,11,15,20]
integer,parameter,dimension(0:6)::ea_fn_idx =[2,6,8,10,14,16,19]
integer,parameter,dimension(0:2)::os_fn_idx =  [5,13,18]
integer,parameter,dimension(0:2)::oa_fn_idx =[4,12,21]
integer,parameter,dimension(0:23)::fn_form  = [0,0,1,0,3,2,1,0,1,0,1,0,3,2,1,0,1,0,2,1,0,3,0,0]
integer,allocatable::max_id(:)
!stride_arr
integer,allocatable::stride_arr_quad_es(:,:),stride_arr_quad_os(:,:),stride_arr_quad_ea(:,:),stride_arr_quad_oa(:,:)

integer,allocatable::stride_arr_qubic_es_gt5(:,:),stride_arr_qubic_es_lt5(:,:),stride_arr_qubic_es_lt2(:,:),stride_arr_qubic_ea_gt5(:,:),stride_arr_qubic_ea_lt5(:,:),stride_arr_qubic_ea_lt2(:,:)
integer,allocatable::stride_arr_qubic_os_gt5(:,:),stride_arr_qubic_os_lt5(:,:),stride_arr_qubic_os_lt2(:,:),stride_arr_qubic_oa_gt5(:,:),stride_arr_qubic_oa_lt5(:,:),stride_arr_qubic_oa_lt2(:,:)
integer,allocatable::stride_arr_quartic_es_gt5(:,:),stride_arr_quartic_es_lt5(:,:),stride_arr_quartic_es_lt2(:,:),stride_arr_quartic_ea_gt5(:,:),stride_arr_quartic_ea_lt5(:,:),stride_arr_quartic_ea_lt2(:,:)
integer,allocatable::stride_arr_quartic_os_gt5(:,:),stride_arr_quartic_os_lt5(:,:),stride_arr_quartic_os_lt2(:,:) ,stride_arr_quartic_oa_gt5(:,:),stride_arr_quartic_oa_lt5(:,:),stride_arr_quartic_oa_lt2(:,:) 
integer,allocatable::stride_arr_pot(:,:)
integer,allocatable::stride_arr_b_es(:,:),stride_arr_b_os(:,:),stride_arr_b_ea(:,:),stride_arr_b_oa(:,:)
integer,allocatable::stride_arr_a_es(:,:),stride_arr_a_os(:,:),stride_arr_a_ea(:,:),stride_arr_a_oa(:,:)
integer,allocatable::stride_arr_d_es(:,:),stride_arr_d_os(:,:),stride_arr_d_ea(:,:),stride_arr_d_oa(:,:)
!fn_order
integer::fn_order_fij_es,fn_order_fij_os,fn_order_fij_ea,fn_order_fij_oa
integer::fn_order_fijk_es_gt5,fn_order_fijk_es_lt5,fn_order_fijk_es_lt2,fn_order_fijk_ea_gt5,fn_order_fijk_ea_lt5,fn_order_fijk_ea_lt2
integer::fn_order_fijk_os_gt5,fn_order_fijk_os_lt5,fn_order_fijk_os_lt2,fn_order_fijk_oa_gt5,fn_order_fijk_oa_lt5,fn_order_fijk_oa_lt2
integer::fn_order_fijkl_es_gt5,fn_order_fijkl_es_lt5,fn_order_fijkl_es_lt2,fn_order_fijkl_ea_gt5,fn_order_fijkl_ea_lt5,fn_order_fijkl_ea_lt2
integer::fn_order_fijkl_os_gt5,fn_order_fijkl_os_lt5,fn_order_fijkl_os_lt2,fn_order_fijkl_oa_gt5,fn_order_fijkl_oa_lt5,fn_order_fijkl_oa_lt2
integer:: fn_order_fijkl_es_ltpt1,fn_order_fijkl_ea_ltpt1,fn_order_fijkl_os_ltpt1,fn_order_fijkl_oa_ltpt1
integer ::fn_order_pot
integer ::fn_order_bonds_ea,fn_order_angs_ea,fn_order_dihs_ea,fn_order_bonds_oa,fn_order_angs_oa,fn_order_dihs_oa
integer ::fn_order_bonds_es,fn_order_angs_es,fn_order_dihs_es,fn_order_bonds_os,fn_order_angs_os,fn_order_dihs_os
!n_params
integer::n_fcs_quad_es,n_fcs_quad_os,n_fcs_quad_ea,n_fcs_quad_oa
integer::n_fcs_qubic_ea_gt5,n_fcs_qubic_ea_lt5,n_fcs_qubic_ea_lt2,n_fcs_qubic_es_gt5,n_fcs_qubic_es_lt5,n_fcs_qubic_es_lt2
integer::n_fcs_qubic_oa_gt5,n_fcs_qubic_oa_lt5,n_fcs_qubic_oa_lt2,n_fcs_qubic_os_gt5,n_fcs_qubic_os_lt5,n_fcs_qubic_os_lt2
integer::n_fcs_quartic_es_gt5,n_fcs_quartic_es_lt5,n_fcs_quartic_es_lt2,n_fcs_quartic_ea_gt5,n_fcs_quartic_ea_lt5,n_fcs_quartic_ea_lt2
integer::n_fcs_quartic_os_gt5,n_fcs_quartic_os_lt5,n_fcs_quartic_os_lt2,n_fcs_quartic_oa_gt5,n_fcs_quartic_oa_lt5,n_fcs_quartic_oa_lt2
integer::n_fcs_quartic_es_ltpt1,n_fcs_quartic_ea_ltpt1,n_fcs_quartic_os_ltpt1,n_fcs_quartic_oa_ltpt1
!Ncoeff
integer :: Ncoeff_pot(0:4)
integer ::Ncoeff_bonds_es(0:4),Ncoeff_bonds_os(0:4),Ncoeff_angs_es(0:4),Ncoeff_angs_os(0:4),Ncoeff_dihs_es(0:4),Ncoeff_dihs_os(0:4) 
integer ::Ncoeff_bonds_ea(0:4),Ncoeff_bonds_oa(0:4),Ncoeff_angs_ea(0:4),Ncoeff_angs_oa(0:4),Ncoeff_dihs_ea(0:4),Ncoeff_dihs_oa(0:4) 
integer ::Ncoeff_fij_es(0:4),Ncoeff_fij_os(0:4),Ncoeff_fij_ea(0:4),Ncoeff_fij_oa(0:4)
integer ::Ncoeff_fijk_es_gt5(0:4),Ncoeff_fijk_es_lt5(0:4),Ncoeff_fijk_es_lt2(0:4),Ncoeff_fijk_ea_gt5(0:4),Ncoeff_fijk_ea_lt5(0:4),Ncoeff_fijk_ea_lt2(0:4)
integer ::Ncoeff_fijk_os_gt5(0:4),Ncoeff_fijk_os_lt5(0:4),Ncoeff_fijk_os_lt2(0:4),Ncoeff_fijk_oa_gt5(0:4),Ncoeff_fijk_oa_lt5(0:4),Ncoeff_fijk_oa_lt2(0:4)
integer ::Ncoeff_fijkl_es_gt5(0:4),Ncoeff_fijkl_es_lt5(0:4),Ncoeff_fijkl_es_lt2(0:4),Ncoeff_fijkl_ea_gt5(0:4),Ncoeff_fijkl_ea_lt5(0:4),Ncoeff_fijkl_ea_lt2(0:4)
integer ::Ncoeff_fijkl_os_gt5(0:4),Ncoeff_fijkl_os_lt5(0:4),Ncoeff_fijkl_os_lt2(0:4) ,Ncoeff_fijkl_oa_gt5(0:4),Ncoeff_fijkl_oa_lt5(0:4),Ncoeff_fijkl_oa_lt2(0:4) 




integer,allocatable::fc_idx_quad_es(:,:),fc_idx_quad_os(:,:),fc_idx_quad_ea(:,:),fc_idx_quad_oa(:,:)
integer,allocatable::fc_idx_qubic_es_gt5(:,:),fc_idx_qubic_es_lt5(:,:),fc_idx_qubic_es_lt2(:,:),fc_idx_qubic_ea_gt5(:,:),fc_idx_qubic_ea_lt5(:,:),fc_idx_qubic_ea_lt2(:,:)
integer,allocatable::fc_idx_qubic_os_gt5(:,:),fc_idx_qubic_os_lt5(:,:),fc_idx_qubic_os_lt2(:,:),fc_idx_qubic_oa_gt5(:,:),fc_idx_qubic_oa_lt5(:,:),fc_idx_qubic_oa_lt2(:,:)
integer,allocatable::fc_idx_quartic_es_gt5(:,:),fc_idx_quartic_es_lt5(:,:),fc_idx_quartic_es_lt2(:,:),fc_idx_quartic_ea_gt5(:,:),fc_idx_quartic_ea_lt5(:,:),fc_idx_quartic_ea_lt2(:,:)
integer,allocatable::fc_idx_quartic_os_gt5(:,:),fc_idx_quartic_os_lt5(:,:),fc_idx_quartic_os_lt2(:,:) , fc_idx_quartic_oa_gt5(:,:),fc_idx_quartic_oa_lt5(:,:),fc_idx_quartic_oa_lt2(:,:) 
integer,allocatable::fc_idx_quartic_es_ltpt1(:,:),fc_idx_quartic_ea_ltpt1(:,:),fc_idx_quartic_os_ltpt1(:,:),fc_idx_quartic_oa_ltpt1(:,:)
!pot coeff
real(8),allocatable::Vt_coeff(:)
!fitted internal coordinates parameter

real(8),allocatable :: fitted_bonds_coeff_es(:,:),fitted_bonds_coeff_os(:,:),fitted_bonds_coeff_ea(:,:),fitted_bonds_coeff_oa(:,:)  
real(8),allocatable :: fitted_angs_coeff_es(:,:),fitted_angs_coeff_os(:,:) , fitted_angs_coeff_ea(:,:),fitted_angs_coeff_oa(:,:)
real(8),allocatable :: fitted_dihs_coeff_es(:,:),fitted_dihs_coeff_os(:,:) , fitted_dihs_coeff_ea(:,:),fitted_dihs_coeff_oa(:,:) 
!force const 
real(8),allocatable::fitted_fc_coeff_ij_es(:,:),fitted_fc_coeff_ij_os(:,:),fitted_fc_coeff_ij_ea(:,:),fitted_fc_coeff_ij_oa(:,:)
real(8),allocatable::fitted_fc_coeff_ijk_es_gt5(:,:),fitted_fc_coeff_ijk_es_lt5(:,:),fitted_fc_coeff_ijk_es_lt2(:,:),fitted_fc_coeff_ijk_ea_gt5(:,:),fitted_fc_coeff_ijk_ea_lt5(:,:),fitted_fc_coeff_ijk_ea_lt2(:,:)
real(8),allocatable::fitted_fc_coeff_ijk_os_gt5(:,:),fitted_fc_coeff_ijk_os_lt5(:,:),fitted_fc_coeff_ijk_os_lt2(:,:),fitted_fc_coeff_ijk_oa_gt5(:,:),fitted_fc_coeff_ijk_oa_lt5(:,:),fitted_fc_coeff_ijk_oa_lt2(:,:)
real(8),allocatable::fitted_fc_coeff_ijkl_es_gt5(:,:),fitted_fc_coeff_ijkl_es_lt5(:,:),fitted_fc_coeff_ijkl_es_lt2(:,:),fitted_fc_coeff_ijkl_ea_gt5(:,:),fitted_fc_coeff_ijkl_ea_lt5(:,:),fitted_fc_coeff_ijkl_ea_lt2(:,:)
real(8),allocatable::fitted_fc_coeff_ijkl_os_gt5(:,:),fitted_fc_coeff_ijkl_os_lt5(:,:),fitted_fc_coeff_ijkl_os_lt2(:,:), fitted_fc_coeff_ijkl_oa_gt5(:,:), fitted_fc_coeff_ijkl_oa_lt5(:,:),fitted_fc_coeff_ijkl_oa_lt2(:,:) 
real(8),allocatable::fitted_fc_coeff_ijkl_es_ltpt1(:),fitted_fc_coeff_ijkl_ea_ltpt1(:),fitted_fc_coeff_ijkl_os_ltpt1(:),fitted_fc_coeff_ijkl_oa_ltpt1(:)
real(8),allocatable:: dim_scal_factor(:),dim_scal_factor_inv(:)
contains

subroutine initializePotvars()
	implicit none 
	integer::i,ii,ctr,id1,id2,j,k,l
	Natoms = 10;Nbonds = 9;N_int =24;Nang = 8;Ncarts = 30
	allocate(dim_scal_factor(0:N_int-1),dim_scal_factor_inv(0:N_int-1))
	allocate(sym_mat(0:N_int-1,0:N_int-1),sym_mat_transpose(0:N_int-1,0:N_int-1))
	sym_mat_transpose = 0.d0
	five_pt_1st_deriv = [1.d0,-8.d0,0.d0,8.d0,-1.d0]
	five_pt_1st_deriv = five_pt_1st_deriv/12.d0
	!sym_idx      = [3,4,5,6,10,11,12,13,16,17,18,19] 
	sym_mat      = 0.d0
	open(unit =13,file ="pot.param",status="old")
	read(13,*)potParamDir
	!print*,'mass=',mass
	!do i = 0,Natoms-1
	!	print*,i
	!	cartmass(3*i:3*i+2) = mass(i)
	!enddo
	!Read_scale_factor
	open(unit =11,file ="scale_factor_eg_avtz.dat",status="old") 
	read(11,*)(dim_scal_factor(j),j=0,N_int-1)
	do i = 0,N_int-1
		dim_scal_factor_inv(i) = 1.d0 / dim_scal_factor(i)
	enddo
	fn_order_fijkl_es_ltpt1 = 0;
	fn_order_fijkl_ea_ltpt1 = 0;
	fn_order_fijkl_os_ltpt1 = 0;
	fn_order_fijkl_oa_ltpt1 = 0;
	!****************************************************************************************#
	!Parameters for Minimum energy path along OCCF and HOCC diherdrals
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_pot,Ncoeff_pot)
	allocate(stride_arr_pot(0:fn_order_pot-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_pot,stride_arr_pot)
	!********************************For bonds********************************#
	!Even symmetric expansion
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_bonds_es,Ncoeff_bonds_es)
	allocate(stride_arr_b_es(0:fn_order_bonds_es-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_bonds_es,stride_arr_b_es)
	!Even asymmetric expansion
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_bonds_ea,Ncoeff_bonds_ea)
	allocate(stride_arr_b_ea(0:fn_order_bonds_ea-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_bonds_ea,stride_arr_b_ea)
	!Odd symmetric expansion 
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_bonds_os,Ncoeff_bonds_os)
	allocate(stride_arr_b_os(0:fn_order_bonds_os-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_bonds_os,stride_arr_b_os)
	!Odd asymmetric expansion 
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_bonds_oa,Ncoeff_bonds_oa)
	allocate(stride_arr_b_oa(0:fn_order_bonds_oa-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_bonds_oa,stride_arr_b_oa)
	!********************************For Angles*******************************#
	!Even symmetric  expansion
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_angs_es,Ncoeff_angs_es)
	allocate(stride_arr_a_es(0:fn_order_angs_es-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_angs_es,stride_arr_a_es)
	!Even asymmetric  expansion
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_angs_ea,Ncoeff_angs_ea)
	allocate(stride_arr_a_ea(0:fn_order_angs_ea-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_angs_ea,stride_arr_a_ea)
	!Odd  symmetric expansion 
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_angs_os,Ncoeff_angs_os )
	allocate(stride_arr_a_os(0:fn_order_angs_os-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_angs_os,stride_arr_a_os)
	!Odd asymmetric  expansion 
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_angs_oa,Ncoeff_angs_oa )
	allocate(stride_arr_a_oa(0:fn_order_angs_oa-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_angs_oa,stride_arr_a_oa)
	!********************************For Dihedrals*******************************#
	!Even symmetric  expansion
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_dihs_es,Ncoeff_dihs_es)
	allocate(stride_arr_d_es(0:fn_order_dihs_es-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_dihs_es,stride_arr_d_es)
	!Even asymmetric  expansion
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_dihs_ea,Ncoeff_dihs_ea)
	allocate(stride_arr_d_ea(0:fn_order_dihs_ea-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_dihs_ea,stride_arr_d_ea)
	!Odd  symmetric expansion 
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_dihs_os,Ncoeff_dihs_os)
	allocate(stride_arr_d_os(0:fn_order_dihs_os-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_dihs_os,stride_arr_d_os)
	!Odd  asymmetric expansion 
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_dihs_oa,Ncoeff_dihs_oa)
	allocate(stride_arr_d_oa(0:fn_order_dihs_oa-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_dihs_oa,stride_arr_d_oa)
	!**************************************************************************#
	!Quadratic force constants
	!Even symmetric  
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_fij_es,Ncoeff_fij_es)
	allocate(stride_arr_quad_es(0:fn_order_fij_es-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_fij_es,stride_arr_quad_es)
	!Even asymmetric  
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_fij_ea,Ncoeff_fij_ea)
	allocate(stride_arr_quad_ea(0:fn_order_fij_ea-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_fij_ea,stride_arr_quad_ea)
	!Odd  symmetric   
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_fij_os,Ncoeff_fij_os)
	allocate(stride_arr_quad_os(0:fn_order_fij_os-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_fij_os,stride_arr_quad_os)
	!Odd asymmetric  
	Mmax = 8;Nmax = 6;Lmax = 6
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_fij_oa,Ncoeff_fij_oa)
	allocate(stride_arr_quad_oa(0:fn_order_fij_oa-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_fij_oa,stride_arr_quad_oa)
	!****************************************************************************************************************************************#
	!Cubic force constants
	!Even  symmetric  
	!One Body 
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_fijk_es_gt5,Ncoeff_fijk_es_gt5)
	allocate(stride_arr_qubic_es_gt5(0:fn_order_fijk_es_gt5-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_fijk_es_gt5,stride_arr_qubic_es_gt5)
	!Two Body 
	Mmax = 6;Nmax = 5;Lmax = 5
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_fijk_es_lt5,Ncoeff_fijk_es_lt5)
	allocate(stride_arr_qubic_es_lt5(0:fn_order_fijk_es_lt5-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_fijk_es_lt5,stride_arr_qubic_es_lt5)
	!Three Body
	Mmax = 5;Nmax = 4;Lmax = 4
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_fijk_es_lt2,Ncoeff_fijk_es_lt2)
	allocate(stride_arr_qubic_es_lt2(0:fn_order_fijk_es_lt2-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_fijk_es_lt2,stride_arr_qubic_es_lt2)
	!Even asymmetric  
	!One Body 
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_fijk_ea_gt5,Ncoeff_fijk_ea_gt5)
	allocate(stride_arr_qubic_ea_gt5(0:fn_order_fijk_ea_gt5-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_fijk_ea_gt5,stride_arr_qubic_ea_gt5)
	!Two Body 
	Mmax = 6;Nmax = 5;Lmax = 5
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_fijk_ea_lt5,Ncoeff_fijk_ea_lt5)
	allocate(stride_arr_qubic_ea_lt5(0:fn_order_fijk_ea_lt5-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_fijk_ea_lt5,stride_arr_qubic_ea_lt5)
	!Three Body
	Mmax = 5;Nmax = 4;Lmax = 4
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_fijk_ea_lt2,Ncoeff_fijk_ea_lt2)
	allocate(stride_arr_qubic_ea_lt2(0:fn_order_fijk_ea_lt2-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_fijk_ea_lt2,stride_arr_qubic_ea_lt2)
	!Odd symmetric  
	!One Body 
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_fijk_os_gt5,Ncoeff_fijk_os_gt5)
	allocate(stride_arr_qubic_os_gt5(0:fn_order_fijk_os_gt5-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_fijk_os_gt5,stride_arr_qubic_os_gt5)
	!Two Body 
	Mmax = 6;Nmax = 5;Lmax = 5
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_fijk_os_lt5,Ncoeff_fijk_os_lt5)
	allocate(stride_arr_qubic_os_lt5(0:fn_order_fijk_os_lt5-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_fijk_os_lt5,stride_arr_qubic_os_lt5)
	!Three Body
	Mmax = 5;Nmax = 4;Lmax = 4
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_fijk_os_lt2,Ncoeff_fijk_os_lt2)
	allocate(stride_arr_qubic_os_lt2(0:fn_order_fijk_os_lt2-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_fijk_os_lt2,stride_arr_qubic_os_lt2)
	!Odd asymmetric  
	!One Body 
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_fijk_oa_gt5,Ncoeff_fijk_oa_gt5)
	allocate(stride_arr_qubic_oa_gt5(0:fn_order_fijk_oa_gt5-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_fijk_oa_gt5,stride_arr_qubic_oa_gt5)
	!Two Body 
	Mmax = 6;Nmax = 5;Lmax = 5
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_fijk_oa_lt5,Ncoeff_fijk_oa_lt5)
	allocate(stride_arr_qubic_oa_lt5(0:fn_order_fijk_oa_lt5-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_fijk_oa_lt5,stride_arr_qubic_oa_lt5)
	!Three Body
	Mmax = 5;Nmax = 4;Lmax = 4
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_fijk_oa_lt2,Ncoeff_fijk_oa_lt2)
	allocate(stride_arr_qubic_oa_lt2(0:fn_order_fijk_oa_lt2-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_fijk_oa_lt2,stride_arr_qubic_oa_lt2)
	!****************************************************************************************************************************************#
	!Quartic force constants
	!Even symmetric  
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_fijkl_es_gt5,Ncoeff_fijkl_es_gt5)
	allocate(stride_arr_quartic_es_gt5(0:fn_order_fijkl_es_gt5-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_fijkl_es_gt5,stride_arr_quartic_es_gt5)
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_fijkl_es_lt5,Ncoeff_fijkl_es_lt5)
	allocate(stride_arr_quartic_es_lt5(0:fn_order_fijkl_es_lt5-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_fijkl_es_lt5,stride_arr_quartic_es_lt5)
	Mmax = 5;Nmax = 4;Lmax = 4
	call get_fn_order_es(Mmax,Nmax,Lmax,fn_order_fijkl_es_lt2,Ncoeff_fijkl_es_lt2)
	allocate(stride_arr_quartic_es_lt2(0:fn_order_fijkl_es_lt2-1,0:2))
	call get_stride_arr_es(Mmax,Nmax,Lmax,fn_order_fijkl_es_lt2,stride_arr_quartic_es_lt2)
	!Quartic force constants
	!Even asymmetric  
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_fijkl_ea_gt5,Ncoeff_fijkl_ea_gt5)
	allocate(stride_arr_quartic_ea_gt5(0:fn_order_fijkl_ea_gt5-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_fijkl_ea_gt5,stride_arr_quartic_ea_gt5)
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_fijkl_ea_lt5,Ncoeff_fijkl_ea_lt5)
	allocate(stride_arr_quartic_ea_lt5(0:fn_order_fijkl_ea_lt5-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_fijkl_ea_lt5,stride_arr_quartic_ea_lt5)
	Mmax = 5;Nmax = 4;Lmax = 4
	call get_fn_order_ea(Mmax,Nmax,Lmax,fn_order_fijkl_ea_lt2,Ncoeff_fijkl_ea_lt2)
	allocate(stride_arr_quartic_ea_lt2(0:fn_order_fijkl_ea_lt2-1,0:2))
	call get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order_fijkl_ea_lt2,stride_arr_quartic_ea_lt2)
	!Odd symmetric  
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_fijkl_os_gt5,Ncoeff_fijkl_os_gt5)
	allocate(stride_arr_quartic_os_gt5(0:fn_order_fijkl_os_gt5-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_fijkl_os_gt5,stride_arr_quartic_os_gt5)
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_fijkl_os_lt5,Ncoeff_fijkl_os_lt5)
	allocate(stride_arr_quartic_os_lt5(0:fn_order_fijkl_os_lt5-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_fijkl_os_lt5,stride_arr_quartic_os_lt5)
	Mmax = 5;Nmax = 4;Lmax = 4
	call get_fn_order_os(Mmax,Nmax,Lmax,fn_order_fijkl_os_lt2,Ncoeff_fijkl_os_lt2)
	allocate(stride_arr_quartic_os_lt2(0:fn_order_fijkl_os_lt2-1,0:2))
	call get_stride_arr_os(Mmax,Nmax,Lmax,fn_order_fijkl_os_lt2,stride_arr_quartic_os_lt2)
	!Odd asymmetric  
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_fijkl_oa_gt5,Ncoeff_fijkl_oa_gt5)
	allocate(stride_arr_quartic_oa_gt5(0:fn_order_fijkl_oa_gt5-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_fijkl_oa_gt5,stride_arr_quartic_oa_gt5)
	Mmax = 7;Nmax = 5;Lmax = 5
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_fijkl_oa_lt5,Ncoeff_fijkl_oa_lt5)
	allocate(stride_arr_quartic_oa_lt5(0:fn_order_fijkl_oa_lt5-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_fijkl_oa_lt5,stride_arr_quartic_oa_lt5)
	Mmax = 5;Nmax = 4;Lmax = 4
	call get_fn_order_oa(Mmax,Nmax,Lmax,fn_order_fijkl_oa_lt2,Ncoeff_fijkl_oa_lt2)
	allocate(stride_arr_quartic_oa_lt2(0:fn_order_fijkl_oa_lt2-1,0:2))
	call get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order_fijkl_oa_lt2,stride_arr_quartic_oa_lt2)
	!****************************************************************************************************************************************#
	open(unit=101,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quad_8_6_6_es_fortran.dat',status='old')
	read(101,*)n_fcs_quad_es
	allocate(fitted_fc_coeff_ij_es(0:fn_order_fij_es-1,0:n_fcs_quad_es-1),fc_idx_quad_es(0:n_fcs_quad_es-1,0:3))
	open(unit=102,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_gt5_7_5_5_es_fortran.dat',status='old')
	read(102,*)n_fcs_qubic_es_gt5
	allocate(fitted_fc_coeff_ijk_es_gt5(0:fn_order_fijk_es_gt5-1,0:n_fcs_qubic_es_gt5-1),fc_idx_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1,0:3))
	open(unit=103,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_lt5_7_5_5_es_fortran.dat',status='old')
	read(103,*)n_fcs_qubic_es_lt5
	allocate(fitted_fc_coeff_ijk_es_lt5(0:fn_order_fijk_es_lt5-1,0:n_fcs_qubic_es_lt5-1),fc_idx_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1,0:3))
	open(unit=104,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_lt2_5_4_4_es_fortran.dat',status='old')
	read(104,*)n_fcs_qubic_es_lt2
	allocate(fitted_fc_coeff_ijk_es_lt2(0:fn_order_fijk_es_lt2-1,0:n_fcs_qubic_es_lt2-1),fc_idx_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1,0:3))
	open(unit=105,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_gt5_7_5_5_es_fortran.dat',status='old')
	read(105,*)n_fcs_quartic_es_gt5
	allocate(fitted_fc_coeff_ijkl_es_gt5(0:fn_order_fijkl_es_gt5-1,0:n_fcs_quartic_es_gt5-1),fc_idx_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1,0:3))
	open(unit=106,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_lt5_7_5_5_es_fortran.dat',status='old')
	read(106,*)n_fcs_quartic_es_lt5
	allocate(fitted_fc_coeff_ijkl_es_lt5(0:fn_order_fijkl_es_lt5-1,0:n_fcs_quartic_es_lt5-1),fc_idx_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1,0:3))
	open(unit=107,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_lt2_5_4_4_es_fortran.dat',status='old')
	read(107,*)n_fcs_quartic_es_lt2
	allocate(fitted_fc_coeff_ijkl_es_lt2(0:fn_order_fijkl_es_lt2-1,0:n_fcs_quartic_es_lt2-1),fc_idx_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1,0:3))
	open(unit=108,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_ltpt1_es_fortran.dat',status='old')
	read(108,*)n_fcs_quartic_es_ltpt1
	allocate(fitted_fc_coeff_ijkl_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_idx_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1,0:3))
	open(unit=109,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quad_8_6_6_ea_fortran.dat',status='old')
	read(109,*)n_fcs_quad_ea
	allocate(fitted_fc_coeff_ij_ea(0:fn_order_fij_ea-1,0:n_fcs_quad_ea-1),fc_idx_quad_ea(0:n_fcs_quad_ea-1,0:3))
	open(unit=110,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_gt5_7_5_5_ea_fortran.dat',status='old')
	read(110,*)n_fcs_qubic_ea_gt5
	allocate(fitted_fc_coeff_ijk_ea_gt5(0:fn_order_fijk_ea_gt5-1,0:n_fcs_qubic_ea_gt5-1),fc_idx_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1,0:3))
	open(unit=111,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_lt5_7_5_5_ea_fortran.dat',status='old')
	read(111,*)n_fcs_qubic_ea_lt5
	allocate(fitted_fc_coeff_ijk_ea_lt5(0:fn_order_fijk_ea_lt5-1,0:n_fcs_qubic_ea_lt5-1),fc_idx_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1,0:3))
	open(unit=112,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_lt2_5_4_4_ea_fortran.dat',status='old')
	read(112,*)n_fcs_qubic_ea_lt2
	allocate(fitted_fc_coeff_ijk_ea_lt2(0:fn_order_fijk_ea_lt2-1,0:n_fcs_qubic_ea_lt2-1),fc_idx_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1,0:3))
	open(unit=113,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_gt5_7_5_5_ea_fortran.dat',status='old')
	read(113,*)n_fcs_quartic_ea_gt5
	allocate(fitted_fc_coeff_ijkl_ea_gt5(0:fn_order_fijkl_ea_gt5-1,0:n_fcs_quartic_ea_gt5-1),fc_idx_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1,0:3))
	open(unit=114,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_lt5_7_5_5_ea_fortran.dat',status='old')
	read(114,*)n_fcs_quartic_ea_lt5
	allocate(fitted_fc_coeff_ijkl_ea_lt5(0:fn_order_fijkl_ea_lt5-1,0:n_fcs_quartic_ea_lt5-1),fc_idx_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1,0:3))
	open(unit=115,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_lt2_5_4_4_ea_fortran.dat',status='old')
	read(115,*)n_fcs_quartic_ea_lt2
	allocate(fitted_fc_coeff_ijkl_ea_lt2(0:fn_order_fijkl_ea_lt2-1,0:n_fcs_quartic_ea_lt2-1),fc_idx_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1,0:3))
	open(unit=116,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_ltpt1_ea_fortran.dat',status='old')
	read(116,*)n_fcs_quartic_ea_ltpt1
	allocate(fitted_fc_coeff_ijkl_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1),fc_idx_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1,0:3))
	call Read_fc_coeff(101,N_int,fn_order_fij_es,4,n_fcs_quad_es,fc_idx_quad_es,fitted_fc_coeff_ij_es)
	call Read_fc_coeff(102,N_int,fn_order_fijk_es_gt5,4,n_fcs_qubic_es_gt5,fc_idx_qubic_es_gt5,fitted_fc_coeff_ijk_es_gt5)
	call Read_fc_coeff(103,N_int,fn_order_fijk_es_lt5,4,n_fcs_qubic_es_lt5,fc_idx_qubic_es_lt5,fitted_fc_coeff_ijk_es_lt5)
	call Read_fc_coeff(104,N_int,fn_order_fijk_es_lt2,4,n_fcs_qubic_es_lt2,fc_idx_qubic_es_lt2,fitted_fc_coeff_ijk_es_lt2)
	call Read_fc_coeff(105,N_int,fn_order_fijkl_es_gt5,4,n_fcs_quartic_es_gt5,fc_idx_quartic_es_gt5,fitted_fc_coeff_ijkl_es_gt5)
	call Read_fc_coeff(106,N_int,fn_order_fijkl_es_lt5,4,n_fcs_quartic_es_lt5,fc_idx_quartic_es_lt5,fitted_fc_coeff_ijkl_es_lt5)
	call Read_fc_coeff(107,N_int,fn_order_fijkl_es_lt2,4,n_fcs_quartic_es_lt2,fc_idx_quartic_es_lt2,fitted_fc_coeff_ijkl_es_lt2)
	call Read_fc_coeff_ltpt1(108,N_int,4,n_fcs_quartic_es_ltpt1,fc_idx_quartic_es_ltpt1,fitted_fc_coeff_ijkl_es_ltpt1)
	call Read_fc_coeff(109,N_int,fn_order_fij_ea,4,n_fcs_quad_ea,fc_idx_quad_ea,fitted_fc_coeff_ij_ea)
	call Read_fc_coeff(110,N_int,fn_order_fijk_ea_gt5,4,n_fcs_qubic_ea_gt5,fc_idx_qubic_ea_gt5,fitted_fc_coeff_ijk_ea_gt5)
	call Read_fc_coeff(111,N_int,fn_order_fijk_ea_lt5,4,n_fcs_qubic_ea_lt5,fc_idx_qubic_ea_lt5,fitted_fc_coeff_ijk_ea_lt5)
	call Read_fc_coeff(112,N_int,fn_order_fijk_ea_lt2,4,n_fcs_qubic_ea_lt2,fc_idx_qubic_ea_lt2,fitted_fc_coeff_ijk_ea_lt2)
	call Read_fc_coeff(113,N_int,fn_order_fijkl_ea_gt5,4,n_fcs_quartic_ea_gt5,fc_idx_quartic_ea_gt5,fitted_fc_coeff_ijkl_ea_gt5)
	call Read_fc_coeff(114,N_int,fn_order_fijkl_ea_lt5,4,n_fcs_quartic_ea_lt5,fc_idx_quartic_ea_lt5,fitted_fc_coeff_ijkl_ea_lt5)
	call Read_fc_coeff(115,N_int,fn_order_fijkl_ea_lt2,4,n_fcs_quartic_ea_lt2,fc_idx_quartic_ea_lt2,fitted_fc_coeff_ijkl_ea_lt2)
	call Read_fc_coeff_ltpt1(116,N_int,4,n_fcs_quartic_ea_ltpt1,fc_idx_quartic_ea_ltpt1,fitted_fc_coeff_ijkl_ea_ltpt1)

	open(unit=201,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quad_8_6_6_os_fortran.dat',status='old')
	read(201,*)n_fcs_quad_os
	allocate(fitted_fc_coeff_ij_os(0:fn_order_fij_os-1,0:n_fcs_quad_os-1),fc_idx_quad_os(0:n_fcs_quad_os-1,0:3))

	open(unit=202,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_gt5_7_5_5_os_fortran.dat',status='old')
	read(202,*)n_fcs_qubic_os_gt5
	allocate(fitted_fc_coeff_ijk_os_gt5(0:fn_order_fijk_os_gt5-1,0:n_fcs_qubic_os_gt5-1),fc_idx_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1,0:3))

	open(unit=203,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_lt5_7_5_5_os_fortran.dat',status='old')
	read(203,*)n_fcs_qubic_os_lt5
	allocate(fitted_fc_coeff_ijk_os_lt5(0:fn_order_fijk_os_lt5-1,0:n_fcs_qubic_os_lt5-1),fc_idx_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1,0:3))

	open(unit=204,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_lt2_5_4_4_os_fortran.dat',status='old')
	read(204,*)n_fcs_qubic_os_lt2
	allocate(fitted_fc_coeff_ijk_os_lt2(0:fn_order_fijk_os_lt2-1,0:n_fcs_qubic_os_lt2-1),fc_idx_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1,0:3))
	open(unit=205,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_gt5_7_5_5_os_fortran.dat',status='old')
	read(205,*)n_fcs_quartic_os_gt5
	allocate(fitted_fc_coeff_ijkl_os_gt5(0:fn_order_fijkl_os_gt5-1,0:n_fcs_quartic_os_gt5-1),fc_idx_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1,0:3))
	open(unit=206,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_lt5_7_5_5_os_fortran.dat',status='old')
	read(206,*)n_fcs_quartic_os_lt5
	allocate(fitted_fc_coeff_ijkl_os_lt5(0:fn_order_fijkl_os_lt5-1,0:n_fcs_quartic_os_lt5-1),fc_idx_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1,0:3))
	open(unit=207,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_lt2_5_4_4_os_fortran.dat',status='old')
	read(207,*)n_fcs_quartic_os_lt2
	allocate(fitted_fc_coeff_ijkl_os_lt2(0:fn_order_fijkl_os_lt2-1,0:n_fcs_quartic_os_lt2-1),fc_idx_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1,0:3))

	open(unit=208,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_ltpt1_os_fortran.dat',status='old')
	read(208,*)n_fcs_quartic_os_ltpt1
	allocate(fitted_fc_coeff_ijkl_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_idx_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1,0:3))

	open(unit=209,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quad_8_6_6_oa_fortran.dat',status='old')
	read(209,*)n_fcs_quad_oa
	allocate(fitted_fc_coeff_ij_oa(0:fn_order_fij_oa-1,0:n_fcs_quad_oa-1),fc_idx_quad_oa(0:n_fcs_quad_oa-1,0:3))

	open(unit=210,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_gt5_7_5_5_oa_fortran.dat',status='old')
	read(210,*)n_fcs_qubic_oa_gt5
	allocate(fitted_fc_coeff_ijk_oa_gt5(0:fn_order_fijk_oa_gt5-1,0:n_fcs_qubic_oa_gt5-1),fc_idx_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1,0:3))

	open(unit=211,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_lt5_7_5_5_oa_fortran.dat',status='old')
	read(211,*)n_fcs_qubic_oa_lt5
	allocate(fitted_fc_coeff_ijk_oa_lt5(0:fn_order_fijk_oa_lt5-1,0:n_fcs_qubic_oa_lt5-1),fc_idx_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1,0:3))

	open(unit=212,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_qubic_lt2_5_4_4_oa_fortran.dat',status='old')
	read(212,*)n_fcs_qubic_oa_lt2
	allocate(fitted_fc_coeff_ijk_oa_lt2(0:fn_order_fijk_oa_lt2-1,0:n_fcs_qubic_oa_lt2-1),fc_idx_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1,0:3))
	open(unit=213,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_gt5_7_5_5_oa_fortran.dat',status='old')
	read(213,*)n_fcs_quartic_oa_gt5
	allocate(fitted_fc_coeff_ijkl_oa_gt5(0:fn_order_fijkl_oa_gt5-1,0:n_fcs_quartic_oa_gt5-1),fc_idx_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1,0:3))
	open(unit=214,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_lt5_7_5_5_oa_fortran.dat',status='old')
	read(214,*)n_fcs_quartic_oa_lt5
	allocate(fitted_fc_coeff_ijkl_oa_lt5(0:fn_order_fijkl_oa_lt5-1,0:n_fcs_quartic_oa_lt5-1),fc_idx_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1,0:3))
	open(unit=215,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_lt2_5_4_4_oa_fortran.dat',status='old')
	read(215,*)n_fcs_quartic_oa_lt2
	allocate(fitted_fc_coeff_ijkl_oa_lt2(0:fn_order_fijkl_oa_lt2-1,0:n_fcs_quartic_oa_lt2-1),fc_idx_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1,0:3))

	open(unit=216,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/fitted_param_quartic_ltpt1_oa_fortran.dat',status='old')
	read(216,*)n_fcs_quartic_oa_ltpt1
	allocate(fitted_fc_coeff_ijkl_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1),fc_idx_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1,0:3))

	call Read_fc_coeff(201,N_int,fn_order_fij_os,4,n_fcs_quad_os,fc_idx_quad_os,fitted_fc_coeff_ij_os)
	call Read_fc_coeff(202,N_int,fn_order_fijk_os_gt5,4,n_fcs_qubic_os_gt5,fc_idx_qubic_os_gt5,fitted_fc_coeff_ijk_os_gt5)
	call Read_fc_coeff(203,N_int,fn_order_fijk_os_lt5,4,n_fcs_qubic_os_lt5,fc_idx_qubic_os_lt5,fitted_fc_coeff_ijk_os_lt5)
	call Read_fc_coeff(204,N_int,fn_order_fijk_os_lt2,4,n_fcs_qubic_os_lt2,fc_idx_qubic_os_lt2,fitted_fc_coeff_ijk_os_lt2)
	call Read_fc_coeff(205,N_int,fn_order_fijkl_os_gt5,4,n_fcs_quartic_os_gt5,fc_idx_quartic_os_gt5,fitted_fc_coeff_ijkl_os_gt5)
	call Read_fc_coeff(206,N_int,fn_order_fijkl_os_lt5,4,n_fcs_quartic_os_lt5,fc_idx_quartic_os_lt5,fitted_fc_coeff_ijkl_os_lt5)
	call Read_fc_coeff(207,N_int,fn_order_fijkl_os_lt2,4,n_fcs_quartic_os_lt2,fc_idx_quartic_os_lt2,fitted_fc_coeff_ijkl_os_lt2)
	call Read_fc_coeff_ltpt1(208,N_int,4,n_fcs_quartic_os_ltpt1,fc_idx_quartic_os_ltpt1,fitted_fc_coeff_ijkl_os_ltpt1)
	call Read_fc_coeff(209,N_int,fn_order_fij_oa,4,n_fcs_quad_oa,fc_idx_quad_oa,fitted_fc_coeff_ij_oa)
	call Read_fc_coeff(210,N_int,fn_order_fijk_oa_gt5,4,n_fcs_qubic_oa_gt5,fc_idx_qubic_oa_gt5,fitted_fc_coeff_ijk_oa_gt5)
	call Read_fc_coeff(211,N_int,fn_order_fijk_oa_lt5,4,n_fcs_qubic_oa_lt5,fc_idx_qubic_oa_lt5,fitted_fc_coeff_ijk_oa_lt5)
	call Read_fc_coeff(212,N_int,fn_order_fijk_oa_lt2,4,n_fcs_qubic_oa_lt2,fc_idx_qubic_oa_lt2,fitted_fc_coeff_ijk_oa_lt2)
	call Read_fc_coeff(213,N_int,fn_order_fijkl_oa_gt5,4,n_fcs_quartic_oa_gt5,fc_idx_quartic_oa_gt5,fitted_fc_coeff_ijkl_oa_gt5)
	call Read_fc_coeff(214,N_int,fn_order_fijkl_oa_lt5,4,n_fcs_quartic_oa_lt5,fc_idx_quartic_oa_lt5,fitted_fc_coeff_ijkl_oa_lt5)
	call Read_fc_coeff(215,N_int,fn_order_fijkl_oa_lt2,4,n_fcs_quartic_oa_lt2,fc_idx_quartic_oa_lt2,fitted_fc_coeff_ijkl_oa_lt2)
	call Read_fc_coeff_ltpt1(216,N_int,4,n_fcs_quartic_oa_ltpt1,fc_idx_quartic_oa_ltpt1,fitted_fc_coeff_ijkl_oa_ltpt1)
	!****************************************************************************************#
	!MEP Path fitted parameters
	allocate(Vt_coeff(0:fn_order_pot-1))
	call Read_Vt_coeff(trim(adjustl(potParamDir))//"/fortran_param_files_avtz/pot_par.dat",fn_order_pot,Vt_coeff )
	!Vt_coeff = In Hartee
	!Symmetrized internal coordinate fit
	open(unit=302,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/eq_par_8_6_6_es_fortran.dat',status='old')
	read(302,*)n_eb,n_ea,n_ed
    allocate(fitted_bonds_coeff_es(0:fn_order_bonds_es-1,0:n_eb-1),fitted_angs_coeff_es(0:fn_order_angs_es-1,0:n_ea-1),fitted_dihs_coeff_es(0:fn_order_dihs_es-1,0:n_ed-1))
    call   Read_eq_int_coeff(302,n_eb,n_ea,n_ed,Nbonds,Nang,N_int,fn_order_bonds_es,fn_order_angs_es,fn_order_dihs_es,es_fn_idx,fitted_bonds_coeff_es,fitted_angs_coeff_es,fitted_dihs_coeff_es)
	!Symmetrized internal coordinate fit
	open(unit=303,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/eq_par_8_6_6_ea_fortran.dat',status='old')
	read(303,*)n_eb,n_ea,n_ed
    allocate(fitted_bonds_coeff_ea(0:fn_order_bonds_ea-1,0:n_eb-1),fitted_angs_coeff_ea(0:fn_order_angs_ea-1,0:n_ea-1),fitted_dihs_coeff_ea(0:fn_order_dihs_ea-1,0:n_ed-1))
    call   Read_eq_int_coeff(303,n_eb,n_ea,n_ed,Nbonds,Nang,N_int,fn_order_bonds_ea,fn_order_angs_ea,fn_order_dihs_ea,ea_fn_idx,fitted_bonds_coeff_ea,fitted_angs_coeff_ea,fitted_dihs_coeff_ea)
	open(unit=304,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/eq_par_8_6_6_os_fortran.dat',status='old')
	read(304,*)n_ob,n_oa,n_od
    allocate(fitted_bonds_coeff_os(0:fn_order_bonds_os-1,0:n_eb-1),fitted_angs_coeff_os(0:fn_order_angs_os-1,0:n_ea-1),fitted_dihs_coeff_os(0:fn_order_dihs_os-1,0:n_ed-1))
    call  Read_eq_int_coeff(304,n_ob,n_oa,n_od,Nbonds,Nang,N_int,fn_order_bonds_os,fn_order_angs_os,fn_order_dihs_os,os_fn_idx,fitted_bonds_coeff_os,fitted_angs_coeff_os,fitted_dihs_coeff_os)
	open(unit=305,file=trim(adjustl(potParamDir))//'/fortran_param_files_avtz/eq_par_8_6_6_oa_fortran.dat',status='old')
	read(305,*)n_ob,n_oa,n_od
    allocate(fitted_bonds_coeff_oa(0:fn_order_bonds_oa-1,0:n_eb-1),fitted_angs_coeff_oa(0:fn_order_angs_oa-1,0:n_ea-1),fitted_dihs_coeff_oa(0:fn_order_dihs_oa-1,0:n_ed-1))
    call  Read_eq_int_coeff(305,n_ob,n_oa,n_od,Nbonds,Nang,N_int,fn_order_bonds_oa,fn_order_angs_oa,fn_order_dihs_oa,oa_fn_idx,fitted_bonds_coeff_oa,fitted_angs_coeff_oa,fitted_dihs_coeff_oa)
	!Fitted_b_coeff in Hartee/angstrom2(or degree2 or angstrom*degree)  
	allocate(atom_ind(0:Natoms-1,0:3))
	call Read_atom_ind("atom_ind.dat",Natoms,atom_ind)
	n_fix_tors = 3 
	allocate(max_id(0:n_fix_tors-1))	
	max_id(0) = 8;max_id(1) = 6;max_id(2) = 6
	sym_mat = 0.d0
	!sym_pair =  reshape((/ 0,-1,-1,-1,1,2,-1,-1,3,4,5,6,7,8,-1,-1,9,10,-1,-1,11,12,13,14,15,16,-1,-1,18,19,20,21/),shape(sym_pair))
	do i = 0,size(tors_idx)-1
		sym_mat(tors_idx(i),tors_idx(i)) = 1.d0
	enddo
	do ii = 0,size(sym_pair(:,0))-1
	    if (sym_pair(ii,1)==-1) then
	       sym_mat(sym_pair(ii,0),sym_pair(ii,0))  = 1.d0
	    else if (sym_pair(ii,2)==-1) then
	       i = sym_pair(ii,0);j=sym_pair(ii,1)
	       sym_mat(i,i) =  sqrt_2_inv  
	       sym_mat(i,j) =  sqrt_2_inv
	       sym_mat(j,i) = -sqrt_2_inv
	       sym_mat(j,j) =  sqrt_2_inv
	    else if (sym_pair(ii,3)/=-1) then
	       i = sym_pair(ii,0);j=sym_pair(ii,1);k = sym_pair(ii,2);l=sym_pair(ii,3)
	       sym_mat(i,i) =  0.5d0; sym_mat(i,j)  =  0.5d0; sym_mat(i,k)  =  0.5d0; sym_mat(i,l)  =  0.5d0;
	       sym_mat(j,i) =  0.5d0; sym_mat(j,j)  = -0.5d0; sym_mat(j,k)  =  0.5d0; sym_mat(j,l)  = -0.5d0;
	       sym_mat(k,i) =  0.5d0; sym_mat(k,j)  = -0.5d0; sym_mat(k,k)  = -0.5d0; sym_mat(k,l)  =  0.5d0;
	       sym_mat(l,i) =  0.5d0; sym_mat(l,j)  =  0.5d0; sym_mat(l,k)  = -0.5d0; sym_mat(l,l)  = -0.5d0;
		endif
	enddo
	sym_mat_transpose = transpose(sym_mat)
end subroutine 
!****************************************************************************************#
! Reads the fitted force constant coefficients into  pot_coeff array.
! pot_coeff is an 5 dimensional array, first four corresponds to internal coordinates,
! and the fifth corresponds to order of the term.
subroutine Read_fc_coeff(f_n,N_int,fn_order,order,n_params,fc_idx,fitted_coeff)
    implicit none
    integer::i,j,k,ctr,idum
    integer,intent(in)::f_n,n_params,N_int,fn_order,order
    integer,intent(out),allocatable::fc_idx(:,:)
    real(8),intent(out),allocatable::fitted_coeff(:,:)
    !Input file has 5 columns. First four columns contains the index of the internal coordinates. 
    !and the last column has the corresponding fitted parameters. fc_idx stores the unique indexes 
    !and the corresponding row of fitted_coeff contains fitted coefficients for that index values
    allocate(fc_idx(0:n_params-1,0:order-1),fitted_coeff(0:fn_order-1,0:n_params-1))
    !n_params: no of coefficients in potential expansion
    fc_idx = 0;fitted_coeff = 0.d0
    ctr = 0
    do i = 1,fn_order * n_params,fn_order
        do k = 0,fn_order-1
            read(f_n,*)(fc_idx(ctr,j),j=0,order-1),idum,fitted_coeff(k,ctr)
        enddo
        ctr = ctr + 1
    enddo
    fc_idx  = fc_idx-1
    close(f_n)
endsubroutine
!****************************************************************************************#
! Reads the fitted force constant coefficients into  pot_coeff array.
! pot_coeff is an 5 dimensional array, first four corresponds to internal coordinates,
! and the fifth corresponds to order of the term.
subroutine Read_fc_coeff_ltpt1(f_n,N_int,order,n_params,fc_idx,fitted_coeff)
    implicit none
    integer::i,j,k,ctr,idum
    integer,intent(in)::f_n,n_params,N_int,order
    integer,intent(out),allocatable::fc_idx(:,:)
    real(8),intent(out),allocatable::fitted_coeff(:)
    !Input file has 5 columns. First four columns contains the index of the internal coordinates. 
    !and the last column has the corresponding fitted parameters. fc_idx stores the unique indexes 
    !and the corresponding row of fitted_coeff contains fitted coefficients for that index values
    allocate(fc_idx(0:n_params-1,0:order-1),fitted_coeff(0:n_params-1))
    !n_params: no of coefficients in potential expansion
    fc_idx = 0;fitted_coeff = 0.d0
    do i = 0,n_params-1
        read(f_n,*)(fc_idx(i,j),j=0,order-1),idum,fitted_coeff(i)
    enddo
    fc_idx  = fc_idx-1
    close(f_n)
endsubroutine
!****************************************************************************************#
subroutine  Read_Vt_coeff(f_name,fn_order,Vt_coeff)
    implicit none
    integer::i,idum1,idum2,idx
    real(8)::val
    character(len=*)::f_name
    integer,intent(in)::fn_order
    real(8),intent(out),dimension(0:fn_order-1)::Vt_coeff
    open(unit=2222,file=f_name,status="old")
    do i = 0,fn_order-1
        read(2222,*)idx,idum1,idum2,val
        if (i/=idx) then
            print*,"Error reading file ",f_name
            exit
        else
            Vt_coeff(idx)    = val
        endif
    enddo
	close(2222)
endsubroutine Read_Vt_coeff
!****************************************************************************************#
subroutine Read_eq_int_coeff(f_n,nb,na,nd,Nbonds,Nang,N_int,fn_order_bonds,fn_order_angs,fn_order_dihs,coordidx,fitted_bond_coeff,fitted_ang_coeff,fitted_dih_coeff)
    implicit none
    integer::i,j,idum,a0,a1,idx1,idx2,odd_b,odd_a,odd_d,b_ctr,a_ctr,d_ctr,n_params
    real(8)::val
    integer,intent(in) :: f_n,nb,na,nd,Nbonds,Nang,N_int,fn_order_bonds,fn_order_angs,fn_order_dihs,coordidx(*)
    !Dimension coordidx(*)
    real(8),intent(out),allocatable::fitted_bond_coeff(:,:),fitted_ang_coeff(:,:),fitted_dih_coeff(:,:)
	allocate(fitted_bond_coeff(0:fn_order_bonds-1,0:nb-1))
	allocate( fitted_ang_coeff(0:fn_order_angs-1,0:na-1))
	allocate( fitted_dih_coeff(0:fn_order_dihs-1,0:nd-1))
    fitted_bond_coeff = 0.d0;fitted_ang_coeff=0.d0;fitted_dih_coeff=0.d0
	n_params = nb + na  + nd 
	a_ctr = 0;b_ctr =0 ;d_ctr=0
    do i = 0,n_params-1
        read(f_n,*)a0,idum,val
        if (a0<Nbonds) then
        fitted_bond_coeff(0,b_ctr) = val
            do j = 0,fn_order_bonds-2
                read(f_n,*)a1,idum,val
                if (a1==a0) then
                    fitted_bond_coeff(j+1,b_ctr) = val
                else
                    exit
                endif
            enddo
            b_ctr = b_ctr + 1
        else  if ((a0>=Nbonds).and.(a0<Nbonds+Nang)) then
            fitted_ang_coeff(0,a_ctr) = val
            do j = 0,fn_order_angs-2
                read(f_n,*)a1,idum,val
                if (a1==a0) then
                    fitted_ang_coeff(j+1,a_ctr) = val
                else
                    exit
                endif
            enddo
            a_ctr = a_ctr + 1
        else
            fitted_dih_coeff(0,d_ctr) = val
            do j = 0,fn_order_dihs-2
                read(f_n,*)a1,idum,val
                if (a1==a0) then
                    fitted_dih_coeff(j+1,d_ctr) = val
                else
                    exit
                endif
            enddo
            d_ctr = d_ctr + 1
        endif
    enddo
	close(f_n)
end subroutine Read_eq_int_coeff
!****************************************************************************************#
subroutine Read_atom_ind(f_name,Natoms,atom_ind)
    implicit none
    integer::i,j
    character(len=*)::f_name
    integer,intent(in)::Natoms
    integer,intent(out),dimension(0:Natoms-1,0:3)::atom_ind
    open(unit=12,file=f_name,status="old")
    atom_ind = 0
    do i = 0,Natoms-1
        read(12,*)(atom_ind(i,j),j=0,3)
    enddo
endsubroutine Read_atom_ind
!******************************************************************************************************************************#
subroutine get_stride_arr_es(Mmax,Nmax,Lmax,fn_order,stride_arr)
	implicit none
    integer,intent(in):: Mmax,Nmax,Lmax,fn_order
    integer,intent(Out)::stride_arr(0:fn_order-1,0:2)
    integer::ctr,i,j,k

    ctr = 0
    do i = 0,Mmax
        do j = 0,Nmax
            stride_arr(ctr,:) = [i,j,j]
            ctr = ctr +  1
		enddo
	enddo
    do i = 0,Mmax
        do j = 0,Nmax
            do k = 0,j-1
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
    do i = 0,Mmax
        do j = 1,Nmax
            stride_arr(ctr,:) = [i,j,j]
            ctr = ctr +  1
		enddo
	enddo
    do i = 0,Mmax
        do j = 1,Nmax
            do k = 1,j-1
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
    do i = 1,Mmax
        do j = 0,Nmax
            do k = 1,Lmax
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine get_stride_arr_ea(Mmax,Nmax,Lmax,fn_order,stride_arr)
	implicit none
    integer,intent(in):: Mmax,Nmax,Lmax,fn_order
    integer,intent(Out)::stride_arr(0:fn_order-1,0:2)
    integer::ctr,i,j,k
    ctr = 0;
    do i = 0,Mmax
        do j = 1,Nmax
            do k = 0,j-1
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
    do i = 0,Mmax
        do j = 2,Nmax
            do k = 1,j-1
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
    do i = 1,Mmax
        do j = 0,Nmax
            do k = 1,Lmax
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine  get_stride_arr_os(Mmax,Nmax,Lmax,fn_order,stride_arr)
	implicit none
    integer,intent(in):: Mmax,Nmax,Lmax,fn_order
    integer,intent(Out)::stride_arr(0:fn_order-1,0:2)
    integer::ctr,i,j,k
    ctr=0
    do i = 1,Mmax
        do j = 0,Nmax
                stride_arr(ctr,:) = [i,j,j]
                ctr = ctr +  1
		enddo
	enddo
    do i = 1,Mmax
        do j = 0,Nmax
            do k = 0,j-1
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
    do i = 1,Mmax
        do j = 1,Nmax
                stride_arr(ctr,:) = [i,j,j]
                ctr = ctr +  1
		enddo
	enddo
    do i = 1,Mmax
        do j = 1,Nmax
            do k = 1,j-1
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo

    do i = 0,Mmax
        do j = 0,Nmax
            do k = 1,Lmax
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine  get_stride_arr_oa(Mmax,Nmax,Lmax,fn_order,stride_arr)
	implicit none
    integer,intent(in):: Mmax,Nmax,Lmax,fn_order
    integer,intent(Out)::stride_arr(0:fn_order-1,0:2)
    integer::ctr,i,j,k
    ctr=0
    do i = 1,Mmax
        do j = 1,Nmax
            do k = 0,j-1
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
    do i = 1,Mmax
        do j = 2,Nmax
            do k = 1,j-1
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
    do i = 0,Mmax
        do j = 0,Nmax
            do k = 1,Lmax
                stride_arr(ctr,:) = [i,j,k]
                ctr = ctr +  1
			enddo
		enddo
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine get_fn_order_es(Mmax,Nmax,Lmax,fn_order,Ncoeff)
    implicit none
    integer,intent(in)::Mmax,Nmax,Lmax
    integer,intent(out)::fn_order,Ncoeff(0:4)
	integer::n1,n2,n3,n4,n5,Nterms(0:4) 
    integer::i,j,k,ctr1,ctr2,ctr3,ctr4,ctr5
	Nterms = 0;Ncoeff = 0
    fn_order = 0;ctr1 = 0;ctr2 = 0;ctr3 = 0;ctr4 = 0;ctr5 = 0
    do i = 0,Mmax
        do j = 0,Nmax
            ctr1 = ctr1 +  1
		enddo
	enddo
    do i = 0,Mmax
        do j = 0,Nmax
            do k = 0,j-1
                ctr2 = ctr2 +  1
			enddo
		enddo
	enddo
    do i = 0,Mmax
        do j = 1,Nmax
            ctr3 = ctr3 +  1
		enddo
	enddo
    do i = 0,Mmax
        do j = 1,Nmax
            do k = 1,j-1
                ctr4 = ctr4 +  1
			enddo
		enddo
	enddo
    do i = 1,Mmax
        do j = 0,Nmax
            do k = 1,Lmax
                ctr5 = ctr5 +  1
			enddo
		enddo
	enddo
    Nterms(:) = [ctr1,ctr2,ctr3,ctr4,ctr5]
    fn_order = sum(Nterms)
    n1 = ctr1;n2 = n1+ctr2;n3=n2+ctr3;n4=n3+ctr4;n5=n4+ctr5
    Ncoeff(:) = [n1,n2,n3,n4,n5]
endsubroutine 
!******************************************************************************************************************************#
subroutine get_fn_order_ea(Mmax,Nmax,Lmax,fn_order,Ncoeff)
    implicit none
    integer,intent(in)::Mmax,Nmax,Lmax
    integer,intent(out)::fn_order,Ncoeff(0:4)
	integer::n1,n2,n3,n4,n5,Nterms(0:4) 
    integer::i,j,k,ctr1,ctr2,ctr3,ctr4,ctr5
	Nterms = 0;Ncoeff = 0
    fn_order = 0
    ctr1 = 0; ctr2 = 0; ctr3 = 0;
    do i = 0,Mmax
        do j = 1,Nmax
            do k = 0,j-1
                ctr1 = ctr1 +  1
			enddo
		enddo
	enddo
    do i = 0,Mmax
        do j = 2,Nmax
            do k = 1,j-1
                ctr2 = ctr2 +  1
			enddo
		enddo
	enddo
    do i = 1,Mmax
        do j = 0,Nmax
            do k = 1,Lmax
                ctr3 = ctr3 +  1
			enddo
		enddo
	enddo
    Nterms(:) = [ctr1,ctr2,ctr3,0,0]
    fn_order = sum(Nterms)
    n1 = ctr1;n2 = n1+ctr2;n3=n2+ctr3
    Ncoeff(:) = [n1,n2,n3,0,0]
endsubroutine 

!******************************************************************************************************************************#
subroutine get_fn_order_os(Mmax,Nmax,Lmax,fn_order,Ncoeff)
    implicit none
    integer,intent(in)::Mmax,Nmax,Lmax
    integer,intent(out)::fn_order,Ncoeff(0:4)
	integer::n1,n2,n3,n4,n5,Nterms(0:4) 
    integer::i,j,k,ctr1,ctr2,ctr3,ctr4,ctr5
	Nterms = 0;Ncoeff = 0
    fn_order = 0;ctr1 = 0;ctr2 = 0;ctr3 = 0;ctr4 = 0;ctr5 = 0
    do i = 1,Mmax
        do j = 0,Nmax
                ctr1 = ctr1 +  1
		enddo
	enddo
    do i = 1,Mmax
        do j = 0,Nmax
            do k = 0,j-1
                ctr2 = ctr2 +  1
			enddo
		enddo
	enddo
    do i = 1,Mmax
        do j = 1,Nmax
                ctr3 = ctr3 +  1
		enddo
	enddo
    do i = 1,Mmax
        do j = 1,Nmax
            do k = 1,j-1
                ctr4 = ctr4 +  1
			enddo
		enddo
	enddo
    do i = 0,Mmax
        do j = 0,Nmax
            do k = 1,Lmax
                ctr5 = ctr5 +  1
			enddo
		enddo
	enddo
    Nterms(:) = [ctr1,ctr2,ctr3,ctr4,ctr5]
    fn_order = sum(Nterms)
    n1 = ctr1;n2 = n1+ctr2;n3=n2+ctr3;n4=n3+ctr4;n5=n4+ctr5
    Ncoeff(:) = [n1,n2,n3,n4,n5]
endsubroutine 

!******************************************************************************************************************************#
subroutine get_fn_order_oa(Mmax,Nmax,Lmax,fn_order,Ncoeff)
    implicit none
    integer,intent(in)::Mmax,Nmax,Lmax
    integer,intent(out)::fn_order,Ncoeff(0:4)
	integer::n1,n2,n3,n4,n5,Nterms(0:4) 
    integer::i,j,k,ctr1,ctr2,ctr3,ctr4,ctr5
	Nterms = 0;Ncoeff = 0
    fn_order = 0
    ctr1 = 0; ctr2 = 0; ctr3 = 0;
    do i = 1,Mmax
        do j = 1,Nmax
            do k = 0,j-1
                ctr1 = ctr1 +  1
			enddo
		enddo
	enddo
    do i = 1,Mmax
        do j = 2,Nmax
            do k = 1,j-1
                ctr2 = ctr2 +  1
			enddo
		enddo
	enddo
    do i = 0,Mmax
        do j = 0,Nmax
            do k = 1,Lmax
                ctr3 = ctr3 +  1
			enddo
		enddo
	enddo
    Nterms(:) = [ctr1,ctr2,ctr3,0,0]
    fn_order = sum(Nterms)
    n1 = ctr1;n2 = n1+ctr2;n3=n2+ctr3
    Ncoeff(:) = [n1,n2,n3,0,0]
endsubroutine 
!******************************************************************************************************************************#
end module
