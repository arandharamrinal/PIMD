module eg_pot_force
use POTVARS
use allvars, only : cartmass, total_mass_inv
use constants
contains
!Change it 
!****************************************************************************************#
subroutine Read_ffext(f_name,Ncarts,X,dVdx,E)
    implicit none
    integer::i,j
    character(len=50)::f_name
    integer,intent(in)::Ncarts
    real(8),intent(out),dimension(0:Ncarts-1)::X,dVdx
    real(8),intent(out)::E
    open(unit=12,file=f_name,status="old")
    X = 0.d0 ;dVdx =0.d0;E = 0.d0
    do i = 0,3
        read(12,*)
    enddo
    !Read Cartesian coordinates (assuming 4th line)
    read(12,*)(X(j),j=0,Ncarts-1)
    do i = 0,3
        read(12,*)
    enddo
    read(12,*)E
    read(12,*)
    read(12,*)(dVdX(j),j=0,Ncarts-1)
    X = X * BohrToAng
endsubroutine Read_ffext
!****************************************************************************************#
!Function that converts cartesian to internal coordinates
subroutine  cart_to_internal(cart_coord,int_arr)
    implicit none
    integer::i,a_ctr,b_ctr,d_ctr,atm1,atm2,atm3,atm4
    real(8)::v_dot,v1_dot,v2_dot
    real(8)::nr1,nr2,nr3,ang,dih
    real(8),dimension(0:2)::r1,r2,r3,r4
    real(8),dimension(0:2)::a1,a2,v1,v2,v3,u1,u2,u3
    real(8),intent(in)::cart_coord(0:Ncarts-1)
    real(8),intent(out)::int_arr(0:N_int-1)
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
	int_arr = 0.d0
    do i = 1,Natoms-1
            atm1 = atom_ind(i,0)-1
            atm2 = atom_ind(i,1)-1
            atm3 = atom_ind(i,2)-1
            atm4 = atom_ind(i,3)-1
            if (i==1) then
                    !compute only bond length.
                    r1 = cart_coord(3*atm1:3*atm1+2)
                    r2 = cart_coord(3*atm2:3*atm2+2)
					v1 = r2-r1
    				nr1  = norm2(v1)
                    int_arr(b_ctr) = nr1
                    b_ctr =b_ctr + 1
            else if (i == 2) then
                    !compute only the bond length and bond angle
                    r1 = cart_coord(3*atm1:3*atm1+2)
                    r2 = cart_coord(3*atm2:3*atm2+2)
                    r3 = cart_coord(3*atm3:3*atm3+2)
					v1 = r2-r1
    				nr1  = norm2(v1)
					v2   = r3-r2
					nr2  = norm2(v2)
					u1 = v1/nr1
					u2 = v2/nr2
                    int_arr(b_ctr) = nr1
                    b_ctr =b_ctr +1
    				v_dot = dot_product(-u1,u2)
    				ang = dacos(v_dot) * RadianToDegree
    				ang = mod(ang,360.d0)
                    int_arr(Nbonds+a_ctr) = ang
                    a_ctr =a_ctr +1
            else
                    !compute all the three
                    r1 = cart_coord(3*atm1:3*atm1+2)
                    r2 = cart_coord(3*atm2:3*atm2+2)
                    r3 = cart_coord(3*atm3:3*atm3+2)
                    r4 = cart_coord(3*atm4:3*atm4+2)
					v1 = r2-r1
    				nr1  = norm2(v1)
					v2   = r3-r2
					nr2  = norm2(v2)
					v3   = r4-r3
					u1 = v1/nr1
					u2 = v2/nr2
                    int_arr(b_ctr) = nr1
                    b_ctr =b_ctr +1
    				v_dot = dot_product(-u1,u2)
    				ang = dacos(v_dot) * RadianToDegree
    				ang = mod(ang,360.d0)
                    int_arr(Nbonds + a_ctr) = ang
                    a_ctr =a_ctr +1
    				a1(0) = v1(1)*v2(2) - v2(1)*v1(2)
    				a1(1) = v2(0)*v1(2) - v2(2)*v1(0)
    				a1(2) = v1(0)*v2(1) - v1(1)*v2(0)
    				a2(0) = v2(1)*v3(2) - v3(1)*v2(2)
    				a2(1) = v3(0)*v2(2) - v3(2)*v2(0)
    				a2(2) = v2(0)*v3(1) - v2(1)*v3(0)
    				v1_dot = norm2(v2) * dot_product(v1,a2)
    				v2_dot = dot_product(a1,a2)
    				dih = datan2(v1_dot,v2_dot) * RadianToDegree
    				dih = mod(dih,360.d0)
    				if (dih>180.0) then
    				    dih = dih - 360.d0
    				endif
                    int_arr(Nbonds+Nang+d_ctr) = dih
                    d_ctr = d_ctr + 1
            endif
    enddo
endsubroutine
!****************************************************************************************#
!Build the rotation matrix 
subroutine get_rotation_matrix(M,theta,u)
    implicit none
    real(8),intent(out),dimension(0:2,0:2)::M
    real(8),intent(in),dimension(0:2) ::u
    real(8),intent(in)::theta
    real(8)::cth,sth

    !u : axis of rotation
    !theta : angle of rotation
    !Given "u" and "theta", M rotates a vector about "u" by an angle "theta".
    cth = dcos(theta)
    sth = dsin(theta)
    M(0,0) = u(0)*u(0) * (1.0 - cth) + cth
    M(0,1) = u(0)*u(1) * (1.0 - cth) - u(2) * sth
    M(0,2) = u(0)*u(2) * (1.0 - cth) + u(1) * sth
    M(1,0) = u(0)*u(1) * (1.0 - cth) + u(2) * sth
    M(1,1) = u(1)*u(1) * (1.0 - cth) + cth
    M(1,2) = u(1)*u(2) * (1.0 - cth) - u(0) * sth
    M(2,0) = u(0)*u(2) * (1.0 - cth) - u(1) * sth
    M(2,1) = u(1)*u(2) * (1.0 - cth) + u(0) * sth
    M(2,2) = u(2)*u(2) * (1.0 - cth) + cth
end subroutine get_rotation_matrix
!*************************************************************************************!
subroutine get_norm(Ndim,v,norm)
    implicit none
    real(8),intent(in),dimension(0:2)::v
    integer,intent(in):: Ndim
    real(8)::sum1
    integer::i
    real(8),intent(out)::norm
    sum1 = 0.0
    do i =0,Ndim-1
        sum1 = sum1 + v(i)**2
    enddo
    norm = dsqrt(sum1)
endsubroutine get_norm
!*************************************************************************************!
subroutine get_cross_prod(v1,v2,v_cross,Ndim)
    implicit none
    integer,intent(in)::Ndim
    real(8),intent(in),dimension(0:2) :: v1,v2
    real(8),intent(out),dimension(0:2)::v_cross

    v_cross=0.0
    v_cross(0) = v1(1)*v2(2) - v1(2)*v2(1)
    v_cross(1) = v1(2)*v2(0) - v1(0)*v2(2)
    v_cross(2) = v1(0)*v2(1) - v1(1)*v2(0)
endsubroutine get_cross_prod
!*************************************************************************************!
subroutine internal_calculator(u1,u2,r2,a2,d2,q4)
    implicit none
    real(8),intent(in),dimension(0:2)::u1,u2
    real(8),intent(in)::a2,d2,r2
    real(8),dimension(0:2)::ua1,ua2,n_v
    real(8)::ang2,dih2
    real(8),intent(out),dimension(0:2)::q4
    real(8)::pi,norm
    real(8),dimension(0:2,0:2)::M
    integer::i
    ua1 = 0.d0;ua2 = 0.d0;n_v = 0.d0
    ang2=a2;dih2=d2
    pi = 4.d0*atan(1.d0)
    ang2    =  pi- ang2 * pi / 180.d0
    dih2    =  dih2 * pi / 180.d0
    ua1 = u1;ua2=u2
    !normalize u2
    call get_norm(3,ua2,norm)
    ua2   = ua2 / norm

    !Compute the normal vector to the plane of the first 3 atoms.
    call get_cross_prod(ua1,ua2,n_v,3)
    call get_norm(3,n_v,norm)
    n_v  = n_v / norm
    !place the fourth atom(D) along the bond between atoms 2 and 3.
    q4   = ua2 * r2

    !Get rotaion matrix for rotation around n_v
    call get_rotation_matrix(M,ang2,n_v)
    q4   = matmul(M,q4)

    !Get rotation matrix for rotation around u2
    call get_rotation_matrix(M,dih2,ua2)
    q4   = matmul(M,q4)
end subroutine
!**********************************************************************************
subroutine get_com(X,R)
    implicit none
    integer::i,j
    real(8),intent(in),dimension(0:Ncarts-1)::X
    real(8),intent(out),dimension(0:2)::R
    R=0.d0

    !Takes cartesian coordinates(3*Natoms dimensional) and cartmasses as inputs and computes the center of cartmass"""
    do i = 0,Natoms-1
        do j = 0,2
            R(j) = R(j) + cartmass(3*i) * X(3*i+j)
        enddo
    enddo
    R =R*total_mass_inv
end subroutine
!**********************************************************************************
subroutine shift_cm_to_origin(X)
    implicit none
    !Shifts the origin of the coordinate system to the CoM"""
    integer::i
    real(8),intent(inout),dimension(0:Ncarts-1)::X
    real(8),dimension(0:2)::R
    R = 0.d0
    call  get_com(X,R)
    do i = 0,Natoms-1
        X(3*i)   = X(3*i)   -R(0)
        X(3*i+1) = X(3*i+1) -R(1)
        X(3*i+2) = X(3*i+2) -R(2)
    enddo
end subroutine
!**********************************************************************************
subroutine  get_moi_tensor(X,M)
    !Computes the Moment of inertia tensor"""
	use POTVARS, only: cartmass,Natoms,Ncarts
    implicit none
    integer::i
    real(8),intent(in),dimension(0:Ncarts-1)::X
    real(8),intent(out),dimension(0:2,0:2)::M
    M = 0.d0
    do i = 0,Natoms-1
        M(0,0) = M(0,0) + cartmass(3*i) * (X(3*i+1)**2 + X(3*i+2)**2)
        M(1,1) = M(1,1) + cartmass(3*i) * (X(3*i)**2    + X(3*i+2)**2)
        M(2,2) = M(2,2) + cartmass(3*i) * (X(3*i+1)**2 + X(3*i)**2  )
        M(0,1) = M(0,1) - cartmass(3*i) *  X(3*i)   * X(3*i+1)
        M(0,2) = M(0,2) - cartmass(3*i) *  X(3*i)   * X(3*i+2)
        M(1,2) = M(1,2) - cartmass(3*i) *  X(3*i+1) * X(3*i+2)
    enddo
    M(1,0) = M(0,1)
    M(2,0) = M(0,2)
    M(2,1) = M(1,2)
end subroutine
!**********************************************************************************
subroutine  transform_to_pa(X)
    !Transforms the coordinates(X) so that MoI is diagonal for the new(output) geometry """
    implicit none
    integer::i
    real(8)::det
    real(8),intent(inout),dimension(0:Ncarts-1)::X
    real(8),dimension(0:2,0:2)::MoI,PA
    real(8),dimension(0:2) :: W
    MoI = 0.d0
    call get_moi_tensor(X,MoI)
    call HOUSEDIAG(MoI,3,3,W)
    PA = MoI;
    det = 0.d0
    det =  PA(0,0) * (PA(1,1)*PA(2,2)-PA(1,2)*PA(2,1)) + PA(0,1) * (PA(1,2)*PA(2,0)-PA(1,0)*PA(2,2))
    det = det + PA(0,2) *(PA(0,0)*PA(2,1)-PA(1,1)*PA(2,0))
    if (det<0.d0) then
        PA(:,0) = -1.0 * PA(:,0)
    endif
    do i = 0,Natoms-1
        X(3*i:3*i+2) = matmul(transpose(PA),X(3*i:3*i+2))
    enddo
end subroutine
!**********************************************************************************

subroutine getVsav(S,fc_idx,fc,n_params,Vsav)
	implicit none
	integer,intent(in):: n_params,fc_idx(0:n_params-1,0:3)
	real(8),intent(in):: S(0:N_int-1),fc(0:n_params-1) 
	real(8),intent(out):: Vsav
	integer:: ii,i,j,k,l
    Vsav = 0.d0
	!print*,"SS=",S
	do ii = 0,n_params-1
        i = fc_idx(ii,0);j = fc_idx(ii,1);k = fc_idx(ii,2);l = fc_idx(ii,3);
        if (k==-1) then   
            Vsav = Vsav + fc(ii) * S(i) * S(j)
        elseif (l==-1) then
            Vsav = Vsav + fc(ii) * S(i) * S(j) * S(k)
        else
            Vsav = Vsav + fc(ii) * S(i) * S(j) * S(k) * S(l)
		endif
	enddo 
endsubroutine getVsav
!****************************************************************************************#
!compute force constants for a given value of phi1 and phi2
!****************************************************************************************#
subroutine get_all_vsav(S,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2 ,fc_quartic_ea_ltpt1,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,Vsav) 
		implicit none 
		real(8),intent(in):: S(0:N_int-1)
		real(8),intent(in):: fc_quad_es(0:n_fcs_quad_es-1),fc_quad_os(0:n_fcs_quad_os-1)
		real(8),intent(in):: fc_quad_ea(0:n_fcs_quad_ea-1),fc_quad_oa(0:n_fcs_quad_oa-1)
		real(8),intent(in):: fc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1),fc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1),fc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1)
		real(8),intent(in):: fc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1),fc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1),fc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1)
		real(8),intent(in):: fc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1),fc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1),fc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1) 
		real(8),intent(in):: fc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1),fc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1),fc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1) 
		real(8),intent(in):: fc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1),fc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1),fc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1) 
		real(8),intent(in):: fc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1),fc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1),fc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1) 
		real(8),intent(in):: fc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1),fc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1),fc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1) 
		real(8),intent(in):: fc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1),fc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1),fc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1) 
		real(8),intent(in)::fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
		real(8),intent(in)::fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1) 
		real(8),intent(out) :: Vsav
		real(8) :: Vsav_quad_ea,Vsav_quad_oa,Vsav_quad_es,Vsav_quad_os
		real(8) :: Vsav_qubic_gt5_es,Vsav_qubic_lt5_es,Vsav_qubic_lt2_es
		real(8) :: Vsav_qubic_gt5_ea,Vsav_qubic_lt5_ea,Vsav_qubic_lt2_ea
		real(8) :: Vsav_qubic_gt5_os,Vsav_qubic_lt5_os,Vsav_qubic_lt2_os
		real(8) :: Vsav_qubic_gt5_oa,Vsav_qubic_lt5_oa,Vsav_qubic_lt2_oa
		real(8) :: Vsav_quartic_gt5_es,Vsav_quartic_lt5_es,Vsav_quartic_lt2_es
		real(8) :: Vsav_quartic_gt5_ea,Vsav_quartic_lt5_ea,Vsav_quartic_lt2_ea
		real(8) :: Vsav_quartic_gt5_os,Vsav_quartic_lt5_os,Vsav_quartic_lt2_os
		real(8) :: Vsav_quartic_gt5_oa,Vsav_quartic_lt5_oa,Vsav_quartic_lt2_oa
		real(8) :: Vsav_quartic_ltpt1_os,Vsav_quartic_ltpt1_oa 
		real(8) :: Vsav_quartic_ltpt1_es,Vsav_quartic_ltpt1_ea 
        Vsav = 0.d0
        !quadratic 
        call getVsav(S,fc_idx_quad_es,fc_quad_es,n_fcs_quad_es,Vsav_quad_es)
        call getVsav(S,fc_idx_quad_ea,fc_quad_ea,n_fcs_quad_ea,Vsav_quad_ea)
        call getVsav(S,fc_idx_quad_os,fc_quad_os,n_fcs_quad_os,Vsav_quad_os)
        call getVsav(S,fc_idx_quad_oa,fc_quad_oa,n_fcs_quad_oa,Vsav_quad_oa)
        !Cubic even
        call getVsav(S,fc_idx_qubic_es_gt5,fc_qubic_es_gt5,n_fcs_qubic_es_gt5,Vsav_qubic_gt5_es)
        call getVsav(S,fc_idx_qubic_es_lt5,fc_qubic_es_lt5,n_fcs_qubic_es_lt5,Vsav_qubic_lt5_es)
        call getVsav(S,fc_idx_qubic_es_lt2,fc_qubic_es_lt2,n_fcs_qubic_es_lt2,Vsav_qubic_lt2_es)
        call getVsav(S,fc_idx_qubic_ea_gt5,fc_qubic_ea_gt5,n_fcs_qubic_ea_gt5,Vsav_qubic_gt5_ea)
        call getVsav(S,fc_idx_qubic_ea_lt5,fc_qubic_ea_lt5,n_fcs_qubic_ea_lt5,Vsav_qubic_lt5_ea)
        call getVsav(S,fc_idx_qubic_ea_lt2,fc_qubic_ea_lt2,n_fcs_qubic_ea_lt2,Vsav_qubic_lt2_ea)
        !Cubic Odd
        call getVsav(S,fc_idx_qubic_os_gt5,fc_qubic_os_gt5,n_fcs_qubic_os_gt5,Vsav_qubic_gt5_os)
        call getVsav(S,fc_idx_qubic_os_lt5,fc_qubic_os_lt5,n_fcs_qubic_os_lt5,Vsav_qubic_lt5_os)
        call getVsav(S,fc_idx_qubic_os_lt2,fc_qubic_os_lt2,n_fcs_qubic_os_lt2,Vsav_qubic_lt2_os)
        call getVsav(S,fc_idx_qubic_oa_gt5,fc_qubic_oa_gt5,n_fcs_qubic_oa_gt5,Vsav_qubic_gt5_oa)
        call getVsav(S,fc_idx_qubic_oa_lt5,fc_qubic_oa_lt5,n_fcs_qubic_oa_lt5,Vsav_qubic_lt5_oa)
        call getVsav(S,fc_idx_qubic_oa_lt2,fc_qubic_oa_lt2,n_fcs_qubic_oa_lt2,Vsav_qubic_lt2_oa)
        !Quartic even   
        call getVsav(S,fc_idx_quartic_es_gt5,fc_quartic_es_gt5,n_fcs_quartic_es_gt5,Vsav_quartic_gt5_es)
        call getVsav(S,fc_idx_quartic_es_lt5,fc_quartic_es_lt5,n_fcs_quartic_es_lt5,Vsav_quartic_lt5_es)
        call getVsav(S,fc_idx_quartic_es_lt2,fc_quartic_es_lt2,n_fcs_quartic_es_lt2,Vsav_quartic_lt2_es)
        call getVsav(S,fc_idx_quartic_es_ltpt1,fc_quartic_es_ltpt1,n_fcs_quartic_es_ltpt1,Vsav_quartic_ltpt1_es)
        call getVsav(S,fc_idx_quartic_ea_gt5,fc_quartic_ea_gt5,n_fcs_quartic_ea_gt5,Vsav_quartic_gt5_ea)
        call getVsav(S,fc_idx_quartic_ea_lt5,fc_quartic_ea_lt5,n_fcs_quartic_ea_lt5,Vsav_quartic_lt5_ea)
        call getVsav(S,fc_idx_quartic_ea_lt2,fc_quartic_ea_lt2,n_fcs_quartic_ea_lt2,Vsav_quartic_lt2_ea)
        call getVsav(S,fc_idx_quartic_ea_ltpt1,fc_quartic_ea_ltpt1,n_fcs_quartic_ea_ltpt1,Vsav_quartic_ltpt1_ea)
        !Quartic odd    
        call getVsav(S,fc_idx_quartic_os_gt5,fc_quartic_os_gt5,n_fcs_quartic_os_gt5,Vsav_quartic_gt5_os)
        call getVsav(S,fc_idx_quartic_os_lt5,fc_quartic_os_lt5,n_fcs_quartic_os_lt5,Vsav_quartic_lt5_os)
        call getVsav(S,fc_idx_quartic_os_lt2,fc_quartic_os_lt2,n_fcs_quartic_os_lt2,Vsav_quartic_lt2_os)
        call getVsav(S,fc_idx_quartic_os_ltpt1,fc_quartic_os_ltpt1,n_fcs_quartic_os_ltpt1,Vsav_quartic_ltpt1_os)
        call getVsav(S,fc_idx_quartic_oa_gt5,fc_quartic_oa_gt5,n_fcs_quartic_oa_gt5,Vsav_quartic_gt5_oa)
        call getVsav(S,fc_idx_quartic_oa_lt5,fc_quartic_oa_lt5,n_fcs_quartic_oa_lt5,Vsav_quartic_lt5_oa)
        call getVsav(S,fc_idx_quartic_oa_lt2,fc_quartic_oa_lt2,n_fcs_quartic_oa_lt2,Vsav_quartic_lt2_oa)
        call getVsav(S,fc_idx_quartic_oa_ltpt1,fc_quartic_oa_ltpt1,n_fcs_quartic_oa_ltpt1,Vsav_quartic_ltpt1_oa)
        !Sum up all 
		Vsav = Vsav + Vsav_quad_ea + Vsav_quad_oa + Vsav_quad_es + Vsav_quad_os
		Vsav = Vsav + Vsav_qubic_gt5_es + Vsav_qubic_lt5_es + Vsav_qubic_lt2_es
		Vsav = Vsav + Vsav_qubic_gt5_ea + Vsav_qubic_lt5_ea + Vsav_qubic_lt2_ea
		Vsav = Vsav + Vsav_qubic_gt5_os + Vsav_qubic_lt5_os + Vsav_qubic_lt2_os
		Vsav = Vsav + Vsav_qubic_gt5_oa + Vsav_qubic_lt5_oa + Vsav_qubic_lt2_oa
		Vsav = Vsav + Vsav_quartic_gt5_es + Vsav_quartic_lt5_es + Vsav_quartic_lt2_es
		Vsav = Vsav + Vsav_quartic_gt5_ea + Vsav_quartic_lt5_ea + Vsav_quartic_lt2_ea
		Vsav = Vsav + Vsav_quartic_lt5_os + Vsav_quartic_lt2_os
		Vsav = Vsav + Vsav_quartic_lt5_oa + Vsav_quartic_lt2_oa
endsubroutine  

!****************************************************************************************#
subroutine  compute_pot(S,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1 ,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,pot,Vt_out) 
	implicit none 
	real(8),intent(in):: S(0:N_int-1),sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(in):: fc_quad_es(0:n_fcs_quad_es-1),fc_quad_os(0:n_fcs_quad_os-1)
	real(8),intent(in):: fc_quad_ea(0:n_fcs_quad_ea-1),fc_quad_oa(0:n_fcs_quad_oa-1)
	real(8),intent(in):: fc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1),fc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1),fc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1)
	real(8),intent(in):: fc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1),fc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1),fc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1)
	real(8),intent(in):: fc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1),fc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1),fc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1) 
	real(8),intent(in):: fc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1),fc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1),fc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1) 
	real(8),intent(in):: fc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1),fc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1),fc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1) 
	real(8),intent(in):: fc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1),fc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1),fc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1) 
	real(8),intent(in):: fc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1),fc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1),fc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1) 
	real(8),intent(in):: fc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1),fc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1),fc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1) 
	real(8),intent(in)::fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
	real(8),intent(in)::fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1) 
	real(8),intent(out)::pot,Vt_out
	real(8)::Vt=0.d0,Vsav=0.d0
    pot = 0.d0
	call getVt(sinth,costh,Vt_out)
	call get_all_vsav(S,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2 ,fc_quartic_ea_ltpt1,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,Vsav) 
    pot  = Vt_out + Vsav 
endsubroutine 
!****************************************************************************************#
subroutine  dvdsQuad(fc_idx,fc,dS_dtors,dfc,S,icoord,n_params,dVds)
	implicit none 
	integer,intent(in)::icoord,n_params
	real(8),intent(in)::fc(0:n_params-1),dfc(0:n_params-1,0:n_fix_tors-1),dS_dtors(0:N_int-1,0:n_fix_tors-1),S(0:N_int-1)
	integer,intent(in) :: fc_idx(0:n_params-1,0:3)
	real(8),intent(inout):: dVds(0:N_int-1)
	integer:: ii,i,itors,j,k,l
	do ii =0,n_params-1
        i    = fc_idx(ii,0);j = fc_idx(ii,1);k = fc_idx(ii,2);l = fc_idx(ii,3);
		dVds(i) = dVds(i) + fc(ii) * S(j)
		dVds(j) = dVds(j) + fc(ii) * S(i)
		dVds(tors_idx(0)) = dVds(tors_idx(0)) + dfc(ii,0) * S(i) * S(j) + fc(ii) * dS_dtors(i,0) * S(j) + fc(ii) * dS_dtors(j,0) * S(i)
		dVds(tors_idx(1)) = dVds(tors_idx(1)) + dfc(ii,1) * S(i) * S(j) + fc(ii) * dS_dtors(i,1) * S(j) + fc(ii) * dS_dtors(j,1) * S(i)
	enddo
endsubroutine 
!****************************************************************************************#
subroutine  dvdsQubic(fc_idx,fc,dS_dtors,dfc,S,icoord,n_params,dVds)
	implicit none 
	integer,intent(in)::icoord,n_params
	real(8),intent(in)::fc(0:n_params-1),dfc(0:n_params-1,0:n_fix_tors-1),dS_dtors(0:N_int-1,0:n_fix_tors-1),S(0:N_int-1)
	integer,intent(in) :: fc_idx(0:n_params-1,0:3)
	real(8),intent(inout):: dVds(0:N_int-1)
	integer:: ii,i,itors,j,k,l
	do ii =0,n_params-1
        i    = fc_idx(ii,0);j = fc_idx(ii,1);k = fc_idx(ii,2);l = fc_idx(ii,3);
		dVds(i) = dVds(i) + fc(ii) * S(j) * S(k)
		dVds(j) = dVds(j) + fc(ii) * S(i) * S(k)
		dVds(k) = dVds(k) + fc(ii) * S(i) * S(j)
        dVds(tors_idx(0)) = dVds(tors_idx(0)) + dfc(ii,0) * S(i) * S(j) * S(k) + fc(ii) * dS_dtors(i,0) * S(j) * S(k)  + fc(ii) * dS_dtors(j,0) * S(i) * S(k) 
        dVds(tors_idx(0)) = dVds(tors_idx(0)) + fc(ii) * dS_dtors(k,0) * S(i) * S(j)
        dVds(tors_idx(1)) = dVds(tors_idx(1)) + dfc(ii,1) * S(i) * S(j) * S(k) + fc(ii) * dS_dtors(i,1) * S(j) * S(k)  + fc(ii) * dS_dtors(j,1) * S(i) * S(k) 
        dVds(tors_idx(1)) = dVds(tors_idx(1)) + fc(ii) * dS_dtors(k,1) * S(i) * S(j)
	enddo
endsubroutine 
!****************************************************************************************#
subroutine  dvdsQuartic(fc_idx,fc,dS_dtors,dfc,S,icoord,n_params,dVds)
	implicit none 
	integer,intent(in)::icoord,n_params
	real(8),intent(in)::fc(0:n_params-1),dfc(0:n_params-1,0:n_fix_tors-1),dS_dtors(0:N_int-1,0:n_fix_tors-1),S(0:N_int-1)
	integer,intent(in) :: fc_idx(0:n_params-1,0:3)
	real(8),intent(inout):: dVds(0:N_int-1)
	integer:: ii,i,itors,j,k,l
	do ii =0,n_params-1
        i    = fc_idx(ii,0);j = fc_idx(ii,1);k = fc_idx(ii,2);l = fc_idx(ii,3);
		dVds(i) = dVds(i) + fc(ii) * S(j) * S(k) * S(l)
		dVds(j) = dVds(j) + fc(ii) * S(i) * S(k) * S(l)
		dVds(k) = dVds(k) + fc(ii) * S(i) * S(j) * S(l)
		dVds(l) = dVds(l) + fc(ii) * S(i) * S(j) * S(k)
        dVds(tors_idx(0)) =dVds(tors_idx(0)) + dfc(ii,0) * S(i) * S(j) * S(k) * S(l)

        dVds(tors_idx(0)) =dVds(tors_idx(0)) + fc(ii) * S(i) * dS_dtors(j,0) * S(k) * S(l) + fc(ii) * S(j) * dS_dtors(i,0) * S(k) * S(l) 
        dVds(tors_idx(0)) =dVds(tors_idx(0)) + fc(ii) * dS_dtors(k,0) * S(i) * S(j) * S(l) + fc(ii) * dS_dtors(l,0) * S(i) * S(j) * S(k)
        dVds(tors_idx(1)) =dVds(tors_idx(1)) + dfc(ii,1) * S(i) * S(j) * S(k) * S(l)
        dVds(tors_idx(1)) =dVds(tors_idx(1)) + fc(ii) * S(i) * dS_dtors(j,1) * S(k) * S(l) + fc(ii) * S(j) * dS_dtors(i,1) * S(k) * S(l) 
        dVds(tors_idx(1)) =dVds(tors_idx(1)) + fc(ii) * dS_dtors(k,1) * S(i) * S(j) * S(l) + fc(ii) * dS_dtors(l,1) * S(i) * S(j) * S(k)
	enddo
endsubroutine 
!****************************************************************************************#
!compute internal gradients.
subroutine getDvds(S,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2 ,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,dfc_quad_es,dfc_quad_os,dfc_quad_ea,dfc_quad_oa,dfc_qubic_es_gt5,dfc_qubic_es_lt5,dfc_qubic_es_lt2,dfc_qubic_ea_gt5,dfc_qubic_ea_lt5,dfc_qubic_ea_lt2,dfc_qubic_os_gt5,dfc_qubic_os_lt5,dfc_qubic_os_lt2,dfc_qubic_oa_gt5,dfc_qubic_oa_lt5,dfc_qubic_oa_lt2,dfc_quartic_es_gt5,dfc_quartic_es_lt5,dfc_quartic_es_lt2,dfc_quartic_es_ltpt1,dfc_quartic_ea_gt5,dfc_quartic_ea_lt5,dfc_quartic_ea_lt2,dfc_quartic_ea_ltpt1, dfc_quartic_os_gt5,dfc_quartic_os_lt5,dfc_quartic_os_ltpt1,dfc_quartic_os_lt2,dfc_quartic_oa_gt5,dfc_quartic_oa_lt5,dfc_quartic_oa_lt2,dfc_quartic_oa_ltpt1,dVds) 
	implicit none 
	integer:: icoord
	real(8) :: dfcds(0:n_fix_tors-1),dVt(0:n_fix_tors-1),dS_dtors(0:N_int-1,0:n_fix_tors-1)
	integer:: t_id,i,idxs,b_ctr,a_ctr,d_ctr
	real(8)::dVdstmp
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(out)::dVds(0:N_int-1)
	real(8),intent(in) :: S(0:N_int-1)
	real(8),intent(in):: fc_quad_es(0:n_fcs_quad_es-1),fc_quad_os(0:n_fcs_quad_os-1)
	real(8),intent(in):: fc_quad_ea(0:n_fcs_quad_ea-1),fc_quad_oa(0:n_fcs_quad_oa-1)
	real(8),intent(in):: fc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1),fc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1),fc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1)
	real(8),intent(in):: fc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1),fc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1),fc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1)
	real(8),intent(in):: fc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1),fc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1),fc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1) 
	real(8),intent(in):: fc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1),fc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1),fc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1) 
	real(8),intent(in):: fc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1),fc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1),fc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1) 
	real(8),intent(in):: fc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1),fc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1),fc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1) 
	real(8),intent(in):: fc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1),fc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1),fc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1) 
	real(8),intent(in):: fc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1),fc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1),fc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1) 
	real(8),intent(in)::fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
	real(8),intent(in)::fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1) 
	real(8),intent(in):: dfc_quad_es(0:n_fcs_quad_es-1,0:n_fix_tors-1),dfc_quad_os(0:n_fcs_quad_os-1,0:n_fix_tors-1),dfc_quad_ea(0:n_fcs_quad_ea-1,0:n_fix_tors-1),dfc_quad_oa(0:n_fcs_quad_oa-1,0:n_fix_tors-1)
	real(8),intent(in):: dfc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1,0:n_fix_tors-1),dfc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1,0:n_fix_tors-1),dfc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1,0:n_fix_tors-1) 
	real(8),intent(in):: dfc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1,0:n_fix_tors-1),dfc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1,0:n_fix_tors-1), dfc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1,0:n_fix_tors-1)
	real(8),intent(in):: dfc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1,0:n_fix_tors-1),dfc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1,0:n_fix_tors-1),dfc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1,0:n_fix_tors-1) 
	real(8),intent(in):: dfc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1,0:n_fix_tors-1),dfc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1,0:n_fix_tors-1),dfc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1,0:n_fix_tors-1) 
	real(8),intent(in):: dfc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1,0:n_fix_tors-1),dfc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1,0:n_fix_tors-1),dfc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1,0:n_fix_tors-1) 
	real(8),intent(in):: dfc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1,0:n_fix_tors-1),dfc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1,0:n_fix_tors-1),dfc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1,0:n_fix_tors-1) 
	real(8),intent(in):: dfc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1,0:n_fix_tors-1),dfc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1,0:n_fix_tors-1),dfc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1,0:n_fix_tors-1) 
	real(8),intent(in):: dfc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1,0:n_fix_tors-1),dfc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1,0:n_fix_tors-1),dfc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1,0:n_fix_tors-1) 
	real(8),intent(in)::dfc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1,0:n_fix_tors-1),dfc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1,0:n_fix_tors-1)
	real(8),intent(in)::dfc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1,0:n_fix_tors-1),dfc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1,0:n_fix_tors-1)
    !""" i'th element of dVds contains the derivatives of potential w.r.t to the i'th internal"""
    !Compute the derivatives of the i'th internal w.r.t the fixed coordinates
	dS_dtors = 0.d0
	dVds = 0.d0
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    do i = 0,size(os_fn_idx)-1
        idxs = os_fn_idx(i)
        if (idxs<Nbonds) then
			dfcds = 0.d0
			call DTrigseries_os(sinth,costh,Ncoeff_bonds_os,fitted_bonds_coeff_os(:,b_ctr),stride_arr_b_os,fn_order_bonds_os,dfcds)
            dS_dtors(idxs,:) = -dfcds(:) 
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
			dfcds = 0.d0
			call DTrigseries_os(sinth,costh,Ncoeff_angs_os,fitted_angs_coeff_os(:,a_ctr),stride_arr_a_os,fn_order_angs_os,dfcds)
            dS_dtors(idxs,:) = -dfcds(:) 
            a_ctr = a_ctr + 1 
        else
			dfcds = 0.d0
			call DTrigseries_os(sinth,costh,Ncoeff_dihs_os,fitted_dihs_coeff_os(:,d_ctr),stride_arr_d_os,fn_order_dihs_os,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:)
            d_ctr = d_ctr + 1 
		endif
	enddo
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    do i = 0,size(oa_fn_idx)-1
        idxs = oa_fn_idx(i)
        if (idxs<Nbonds) then
			dfcds = 0.d0
			call DTrigseries_oa(sinth,costh,Ncoeff_bonds_oa,fitted_bonds_coeff_oa(:,b_ctr),stride_arr_b_oa,fn_order_bonds_oa,dfcds)
            dS_dtors(idxs,:) = -dfcds(:) 
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
			dfcds = 0.d0
			call DTrigseries_oa(sinth,costh,Ncoeff_angs_oa,fitted_angs_coeff_oa(:,a_ctr),stride_arr_a_oa,fn_order_angs_oa,dfcds)
            dS_dtors(idxs,:) = -dfcds(:) 
            a_ctr = a_ctr + 1 
        else
			dfcds = 0.d0
			call DTrigseries_oa(sinth,costh,Ncoeff_dihs_oa,fitted_dihs_coeff_oa(:,d_ctr),stride_arr_d_oa,fn_order_dihs_oa,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:)
            d_ctr = d_ctr + 1 
		endif
	enddo
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    do i = 0,size(es_fn_idx)-1
        idxs = es_fn_idx(i)
        if (idxs<Nbonds) then
			dfcds = 0.d0
			call DTrigseries_es(sinth,costh,Ncoeff_bonds_es,fitted_bonds_coeff_es(:,b_ctr),stride_arr_b_es,fn_order_bonds_es,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:)
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
			dfcds = 0.d0
			call DTrigseries_es(sinth,costh,Ncoeff_angs_es,fitted_angs_coeff_es(:,a_ctr),stride_arr_a_es,fn_order_angs_es,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:) 
            a_ctr = a_ctr + 1 
        else
			dfcds = 0.d0
			call DTrigseries_es(sinth,costh,Ncoeff_dihs_es,fitted_dihs_coeff_es(:,d_ctr),stride_arr_d_es,fn_order_dihs_es,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:) 
            d_ctr = d_ctr + 1 
		endif
	enddo
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    do i = 0,size(ea_fn_idx)-1
        idxs = ea_fn_idx(i)
        if (idxs<Nbonds) then
			dfcds = 0.d0
			call DTrigseries_ea(sinth,costh,Ncoeff_bonds_ea,fitted_bonds_coeff_ea(:,b_ctr),stride_arr_b_ea,fn_order_bonds_ea,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:)
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
			dfcds = 0.d0
			call DTrigseries_ea(sinth,costh,Ncoeff_angs_ea,fitted_angs_coeff_ea(:,a_ctr),stride_arr_a_ea,fn_order_angs_ea,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:) 
            a_ctr = a_ctr + 1 
        else
			dfcds = 0.d0
			call DTrigseries_ea(sinth,costh,Ncoeff_dihs_ea,fitted_dihs_coeff_ea(:,d_ctr),stride_arr_d_ea,fn_order_dihs_ea,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:) 
            d_ctr = d_ctr + 1 
		endif
	enddo
	do i = 0,N_int-1 
		dS_dtors(i,:) = dS_dtors(i,:) * dim_scal_factor_inv(i)
		if (i>= Nbonds) then
				dS_dtors(i,:) = dS_dtors(i,:) * degreeToRadian
		endif
	enddo
    !Compute the derivative of the potential along the minimum path(zero'th order term in the expansion) w.r.t the fixed coordinates
    call getDVt(sinth,costh,dVt)
    dVds(tors_idx) =  dVds(tors_idx) +  dVt(:)
    !Compute derivative of the 2nd order, 3rd order and fourth order terms w.r.t to all internal coordinates
    !!For quadratic even
	!ES
	call dvdsQuad(fc_idx_quad_es,fc_quad_es,dS_dtors,dfc_quad_es,S,icoord,n_fcs_quad_es,dVds)
	!EA
	call dvdsQuad(fc_idx_quad_ea,fc_quad_ea,dS_dtors,dfc_quad_ea,S,icoord,n_fcs_quad_ea,dVds)
	!OS
    call dvdsQuad(fc_idx_quad_os,fc_quad_os,dS_dtors,dfc_quad_os,S,icoord,n_fcs_quad_os,dVds)
	!OA
    call dvdsQuad(fc_idx_quad_oa,fc_quad_oa,dS_dtors,dfc_quad_oa,S,icoord,n_fcs_quad_oa,dVds)

    !Cubic:
	!ES
    call dvdsQubic(fc_idx_qubic_es_gt5,fc_qubic_es_gt5,dS_dtors,dfc_qubic_es_gt5,S,icoord,n_fcs_qubic_es_gt5,dVds)
    call dvdsQubic(fc_idx_qubic_es_lt5,fc_qubic_es_lt5,dS_dtors,dfc_qubic_es_lt5,S,icoord,n_fcs_qubic_es_lt5,dVds)
    call dvdsQubic(fc_idx_qubic_es_lt2,fc_qubic_es_lt2,dS_dtors,dfc_qubic_es_lt2,S,icoord,n_fcs_qubic_es_lt2,dVds)

	!EA
    call dvdsQubic(fc_idx_qubic_ea_gt5,fc_qubic_ea_gt5,dS_dtors,dfc_qubic_ea_gt5,S,icoord,n_fcs_qubic_ea_gt5,dVds)
    call dvdsQubic(fc_idx_qubic_ea_lt5,fc_qubic_ea_lt5,dS_dtors,dfc_qubic_ea_lt5,S,icoord,n_fcs_qubic_ea_lt5,dVds)
    call dvdsQubic(fc_idx_qubic_ea_lt2,fc_qubic_ea_lt2,dS_dtors,dfc_qubic_ea_lt2,S,icoord,n_fcs_qubic_ea_lt2,dVds)

	!OS
    call dvdsQubic(fc_idx_qubic_os_gt5,fc_qubic_os_gt5,dS_dtors,dfc_qubic_os_gt5,S,icoord,n_fcs_qubic_os_gt5,dVds)
    call dvdsQubic(fc_idx_qubic_os_lt5,fc_qubic_os_lt5,dS_dtors,dfc_qubic_os_lt5,S,icoord,n_fcs_qubic_os_lt5,dVds)
    call dvdsQubic(fc_idx_qubic_os_lt2,fc_qubic_os_lt2,dS_dtors,dfc_qubic_os_lt2,S,icoord,n_fcs_qubic_os_lt2,dVds)

	!EA
    call dvdsQubic(fc_idx_qubic_oa_gt5,fc_qubic_oa_gt5,dS_dtors,dfc_qubic_oa_gt5,S,icoord,n_fcs_qubic_oa_gt5,dVds)
    call dvdsQubic(fc_idx_qubic_oa_lt5,fc_qubic_oa_lt5,dS_dtors,dfc_qubic_oa_lt5,S,icoord,n_fcs_qubic_oa_lt5,dVds)
    call dvdsQubic(fc_idx_qubic_oa_lt2,fc_qubic_oa_lt2,dS_dtors,dfc_qubic_oa_lt2,S,icoord,n_fcs_qubic_oa_lt2,dVds)
    !Quartic:    
	!ES
    call dvdsQuartic(fc_idx_quartic_es_gt5,fc_quartic_es_gt5,dS_dtors,dfc_quartic_es_gt5,S,icoord,n_fcs_quartic_es_gt5,dVds)
    call dvdsQuartic(fc_idx_quartic_es_lt5,fc_quartic_es_lt5,dS_dtors,dfc_quartic_es_lt5,S,icoord,n_fcs_quartic_es_lt5,dVds)
    call dvdsQuartic(fc_idx_quartic_es_lt2,fc_quartic_es_lt2,dS_dtors,dfc_quartic_es_lt2,S,icoord,n_fcs_quartic_es_lt2,dVds)
    call dvdsQuartic(fc_idx_quartic_es_ltpt1,fc_quartic_es_ltpt1,dS_dtors,dfc_quartic_es_ltpt1,S,icoord,n_fcs_quartic_es_ltpt1,dVds)

	!EA
    call dvdsQuartic(fc_idx_quartic_ea_gt5,fc_quartic_ea_gt5,dS_dtors,dfc_quartic_ea_gt5,S,icoord,n_fcs_quartic_ea_gt5,dVds)
    call dvdsQuartic(fc_idx_quartic_ea_lt5,fc_quartic_ea_lt5,dS_dtors,dfc_quartic_ea_lt5,S,icoord,n_fcs_quartic_ea_lt5,dVds)
    call dvdsQuartic(fc_idx_quartic_ea_lt2,fc_quartic_ea_lt2,dS_dtors,dfc_quartic_ea_lt2,S,icoord,n_fcs_quartic_ea_lt2,dVds)
    call dvdsQuartic(fc_idx_quartic_ea_ltpt1,fc_quartic_ea_ltpt1,dS_dtors,dfc_quartic_ea_ltpt1,S,icoord,n_fcs_quartic_ea_ltpt1,dVds)

	!OS
    call dvdsQuartic(fc_idx_quartic_os_gt5,fc_quartic_os_gt5,dS_dtors,dfc_quartic_os_gt5,S,icoord,n_fcs_quartic_os_gt5,dVds)
    call dvdsQuartic(fc_idx_quartic_os_lt5,fc_quartic_os_lt5,dS_dtors,dfc_quartic_os_lt5,S,icoord,n_fcs_quartic_os_lt5,dVds)
    call dvdsQuartic(fc_idx_quartic_os_lt2,fc_quartic_os_lt2,dS_dtors,dfc_quartic_os_lt2,S,icoord,n_fcs_quartic_os_lt2,dVds)
    call dvdsQuartic(fc_idx_quartic_os_ltpt1,fc_quartic_os_ltpt1,dS_dtors,dfc_quartic_os_ltpt1,S,icoord,n_fcs_quartic_os_ltpt1,dVds)

	!OA
    call dvdsQuartic(fc_idx_quartic_oa_gt5,fc_quartic_oa_gt5,dS_dtors,dfc_quartic_oa_gt5,S,icoord,n_fcs_quartic_oa_gt5,dVds)
    call dvdsQuartic(fc_idx_quartic_oa_lt5,fc_quartic_oa_lt5,dS_dtors,dfc_quartic_oa_lt5,S,icoord,n_fcs_quartic_oa_lt5,dVds)
    call dvdsQuartic(fc_idx_quartic_oa_lt2,fc_quartic_oa_lt2,dS_dtors,dfc_quartic_oa_lt2,S,icoord,n_fcs_quartic_oa_lt2,dVds)
    call dvdsQuartic(fc_idx_quartic_oa_ltpt1,fc_quartic_oa_ltpt1,dS_dtors,dfc_quartic_oa_ltpt1,S,icoord,n_fcs_quartic_oa_ltpt1,dVds)
endsubroutine
!****************************************************************************************#
!****************************************************************************************#
subroutine get_ccc_deriv(i,j,k,sinth,costh,coeff,dfc_dtors)
    implicit none
    integer,intent(in):: i,j,k
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    dfc_dtors = 0.d0
    dfc_dtors(0) =  - dble(i) * coeff * sinth(0,i) * costh(1,j) * costh(2,k)
    dfc_dtors(1) =  - dble(j) * coeff * costh(0,i) * sinth(1,j) * costh(2,k) 
    dfc_dtors(2) =  - dble(k) * coeff * costh(0,i) * costh(1,j) * sinth(2,k) 
end subroutine
!****************************************************************************************#
subroutine get_ccs_deriv(i,j,k,sinth,costh,coeff,dfc_dtors)
    implicit none
    integer,intent(in):: i,j,k
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    dfc_dtors = 0.d0
    dfc_dtors(0) =  - dble(i) * coeff * sinth(0,i) * costh(1,j) * sinth(2,k)
    dfc_dtors(1) =  - dble(j) * coeff * costh(0,i) * sinth(1,j) * sinth(2,k) 
    dfc_dtors(2) =  + dble(k) * coeff * costh(0,i) * costh(1,j) * costh(2,k) 
end subroutine
!****************************************************************************************#
subroutine  get_sss_deriv(i,j,k,sinth,costh,coeff,dfc_dtors)
    implicit none
    integer,intent(in):: i,j,k
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    dfc_dtors = 0.d0
    dfc_dtors(0) = + dble(i) * coeff * costh(0,i) * sinth(1,j) * sinth(2,k)
    dfc_dtors(1) = + dble(j) * coeff * sinth(0,i) * costh(1,j) * sinth(2,k)
    dfc_dtors(2) = + dble(k) * coeff * sinth(0,i) * sinth(1,j) * costh(2,k)
end subroutine
!****************************************************************************************#
subroutine  get_ssc_deriv(i,j,k,sinth,costh,coeff,dfc_dtors)
    implicit none
    integer,intent(in):: i,j,k
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    dfc_dtors = 0.d0
    dfc_dtors(0) = + dble(i) * coeff * costh(0,i) * sinth(1,j) * costh(2,k)
    dfc_dtors(1) = + dble(j) * coeff * sinth(0,i) * costh(1,j) * costh(2,k)
    dfc_dtors(2) = - dble(k) * coeff * sinth(0,i) * sinth(1,j) * sinth(2,k)
end subroutine
!****************************************************************************************#
subroutine get_csc_deriv(i,j,k,sinth,costh,coeff,dfc_dtors)
    implicit none
    integer,intent(in):: i,j,k
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    dfc_dtors = 0.d0
    dfc_dtors(0) =  - dble(i) * coeff * sinth(0,i) * sinth(1,j) * costh(2,k)
    dfc_dtors(1) =  + dble(j) * coeff * costh(0,i) * costh(1,j) * costh(2,k)
    dfc_dtors(2) =  - dble(k) * coeff * costh(0,i) * sinth(1,j) * sinth(2,k)
end subroutine
!****************************************************************************************#
subroutine get_css_deriv(i,j,k,sinth,costh,coeff,dfc_dtors)
    implicit none
    integer,intent(in):: i,j,k
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    dfc_dtors = 0.d0
    dfc_dtors(0) =  - dble(i) * coeff * sinth(0,i) * sinth(1,j) * sinth(2,k)
    dfc_dtors(1) =  + dble(j) * coeff * costh(0,i) * costh(1,j) * sinth(2,k)
    dfc_dtors(2) =  + dble(k) * coeff * costh(0,i) * sinth(1,j) * costh(2,k)
end subroutine
!****************************************************************************************#
subroutine get_scs_deriv(i,j,k,sinth,costh,coeff,dfc_dtors)
    implicit none
    integer,intent(in):: i,j,k
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    dfc_dtors = 0.d0
    dfc_dtors(0) = + dble(i) * coeff * costh(0,i) * costh(1,j) * sinth(2,k)
    dfc_dtors(1) = - dble(j) * coeff * sinth(0,i) * sinth(1,j) * sinth(2,k)
    dfc_dtors(2) = + dble(k) * coeff * sinth(0,i) * costh(1,j) * costh(2,k)
end subroutine
!****************************************************************************************#
subroutine get_scc_deriv(i,j,k,sinth,costh,coeff,dfc_dtors)
    implicit none
    integer,intent(in):: i,j,k
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    dfc_dtors = 0.d0
    dfc_dtors(0) = + dble(i) * coeff * costh(0,i) * costh(1,j) * costh(2,k)
    dfc_dtors(1) = - dble(j) * coeff * sinth(0,i) * sinth(1,j) * costh(2,k)
    dfc_dtors(2) = - dble(k) * coeff * sinth(0,i) * costh(1,j) * sinth(2,k)
end subroutine
!****************************************************************************************#
subroutine DTrigseries_es(sinth,costh,Ncoeff,params,stride_arr,fn_order,dfc_dtors)
    implicit none
    integer,intent(in):: fn_order
    integer,intent(in)::Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),params(0:fn_order-1)
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    integer:: i,j,k,ii
    real(8)::dfdt(0:n_fix_tors-1)
	dfdt=0.d0
    dfc_dtors = 0.d0
    do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
			call get_ccc_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt 
        elseif ((ii>=Ncoeff(0)).and.(ii<Ncoeff(1))) then
			call get_ccc_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_ccc_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
        elseif ((ii>=Ncoeff(1)).and.(ii <Ncoeff(2))) then
			call get_css_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
        elseif ((ii>=Ncoeff(2)).and.(ii <Ncoeff(3))) then
			call get_css_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_css_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
        else
			call get_ssc_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_scs_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
		endif
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine DTrigseries_ea(sinth,costh,Ncoeff,params,stride_arr,fn_order,dfc_dtors)
    implicit none
    integer,intent(in):: fn_order
    integer,intent(in)::Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),params(0:fn_order-1)
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    integer:: i,j,k,ii
    real(8)::dfdt(0:n_fix_tors-1)
	dfdt=0.d0
    dfc_dtors = 0.d0
    do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
			call get_ccc_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_ccc_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors - dfdt
        elseif ((ii>=Ncoeff(0)).and.(ii<Ncoeff(1))) then
			call get_css_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_css_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors - dfdt
        else
			call get_scc_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_scs_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors - dfdt
		endif
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine DTrigseries_os(sinth,costh,Ncoeff,params,stride_arr,fn_order,dfc_dtors)
    implicit none
    integer,intent(in):: fn_order
    integer,intent(in)::Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),params(0:fn_order-1)
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    integer:: i,j,k,ii
    real(8)::dfdt(0:n_fix_tors-1)
	dfdt=0.d0
    dfc_dtors = 0.d0
    do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
			call get_scc_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
			call get_scc_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_scc_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
        elseif ((ii>=Ncoeff(1)).and.(ii <Ncoeff(2))) then
			call get_sss_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
        elseif ((ii>=Ncoeff(2)).and.(ii <Ncoeff(3))) then
			call get_sss_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_sss_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
        else
			call get_csc_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_ccs_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
		endif
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine DTrigseries_oa(sinth,costh,Ncoeff,params,stride_arr,fn_order,dfc_dtors)
    implicit none
    integer,intent(in):: fn_order
    integer,intent(in)::Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),params(0:fn_order-1)
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    integer:: i,j,k,ii
    real(8)::dfdt(0:n_fix_tors-1)
	dfdt=0.d0
    dfc_dtors = 0.d0
    do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
			call get_scc_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_scc_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors - dfdt
        elseif ((ii>=Ncoeff(0)).and.(ii<Ncoeff(1))) then
			call get_sss_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_sss_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors - dfdt
        else
			call get_csc_deriv(i,k,j,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
			call get_ccs_deriv(i,j,k,sinth,costh,params(ii),dfdt)
            dfc_dtors = dfc_dtors + dfdt
		endif
	enddo
endsubroutine 
!****************************************************************************************#
subroutine Trigseries_es(sinth,costh,Ncoeff,params,stride_arr,fn_order,fit_val)
    implicit none
    integer,intent(in):: fn_order
    integer,intent(in)::Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),params(0:fn_order-1)
    real(8),intent(out)::fit_val
    integer:: i,j,k,ii
    fit_val = 0.d0
    do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
            fit_val = fit_val + params(ii)  *  costh(0,i) * costh(1,j) * costh(2,k)
        elseif ((ii>=Ncoeff(0)).and.(ii<Ncoeff(1))) then
            fit_val = fit_val + params(ii)  *  costh(0,i) * (costh(1,j) * costh(2,k) + costh(1,k) * costh(2,j))
        elseif ((ii>=Ncoeff(1)).and.(ii <Ncoeff(2))) then
            fit_val = fit_val + params(ii)  * costh(0,i) * sinth(1,j) * sinth(2,k)
        elseif ((ii>=Ncoeff(2)).and.(ii <Ncoeff(3))) then
            fit_val = fit_val + params(ii)  * costh(0,i) * (sinth(1,j) * sinth(2,k) + sinth(1,k) * sinth(2,j))
        else
            fit_val = fit_val + params(ii)  * sinth(0,i) * (sinth(1,k) * costh(2,j) + costh(1,j) * sinth(2,k))
		endif
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine Trigseries_ea(sinth,costh,Ncoeff,params,stride_arr,fn_order,fit_val)
    implicit none
    integer,intent(in):: fn_order
    integer,intent(in)::Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),params(0:fn_order-1)
    real(8),intent(out)::fit_val
    integer:: i,j,k,ii
    fit_val = 0.d0
    do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
            fit_val = fit_val + params(ii)  *   costh(0,i) *(costh(1,j) * costh(2,k)-costh(1,k) * costh(2,j))
        elseif ((ii>=Ncoeff(0)).and.(ii<Ncoeff(1))) then
            fit_val = fit_val + params(ii)  *   costh(0,i) * (sinth(1,j) * sinth(2,k) - sinth(2,j) * sinth(1,k))
        else
            fit_val = fit_val + params(ii)  *   sinth(0,i) * (sinth(1,k) * costh(2,j) - sinth(2,k) * costh(1,j))
		endif
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine Trigseries_os(sinth,costh,Ncoeff,params,stride_arr,fn_order,fit_val)
    implicit none
    integer,intent(in):: fn_order
    integer,intent(in)::Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),params(0:fn_order-1)
    real(8),intent(out)::fit_val
    integer:: i,j,k,ii
    fit_val = 0.d0
    do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
            fit_val = fit_val + params(ii)  *   sinth(0,i)  * costh(1,j) * costh(2,k)
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
            fit_val = fit_val + params(ii)  *   sinth(0,i)  * (costh(1,j) * costh(2,k) + costh(1,k) * costh(2,j))
        elseif ((ii>=Ncoeff(1)).and.(ii <Ncoeff(2))) then
            fit_val = fit_val + params(ii)  *   sinth(0,i)  * sinth(1,j) * sinth(2,k)
        elseif ((ii>=Ncoeff(2)).and.(ii <Ncoeff(3))) then
            fit_val = fit_val + params(ii)  *   sinth(0,i)  * (sinth(1,j) * sinth(2,k) + sinth(1,k) * sinth(2,j))
        else
            fit_val = fit_val + params(ii)  *   costh(0,i)  * (sinth(1,k) * costh(2,j) + costh(1,j) * sinth(2,k))
		endif
	enddo
endsubroutine 
!******************************************************************************************************************************#
subroutine Trigseries_oa(sinth,costh,Ncoeff,params,stride_arr,fn_order,fit_val)
    implicit none
    integer,intent(in):: fn_order
    integer,intent(in)::Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),params(0:fn_order-1)
    real(8),intent(out)::fit_val
    integer:: i,j,k,ii
    fit_val = 0.d0
    do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
            fit_val = fit_val +  params(ii)  *    sinth(0,i) * ( costh(1,j) * costh(2,k) - costh(1,k) * costh(2,j))
        elseif ((ii>=Ncoeff(0)).and.(ii<Ncoeff(1))) then
            fit_val = fit_val +  params(ii)  *    sinth(0,i) * ( sinth(1,j) * sinth(2,k) - sinth(2,j) * sinth(1,k))
        else
            fit_val = fit_val +  params(ii)  *    costh(0,i) * ( sinth(1,k) * costh(2,j) - sinth(2,k) * costh(1,j))
		endif
	enddo
endsubroutine 
!****************************************************************************************#
subroutine  getVt(sinth,costh,Vt)
	implicit none 
	real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(out)::Vt
    Vt = 0.0
    call Trigseries_es(sinth,costh,Ncoeff_pot,Vt_coeff,stride_arr_pot,fn_order_pot,Vt)
endsubroutine 
!****************************************************************************************#
subroutine  getDVt(sinth,costh,dVt)
	implicit none
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
    real(8),intent(out)::dVt(0:n_fix_tors-1)
	integer:: jj,i,j
    dVt = 0.d0;
	call DTrigseries_es(sinth,costh,Ncoeff_pot,Vt_coeff,stride_arr_pot,fn_order_pot,dVt)
endsubroutine 
!Compute Vsav
!****************************************************************************************#
!given q6 and order of expansion compute the force constants:
subroutine get_fc_es(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,fc)
	implicit none 
    integer::ii
	integer,intent(in)::Ncoeff(0:4),fn_order,n_params,stride_arr(0:fn_order-1,0:n_fix_tors-1) 
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::fc(0:n_params-1) 
    fc = 0.d0
    do ii=0,n_params-1
        call Trigseries_es(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,fc(ii)) 
    enddo 
endsubroutine 
!***************************************************************************************#

subroutine get_fc_ea(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,fc)
	implicit none 
    integer::ii
	integer,intent(in)::Ncoeff(0:4),fn_order,n_params,stride_arr(0:fn_order-1,0:n_fix_tors-1) 
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::fc(0:n_params-1) 
    fc = 0.d0
        do ii=0,n_params-1
            call Trigseries_ea(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,fc(ii)) 
        enddo 
endsubroutine 
!***************************************************************************************#

subroutine get_fc_os(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,fc)
	implicit none 
    integer::ii
	integer,intent(in)::Ncoeff(0:4),fn_order,n_params,stride_arr(0:fn_order-1,0:n_fix_tors-1) 
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::fc(0:n_params-1) 
    fc = 0.d0
        do ii=0,n_params-1
            call Trigseries_os(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,fc(ii)) 
        enddo 
endsubroutine 

!***************************************************************************************#
subroutine get_fc_oa(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,fc)
	implicit none 
    integer::ii
	integer,intent(in)::Ncoeff(0:4),fn_order,n_params,stride_arr(0:fn_order-1,0:n_fix_tors-1) 
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::fc(0:n_params-1) 
    fc = 0.d0
        do ii=0,n_params-1
            call Trigseries_oa(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,fc(ii)) 
        enddo 
endsubroutine 

!***************************************************************************************#
!given  function order compute the force constants:
subroutine get_dfc_es(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,dfc)
	implicit none 
    integer::i1,j1,ii,jj
	real(8):: th1,th2
	integer,intent(in):: Ncoeff(0:n_fix_tors-1),fn_order,n_params,stride_arr(0:fn_order-1,0:n_fix_tors-1) 
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::dfc(0:n_params-1,0:n_fix_tors-1)
	real(8)::dfc_dtors(0:n_fix_tors-1)
	dfc = 0.d0;dfc_dtors=0.d0	
    do ii = 0,n_params-1
	    call DTrigseries_es(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,dfc_dtors)
	    dfc(ii,:) = dfc_dtors(:)
	enddo
endsubroutine 
!***************************************************************************************#
!given  function order compute the force constants:
subroutine get_dfc_ea(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,dfc)
	implicit none 
    integer::i1,j1,ii,jj
	real(8):: th1,th2
	integer,intent(in):: Ncoeff(0:4),fn_order,n_params,stride_arr(0:fn_order-1,0:n_fix_tors-1) 
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::dfc(0:n_params-1,0:n_fix_tors-1)
	real(8)::dfc_dtors(0:n_fix_tors-1)
	dfc = 0.d0;dfc_dtors=0.d0	
    do ii = 0,n_params-1
	    call DTrigseries_ea(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,dfc_dtors)
	    dfc(ii,:) = dfc_dtors(:)
	enddo
endsubroutine 
!***************************************************************************************#
!given  function order compute the force constants:
subroutine get_dfc_os(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,dfc)
	implicit none 
    integer::i1,j1,ii,jj
	real(8):: th1,th2
	integer,intent(in):: Ncoeff(0:4),fn_order,n_params,stride_arr(0:fn_order-1,0:n_fix_tors-1) 
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::dfc(0:n_params-1,0:n_fix_tors-1)
	real(8)::dfc_dtors(0:n_fix_tors-1)
	dfc = 0.d0;dfc_dtors=0.d0	
    do ii = 0,n_params-1
	    call DTrigseries_os(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,dfc_dtors)
	    dfc(ii,:) = dfc_dtors(:)
	enddo
endsubroutine 
!***************************************************************************************#
!given  function order compute the force constants:
subroutine get_dfc_oa(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,dfc)
	implicit none 
    integer::i1,j1,ii,jj
	real(8):: th1,th2
	integer,intent(in):: Ncoeff(0:4),fn_order,n_params,stride_arr(0:fn_order-1,0:n_fix_tors-1) 
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::dfc(0:n_params-1,0:n_fix_tors-1)
	real(8)::dfc_dtors(0:n_fix_tors-1)
	dfc = 0.d0;dfc_dtors=0.d0	
    do ii = 0,n_params-1
	    call DTrigseries_oa(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,dfc_dtors)
	    dfc(ii,:) = dfc_dtors(:)
	enddo
endsubroutine 
!*******************************************************************************************************************#
subroutine get_se(sinth,costh,S_e)
	implicit none
	real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(out):: S_e(0:N_int-1)
	integer:: b_ctr,a_ctr,i,idxs
	real(8) :: dum = 0.d0
    S_e = 0.d0
    b_ctr =0;a_ctr = 0
    do i = 0,size(es_fn_idx)-1
        idxs = es_fn_idx(i)
        if (idxs<Nbonds) then
            call  Trigseries_es(sinth,costh,Ncoeff_bonds_es,fitted_bonds_coeff_es(:,i),stride_arr_b_es,fn_order_bonds_es,S_e(idxs))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            call Trigseries_es(sinth,costh,Ncoeff_angs_es,fitted_angs_coeff_es(:,i-b_ctr),stride_arr_a_es,fn_order_angs_es,S_e(idxs))  
            a_ctr = a_ctr + 1
        else
            call Trigseries_es(sinth,costh,Ncoeff_dihs_es,fitted_dihs_coeff_es(:,i-b_ctr-a_ctr),stride_arr_d_es,fn_order_dihs_es,S_e(idxs))    
		endif
	enddo
    b_ctr =0;a_ctr = 0
    do i = 0,size(ea_fn_idx)-1
        idxs = ea_fn_idx(i)
        if (idxs<Nbonds) then
            call  Trigseries_ea(sinth,costh,Ncoeff_bonds_ea,fitted_bonds_coeff_ea(:,i),stride_arr_b_ea,fn_order_bonds_ea,S_e(idxs))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            call Trigseries_ea(sinth,costh,Ncoeff_angs_ea,fitted_angs_coeff_ea(:,i-b_ctr),stride_arr_a_ea,fn_order_angs_ea,S_e(idxs))  
            a_ctr = a_ctr + 1
        else
            call Trigseries_ea(sinth,costh,Ncoeff_dihs_ea,fitted_dihs_coeff_ea(:,i-b_ctr-a_ctr),stride_arr_d_ea,fn_order_dihs_ea,S_e(idxs))    
		endif
	enddo
    b_ctr =0;a_ctr = 0
    do i = 0,size(os_fn_idx)-1
        idxs = os_fn_idx(i)
        !print(idxs,fitted_bonds_coeff_os[i])
        if (idxs<Nbonds) then
            call Trigseries_os(sinth,costh,Ncoeff_bonds_os,fitted_bonds_coeff_os(:,i),stride_arr_b_os,fn_order_bonds_os,S_e(idxs))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            call Trigseries_os(sinth,costh,Ncoeff_angs_os,fitted_angs_coeff_os(:,i-b_ctr),stride_arr_a_os,fn_order_angs_os,S_e(idxs))  
            a_ctr = a_ctr + 1
        else
            call Trigseries_os(sinth,costh,Ncoeff_dihs_os,fitted_dihs_coeff_os(:,i-b_ctr-a_ctr),stride_arr_d_os,fn_order_dihs_os,S_e(idxs) )    
        endif
    enddo
    b_ctr =0;a_ctr = 0
    do i = 0,size(oa_fn_idx)-1
        idxs = oa_fn_idx(i)
        !print(idxs,fitted_bonds_coeff_oa[i])
        if (idxs<Nbonds) then
            call Trigseries_oa(sinth,costh,Ncoeff_bonds_oa,fitted_bonds_coeff_oa(:,i),stride_arr_b_oa,fn_order_bonds_oa,S_e(idxs))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            call Trigseries_oa(sinth,costh,Ncoeff_angs_oa,fitted_angs_coeff_oa(:,i-b_ctr),stride_arr_a_oa,fn_order_angs_oa,S_e(idxs))  
            a_ctr = a_ctr + 1
        else
            call Trigseries_oa(sinth,costh,Ncoeff_dihs_oa,fitted_dihs_coeff_oa(:,i-b_ctr-a_ctr),stride_arr_d_oa,fn_order_dihs_oa,S_e(idxs) )    
        endif
    enddo
endsubroutine 
!****************************************************************************************************************************************************#
subroutine  get_all_dfc(sinth,costh,dfc_quad_es,dfc_quad_os,dfc_quad_ea,dfc_quad_oa,dfc_qubic_es_gt5,dfc_qubic_es_lt5,dfc_qubic_es_lt2,dfc_qubic_ea_gt5,dfc_qubic_ea_lt5,dfc_qubic_ea_lt2,dfc_qubic_os_gt5,dfc_qubic_os_lt5,dfc_qubic_os_lt2,dfc_qubic_oa_gt5,dfc_qubic_oa_lt5,dfc_qubic_oa_lt2,dfc_quartic_es_gt5,dfc_quartic_es_lt5,dfc_quartic_es_lt2,dfc_quartic_ea_gt5,dfc_quartic_ea_lt5,dfc_quartic_ea_lt2,dfc_quartic_os_gt5,dfc_quartic_os_lt5,dfc_quartic_os_lt2,dfc_quartic_oa_gt5,dfc_quartic_oa_lt5,dfc_quartic_oa_lt2) 
	implicit none
	real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(out):: dfc_quad_es(0:n_fcs_quad_es-1,0:n_fix_tors-1),dfc_quad_os(0:n_fcs_quad_os-1,0:n_fix_tors-1),dfc_quad_ea(0:n_fcs_quad_ea-1,0:n_fix_tors-1),dfc_quad_oa(0:n_fcs_quad_oa-1,0:n_fix_tors-1)
	real(8),intent(out):: dfc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1,0:n_fix_tors-1),dfc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1,0:n_fix_tors-1),dfc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1,0:n_fix_tors-1),dfc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1,0:n_fix_tors-1), dfc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1,0:n_fix_tors-1)
	real(8),intent(out):: dfc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1,0:n_fix_tors-1),dfc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1,0:n_fix_tors-1),dfc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1,0:n_fix_tors-1),dfc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1,0:n_fix_tors-1),dfc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1,0:n_fix_tors-1),dfc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1,0:n_fix_tors-1),dfc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1,0:n_fix_tors-1),dfc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1,0:n_fix_tors-1),dfc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1,0:n_fix_tors-1),dfc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1,0:n_fix_tors-1),dfc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1,0:n_fix_tors-1),dfc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1,0:n_fix_tors-1),dfc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1,0:n_fix_tors-1) 
	!Initialize to zero
	dfc_quad_es = 0.d0;dfc_quad_os = 0.d0; 
	dfc_quad_ea = 0.d0;dfc_quad_oa = 0.d0; 
	dfc_qubic_es_gt5 = 0.d0;dfc_qubic_es_lt5 = 0.d0;dfc_qubic_es_lt2 = 0.d0;
	dfc_qubic_ea_gt5 = 0.d0;dfc_qubic_ea_lt5 = 0.d0;dfc_qubic_ea_lt2 = 0.d0;
	dfc_qubic_os_gt5 = 0.d0;dfc_qubic_os_lt5 = 0.d0;dfc_qubic_os_lt2 = 0.d0;
	dfc_qubic_oa_gt5 = 0.d0;dfc_qubic_oa_lt5 = 0.d0;dfc_qubic_oa_lt2 = 0.d0;
	dfc_quartic_es_gt5 = 0.d0; dfc_quartic_es_lt5 = 0.d0;dfc_quartic_es_lt2 = 0.d0;
	dfc_quartic_ea_gt5 = 0.d0; dfc_quartic_ea_lt5 = 0.d0;dfc_quartic_ea_lt2 = 0.d0;
	dfc_quartic_os_lt5 = 0.d0;dfc_quartic_os_lt2 = 0.d0; 
	dfc_quartic_oa_lt5 = 0.d0;dfc_quartic_oa_lt2 = 0.d0; 
    !For quadratic
    !ES
	call get_dfc_es(sinth,costh,Ncoeff_fij_es,fitted_fc_coeff_ij_es,stride_arr_quad_es,fn_order_fij_es,n_fcs_quad_es,dfc_quad_es)
    !EA
	call get_dfc_ea(sinth,costh,Ncoeff_fij_ea,fitted_fc_coeff_ij_ea,stride_arr_quad_ea,fn_order_fij_ea,n_fcs_quad_ea,dfc_quad_ea)
    !OS
    call get_dfc_os(sinth,costh,Ncoeff_fij_os,fitted_fc_coeff_ij_os,stride_arr_quad_os,fn_order_fij_os,n_fcs_quad_os,dfc_quad_os)
    !OA
    call get_dfc_oa(sinth,costh,Ncoeff_fij_oa,fitted_fc_coeff_ij_oa,stride_arr_quad_oa,fn_order_fij_oa,n_fcs_quad_oa,dfc_quad_oa)
   	 
    !For qubic
    !ES
    call get_dfc_es(sinth,costh,Ncoeff_fijk_es_gt5,fitted_fc_coeff_ijk_es_gt5,stride_arr_qubic_es_gt5,fn_order_fijk_es_gt5,n_fcs_qubic_es_gt5,dfc_qubic_es_gt5)
   	call get_dfc_es(sinth,costh,Ncoeff_fijk_es_lt5,fitted_fc_coeff_ijk_es_lt5,stride_arr_qubic_es_lt5,fn_order_fijk_es_lt5,n_fcs_qubic_es_lt5,dfc_qubic_es_lt5)
    call get_dfc_es(sinth,costh,Ncoeff_fijk_es_lt2,fitted_fc_coeff_ijk_es_lt2,stride_arr_qubic_es_lt2,fn_order_fijk_es_lt2,n_fcs_qubic_es_lt2,dfc_qubic_es_lt2)

    !EA
    call get_dfc_ea(sinth,costh,Ncoeff_fijk_ea_gt5,fitted_fc_coeff_ijk_ea_gt5,stride_arr_qubic_ea_gt5,fn_order_fijk_ea_gt5,n_fcs_qubic_ea_gt5,dfc_qubic_ea_gt5)
   	call get_dfc_ea(sinth,costh,Ncoeff_fijk_ea_lt5,fitted_fc_coeff_ijk_ea_lt5,stride_arr_qubic_ea_lt5,fn_order_fijk_ea_lt5,n_fcs_qubic_ea_lt5,dfc_qubic_ea_lt5)
    call get_dfc_ea(sinth,costh,Ncoeff_fijk_ea_lt2,fitted_fc_coeff_ijk_ea_lt2,stride_arr_qubic_ea_lt2,fn_order_fijk_ea_lt2,n_fcs_qubic_ea_lt2,dfc_qubic_ea_lt2)
    
    !OS
    call get_dfc_os(sinth,costh,Ncoeff_fijk_os_gt5,fitted_fc_coeff_ijk_os_gt5,stride_arr_qubic_os_gt5,fn_order_fijk_os_gt5,n_fcs_qubic_os_gt5,dfc_qubic_os_gt5)
    call get_dfc_os(sinth,costh,Ncoeff_fijk_os_lt5,fitted_fc_coeff_ijk_os_lt5,stride_arr_qubic_os_lt5,fn_order_fijk_os_lt5,n_fcs_qubic_os_lt5,dfc_qubic_os_lt5)
    call get_dfc_os(sinth,costh,Ncoeff_fijk_os_lt2,fitted_fc_coeff_ijk_os_lt2,stride_arr_qubic_os_lt2,fn_order_fijk_os_lt2,n_fcs_qubic_os_lt2,dfc_qubic_os_lt2)
    
    !OA
    call get_dfc_oa(sinth,costh,Ncoeff_fijk_oa_gt5,fitted_fc_coeff_ijk_oa_gt5,stride_arr_qubic_oa_gt5,fn_order_fijk_oa_gt5,n_fcs_qubic_oa_gt5,dfc_qubic_oa_gt5)
    call get_dfc_oa(sinth,costh,Ncoeff_fijk_oa_lt5,fitted_fc_coeff_ijk_oa_lt5,stride_arr_qubic_oa_lt5,fn_order_fijk_oa_lt5,n_fcs_qubic_oa_lt5,dfc_qubic_oa_lt5)
    call get_dfc_oa(sinth,costh,Ncoeff_fijk_oa_lt2,fitted_fc_coeff_ijk_oa_lt2,stride_arr_qubic_oa_lt2,fn_order_fijk_oa_lt2,n_fcs_qubic_oa_lt2,dfc_qubic_oa_lt2)
    
    !For quartic
    !ES
    call get_dfc_es(sinth,costh,Ncoeff_fijkl_es_gt5,fitted_fc_coeff_ijkl_es_gt5,stride_arr_quartic_es_gt5,fn_order_fijkl_es_gt5,n_fcs_quartic_es_gt5,dfc_quartic_es_gt5)
    call get_dfc_es(sinth,costh,Ncoeff_fijkl_es_lt5,fitted_fc_coeff_ijkl_es_lt5,stride_arr_quartic_es_lt5,fn_order_fijkl_es_lt5,n_fcs_quartic_es_lt5,dfc_quartic_es_lt5)
    call get_dfc_es(sinth,costh,Ncoeff_fijkl_es_lt2,fitted_fc_coeff_ijkl_es_lt2,stride_arr_quartic_es_lt2,fn_order_fijkl_es_lt2,n_fcs_quartic_es_lt2,dfc_quartic_es_lt2)
    
    !EA
    call get_dfc_ea(sinth,costh,Ncoeff_fijkl_ea_gt5,fitted_fc_coeff_ijkl_ea_gt5,stride_arr_quartic_ea_gt5,fn_order_fijkl_ea_gt5,n_fcs_quartic_ea_gt5,dfc_quartic_ea_gt5)
    call get_dfc_ea(sinth,costh,Ncoeff_fijkl_ea_lt5,fitted_fc_coeff_ijkl_ea_lt5,stride_arr_quartic_ea_lt5,fn_order_fijkl_ea_lt5,n_fcs_quartic_ea_lt5,dfc_quartic_ea_lt5)
    call get_dfc_ea(sinth,costh,Ncoeff_fijkl_ea_lt2,fitted_fc_coeff_ijkl_ea_lt2,stride_arr_quartic_ea_lt2,fn_order_fijkl_ea_lt2,n_fcs_quartic_ea_lt2,dfc_quartic_ea_lt2)
    
    !OS
    call get_dfc_os(sinth,costh,Ncoeff_fijkl_os_gt5,fitted_fc_coeff_ijkl_os_gt5,stride_arr_quartic_os_gt5,fn_order_fijkl_os_gt5,n_fcs_quartic_os_gt5,dfc_quartic_os_gt5)
    call get_dfc_os(sinth,costh,Ncoeff_fijkl_os_lt5,fitted_fc_coeff_ijkl_os_lt5,stride_arr_quartic_os_lt5,fn_order_fijkl_os_lt5,n_fcs_quartic_os_lt5,dfc_quartic_os_lt5)
    call get_dfc_os(sinth,costh,Ncoeff_fijkl_os_lt2,fitted_fc_coeff_ijkl_os_lt2,stride_arr_quartic_os_lt2,fn_order_fijkl_os_lt2,n_fcs_quartic_os_lt2,dfc_quartic_os_lt2)

    !OA
    call get_dfc_oa(sinth,costh,Ncoeff_fijkl_oa_gt5,fitted_fc_coeff_ijkl_oa_gt5,stride_arr_quartic_oa_gt5,fn_order_fijkl_oa_gt5,n_fcs_quartic_oa_gt5,dfc_quartic_oa_gt5)
    call get_dfc_oa(sinth,costh,Ncoeff_fijkl_oa_lt5,fitted_fc_coeff_ijkl_oa_lt5,stride_arr_quartic_oa_lt5,fn_order_fijkl_oa_lt5,n_fcs_quartic_oa_lt5,dfc_quartic_oa_lt5)
    call get_dfc_oa(sinth,costh,Ncoeff_fijkl_oa_lt2,fitted_fc_coeff_ijkl_oa_lt2,stride_arr_quartic_oa_lt2,fn_order_fijkl_oa_lt2,n_fcs_quartic_oa_lt2,dfc_quartic_oa_lt2)
endsubroutine
!****************************************************************************************************************************************************#
subroutine get_all_fc(sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1) 
	implicit none 
	real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(out):: fc_quad_es(0:n_fcs_quad_es-1),fc_quad_os(0:n_fcs_quad_os-1)
	real(8),intent(out):: fc_quad_ea(0:n_fcs_quad_ea-1),fc_quad_oa(0:n_fcs_quad_oa-1)
	real(8),intent(out)::fc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1),fc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1),fc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1)
	real(8),intent(out)::fc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1),fc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1),fc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1)
	real(8),intent(out):: fc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1),fc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1),fc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1) 
	real(8),intent(out):: fc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1),fc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1),fc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1) 
	real(8),intent(out):: fc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1),fc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1),fc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1) 
	real(8),intent(out):: fc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1),fc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1),fc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1) 
	real(8),intent(out):: fc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1),fc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1),fc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1) 
	real(8),intent(out):: fc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1),fc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1),fc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1) 
	real(8),intent(out):: fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
	real(8),intent(out):: fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1) 
	!Initialize to zero
	fc_quad_ea = 0.d0;fc_quad_oa = 0.d0; 
	fc_quad_es = 0.d0;fc_quad_os = 0.d0; 
	fc_qubic_es_gt5 = 0.d0;fc_qubic_es_lt5 = 0.d0;fc_qubic_es_lt2 = 0.d0;
	fc_qubic_ea_gt5 = 0.d0;fc_qubic_ea_lt5 = 0.d0;fc_qubic_ea_lt2 = 0.d0;
	fc_qubic_os_gt5 = 0.d0;fc_qubic_os_lt5 = 0.d0;fc_qubic_os_lt2 = 0.d0;
	fc_qubic_oa_gt5 = 0.d0;fc_qubic_oa_lt5 = 0.d0;fc_qubic_oa_lt2 = 0.d0;
	fc_quartic_es_gt5 = 0.d0;fc_quartic_es_lt5 = 0.d0;fc_quartic_es_lt2 = 0.d0;
	fc_quartic_ea_gt5 = 0.d0;fc_quartic_ea_lt5 = 0.d0;fc_quartic_ea_lt2 = 0.d0;
	fc_quartic_os_gt5 = 0.d0;fc_quartic_os_lt5 = 0.d0;fc_quartic_os_lt2 = 0.d0; 
	fc_quartic_oa_gt5 = 0.d0;fc_quartic_oa_lt5 = 0.d0;fc_quartic_oa_lt2 = 0.d0; 
    !For quadratic
    !ES
	!print*,'quadratic es'
    call get_fc_es(sinth,costh,Ncoeff_fij_es,fitted_fc_coeff_ij_es,stride_arr_quad_es,fn_order_fij_es,n_fcs_quad_es,fc_quad_es)
	!print*,fc_quad_es
	!print*,'quadratic ea'
    !EA
    call get_fc_ea(sinth,costh,Ncoeff_fij_ea,fitted_fc_coeff_ij_ea,stride_arr_quad_ea,fn_order_fij_ea,n_fcs_quad_ea,fc_quad_ea)
    !OS
    call get_fc_os(sinth,costh,Ncoeff_fij_os,fitted_fc_coeff_ij_os,stride_arr_quad_os,fn_order_fij_os,n_fcs_quad_os,fc_quad_os)
    !OA
    call get_fc_oa(sinth,costh,Ncoeff_fij_oa,fitted_fc_coeff_ij_oa,stride_arr_quad_oa,fn_order_fij_oa,n_fcs_quad_oa,fc_quad_oa)
    
    !For qubic
    !ES
    call get_fc_es(sinth,costh,Ncoeff_fijk_es_gt5,fitted_fc_coeff_ijk_es_gt5,stride_arr_qubic_es_gt5,fn_order_fijk_es_gt5,n_fcs_qubic_es_gt5,fc_qubic_es_gt5)
   	call get_fc_es(sinth,costh,Ncoeff_fijk_es_lt5,fitted_fc_coeff_ijk_es_lt5,stride_arr_qubic_es_lt5,fn_order_fijk_es_lt5,n_fcs_qubic_es_lt5,fc_qubic_es_lt5)
    call get_fc_es(sinth,costh,Ncoeff_fijk_es_lt2,fitted_fc_coeff_ijk_es_lt2,stride_arr_qubic_es_lt2,fn_order_fijk_es_lt2,n_fcs_qubic_es_lt2,fc_qubic_es_lt2)
    
    !EA
    call get_fc_ea(sinth,costh,Ncoeff_fijk_ea_gt5,fitted_fc_coeff_ijk_ea_gt5,stride_arr_qubic_ea_gt5,fn_order_fijk_ea_gt5,n_fcs_qubic_ea_gt5,fc_qubic_ea_gt5)
   	call get_fc_ea(sinth,costh,Ncoeff_fijk_ea_lt5,fitted_fc_coeff_ijk_ea_lt5,stride_arr_qubic_ea_lt5,fn_order_fijk_ea_lt5,n_fcs_qubic_ea_lt5,fc_qubic_ea_lt5)
    call get_fc_ea(sinth,costh,Ncoeff_fijk_ea_lt2,fitted_fc_coeff_ijk_ea_lt2,stride_arr_qubic_ea_lt2,fn_order_fijk_ea_lt2,n_fcs_qubic_ea_lt2,fc_qubic_ea_lt2)
    
    !OS
    call get_fc_os(sinth,costh,Ncoeff_fijk_os_gt5,fitted_fc_coeff_ijk_os_gt5,stride_arr_qubic_os_gt5,fn_order_fijk_os_gt5,n_fcs_qubic_os_gt5,fc_qubic_os_gt5)
    call get_fc_os(sinth,costh,Ncoeff_fijk_os_lt5,fitted_fc_coeff_ijk_os_lt5,stride_arr_qubic_os_lt5,fn_order_fijk_os_lt5,n_fcs_qubic_os_lt5,fc_qubic_os_lt5)
    call get_fc_os(sinth,costh,Ncoeff_fijk_os_lt2,fitted_fc_coeff_ijk_os_lt2,stride_arr_qubic_os_lt2,fn_order_fijk_os_lt2,n_fcs_qubic_os_lt2,fc_qubic_os_lt2)
    
    !OA
    call get_fc_oa(sinth,costh,Ncoeff_fijk_oa_gt5,fitted_fc_coeff_ijk_oa_gt5,stride_arr_qubic_oa_gt5,fn_order_fijk_oa_gt5,n_fcs_qubic_oa_gt5,fc_qubic_oa_gt5)
    call get_fc_oa(sinth,costh,Ncoeff_fijk_oa_lt5,fitted_fc_coeff_ijk_oa_lt5,stride_arr_qubic_oa_lt5,fn_order_fijk_oa_lt5,n_fcs_qubic_oa_lt5,fc_qubic_oa_lt5)
    call get_fc_oa(sinth,costh,Ncoeff_fijk_oa_lt2,fitted_fc_coeff_ijk_oa_lt2,stride_arr_qubic_oa_lt2,fn_order_fijk_oa_lt2,n_fcs_qubic_oa_lt2,fc_qubic_oa_lt2)
    
    !For quartic
    !ES
    call get_fc_es(sinth,costh,Ncoeff_fijkl_es_gt5,fitted_fc_coeff_ijkl_es_gt5,stride_arr_quartic_es_gt5,fn_order_fijkl_es_gt5,n_fcs_quartic_es_gt5,fc_quartic_es_gt5)
    call get_fc_es(sinth,costh,Ncoeff_fijkl_es_lt5,fitted_fc_coeff_ijkl_es_lt5,stride_arr_quartic_es_lt5,fn_order_fijkl_es_lt5,n_fcs_quartic_es_lt5,fc_quartic_es_lt5)
    call get_fc_es(sinth,costh,Ncoeff_fijkl_es_lt2,fitted_fc_coeff_ijkl_es_lt2,stride_arr_quartic_es_lt2,fn_order_fijkl_es_lt2,n_fcs_quartic_es_lt2,fc_quartic_es_lt2)
	fc_quartic_es_ltpt1 = fitted_fc_coeff_ijkl_es_ltpt1
    
    !EA
    call get_fc_ea(sinth,costh,Ncoeff_fijkl_ea_gt5,fitted_fc_coeff_ijkl_ea_gt5,stride_arr_quartic_ea_gt5,fn_order_fijkl_ea_gt5,n_fcs_quartic_ea_gt5,fc_quartic_ea_gt5)
    call get_fc_ea(sinth,costh,Ncoeff_fijkl_ea_lt5,fitted_fc_coeff_ijkl_ea_lt5,stride_arr_quartic_ea_lt5,fn_order_fijkl_ea_lt5,n_fcs_quartic_ea_lt5,fc_quartic_ea_lt5)
    call get_fc_ea(sinth,costh,Ncoeff_fijkl_ea_lt2,fitted_fc_coeff_ijkl_ea_lt2,stride_arr_quartic_ea_lt2,fn_order_fijkl_ea_lt2,n_fcs_quartic_ea_lt2,fc_quartic_ea_lt2)
	fc_quartic_ea_ltpt1 = fitted_fc_coeff_ijkl_ea_ltpt1
    
    !OS
    call get_fc_os(sinth,costh,Ncoeff_fijkl_os_gt5,fitted_fc_coeff_ijkl_os_gt5,stride_arr_quartic_os_gt5,fn_order_fijkl_os_gt5,n_fcs_quartic_os_gt5,fc_quartic_os_gt5)
    call get_fc_os(sinth,costh,Ncoeff_fijkl_os_lt5,fitted_fc_coeff_ijkl_os_lt5,stride_arr_quartic_os_lt5,fn_order_fijkl_os_lt5,n_fcs_quartic_os_lt5,fc_quartic_os_lt5)
    call get_fc_os(sinth,costh,Ncoeff_fijkl_os_lt2,fitted_fc_coeff_ijkl_os_lt2,stride_arr_quartic_os_lt2,fn_order_fijkl_os_lt2,n_fcs_quartic_os_lt2,fc_quartic_os_lt2)
	fc_quartic_os_ltpt1 = fitted_fc_coeff_ijkl_os_ltpt1
    !OA
    call get_fc_oa(sinth,costh,Ncoeff_fijkl_oa_gt5,fitted_fc_coeff_ijkl_oa_gt5,stride_arr_quartic_oa_gt5,fn_order_fijkl_oa_gt5,n_fcs_quartic_oa_gt5,fc_quartic_oa_gt5)
    call get_fc_oa(sinth,costh,Ncoeff_fijkl_oa_lt5,fitted_fc_coeff_ijkl_oa_lt5,stride_arr_quartic_oa_lt5,fn_order_fijkl_oa_lt5,n_fcs_quartic_oa_lt5,fc_quartic_oa_lt5)
    call get_fc_oa(sinth,costh,Ncoeff_fijkl_oa_lt2,fitted_fc_coeff_ijkl_oa_lt2,stride_arr_quartic_oa_lt2,fn_order_fijkl_oa_lt2,n_fcs_quartic_oa_lt2,fc_quartic_oa_lt2)
	fc_quartic_oa_ltpt1 = fitted_fc_coeff_ijkl_oa_ltpt1
endsubroutine

!!****************************************************************************************#
!Compute dsdx
subroutine  get_dsdx(X,dsdx)
    implicit none
    integer::i,j,k,istep
    real(8)::d1,d2,d3,d4
    real(8)::Xtmp(0:Ncarts-1),S(0:N_int-1),Stmp(0:N_int-1),dS_ar(0:4,0:N_int-1)
    real(8),intent(in)::X(0:Ncarts-1)
    real(8),intent(out)::dsdx(0:Ncarts-1,0:N_int-1)
    dsdx = 0.d0;Xtmp = 0.d0;S = 0.d0;Stmp = 0.d0
    d1 =0.d0;d2 = 0.d0;d3=0.d0;d4=0.d0
    !call cart_to_internal(X,S)
    do i =0,Ncarts-1
        dS_ar = 0.d0
        Xtmp =X
        do istep = -2,2
            Xtmp(i) = X(i) + dble(istep) * dx
            call cart_to_internal(Xtmp,Stmp)
            dS_ar(istep+2,:)   = Stmp(:)
        enddo
        do j = Nbonds+Nang,N_int-1
            if (abs(abs(dS_ar(2,j))-180.d0)<5.d0) then
                d1 = dS_ar(1,j)-dS_ar(0,j)
                d2 = dS_ar(2,j)-dS_ar(1,j)
                d3 = dS_ar(3,j)-dS_ar(2,j)
                d4 = dS_ar(4,j)-dS_ar(3,j)
                if (abs(d1)>10.d0) then
                   if (d1>0.d0) then
                        do k = 1,4
                           dS_ar(k,j) = dS_ar(k,j) - 360.d0
                        enddo
                   else
                        do k = 1,4
                           dS_ar(k,j) = dS_ar(k,j) + 360.d0
                        enddo
                   endif
                else if (abs(d2)>10.d0) then
                   if (d2>0.d0) then
                        do k = 2,4
                           dS_ar(k,j) = dS_ar(k,j) - 360.d0
                        enddo
                   else
                        do k = 2,4
                           dS_ar(k,j) = dS_ar(k,j) + 360.d0
                        enddo
                   endif

                else if (abs(d3)>10.d0) then
                   if (d3>0.d0) then
                        do k = 3,4
                           dS_ar(k,j) = dS_ar(k,j) - 360.d0
                        enddo
                   else
                        do k = 3,4
                           dS_ar(k,j) = dS_ar(k,j) + 360.d0
                        enddo
                   endif

                else if (abs(d4)>10.d0) then
                   k = 4
                   if (d4>0.d0) then
                           dS_ar(k,j) = dS_ar(k,j) - 360.d0
                   else
                           dS_ar(k,j) = dS_ar(k,j) + 360.d0
                   endif

                endif
            endif
        enddo
        dsdx(i,:) = dx_inv * matmul(five_pt_1st_deriv(:),dS_ar(:,:))
    enddo
endsubroutine 
!****************************************************************************************#
subroutine get_csarr(fix_tors,max_id,sinth,costh)
    implicit  none
    integer:: i,j,k
    integer,intent(in)  :: max_id(0:n_fix_tors-1)
    real(8),intent(in)  :: fix_tors(0:n_fix_tors-1)
    real(8),intent(out) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
    sinth = 0.d0; costh =0.d0
    sinth(:,0) = 0.d0; costh(:,0) = 1.d0
    do i = 0,n_fix_tors-1
    	do j = 1,max_id(i)
            	sinth(i,j)  = dsin(dble(j)*fix_tors(i))
            	costh(i,j)  = dcos(dble(j)*fix_tors(i))
        enddo
    enddo

endsubroutine
!****************************************************************************************#
!Convert internal gradients to cartesian.

subroutine int_to_cart_grad(dsdx,dVds,dVdx)
	implicit none 
	integer::i 
	real(8),intent(in)::dsdx(0:Ncarts-1,0:N_int-1),dVds(0:N_int-1)
	real(8),intent(out)::dVdx(0:Ncarts-1)
    dVdx =0.d0 
    do i = 0,Ncarts-1
       dVdx(i) = dot_product(dVds(:),dsdx(i,:))
	enddo
endsubroutine 
!****************************************************************************************#
subroutine scale_coord(S)
implicit none 
integer::i
real(8),intent(inout)::S(0:N_int-1)
do i = 0,N_int-1
	S(i) = S(i) * dim_scal_factor_inv(i)
	if (i>=Nbonds) then
		S(i) = S(i) * degreeToRadian
	endif
enddo
endsubroutine 
!****************************************************************************************#
subroutine get_pot(X,pot,Vt_out)
	implicit none 
	real(8),intent(in):: X(0:Ncarts-1)
	real(8) :: Xtmp(0:Ncarts-1)
	real(8),intent(out)::pot,Vt_out
	real(8) :: S(0:N_int-1) ,sym_S(0:N_int-1) ,S_eq(0:N_int-1) ,S_diff(0:N_int-1) ,sym_diff(0:N_int-1) ,sym_seq(0:N_int-1),dVdsym(0:N_int-1), dVds(0:N_int-1), dsdx(0:Ncarts-1,0:N_int-1) 
    real(8) :: fix_tors(0:n_fix_tors-1),sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8):: fc_quad_es(0:n_fcs_quad_es-1),fc_quad_os(0:n_fcs_quad_os-1)
	real(8):: fc_quad_ea(0:n_fcs_quad_ea-1),fc_quad_oa(0:n_fcs_quad_oa-1)
	real(8):: fc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1),fc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1),fc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1)
	real(8):: fc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1),fc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1),fc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1)
	real(8):: fc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1),fc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1),fc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1) 
	real(8):: fc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1),fc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1),fc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1) 
	real(8):: fc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1),fc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1),fc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1) 
	real(8):: fc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1),fc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1),fc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1) 
	real(8):: fc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1),fc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1),fc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1) 
	real(8):: fc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1),fc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1),fc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1) 
	real(8)::fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1)
	real(8)::fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
	integer::i
	fc_quad_ea = 0.d0;fc_quad_oa = 0.d0; 
	fc_quad_es = 0.d0;fc_quad_os = 0.d0; 
	fc_qubic_es_gt5 = 0.d0;fc_qubic_es_lt5 = 0.d0;fc_qubic_es_lt2 = 0.d0;
	fc_qubic_ea_gt5 = 0.d0;fc_qubic_ea_lt5 = 0.d0;fc_qubic_ea_lt2 = 0.d0;
	fc_qubic_os_gt5 = 0.d0;fc_qubic_os_lt5 = 0.d0;fc_qubic_os_lt2 = 0.d0;
	fc_qubic_oa_gt5 = 0.d0;fc_qubic_oa_lt5 = 0.d0;fc_qubic_oa_lt2 = 0.d0;
	fc_quartic_es_gt5 = 0.d0;fc_quartic_es_lt5 = 0.d0;fc_quartic_es_lt2 = 0.d0;
	fc_quartic_ea_gt5 = 0.d0;fc_quartic_ea_lt5 = 0.d0;fc_quartic_ea_lt2 = 0.d0;
	fc_quartic_os_lt5 = 0.d0;fc_quartic_os_lt2 = 0.d0; 
	fc_quartic_oa_lt5 = 0.d0;fc_quartic_oa_lt2 = 0.d0; 
	pot =0.d0;Vt_out=0.d0;sinth = 0.d0;costh=0.d0 
	!Convert Bohr to Angstrom
	Xtmp =  X * bohrToAng 
	!Convert cartesian to internal
	call cart_to_internal(Xtmp,S)
	!Fixed torsions array
	fix_tors = S(tors_idx) * degreeToRadian
	!Get sin and cosine arrays
	call get_csarr(fix_tors,max_id,sinth,costh)
	!Get equlibrium S
	call get_se(sinth,costh,sym_Seq)
	!Convert internal to symmetric internal
	sym_Seq(11) = sym_Seq(11) +360.0
	sym_Seq(18) = sym_Seq(18) +360.0
	S_eq  =  matmul(sym_mat_transpose,sym_Seq)
	do i = Nbonds+Nang,N_int-1
        !Check if S(i) is close to 180 degree
        !if (abs(abs(S(i))-180.d0)<90.d0) then
            !If S(i) and S_eq(i) have opposite signs
            if (S(i)*S_eq(i)<0.d0) then
               if (S_eq(i)>0.d0) then
                  S(i) = S(i) + 360.d0
               else
                  S(i) = S(i) - 360.d0
               endif
            endif
        !endif
    enddo
    s_diff = S-S_eq
    Sym_diff = matmul(sym_mat,S_diff)
	!Convert to scaled coordinates
	call scale_coord(Sym_diff)

	!compute fitted force constants(fc)
	call get_all_fc(sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1) 

	!Compute derivative of fitted fc (dfc)
		call compute_pot(sym_diff,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1 ,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,pot,Vt_out) 
	pot = pot * CminvToHartee
    Vt_out = Vt_out * CminvToHartee
endsubroutine 
!****************************************************************************************#
subroutine get_force(X,dVdx,pot,flag)
	implicit none 
	real(8),intent(in) :: X(0:Ncarts-1)
	real(8) :: Xtmp(0:Ncarts-1)
	integer,intent(in) ::flag
	real(8),intent(out) :: dVdx(0:Ncarts-1),pot
	integer::i
	real(8) :: Vt_out,S(0:N_int-1) ,sym_S(0:N_int-1) ,sym_Seq(0:N_int-1) ,S_eq(0:N_int-1) ,S_diff(0:N_int-1) ,sym_diff(0:N_int-1) ,dVdsym(0:N_int-1), dVds(0:N_int-1), dsdx(0:Ncarts-1,0:N_int-1) 
    real(8) :: fix_tors(0:n_fix_tors-1),sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8):: fc_quad_es(0:n_fcs_quad_es-1),fc_quad_os(0:n_fcs_quad_os-1)
	real(8):: fc_quad_ea(0:n_fcs_quad_ea-1),fc_quad_oa(0:n_fcs_quad_oa-1)
	real(8):: fc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1),fc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1),fc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1)
	real(8):: fc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1),fc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1),fc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1)
	real(8):: fc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1),fc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1),fc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1) 
	real(8):: fc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1),fc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1),fc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1) 
	real(8):: fc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1),fc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1),fc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1) 
	real(8):: fc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1),fc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1),fc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1) 
	real(8):: fc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1),fc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1),fc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1) 
	real(8):: fc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1),fc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1),fc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1) 
	real(8):: dfc_quad_es(0:n_fcs_quad_es-1,0:n_fix_tors-1),dfc_quad_os(0:n_fcs_quad_os-1,0:n_fix_tors-1),dfc_quad_ea(0:n_fcs_quad_ea-1,0:n_fix_tors-1),dfc_quad_oa(0:n_fcs_quad_oa-1,0:n_fix_tors-1)
	real(8):: dfc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1,0:n_fix_tors-1),dfc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1,0:n_fix_tors-1),dfc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1,0:n_fix_tors-1) 
	real(8):: dfc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1,0:n_fix_tors-1),dfc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1,0:n_fix_tors-1), dfc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1,0:n_fix_tors-1)
	real(8):: dfc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1,0:n_fix_tors-1),dfc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1,0:n_fix_tors-1),dfc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1,0:n_fix_tors-1) 
	real(8):: dfc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1,0:n_fix_tors-1),dfc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1,0:n_fix_tors-1),dfc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1,0:n_fix_tors-1) 
	real(8):: dfc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1,0:n_fix_tors-1),dfc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1,0:n_fix_tors-1),dfc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1,0:n_fix_tors-1) 
	real(8):: dfc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1,0:n_fix_tors-1),dfc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1,0:n_fix_tors-1),dfc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1,0:n_fix_tors-1) 
	real(8):: dfc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1,0:n_fix_tors-1),dfc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1,0:n_fix_tors-1),dfc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1,0:n_fix_tors-1) 
	real(8):: dfc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1,0:n_fix_tors-1),dfc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1,0:n_fix_tors-1),dfc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1,0:n_fix_tors-1) 
	real(8)::fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1)
	real(8)::fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
	real(8)::dfc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1,0:n_fix_tors-1),dfc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1,0:n_fix_tors-1)
	real(8)::dfc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1,0:n_fix_tors-1),dfc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1,0:n_fix_tors-1)
	fc_quad_ea = 0.d0;fc_quad_oa = 0.d0; 
	fc_quad_es = 0.d0;fc_quad_os = 0.d0; 
	fc_qubic_es_gt5 = 0.d0;fc_qubic_es_lt5 = 0.d0;fc_qubic_es_lt2 = 0.d0;
	fc_qubic_ea_gt5 = 0.d0;fc_qubic_ea_lt5 = 0.d0;fc_qubic_ea_lt2 = 0.d0;
	fc_qubic_os_gt5 = 0.d0;fc_qubic_os_lt5 = 0.d0;fc_qubic_os_lt2 = 0.d0;
	fc_qubic_oa_gt5 = 0.d0;fc_qubic_oa_lt5 = 0.d0;fc_qubic_oa_lt2 = 0.d0;
	fc_quartic_es_gt5 = 0.d0;fc_quartic_es_lt5 = 0.d0;fc_quartic_es_lt2 = 0.d0;
	fc_quartic_ea_gt5 = 0.d0;fc_quartic_ea_lt5 = 0.d0;fc_quartic_ea_lt2 = 0.d0;
	fc_quartic_os_gt5 = 0.d0;fc_quartic_os_lt5 = 0.d0;fc_quartic_os_lt2 = 0.d0; 
	fc_quartic_oa_gt5 = 0.d0;fc_quartic_oa_lt5 = 0.d0;fc_quartic_oa_lt2 = 0.d0; 
	dfc_quad_es = 0.d0;dfc_quad_os = 0.d0; 
	dfc_quad_ea = 0.d0;dfc_quad_oa = 0.d0; 
	dfc_qubic_es_gt5 = 0.d0;dfc_qubic_es_lt5 = 0.d0;dfc_qubic_es_lt2 = 0.d0;
	dfc_qubic_ea_gt5 = 0.d0;dfc_qubic_ea_lt5 = 0.d0;dfc_qubic_ea_lt2 = 0.d0;
	dfc_qubic_os_gt5 = 0.d0;dfc_qubic_os_lt5 = 0.d0;dfc_qubic_os_lt2 = 0.d0;
	dfc_qubic_oa_gt5 = 0.d0;dfc_qubic_oa_lt5 = 0.d0;dfc_qubic_oa_lt2 = 0.d0;
	dfc_quartic_es_gt5 = 0.d0; dfc_quartic_es_lt5 = 0.d0;dfc_quartic_es_lt2 = 0.d0;
	dfc_quartic_ea_gt5 = 0.d0; dfc_quartic_ea_lt5 = 0.d0;dfc_quartic_ea_lt2 = 0.d0;
	dfc_quartic_os_gt5 = 0.d0;dfc_quartic_os_lt5 = 0.d0;dfc_quartic_os_lt2 = 0.d0; 
	dfc_quartic_oa_gt5 = 0.d0;dfc_quartic_oa_lt5 = 0.d0;dfc_quartic_oa_lt2 = 0.d0; 
	dfc_quartic_os_ltpt1 = 0.d0;dfc_quartic_oa_ltpt1=0.d0;
	dfc_quartic_es_ltpt1 = 0.d0;dfc_quartic_ea_ltpt1=0.d0;
	S_eq = 0.d0; sym_Seq = 0.d0; sym_S = 0.d0
	Vt_out = 0.d0;dVdx = 0.d0;sinth = 0.d0;costh=0.d0;dVdsym=0.d0;dVds=0.d0; 
	!Convert Bohr to Angstrom
	Xtmp =  X * bohrToAng 
	!Convert cartesian to internal
	call cart_to_internal(Xtmp,S)
	!Fixed torsions array
	fix_tors = S(tors_idx) * degreeToRadian
	!Get sin and cosine arrays
	call get_csarr(fix_tors,max_id,sinth,costh)
	!Get equlibrium symmetrized internal
	call get_se(sinth,costh,sym_Seq)
	!Convert eqilibrium symmetric internal to unsymmetrized internal
	sym_Seq(11) = sym_Seq(11) +360.0
	sym_Seq(18) = sym_Seq(18) +360.0
	S_eq  =  matmul(sym_mat_transpose,sym_Seq)
	do i = Nbonds+Nang,N_int-1
        !Check if S(i) is close to 180 degree
        !if (abs(abs(S(i))-180.d0)<90.d0) then
            !If S(i) and S_eq(i) have opposite signs
            if (S(i)*S_eq(i)<0.d0) then
               if (S_eq(i)>0.d0) then
                  S(i) = S(i) + 360.d0
               else
                  S(i) = S(i) - 360.d0
               endif
            endif
        !endif
    enddo
    s_diff = S-S_eq
	!print*,'2=',matmul(sym_mat,S);read(*,*)
    Sym_diff = matmul(sym_mat,s_diff)
	!Convert to scaled coordinates
	call scale_coord(Sym_diff)
	!compute fitted force constants(fc)
	call get_all_fc(sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1) 
	call get_all_dfc(sinth,costh,dfc_quad_es,dfc_quad_os,dfc_quad_ea,dfc_quad_oa,dfc_qubic_es_gt5,dfc_qubic_es_lt5,dfc_qubic_es_lt2,dfc_qubic_ea_gt5,dfc_qubic_ea_lt5,dfc_qubic_ea_lt2,dfc_qubic_os_gt5,dfc_qubic_os_lt5,dfc_qubic_os_lt2,dfc_qubic_oa_gt5,dfc_qubic_oa_lt5,dfc_qubic_oa_lt2,dfc_quartic_es_gt5,dfc_quartic_es_lt5,dfc_quartic_es_lt2,dfc_quartic_ea_gt5,dfc_quartic_ea_lt5,dfc_quartic_ea_lt2, dfc_quartic_os_gt5,dfc_quartic_os_lt5,dfc_quartic_os_lt2,dfc_quartic_oa_gt5,dfc_quartic_oa_lt5,dfc_quartic_oa_lt2) 
	!Compute derivative of fitted fc (dfc)
	
	!Compute Dvdsym(cm-1)
call getDvds(sym_diff,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2 ,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,dfc_quad_es,dfc_quad_os,dfc_quad_ea,dfc_quad_oa,dfc_qubic_es_gt5,dfc_qubic_es_lt5,dfc_qubic_es_lt2,dfc_qubic_ea_gt5,dfc_qubic_ea_lt5,dfc_qubic_ea_lt2,dfc_qubic_os_gt5,dfc_qubic_os_lt5,dfc_qubic_os_lt2,dfc_qubic_oa_gt5,dfc_qubic_oa_lt5,dfc_qubic_oa_lt2,dfc_quartic_es_gt5,dfc_quartic_es_lt5,dfc_quartic_es_lt2,dfc_quartic_es_ltpt1,dfc_quartic_ea_gt5,dfc_quartic_ea_lt5,dfc_quartic_ea_lt2,dfc_quartic_ea_ltpt1, dfc_quartic_os_gt5,dfc_quartic_os_lt5,dfc_quartic_os_ltpt1,dfc_quartic_os_lt2,dfc_quartic_oa_gt5,dfc_quartic_oa_lt5,dfc_quartic_oa_lt2,dfc_quartic_oa_ltpt1,dVdsym) 
	!Convert to cm^-1/degree or cm^-1/angstrom 
	do i = 0,N_int-1
		dVdsym(i) = dVdsym(i) * dim_scal_factor_inv(i)
		!Radian to degree
		if (i>= Nbonds) then
				dVdsym(i) = dVdsym(i) * degreeToRadian
		endif
	enddo
	!Compute dVds (cm^-1)
	dVds = matmul(sym_mat_transpose,dVdsym)
	!Get dsdx (degree/angstrom)
	call get_dsdx(Xtmp,dsdx)
	!get cart_grad
	call int_to_cart_grad(dsdx,dVds,dVdx)
	dVdx = dVdx * CminvToHartee * BohrToAng
	if (flag==1) then
		call compute_pot(sym_diff,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1 ,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,pot,Vt_out) 
	endif
endsubroutine 
!****************************************************************************************#
subroutine numerical_dVdx(X,dVdx_num)
implicit none 
real(8),intent(inout)::X(0:Ncarts-1)
real(8),intent(out)::dVdx_num(0:Ncarts-1)
integer::i,istep,flag=1
real(8)::Xtmp(0:Ncarts-1),pot_arr(0:4),pot=0.d0,Vt_out,dum(0:Ncarts-1)
Xtmp = 0.d0
pot = 0.d0
pot_arr = 0.d0;dum = 0.d0
dVdx_num = 0.d0
!print*,"Xinside=",X
X(:)= X(:) * BohrToAng
!loop over the cartesian coordinates
do i = 0,Ncarts-1
	Xtmp(:) = X(:)
	pot = 0.d0
	!displace five points and save the potential at each point
	do istep = 0,4
		Xtmp(i) = X(i)+ dx*dble(istep-2)
		call get_pot(Xtmp,pot,Vt_out)
		pot_arr(istep) = pot	
	enddo
	dVdx_num(i) = dx_inv * dot_product(five_pt_1st_deriv,pot_arr)  
enddo
endsubroutine
end module 
