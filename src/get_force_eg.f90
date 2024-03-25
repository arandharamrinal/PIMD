module eg_pot_force
use POTVARS
use allvars, only : cartmass, total_mass_inv
use constants
implicit none 
contains
!Change it 
!****************************************************************************************#
subroutine Read_ffext(fn,Ncarts,X,dVdx,E)
    implicit none
    integer::i,j
    integer,intent(in)::fn,Ncarts
    real(8),intent(out),dimension(0:Ncarts-1)::X,dVdx
    real(8),intent(out)::E
    X = 0.d0 ;dVdx =0.d0;E = 0.d0
    do i = 0,3
        read(fn,*)
    enddo
    !Read Cartesian coordinates (assuming 4th line)
    read(fn,*)(X(j),j=0,Ncarts-1)
    do i = 0,2
        read(fn,*)
    enddo
    read(fn,*)E
    !read(fn,*)
    !read(fn,*)(dVdX(j),j=0,Ncarts-1)
    !X = X * BohrToAng
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
	!use POTVARS, only: cartmass,Natoms,Ncarts
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

subroutine getVsav(S,fc_idx,fc,n_params,Vs)
	implicit none
	integer,intent(in):: n_params,fc_idx(0:n_params-1,0:3)
	real(8),intent(in):: S(0:N_int-1),fc(0:n_params-1) 
	real(8),intent(out):: Vs
	integer:: ii,i,j,k,l
    Vs = 0.d0
	do ii = 0,n_params-1
        i = fc_idx(ii,0);j = fc_idx(ii,1);k = fc_idx(ii,2);l = fc_idx(ii,3);
      	if (k==-1) then   
            Vs = Vs + fc(ii) * S(i) * S(j)
        elseif (l==-1) then
            Vs = Vs + fc(ii) * S(i) * S(j) * S(k)
        else
            Vs = Vs + fc(ii) * S(i) * S(j) * S(k) * S(l)
		endif
	enddo 
endsubroutine getVsav
!****************************************************************************************#
!****************************************************************************************#
subroutine  getVsav_qo(S,fc_idx,fc,n_params,Vsav)
    implicit none 
    integer,intent(in):: n_params,fc_idx(0:n_params-1,0:7)
    real(8),intent(in):: S(0:N_int-1),fc(0:n_params-1) 
    real(8),intent(out):: Vsav
    integer:: ii,i,j,k,l,nn,jj,kk,ll
    Vsav = 0.d0
    do nn = 0,n_params-1
        i = fc_idx(nn,0);j = fc_idx(nn,1); k = fc_idx(nn,2);l = fc_idx(nn,3);
        ii = fc_idx(nn,4);jj = fc_idx(nn,5); kk = fc_idx(nn,6);ll = fc_idx(nn,7);
        if (jj==-1) then 
            Vsav = Vsav + fc(nn) * S(i) * S(j) * S(k) * S(l) * S(ii)
        else if (kk==-1) then
            Vsav = Vsav + fc(nn) * S(i) * S(j) * S(k) * S(l) * S(ii) * S(jj)
        else if (ll==-1) then 
            Vsav = Vsav + fc(nn) * S(i) * S(j) * S(k) * S(l) * S(ii) * S(jj) * S(kk)
        else
            Vsav = Vsav + fc(nn) * S(i) * S(j) * S(k) * S(l) * S(ii) * S(jj) * S(kk) * S(ll)
        endif
    enddo
endsubroutine getVsav_qo
!****************************************************************************************#
!****************************************************************************************#
!compute force constants for a given value of phi1 and phi2
!****************************************************************************************#
subroutine get_all_vsav(S,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2 ,fc_quartic_ea_ltpt1,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,Vsav) 
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
		real(8),intent(in):: fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
		real(8),intent(in):: fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1) 
		real(8),intent(in):: fc_qo_es(0:n_fcs_qo_es-1),fc_qo_ea(0:n_fcs_qo_ea-1),fc_qo_os(0:n_fcs_qo_os-1),fc_qo_oa(0:n_fcs_qo_oa-1)
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
		real(8) :: Vsav_qo_es,Vsav_qo_ea,Vsav_qo_os,Vsav_qo_oa
		real(8)  :: Vdisp
		Vsav_quad_ea=0.d0;Vsav_quad_oa=0.d0;Vsav_quad_es=0.d0;Vsav_quad_os=0.d0
		Vsav_qubic_gt5_es=0.d0;Vsav_qubic_lt5_es=0.d0;Vsav_qubic_lt2_es=0.d0 
		Vsav_qubic_gt5_ea=0.d0;Vsav_qubic_lt5_ea=0.d0;Vsav_qubic_lt2_ea=0.d0
		Vsav_qubic_gt5_os=0.d0;Vsav_qubic_lt5_os=0.d0;Vsav_qubic_lt2_os=0.d0
		Vsav_qubic_gt5_oa=0.d0;Vsav_qubic_lt5_oa=0.d0;Vsav_qubic_lt2_oa=0.d0
		Vsav_quartic_gt5_es=0.d0;Vsav_quartic_lt5_es=0.d0;Vsav_quartic_lt2_es=0.d0 
		Vsav_quartic_gt5_ea=0.d0;Vsav_quartic_lt5_ea=0.d0;Vsav_quartic_lt2_ea=0.d0
		Vsav_quartic_gt5_os=0.d0;Vsav_quartic_lt5_os=0.d0;Vsav_quartic_lt2_os=0.d0
		Vsav_quartic_gt5_oa=0.d0;Vsav_quartic_lt5_oa=0.d0;Vsav_quartic_lt2_oa=0.d0
		Vsav_quartic_ltpt1_os=0.d0;Vsav_quartic_ltpt1_oa=0.d0 
		Vsav_quartic_ltpt1_es=0.d0;Vsav_quartic_ltpt1_ea=0.d0 
		Vsav_qo_es=0.d0;Vsav_qo_ea=0.d0;Vsav_qo_os=0.d0;Vsav_qo_oa=0.d0 
        Vdisp = 0.d0
        !!quadratic 
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

        call getVsav_qo(S,fc_idx_qo_es,fc_qo_es,n_fcs_qo_es,Vsav_qo_es)
        call getVsav_qo(S,fc_idx_qo_ea,fc_qo_ea,n_fcs_qo_ea,Vsav_qo_ea)
        call getVsav_qo(S,fc_idx_qo_os,fc_qo_os,n_fcs_qo_os,Vsav_qo_os)
        call getVsav_qo(S,fc_idx_qo_oa,fc_qo_oa,n_fcs_qo_oa,Vsav_qo_oa)

        !Sum up all 
		Vdisp = Vdisp + Vsav_quad_ea + Vsav_quad_oa + Vsav_quad_es + Vsav_quad_os
		Vdisp = Vdisp + Vsav_qubic_gt5_es + Vsav_qubic_lt5_es + Vsav_qubic_lt2_es
		Vdisp = Vdisp + Vsav_qubic_gt5_ea + Vsav_qubic_lt5_ea + Vsav_qubic_lt2_ea
		Vdisp = Vdisp + Vsav_qubic_gt5_os + Vsav_qubic_lt5_os + Vsav_qubic_lt2_os
		Vdisp = Vdisp + Vsav_qubic_gt5_oa + Vsav_qubic_lt5_oa + Vsav_qubic_lt2_oa
		Vdisp = Vdisp + Vsav_quartic_gt5_es + Vsav_quartic_lt5_es + Vsav_quartic_lt2_es
		Vdisp = Vdisp + Vsav_quartic_gt5_ea + Vsav_quartic_lt5_ea + Vsav_quartic_lt2_ea
		Vdisp = Vdisp + Vsav_quartic_gt5_os + Vsav_quartic_lt5_os + Vsav_quartic_lt2_os
		Vdisp = Vdisp + Vsav_quartic_gt5_oa + Vsav_quartic_lt5_oa + Vsav_quartic_lt2_oa
		Vdisp = Vdisp + Vsav_quartic_ltpt1_oa + Vsav_quartic_ltpt1_os
		Vdisp = Vdisp + Vsav_quartic_ltpt1_ea + Vsav_quartic_ltpt1_es
		Vdisp = Vdisp + Vsav_qo_es + Vsav_qo_ea + Vsav_qo_os + Vsav_qo_oa
		Vsav = Vdisp 
endsubroutine  

!****************************************************************************************#
subroutine  computepot(S,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1 ,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,DMES_eqm,egpot,VrsEg) 
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
	real(8),intent(in):: fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
	real(8),intent(in):: fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1) 
	real(8),intent(in):: fc_qo_es(0:n_fcs_qo_es-1),fc_qo_ea(0:n_fcs_qo_ea-1),fc_qo_os(0:n_fcs_qo_os-1),fc_qo_oa(0:n_fcs_qo_oa-1)
	real(8),intent(in):: DMES_eqm(0:fn_order_bonds_es-1)
	real(8),intent(out)::egpot,VrsEg
	real(8)::Vt=0.d0,Vsav=0.d0
    egpot = 0.d0
	VrsEg = dot_product(DMES_eqm,Vt_coeff)
	call get_all_vsav(S,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2 ,fc_quartic_ea_ltpt1,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,Vsav) 
    egpot  = VrsEg + Vsav 
endsubroutine 
!****************************************************************************************#
subroutine  dvdsQuad(fc_idx,fc,dS_dtors,dfc,S,n_params,dVds)
	implicit none 
	integer,intent(in)::n_params
	real(8),intent(in)::fc(0:n_params-1),dfc(0:n_params-1,0:n_fix_tors-1),dS_dtors(0:N_int-1,0:n_fix_tors-1),S(0:N_int-1)
	integer,intent(in) :: fc_idx(0:n_params-1,0:3)
	real(8),intent(inout):: dVds(0:N_int-1)
	integer:: ii,i,itors,j,k,l
	do ii =0,n_params-1
        i    = fc_idx(ii,0);j = fc_idx(ii,1);k = fc_idx(ii,2);l = fc_idx(ii,3);
		dVds(i) = dVds(i) + fc(ii) * S(j)
		dVds(j) = dVds(j) + fc(ii) * S(i)
		do itors = 0,n_fix_tors-1
			dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + dfc(ii,itors) * S(i) * S(j) + fc(ii) * dS_dtors(i,itors) * S(j) + fc(ii) * dS_dtors(j,itors) * S(i)
		enddo
	enddo
endsubroutine 
!****************************************************************************************#
subroutine  dvdsQubic(fc_idx,fc,dS_dtors,dfc,S,n_params,dVds)
	implicit none 
	integer,intent(in)::n_params
	real(8),intent(in)::fc(0:n_params-1),dfc(0:n_params-1,0:n_fix_tors-1),dS_dtors(0:N_int-1,0:n_fix_tors-1),S(0:N_int-1)
	integer,intent(in) :: fc_idx(0:n_params-1,0:3)
	real(8),intent(inout):: dVds(0:N_int-1)
	integer:: ii,i,itors,j,k,l
	do ii =0,n_params-1
        i    = fc_idx(ii,0);j = fc_idx(ii,1);k = fc_idx(ii,2);l = fc_idx(ii,3);
		dVds(i) = dVds(i) + fc(ii) * S(j) * S(k)
		dVds(j) = dVds(j) + fc(ii) * S(i) * S(k)
		dVds(k) = dVds(k) + fc(ii) * S(i) * S(j)
		do itors = 0,n_fix_tors-1
        	dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + dfc(ii,itors) * S(i) * S(j) * S(k) + fc(ii) * dS_dtors(i,itors) * S(j) * S(k)  + fc(ii) * dS_dtors(j,itors) * S(i) * S(k) 
        	dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(ii) * dS_dtors(k,itors) * S(i) * S(j)
		enddo
	enddo
endsubroutine 
!****************************************************************************************#
subroutine  dvdsQuartic(fc_idx,fc,dS_dtors,dfc,S,n_params,dVds)
	implicit none 
	integer,intent(in)::n_params
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
		do itors = 0,n_fix_tors-1
        	dVds(tors_idx(itors)) =dVds(tors_idx(itors)) + dfc(ii,itors) * S(i) * S(j) * S(k) * S(l)
        	dVds(tors_idx(itors)) =dVds(tors_idx(itors)) + fc(ii) * S(i) * dS_dtors(j,itors) * S(k) * S(l) + fc(ii) * S(j) * dS_dtors(i,itors) * S(k) * S(l) 
        	dVds(tors_idx(itors)) =dVds(tors_idx(itors)) + fc(ii) * dS_dtors(k,itors) * S(i) * S(j) * S(l) + fc(ii) * dS_dtors(l,itors) * S(i) * S(j) * S(k)
		enddo
	enddo
endsubroutine 
!****************************************************************************************#
!****************************************************************************************#
subroutine  dvdsQO(fc_idx,fc,dS_dtors,dfc,S,n_params,dVds)
    implicit none
    integer,intent(in)::n_params
    real(8),intent(in)::fc(0:n_params-1),dfc(0:n_params-1,0:n_fix_tors-1),dS_dtors(0:N_int-1,0:n_fix_tors-1),S(0:N_int-1)
    integer,intent(in) :: fc_idx(0:n_params-1,0:7)
    real(8),intent(inout):: dVds(0:N_int-1)
    integer:: ii,i,itors,j,k,l,jj,kk,ll,nn
    do nn =0,n_params-1
         i = fc_idx(nn,0);j = fc_idx(nn,1);k = fc_idx(nn,2);l = fc_idx(nn,3);
         ii = fc_idx(nn,4);jj = fc_idx(nn,5);kk = fc_idx(nn,6);ll = fc_idx(nn,7);
         if (jj==-1) then
                 dVds(i) =  dVds(i) + fc(nn) * S(j) * S(k) * S(l) * S(ii)
                 dVds(j) =  dVds(j) + fc(nn) * S(i) * S(k) * S(l) * S(ii)
                 dVds(k) =  dVds(k) + fc(nn) * S(i) * S(j) * S(l) * S(ii)
                 dVds(l) =  dVds(l) + fc(nn) * S(i) * S(j) * S(k) * S(ii)
                dVds(ii) = dVds(ii) + fc(nn) * S(j) * S(k) * S(l) * S(i)
                 do itors = 0,n_fix_tors-1
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + dfc(nn,itors) * S(i) * S(j) * S(k) * S(l) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(j,itors) * S(i) * S(k) * S(l) *  S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(i,itors) * S(j) * S(k) * S(l) *  S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(k,itors) * S(i) * S(j) * S(l) *  S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(l,itors) * S(i) * S(j) * S(k) *  S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(ii,itors) * S(j) * S(k) * S(l) *  S(i)
                enddo
        else if (kk==-1) then
                 dVds(i) =  dVds(i) + fc(nn) * S(j) * S(k) * S(l) * S(jj) * S(ii)
                 dVds(j) =  dVds(j) + fc(nn) * S(i) * S(k) * S(l) * S(jj) * S(ii)
                 dVds(k) =  dVds(k) + fc(nn) * S(i) * S(j) * S(l) * S(jj) * S(ii)
                 dVds(l) =  dVds(l) + fc(nn) * S(i) * S(j) * S(k) * S(jj) * S(ii)
                dVds(ii) = dVds(ii) + fc(nn) * S(j) * S(k) * S(l) * S(i)  * S(jj)
                dVds(jj) = dVds(jj) + fc(nn) * S(i) * S(k) * S(l) * S(j)  * S(ii)
                 do itors = 0,n_fix_tors-1
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + dfc(nn,itors) * S(i) * S(j) * S(k) * S(l) * S(jj) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(j,itors) * S(i) * S(k) * S(l)  * S(jj) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(i,itors) * S(j) * S(k) * S(l)  * S(jj) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(k,itors) * S(i) * S(j) * S(l)  * S(jj) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(l,itors) * S(i) * S(j) * S(k)  * S(jj) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(jj,itors) * S(i) * S(k) * S(l)  * S(j)  * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(ii,itors) * S(j) * S(k) * S(l)  * S(jj) * S(i)
                enddo
        else if (ll==-1) then
                 dVds(i) =  dVds(i) + fc(nn) * S(j) * S(k) * S(l) * S(jj) * S(kk) * S(ii)
                 dVds(j) =  dVds(j) + fc(nn) * S(i) * S(k) * S(l) * S(jj) * S(kk) * S(ii)
                 dVds(k) =  dVds(k) + fc(nn) * S(i) * S(j) * S(l) * S(jj) * S(kk) * S(ii)
                 dVds(l) =  dVds(l) + fc(nn) * S(i) * S(j) * S(k) * S(jj) * S(kk) * S(ii)
                dVds(ii) = dVds(ii) + fc(nn) * S(j) * S(k) * S(l) * S(i)  * S(jj) * S(kk)
                dVds(jj) = dVds(jj) + fc(nn) * S(i) * S(k) * S(l) * S(j)  * S(kk) * S(ii)
                dVds(kk) = dVds(kk) + fc(nn) * S(i) * S(j) * S(l) * S(k)  * S(jj) * S(ii)
                 do itors = 0,n_fix_tors-1
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + dfc(nn,itors) * S(i) * S(j) * S(k) * S(l) * S(jj) * S(kk) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(j,itors) * S(i) * S(k) * S(l)  * S(jj) * S(kk) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(i,itors) * S(j) * S(k) * S(l)  * S(jj) * S(kk) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(k,itors) * S(i) * S(j) * S(l)  * S(jj) * S(kk) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(l,itors) * S(i) * S(j) * S(k)  * S(jj) * S(kk) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(jj,itors) * S(i) * S(k) * S(l)  * S(j)  * S(kk) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(ii,itors) * S(j) * S(k) * S(l)  * S(jj) * S(kk) * S(i)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(kk,itors) * S(i) * S(j) * S(l)  * S(jj) * S(k)  * S(ii)
                enddo
        else
                 dVds(i) =  dVds(i) + fc(nn) * S(j) * S(k) * S(l) * S(jj) * S(kk) * S(ll) * S(ii)
                 dVds(j) =  dVds(j) + fc(nn) * S(i) * S(k) * S(l) * S(jj) * S(kk) * S(ll) * S(ii)
                 dVds(k) =  dVds(k) + fc(nn) * S(i) * S(j) * S(l) * S(jj) * S(kk) * S(ll) * S(ii)
                 dVds(l) =  dVds(l) + fc(nn) * S(i) * S(j) * S(k) * S(jj) * S(kk) * S(ll) * S(ii)
                dVds(ii) = dVds(ii) + fc(nn) * S(j) * S(k) * S(l) * S(i)  * S(jj) * S(kk) * S(ll)
                dVds(jj) = dVds(jj) + fc(nn) * S(i) * S(k) * S(l) * S(j)  * S(kk) * S(ll) * S(ii)
                dVds(kk) = dVds(kk) + fc(nn) * S(i) * S(j) * S(l) * S(k)  * S(jj) * S(ll) * S(ii)
                dVds(ll) = dVds(ll) + fc(nn) * S(i) * S(j) * S(k) * S(l)  * S(jj) * S(kk) * S(ii)
                 do itors = 0,n_fix_tors-1
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + dfc(nn,itors) * S(i) * S(j) * S(k) * S(l) * S(jj) * S(kk) * S(ll) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(j,itors) * S(i) * S(k) * S(l)  * S(jj) * S(kk) * S(ll) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(i,itors) * S(j) * S(k) * S(l)  * S(jj) * S(kk) * S(ll) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(k,itors) * S(i) * S(j) * S(l)  * S(jj) * S(kk) * S(ll) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) *  dS_dtors(l,itors) * S(i) * S(j) * S(k)  * S(jj) * S(kk) * S(ll) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(jj,itors) * S(i) * S(k) * S(l)  * S(j)  * S(kk) * S(ll) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(ii,itors) * S(j) * S(k) * S(l)  * S(jj) * S(kk) * S(ll) * S(i)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(kk,itors) * S(i) * S(j) * S(l)  * S(jj) * S(k)  * S(ll) * S(ii)
                    dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(nn) * dS_dtors(ll,itors) * S(i) * S(j) * S(k)  * S(jj) * S(kk) * S(l)  * S(ii)
                enddo

        endif
    enddo
endsubroutine
!****************************************************************************************#
!compute internal gradients.
subroutine getDvds(S,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2 ,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,dfc_quad_es,dfc_quad_ea,dfc_quad_os,dfc_quad_oa,dfc_qubic_es_gt5,dfc_qubic_es_lt5,dfc_qubic_es_lt2,dfc_qubic_ea_gt5,dfc_qubic_ea_lt5,dfc_qubic_ea_lt2,dfc_qubic_os_gt5,dfc_qubic_os_lt5,dfc_qubic_os_lt2,dfc_qubic_oa_gt5,dfc_qubic_oa_lt5,dfc_qubic_oa_lt2,dfc_quartic_es_gt5,dfc_quartic_es_lt5,dfc_quartic_es_lt2,dfc_quartic_es_ltpt1,dfc_quartic_ea_gt5,dfc_quartic_ea_lt5,dfc_quartic_ea_lt2,dfc_quartic_ea_ltpt1, dfc_quartic_os_gt5,dfc_quartic_os_lt5,dfc_quartic_os_ltpt1,dfc_quartic_os_lt2,dfc_quartic_oa_gt5,dfc_quartic_oa_lt5,dfc_quartic_oa_lt2,dfc_quartic_oa_ltpt1,dfc_qo_es,dfc_qo_ea,dfc_qo_os,dfc_qo_oa,DDMES_eqm,DDMEA_eqm,DDMOS_eqm,DDMOA_eqm,dVds) 
	implicit none 
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
	real(8),intent(in)::fc_qo_es(0:n_fcs_qo_es-1),fc_qo_ea(0:n_fcs_qo_ea-1)
	real(8),intent(in)::fc_qo_os(0:n_fcs_qo_os-1),fc_qo_oa(0:n_fcs_qo_oa-1) 
	 
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
	real(8),intent(in)::dfc_qo_os(0:n_fcs_qo_os-1,0:n_fix_tors-1),dfc_qo_oa(0:n_fcs_qo_oa-1,0:n_fix_tors-1)
	real(8),intent(in)::dfc_qo_es(0:n_fcs_qo_es-1,0:n_fix_tors-1),dfc_qo_ea(0:n_fcs_qo_ea-1,0:n_fix_tors-1)
	real(8),intent(in):: DDMES_eqm(0:n_fix_tors-1,0:fn_order_bonds_es-1),DDMEA_eqm(0:n_fix_tors-1,0:fn_order_bonds_ea-1),DDMOS_eqm(0:n_fix_tors-1,0:fn_order_bonds_os-1),DDMOA_eqm(0:n_fix_tors-1,0:fn_order_bonds_oa-1)
    !""" i'th element of dVds contains the derivatives of egpotential w.r.t to the i'th internal"""
    !Compute the derivatives of the i'th internal w.r.t the fixed coordinates
	dS_dtors = 0.d0
	dVds = 0.d0
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    do i = 0,size(os_fn_idx)-1
        idxs = os_fn_idx(i)
        if (idxs<Nbonds) then
            dS_dtors(idxs,:) = -matmul(DDMOS_eqm(:,:),fitted_bonds_coeff_os(:,b_ctr))
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            dS_dtors(idxs,:) = -matmul(DDMOS_eqm(:,:),fitted_angs_coeff_os(:,a_ctr)) 
            a_ctr = a_ctr + 1 
        else
            dS_dtors(idxs,:) = -matmul(DDMOS_eqm(:,:),fitted_dihs_coeff_os(:,d_ctr))
            d_ctr = d_ctr + 1 
		endif
	enddo
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    do i = 0,size(oa_fn_idx)-1
        idxs = oa_fn_idx(i)
        if (idxs<Nbonds) then
            dS_dtors(idxs,:) = -matmul(DDMOA_eqm(:,:),fitted_bonds_coeff_oa(:,b_ctr))
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            dS_dtors(idxs,:) = -matmul(DDMOA_eqm(:,:),fitted_angs_coeff_oa(:,a_ctr))
            a_ctr = a_ctr + 1 
        else
            dS_dtors(idxs,:) = -matmul(DDMOA_eqm(:,:),fitted_dihs_coeff_oa(:,d_ctr))
            d_ctr = d_ctr + 1 
		endif
	enddo
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    do i = 0,size(es_fn_idx)-1
        idxs = es_fn_idx(i)
        if (idxs<Nbonds) then
            dS_dtors(idxs,:) = -matmul(DDMES_eqm(:,:),fitted_bonds_coeff_es(:,b_ctr))
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            dS_dtors(idxs,:) = -matmul(DDMES_eqm(:,:),fitted_angs_coeff_es(:,a_ctr))
            a_ctr = a_ctr + 1 
        else
            dS_dtors(idxs,:) = -matmul(DDMES_eqm(:,:),fitted_dihs_coeff_es(:,d_ctr))
            d_ctr = d_ctr + 1 
		endif
	enddo
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    do i = 0,size(ea_fn_idx)-1
        idxs = ea_fn_idx(i)
        if (idxs<Nbonds) then
            dS_dtors(idxs,:) =  -matmul(DDMEA_eqm(:,:),fitted_bonds_coeff_ea(:,b_ctr))
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            dS_dtors(idxs,:) =  -matmul(DDMEA_eqm(:,:),fitted_angs_coeff_ea(:,a_ctr))
            a_ctr = a_ctr + 1 
        else
            dS_dtors(idxs,:) =  -matmul(DDMEA_eqm(:,:),fitted_dihs_coeff_ea(:,d_ctr))
            d_ctr = d_ctr + 1 
		endif
	enddo
	do i = 0,N_int-1 
		dS_dtors(i,:) = dS_dtors(i,:) * dim_scal_factor_inv(i)
		if (i>= Nbonds) then
				dS_dtors(i,:) = dS_dtors(i,:) * degreeToRadian
		endif
	enddo
    !Compute the derivative of the egpotential along the minimum path(zero'th order term in the expansion) w.r.t the fixed coordinates
	call  getDVt(sinth,costh,DDMES_eqm,dVt)
    dVds(tors_idx) =  dVds(tors_idx) +  dVt(:)
    !Compute derivative of the 2nd order, 3rd order and fourth order terms w.r.t to all internal coordinates
    !!For quadratic even
	!ES
	call dvdsQuad(fc_idx_quad_es,fc_quad_es,dS_dtors,dfc_quad_es,S,n_fcs_quad_es,dVds)
	!EA
	call dvdsQuad(fc_idx_quad_ea,fc_quad_ea,dS_dtors,dfc_quad_ea,S,n_fcs_quad_ea,dVds)
	!OS
    call dvdsQuad(fc_idx_quad_os,fc_quad_os,dS_dtors,dfc_quad_os,S,n_fcs_quad_os,dVds)
	!OA
    call dvdsQuad(fc_idx_quad_oa,fc_quad_oa,dS_dtors,dfc_quad_oa,S,n_fcs_quad_oa,dVds)

    !Cubic:
	!ES
    call dvdsQubic(fc_idx_qubic_es_gt5,fc_qubic_es_gt5,dS_dtors,dfc_qubic_es_gt5,S,n_fcs_qubic_es_gt5,dVds)
    call dvdsQubic(fc_idx_qubic_es_lt5,fc_qubic_es_lt5,dS_dtors,dfc_qubic_es_lt5,S,n_fcs_qubic_es_lt5,dVds)
    call dvdsQubic(fc_idx_qubic_es_lt2,fc_qubic_es_lt2,dS_dtors,dfc_qubic_es_lt2,S,n_fcs_qubic_es_lt2,dVds)

	!EA
    call dvdsQubic(fc_idx_qubic_ea_gt5,fc_qubic_ea_gt5,dS_dtors,dfc_qubic_ea_gt5,S,n_fcs_qubic_ea_gt5,dVds)
    call dvdsQubic(fc_idx_qubic_ea_lt5,fc_qubic_ea_lt5,dS_dtors,dfc_qubic_ea_lt5,S,n_fcs_qubic_ea_lt5,dVds)
    call dvdsQubic(fc_idx_qubic_ea_lt2,fc_qubic_ea_lt2,dS_dtors,dfc_qubic_ea_lt2,S,n_fcs_qubic_ea_lt2,dVds)

	!OS
    call dvdsQubic(fc_idx_qubic_os_gt5,fc_qubic_os_gt5,dS_dtors,dfc_qubic_os_gt5,S,n_fcs_qubic_os_gt5,dVds)
    call dvdsQubic(fc_idx_qubic_os_lt5,fc_qubic_os_lt5,dS_dtors,dfc_qubic_os_lt5,S,n_fcs_qubic_os_lt5,dVds)
    call dvdsQubic(fc_idx_qubic_os_lt2,fc_qubic_os_lt2,dS_dtors,dfc_qubic_os_lt2,S,n_fcs_qubic_os_lt2,dVds)

	!EA
    call dvdsQubic(fc_idx_qubic_oa_gt5,fc_qubic_oa_gt5,dS_dtors,dfc_qubic_oa_gt5,S,n_fcs_qubic_oa_gt5,dVds)
    call dvdsQubic(fc_idx_qubic_oa_lt5,fc_qubic_oa_lt5,dS_dtors,dfc_qubic_oa_lt5,S,n_fcs_qubic_oa_lt5,dVds)
    call dvdsQubic(fc_idx_qubic_oa_lt2,fc_qubic_oa_lt2,dS_dtors,dfc_qubic_oa_lt2,S,n_fcs_qubic_oa_lt2,dVds)
    !Quartic:    
	!ES
    call dvdsQuartic(fc_idx_quartic_es_gt5,fc_quartic_es_gt5,dS_dtors,dfc_quartic_es_gt5,S,n_fcs_quartic_es_gt5,dVds)
    call dvdsQuartic(fc_idx_quartic_es_lt5,fc_quartic_es_lt5,dS_dtors,dfc_quartic_es_lt5,S,n_fcs_quartic_es_lt5,dVds)
    call dvdsQuartic(fc_idx_quartic_es_lt2,fc_quartic_es_lt2,dS_dtors,dfc_quartic_es_lt2,S,n_fcs_quartic_es_lt2,dVds)
    call dvdsQuartic(fc_idx_quartic_es_ltpt1,fc_quartic_es_ltpt1,dS_dtors,dfc_quartic_es_ltpt1,S,n_fcs_quartic_es_ltpt1,dVds)

	!EA
    call dvdsQuartic(fc_idx_quartic_ea_gt5,fc_quartic_ea_gt5,dS_dtors,dfc_quartic_ea_gt5,S,n_fcs_quartic_ea_gt5,dVds)
    call dvdsQuartic(fc_idx_quartic_ea_lt5,fc_quartic_ea_lt5,dS_dtors,dfc_quartic_ea_lt5,S,n_fcs_quartic_ea_lt5,dVds)
    call dvdsQuartic(fc_idx_quartic_ea_lt2,fc_quartic_ea_lt2,dS_dtors,dfc_quartic_ea_lt2,S,n_fcs_quartic_ea_lt2,dVds)
    call dvdsQuartic(fc_idx_quartic_ea_ltpt1,fc_quartic_ea_ltpt1,dS_dtors,dfc_quartic_ea_ltpt1,S,n_fcs_quartic_ea_ltpt1,dVds)

	!OS
    call dvdsQuartic(fc_idx_quartic_os_gt5,fc_quartic_os_gt5,dS_dtors,dfc_quartic_os_gt5,S,n_fcs_quartic_os_gt5,dVds)
    call dvdsQuartic(fc_idx_quartic_os_lt5,fc_quartic_os_lt5,dS_dtors,dfc_quartic_os_lt5,S,n_fcs_quartic_os_lt5,dVds)
    call dvdsQuartic(fc_idx_quartic_os_lt2,fc_quartic_os_lt2,dS_dtors,dfc_quartic_os_lt2,S,n_fcs_quartic_os_lt2,dVds)
    call dvdsQuartic(fc_idx_quartic_os_ltpt1,fc_quartic_os_ltpt1,dS_dtors,dfc_quartic_os_ltpt1,S,n_fcs_quartic_os_ltpt1,dVds)

	!OA
    call dvdsQuartic(fc_idx_quartic_oa_gt5,fc_quartic_oa_gt5,dS_dtors,dfc_quartic_oa_gt5,S,n_fcs_quartic_oa_gt5,dVds)
    call dvdsQuartic(fc_idx_quartic_oa_lt5,fc_quartic_oa_lt5,dS_dtors,dfc_quartic_oa_lt5,S,n_fcs_quartic_oa_lt5,dVds)
    call dvdsQuartic(fc_idx_quartic_oa_lt2,fc_quartic_oa_lt2,dS_dtors,dfc_quartic_oa_lt2,S,n_fcs_quartic_oa_lt2,dVds)
    call dvdsQuartic(fc_idx_quartic_oa_ltpt1,fc_quartic_oa_ltpt1,dS_dtors,dfc_quartic_oa_ltpt1,S,n_fcs_quartic_oa_ltpt1,dVds)

    call dvdsQO(fc_idx_qo_es,fc_qo_es,dS_dtors,dfc_qo_es,S,n_fcs_qo_es,dVds)
    call dvdsQO(fc_idx_qo_ea,fc_qo_ea,dS_dtors,dfc_qo_ea,S,n_fcs_qo_ea,dVds)
    call dvdsQO(fc_idx_qo_os,fc_qo_os,dS_dtors,dfc_qo_os,S,n_fcs_qo_os,dVds)
    call dvdsQO(fc_idx_qo_oa,fc_qo_oa,dS_dtors,dfc_qo_oa,S,n_fcs_qo_oa,dVds)

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
subroutine  getVt(sinth,costh,DMES_eqm,Vt)
	implicit none 
	real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(out)::Vt
	real(8),intent(in)::DMES_eqm(0:fn_order_bonds_es-1)
    Vt = 0.d0
	Vt = dot_product(DMES_eqm,Vt_coeff)
endsubroutine 
!****************************************************************************************#
subroutine  getDVt(sinth,costh,DDMES_eqm,dVt)
	implicit none
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(in)::DDMES_eqm(0:n_fix_tors-1,0:fn_order_bonds_es-1)
    real(8),intent(out)::dVt(0:n_fix_tors-1)
	integer:: jj,i,j
    dVt = 0.d0;
	dVt(:) = matmul(DDMES_eqm(:,:),Vt_coeff(:))
endsubroutine 
!Compute Vsav
!****************************************************************************************#
!given q6 and order of expansion compute the force constants:
subroutine get_fc(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,DesMat,fc)
	implicit none 
    integer::ii
	integer,intent(in)::n_params ,fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),DesMat(0:fn_order-1),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::fc(0:n_params-1) 
    fc = 0.d0
    do ii=0,n_params-1
        fc(ii) = dot_product(DesMat,fitted_fc_coeff(:,ii))
    enddo 
endsubroutine 
!***************************************************************************************#
!given  function order compute the force constants:
subroutine get_dfc(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,DDesMat,dfc)
	implicit none 
    integer::ii
	integer,intent(in)::n_params ,fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),DDesMat(0:n_fix_tors-1,0:fn_order-1),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::dfc(0:n_params-1,0:n_fix_tors-1)
	dfc = 0.d0;	
    do ii = 0,n_params-1
	    dfc(ii,0) = dot_product(DDesMat(0,:),fitted_fc_coeff(:,ii))
	    dfc(ii,1) = dot_product(DDesMat(1,:),fitted_fc_coeff(:,ii))
	    dfc(ii,2) = dot_product(DDesMat(2,:),fitted_fc_coeff(:,ii))
	enddo
endsubroutine 
!*******************************************************************************************************************#
subroutine get_se(sinth,costh,DMES_eqm,DMEA_eqm,DMOS_eqm,DMOA_eqm,S_e)
	implicit none
	real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(in)::DMES_eqm(0:fn_order_bonds_es-1),DMEA_eqm(0:fn_order_bonds_ea-1),DMOS_eqm(0:fn_order_bonds_os-1),DMOA_eqm(0:fn_order_bonds_oa-1)
	real(8),intent(out):: S_e(0:N_int-1)
	integer:: b_ctr,a_ctr,i,idxs
	real(8) :: dum = 0.d0
    S_e = 0.d0
    b_ctr =0;a_ctr = 0
    do i = 0,size(es_fn_idx)-1
        idxs = es_fn_idx(i)
        if (idxs<Nbonds) then
			S_e(idxs) = dot_product(DMES_eqm,fitted_bonds_coeff_es(:,i))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
			S_e(idxs) = dot_product(DMES_eqm,fitted_angs_coeff_es(:,i-b_ctr))
            a_ctr = a_ctr + 1
        else
			S_e(idxs) = dot_product(DMES_eqm,fitted_dihs_coeff_es(:,i-b_ctr-a_ctr))
		endif
	enddo
    b_ctr =0;a_ctr = 0
    do i = 0,size(ea_fn_idx)-1
        idxs = ea_fn_idx(i)
        if (idxs<Nbonds) then
			S_e(idxs) = dot_product(DMEA_eqm,fitted_bonds_coeff_ea(:,i))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
			S_e(idxs) = dot_product(DMEA_eqm,fitted_angs_coeff_ea(:,i-b_ctr))
            a_ctr = a_ctr + 1
        else
			S_e(idxs) = dot_product(DMEA_eqm,fitted_dihs_coeff_ea(:,i-b_ctr-a_ctr))
		endif
	enddo
    b_ctr =0;a_ctr = 0
    do i = 0,size(os_fn_idx)-1
        idxs = os_fn_idx(i)
        if (idxs<Nbonds) then
			S_e(idxs) = dot_product(DMOS_eqm,fitted_bonds_coeff_os(:,i))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
			S_e(idxs) = dot_product(DMOS_eqm,fitted_angs_coeff_os(:,i-b_ctr))
            a_ctr = a_ctr + 1
        else
			S_e(idxs) = dot_product(DMOS_eqm,fitted_dihs_coeff_os(:,i-b_ctr-a_ctr))
        endif
    enddo
    b_ctr =0;a_ctr = 0
    do i = 0,size(oa_fn_idx)-1
        idxs = oa_fn_idx(i)
        if (idxs<Nbonds) then
			S_e(idxs) = dot_product(DMOA_eqm,fitted_bonds_coeff_oa(:,i))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
			S_e(idxs) = dot_product(DMOA_eqm,fitted_angs_coeff_oa(:,i-b_ctr))
            a_ctr = a_ctr + 1
        else
			S_e(idxs) = dot_product(DMOA_eqm,fitted_dihs_coeff_oa(:,i-b_ctr-a_ctr))
        endif
    enddo
endsubroutine 
!****************************************************************************************************************************************************#
subroutine  get_all_dfc(sinth,costh,dfc_quad_es,dfc_quad_os,dfc_quad_ea,dfc_quad_oa,dfc_qubic_es_gt5,dfc_qubic_es_lt5,dfc_qubic_es_lt2,dfc_qubic_ea_gt5,dfc_qubic_ea_lt5,dfc_qubic_ea_lt2,dfc_qubic_os_gt5,dfc_qubic_os_lt5,dfc_qubic_os_lt2,dfc_qubic_oa_gt5,dfc_qubic_oa_lt5,dfc_qubic_oa_lt2,dfc_quartic_es_gt5,dfc_quartic_es_lt5,dfc_quartic_es_lt2,dfc_quartic_ea_gt5,dfc_quartic_ea_lt5,dfc_quartic_ea_lt2,dfc_quartic_os_gt5,dfc_quartic_os_lt5,dfc_quartic_os_lt2,dfc_quartic_oa_gt5,dfc_quartic_oa_lt5,dfc_quartic_oa_lt2,dfc_qo_es,dfc_qo_ea,dfc_qo_os,dfc_qo_oa,DDMES_quad,DDMES_gt5,DDMES_qo,DDMES_lt2,DDMEA_quad,DDMEA_gt5,DDMEA_qo,DDMEA_lt2,DDMOS_quad,DDMOS_gt5,DDMOS_qo,DDMOS_lt2,DDMOA_quad,DDMOA_gt5,DDMOA_qo,DDMOA_lt2) 
	implicit none
	real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(in):: DDMES_quad(0:n_fix_tors-1,0:fn_order_quad_es-1),DDMES_gt5(0:n_fix_tors-1,0:fn_order_qubic_es_gt5-1),DDMES_qo(0:n_fix_tors-1,0:fn_order_qo_es-1),DDMES_lt2(0:n_fix_tors-1,0:fn_order_qubic_es_lt2-1)
	real(8),intent(in):: DDMEA_quad(0:n_fix_tors-1,0:fn_order_quad_ea-1),DDMEA_gt5(0:n_fix_tors-1,0:fn_order_qubic_ea_gt5-1),DDMEA_qo(0:n_fix_tors-1,0:fn_order_qo_ea-1),DDMEA_lt2(0:n_fix_tors-1,0:fn_order_qubic_ea_lt2-1)
	real(8),intent(in):: DDMOS_quad(0:n_fix_tors-1,0:fn_order_quad_os-1),DDMOS_gt5(0:n_fix_tors-1,0:fn_order_qubic_os_gt5-1),DDMOS_qo(0:n_fix_tors-1,0:fn_order_qo_os-1),DDMOS_lt2(0:n_fix_tors-1,0:fn_order_qubic_os_lt2-1)
	real(8),intent(in):: DDMOA_quad(0:n_fix_tors-1,0:fn_order_quad_oa-1),DDMOA_gt5(0:n_fix_tors-1,0:fn_order_qubic_oa_gt5-1),DDMOA_qo(0:n_fix_tors-1,0:fn_order_qo_oa-1),DDMOA_lt2(0:n_fix_tors-1,0:fn_order_qubic_oa_lt2-1)
	real(8),intent(out):: dfc_quad_es(0:n_fcs_quad_es-1,0:n_fix_tors-1),dfc_quad_os(0:n_fcs_quad_os-1,0:n_fix_tors-1),dfc_quad_ea(0:n_fcs_quad_ea-1,0:n_fix_tors-1),dfc_quad_oa(0:n_fcs_quad_oa-1,0:n_fix_tors-1)
	real(8),intent(out):: dfc_qubic_es_gt5(0:n_fcs_qubic_es_gt5-1,0:n_fix_tors-1),dfc_qubic_es_lt5(0:n_fcs_qubic_es_lt5-1,0:n_fix_tors-1),dfc_qubic_es_lt2(0:n_fcs_qubic_es_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_qubic_ea_gt5(0:n_fcs_qubic_ea_gt5-1,0:n_fix_tors-1),dfc_qubic_ea_lt5(0:n_fcs_qubic_ea_lt5-1,0:n_fix_tors-1), dfc_qubic_ea_lt2(0:n_fcs_qubic_ea_lt2-1,0:n_fix_tors-1)
	real(8),intent(out):: dfc_qubic_os_gt5(0:n_fcs_qubic_os_gt5-1,0:n_fix_tors-1),dfc_qubic_os_lt5(0:n_fcs_qubic_os_lt5-1,0:n_fix_tors-1),dfc_qubic_os_lt2(0:n_fcs_qubic_os_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_qubic_oa_gt5(0:n_fcs_qubic_oa_gt5-1,0:n_fix_tors-1),dfc_qubic_oa_lt5(0:n_fcs_qubic_oa_lt5-1,0:n_fix_tors-1),dfc_qubic_oa_lt2(0:n_fcs_qubic_oa_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_quartic_es_gt5(0:n_fcs_quartic_es_gt5-1,0:n_fix_tors-1),dfc_quartic_es_lt5(0:n_fcs_quartic_es_lt5-1,0:n_fix_tors-1),dfc_quartic_es_lt2(0:n_fcs_quartic_es_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_quartic_ea_gt5(0:n_fcs_quartic_ea_gt5-1,0:n_fix_tors-1),dfc_quartic_ea_lt5(0:n_fcs_quartic_ea_lt5-1,0:n_fix_tors-1),dfc_quartic_ea_lt2(0:n_fcs_quartic_ea_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_quartic_os_gt5(0:n_fcs_quartic_os_gt5-1,0:n_fix_tors-1),dfc_quartic_os_lt5(0:n_fcs_quartic_os_lt5-1,0:n_fix_tors-1),dfc_quartic_os_lt2(0:n_fcs_quartic_os_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_quartic_oa_gt5(0:n_fcs_quartic_oa_gt5-1,0:n_fix_tors-1),dfc_quartic_oa_lt5(0:n_fcs_quartic_oa_lt5-1,0:n_fix_tors-1),dfc_quartic_oa_lt2(0:n_fcs_quartic_oa_lt2-1,0:n_fix_tors-1) 
	real(8),intent(out):: dfc_qo_es(0:n_fcs_qo_es-1,0:n_fix_tors-1),dfc_qo_ea(0:n_fcs_qo_ea-1,0:n_fix_tors-1),dfc_qo_os(0:n_fcs_qo_os-1,0:n_fix_tors-1),dfc_qo_oa(0:n_fcs_qo_oa-1,0:n_fix_tors-1)
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
	dfc_qo_es = 0.d0;dfc_qo_ea = 0.d0;dfc_qo_os = 0.d0;dfc_qo_oa = 0.d0 
    !For quadratic
    !ES
	call get_dfc(sinth,costh,Ncoeff_quad_es,fitted_fc_coeff_quad_es,stride_arr_quad_es,fn_order_quad_es,n_fcs_quad_es,DDMES_quad,dfc_quad_es)
    !EA
	call get_dfc(sinth,costh,Ncoeff_quad_ea,fitted_fc_coeff_quad_ea,stride_arr_quad_ea,fn_order_quad_ea,n_fcs_quad_ea,DDMEA_quad,dfc_quad_ea)
    !OS
    call get_dfc(sinth,costh,Ncoeff_quad_os,fitted_fc_coeff_quad_os,stride_arr_quad_os,fn_order_quad_os,n_fcs_quad_os,DDMOS_quad,dfc_quad_os)
    !OA
    call get_dfc(sinth,costh,Ncoeff_quad_oa,fitted_fc_coeff_quad_oa,stride_arr_quad_oa,fn_order_quad_oa,n_fcs_quad_oa,DDMOA_quad,dfc_quad_oa)
   	 
    !For qubic
    !ES
    call get_dfc(sinth,costh,Ncoeff_qubic_es_gt5,fitted_fc_coeff_qubic_es_gt5,stride_arr_qubic_es_gt5,fn_order_qubic_es_gt5,n_fcs_qubic_es_gt5,DDMES_gt5,dfc_qubic_es_gt5)
   	call get_dfc(sinth,costh,Ncoeff_qubic_es_lt5,fitted_fc_coeff_qubic_es_lt5,stride_arr_qubic_es_lt5,fn_order_qubic_es_lt5,n_fcs_qubic_es_lt5,DDMES_gt5,dfc_qubic_es_lt5)
    call get_dfc(sinth,costh,Ncoeff_qubic_es_lt2,fitted_fc_coeff_qubic_es_lt2,stride_arr_qubic_es_lt2,fn_order_qubic_es_lt2,n_fcs_qubic_es_lt2,DDMES_lt2,dfc_qubic_es_lt2)

    !EA
    call get_dfc(sinth,costh,Ncoeff_qubic_ea_gt5,fitted_fc_coeff_qubic_ea_gt5,stride_arr_qubic_ea_gt5,fn_order_qubic_ea_gt5,n_fcs_qubic_ea_gt5,DDMEA_gt5,dfc_qubic_ea_gt5)
   	call get_dfc(sinth,costh,Ncoeff_qubic_ea_lt5,fitted_fc_coeff_qubic_ea_lt5,stride_arr_qubic_ea_lt5,fn_order_qubic_ea_lt5,n_fcs_qubic_ea_lt5,DDMEA_gt5,dfc_qubic_ea_lt5)
    call get_dfc(sinth,costh,Ncoeff_qubic_ea_lt2,fitted_fc_coeff_qubic_ea_lt2,stride_arr_qubic_ea_lt2,fn_order_qubic_ea_lt2,n_fcs_qubic_ea_lt2,DDMEA_lt2,dfc_qubic_ea_lt2)
    
    !OS
    call get_dfc(sinth,costh,Ncoeff_qubic_os_gt5,fitted_fc_coeff_qubic_os_gt5,stride_arr_qubic_os_gt5,fn_order_qubic_os_gt5,n_fcs_qubic_os_gt5,DDMOS_gt5,dfc_qubic_os_gt5)
    call get_dfc(sinth,costh,Ncoeff_qubic_os_lt5,fitted_fc_coeff_qubic_os_lt5,stride_arr_qubic_os_lt5,fn_order_qubic_os_lt5,n_fcs_qubic_os_lt5,DDMOS_gt5,dfc_qubic_os_lt5)
    call get_dfc(sinth,costh,Ncoeff_qubic_os_lt2,fitted_fc_coeff_qubic_os_lt2,stride_arr_qubic_os_lt2,fn_order_qubic_os_lt2,n_fcs_qubic_os_lt2,DDMOS_lt2,dfc_qubic_os_lt2)
    
    !OA
    call get_dfc(sinth,costh,Ncoeff_qubic_oa_gt5,fitted_fc_coeff_qubic_oa_gt5,stride_arr_qubic_oa_gt5,fn_order_qubic_oa_gt5,n_fcs_qubic_oa_gt5,DDMOA_gt5,dfc_qubic_oa_gt5)
    call get_dfc(sinth,costh,Ncoeff_qubic_oa_lt5,fitted_fc_coeff_qubic_oa_lt5,stride_arr_qubic_oa_lt5,fn_order_qubic_oa_lt5,n_fcs_qubic_oa_lt5,DDMOA_gt5,dfc_qubic_oa_lt5)
    call get_dfc(sinth,costh,Ncoeff_qubic_oa_lt2,fitted_fc_coeff_qubic_oa_lt2,stride_arr_qubic_oa_lt2,fn_order_qubic_oa_lt2,n_fcs_qubic_oa_lt2,DDMOA_lt2,dfc_qubic_oa_lt2)
    
    !For quartic
    !ES
    call get_dfc(sinth,costh,Ncoeff_quartic_es_gt5,fitted_fc_coeff_quartic_es_gt5,stride_arr_quartic_es_gt5,fn_order_quartic_es_gt5,n_fcs_quartic_es_gt5,DDMES_gt5,dfc_quartic_es_gt5)
    call get_dfc(sinth,costh,Ncoeff_quartic_es_lt5,fitted_fc_coeff_quartic_es_lt5,stride_arr_quartic_es_lt5,fn_order_quartic_es_lt5,n_fcs_quartic_es_lt5,DDMES_gt5,dfc_quartic_es_lt5)
    call get_dfc(sinth,costh,Ncoeff_quartic_es_lt2,fitted_fc_coeff_quartic_es_lt2,stride_arr_quartic_es_lt2,fn_order_quartic_es_lt2,n_fcs_quartic_es_lt2,DDMES_lt2,dfc_quartic_es_lt2)
    
    !EA
    call get_dfc(sinth,costh,Ncoeff_quartic_ea_gt5,fitted_fc_coeff_quartic_ea_gt5,stride_arr_quartic_ea_gt5,fn_order_quartic_ea_gt5,n_fcs_quartic_ea_gt5,DDMEA_gt5,dfc_quartic_ea_gt5)
    call get_dfc(sinth,costh,Ncoeff_quartic_ea_lt5,fitted_fc_coeff_quartic_ea_lt5,stride_arr_quartic_ea_lt5,fn_order_quartic_ea_lt5,n_fcs_quartic_ea_lt5,DDMEA_gt5,dfc_quartic_ea_lt5)
    call get_dfc(sinth,costh,Ncoeff_quartic_ea_lt2,fitted_fc_coeff_quartic_ea_lt2,stride_arr_quartic_ea_lt2,fn_order_quartic_ea_lt2,n_fcs_quartic_ea_lt2,DDMEA_lt2,dfc_quartic_ea_lt2)
    
    !OS
    call get_dfc(sinth,costh,Ncoeff_quartic_os_gt5,fitted_fc_coeff_quartic_os_gt5,stride_arr_quartic_os_gt5,fn_order_quartic_os_gt5,n_fcs_quartic_os_gt5,DDMOS_gt5,dfc_quartic_os_gt5)
    call get_dfc(sinth,costh,Ncoeff_quartic_os_lt5,fitted_fc_coeff_quartic_os_lt5,stride_arr_quartic_os_lt5,fn_order_quartic_os_lt5,n_fcs_quartic_os_lt5,DDMOS_gt5,dfc_quartic_os_lt5)
    call get_dfc(sinth,costh,Ncoeff_quartic_os_lt2,fitted_fc_coeff_quartic_os_lt2,stride_arr_quartic_os_lt2,fn_order_quartic_os_lt2,n_fcs_quartic_os_lt2,DDMOS_lt2,dfc_quartic_os_lt2)

    !OA
    call get_dfc(sinth,costh,Ncoeff_quartic_oa_gt5,fitted_fc_coeff_quartic_oa_gt5,stride_arr_quartic_oa_gt5,fn_order_quartic_oa_gt5,n_fcs_quartic_oa_gt5,DDMOA_gt5,dfc_quartic_oa_gt5)
    call get_dfc(sinth,costh,Ncoeff_quartic_oa_lt5,fitted_fc_coeff_quartic_oa_lt5,stride_arr_quartic_oa_lt5,fn_order_quartic_oa_lt5,n_fcs_quartic_oa_lt5,DDMOA_gt5,dfc_quartic_oa_lt5)
    call get_dfc(sinth,costh,Ncoeff_quartic_oa_lt2,fitted_fc_coeff_quartic_oa_lt2,stride_arr_quartic_oa_lt2,fn_order_quartic_oa_lt2,n_fcs_quartic_oa_lt2,DDMOA_lt2,dfc_quartic_oa_lt2)
	!QO terms
    call get_dfc(sinth,costh,Ncoeff_qo_es,fitted_fc_coeff_qo_es,stride_arr_qo_es,fn_order_qo_es,n_fcs_qo_es,DDMES_qo,dfc_qo_es)
    call get_dfc(sinth,costh,Ncoeff_qo_ea,fitted_fc_coeff_qo_ea,stride_arr_qo_ea,fn_order_qo_ea,n_fcs_qo_ea,DDMEA_qo,dfc_qo_ea)
    call get_dfc(sinth,costh,Ncoeff_qo_os,fitted_fc_coeff_qo_os,stride_arr_qo_os,fn_order_qo_os,n_fcs_qo_os,DDMOS_qo,dfc_qo_os)
    call get_dfc(sinth,costh,Ncoeff_qo_oa,fitted_fc_coeff_qo_oa,stride_arr_qo_oa,fn_order_qo_oa,n_fcs_qo_oa,DDMOA_qo,dfc_qo_oa)
endsubroutine
!****************************************************************************************************************************************************#
subroutine get_all_fc(sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,DMES_quad,DMES_gt5,DMES_qo,DMES_lt2,DMEA_quad,DMEA_gt5,DMEA_qo,DMEA_lt2,DMOS_quad,DMOS_gt5,DMOS_qo,DMOS_lt2,DMOA_quad,DMOA_gt5,DMOA_qo,DMOA_lt2) 
	implicit none 
	real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(in):: DMES_quad(0:fn_order_quad_es-1),DMES_gt5(0:fn_order_qubic_es_gt5-1),DMES_qo(0:fn_order_qo_es-1),DMES_lt2(0:fn_order_qubic_es_lt2-1) 
	real(8),intent(in):: DMEA_quad(0:fn_order_quad_ea-1),DMEA_gt5(0:fn_order_qubic_ea_gt5-1),DMEA_qo(0:fn_order_qo_ea-1),DMEA_lt2(0:fn_order_qubic_ea_lt2-1)
	real(8),intent(in):: DMOS_quad(0:fn_order_quad_os-1),DMOS_gt5(0:fn_order_qubic_os_gt5-1),DMOS_qo(0:fn_order_qo_os-1),DMOS_lt2(0:fn_order_qubic_os_lt2-1)
	real(8),intent(in):: DMOA_quad(0:fn_order_quad_oa-1),DMOA_gt5(0:fn_order_qubic_oa_gt5-1),DMOA_qo(0:fn_order_qo_oa-1),DMOA_lt2(0:fn_order_qubic_oa_lt2-1)
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
	real(8),intent(out):: fc_qo_es(0:n_fcs_qo_es-1),fc_qo_ea(0:n_fcs_qo_ea-1),fc_qo_os(0:n_fcs_qo_os-1),fc_qo_oa(0:n_fcs_qo_oa-1)
	integer:: i
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
	fc_quartic_es_ltpt1 = 0.d0; fc_quartic_ea_ltpt1 = 0.d0;
	fc_quartic_os_ltpt1 = 0.d0; fc_quartic_oa_ltpt1 = 0.d0;
	fc_qo_es = 0.d0;fc_qo_ea = 0.d0;fc_qo_os = 0.d0;fc_qo_oa = 0.d0 
    !For quadratic
    !ES
    call get_fc(sinth,costh,Ncoeff_quad_es,fitted_fc_coeff_quad_es,stride_arr_quad_es,fn_order_quad_es,n_fcs_quad_es,DMES_quad,fc_quad_es)
    !EA                                                                                                                   
    call get_fc(sinth,costh,Ncoeff_quad_ea,fitted_fc_coeff_quad_ea,stride_arr_quad_ea,fn_order_quad_ea,n_fcs_quad_ea,DMEA_quad,fc_quad_ea)
    !OS                                                                                                                   
    call get_fc(sinth,costh,Ncoeff_quad_os,fitted_fc_coeff_quad_os,stride_arr_quad_os,fn_order_quad_os,n_fcs_quad_os,DMOS_quad,fc_quad_os)
    !OA                                                                                                                   
    call get_fc(sinth,costh,Ncoeff_quad_oa,fitted_fc_coeff_quad_oa,stride_arr_quad_oa,fn_order_quad_oa,n_fcs_quad_oa,DMOA_quad,fc_quad_oa)
    
    !For qubic
    !ES
    call get_fc(sinth,costh,Ncoeff_qubic_es_gt5,fitted_fc_coeff_qubic_es_gt5,stride_arr_qubic_es_gt5,fn_order_qubic_es_gt5,n_fcs_qubic_es_gt5,DMES_gt5,fc_qubic_es_gt5)
   	call get_fc(sinth,costh,Ncoeff_qubic_es_lt5,fitted_fc_coeff_qubic_es_lt5,stride_arr_qubic_es_lt5,fn_order_qubic_es_lt5,n_fcs_qubic_es_lt5,DMES_gt5,fc_qubic_es_lt5)
    call get_fc(sinth,costh,Ncoeff_qubic_es_lt2,fitted_fc_coeff_qubic_es_lt2,stride_arr_qubic_es_lt2,fn_order_qubic_es_lt2,n_fcs_qubic_es_lt2,DMES_lt2,fc_qubic_es_lt2)
    
    !EA
    call get_fc(sinth,costh,Ncoeff_qubic_ea_gt5,fitted_fc_coeff_qubic_ea_gt5,stride_arr_qubic_ea_gt5,fn_order_qubic_ea_gt5,n_fcs_qubic_ea_gt5,DMEA_gt5,fc_qubic_ea_gt5)
   	call get_fc(sinth,costh,Ncoeff_qubic_ea_lt5,fitted_fc_coeff_qubic_ea_lt5,stride_arr_qubic_ea_lt5,fn_order_qubic_ea_lt5,n_fcs_qubic_ea_lt5,DMEA_gt5,fc_qubic_ea_lt5)
    call get_fc(sinth,costh,Ncoeff_qubic_ea_lt2,fitted_fc_coeff_qubic_ea_lt2,stride_arr_qubic_ea_lt2,fn_order_qubic_ea_lt2,n_fcs_qubic_ea_lt2,DMEA_lt2,fc_qubic_ea_lt2)
    
    !OS
    call get_fc(sinth,costh,Ncoeff_qubic_os_gt5,fitted_fc_coeff_qubic_os_gt5,stride_arr_qubic_os_gt5,fn_order_qubic_os_gt5,n_fcs_qubic_os_gt5,DMOS_gt5,fc_qubic_os_gt5)
    call get_fc(sinth,costh,Ncoeff_qubic_os_lt5,fitted_fc_coeff_qubic_os_lt5,stride_arr_qubic_os_lt5,fn_order_qubic_os_lt5,n_fcs_qubic_os_lt5,DMOS_gt5,fc_qubic_os_lt5)
    call get_fc(sinth,costh,Ncoeff_qubic_os_lt2,fitted_fc_coeff_qubic_os_lt2,stride_arr_qubic_os_lt2,fn_order_qubic_os_lt2,n_fcs_qubic_os_lt2,DMOS_lt2,fc_qubic_os_lt2)
    
    !OA
    call get_fc(sinth,costh,Ncoeff_qubic_oa_gt5,fitted_fc_coeff_qubic_oa_gt5,stride_arr_qubic_oa_gt5,fn_order_qubic_oa_gt5,n_fcs_qubic_oa_gt5,DMOA_gt5,fc_qubic_oa_gt5)
    call get_fc(sinth,costh,Ncoeff_qubic_oa_lt5,fitted_fc_coeff_qubic_oa_lt5,stride_arr_qubic_oa_lt5,fn_order_qubic_oa_lt5,n_fcs_qubic_oa_lt5,DMOA_gt5,fc_qubic_oa_lt5)
    call get_fc(sinth,costh,Ncoeff_qubic_oa_lt2,fitted_fc_coeff_qubic_oa_lt2,stride_arr_qubic_oa_lt2,fn_order_qubic_oa_lt2,n_fcs_qubic_oa_lt2,DMOA_lt2,fc_qubic_oa_lt2)
    
    !For quartic
    !ES
    call get_fc(sinth,costh,Ncoeff_quartic_es_gt5,fitted_fc_coeff_quartic_es_gt5,stride_arr_quartic_es_gt5,fn_order_quartic_es_gt5,n_fcs_quartic_es_gt5,DMES_gt5,fc_quartic_es_gt5)
    call get_fc(sinth,costh,Ncoeff_quartic_es_lt5,fitted_fc_coeff_quartic_es_lt5,stride_arr_quartic_es_lt5,fn_order_quartic_es_lt5,n_fcs_quartic_es_lt5,DMES_gt5,fc_quartic_es_lt5)
    call get_fc(sinth,costh,Ncoeff_quartic_es_lt2,fitted_fc_coeff_quartic_es_lt2,stride_arr_quartic_es_lt2,fn_order_quartic_es_lt2,n_fcs_quartic_es_lt2,DMES_lt2,fc_quartic_es_lt2)
	fc_quartic_es_ltpt1 = fitted_fc_coeff_quartic_es_ltpt1
    !EA
    call get_fc(sinth,costh,Ncoeff_quartic_ea_gt5,fitted_fc_coeff_quartic_ea_gt5,stride_arr_quartic_ea_gt5,fn_order_quartic_ea_gt5,n_fcs_quartic_ea_gt5,DMEA_gt5,fc_quartic_ea_gt5)
    call get_fc(sinth,costh,Ncoeff_quartic_ea_lt5,fitted_fc_coeff_quartic_ea_lt5,stride_arr_quartic_ea_lt5,fn_order_quartic_ea_lt5,n_fcs_quartic_ea_lt5,DMEA_gt5,fc_quartic_ea_lt5)
    call get_fc(sinth,costh,Ncoeff_quartic_ea_lt2,fitted_fc_coeff_quartic_ea_lt2,stride_arr_quartic_ea_lt2,fn_order_quartic_ea_lt2,n_fcs_quartic_ea_lt2,DMEA_lt2,fc_quartic_ea_lt2)
	fc_quartic_ea_ltpt1 = fitted_fc_coeff_quartic_ea_ltpt1
    
    !OS
    call get_fc(sinth,costh,Ncoeff_quartic_os_gt5,fitted_fc_coeff_quartic_os_gt5,stride_arr_quartic_os_gt5,fn_order_quartic_os_gt5,n_fcs_quartic_os_gt5,DMOS_gt5,fc_quartic_os_gt5)
    call get_fc(sinth,costh,Ncoeff_quartic_os_lt5,fitted_fc_coeff_quartic_os_lt5,stride_arr_quartic_os_lt5,fn_order_quartic_os_lt5,n_fcs_quartic_os_lt5,DMOS_gt5,fc_quartic_os_lt5)
    call get_fc(sinth,costh,Ncoeff_quartic_os_lt2,fitted_fc_coeff_quartic_os_lt2,stride_arr_quartic_os_lt2,fn_order_quartic_os_lt2,n_fcs_quartic_os_lt2,DMOS_lt2,fc_quartic_os_lt2)
	fc_quartic_os_ltpt1 = fitted_fc_coeff_quartic_os_ltpt1
    !OA
    call get_fc(sinth,costh,Ncoeff_quartic_oa_gt5,fitted_fc_coeff_quartic_oa_gt5,stride_arr_quartic_oa_gt5,fn_order_quartic_oa_gt5,n_fcs_quartic_oa_gt5,DMOA_gt5,fc_quartic_oa_gt5)
    call get_fc(sinth,costh,Ncoeff_quartic_oa_lt5,fitted_fc_coeff_quartic_oa_lt5,stride_arr_quartic_oa_lt5,fn_order_quartic_oa_lt5,n_fcs_quartic_oa_lt5,DMOA_gt5,fc_quartic_oa_lt5)
    call get_fc(sinth,costh,Ncoeff_quartic_oa_lt2,fitted_fc_coeff_quartic_oa_lt2,stride_arr_quartic_oa_lt2,fn_order_quartic_oa_lt2,n_fcs_quartic_oa_lt2,DMOA_lt2,fc_quartic_oa_lt2)
	fc_quartic_oa_ltpt1 = fitted_fc_coeff_quartic_oa_ltpt1
  	!QO
	call get_fc(sinth,costh,Ncoeff_qo_es,fitted_fc_coeff_qo_es,stride_arr_qo_es,fn_order_qo_es,n_fcs_qo_es,DMES_qo,fc_qo_es)
	call get_fc(sinth,costh,Ncoeff_qo_ea,fitted_fc_coeff_qo_ea,stride_arr_qo_ea,fn_order_qo_ea,n_fcs_qo_ea,DMEA_qo,fc_qo_ea)
	call get_fc(sinth,costh,Ncoeff_qo_os,fitted_fc_coeff_qo_os,stride_arr_qo_os,fn_order_qo_os,n_fcs_qo_os,DMOS_qo,fc_qo_os)
	call get_fc(sinth,costh,Ncoeff_qo_oa,fitted_fc_coeff_qo_oa,stride_arr_qo_oa,fn_order_qo_oa,n_fcs_qo_oa,DMOA_qo,fc_qo_oa)
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
subroutine getDDesignMatES(sinth,costh,fn_order,Ncoeff,stride_arr,DDesignMat)
    implicit none
	integer::ii,i,j,k,ctr
	integer,intent(in)::fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(out):: DDesignMat(0:n_fix_tors-1,0:fn_order-1) 
	DDesignMat  = 0.d0
	do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
    		DDesignmat(0,ii) =   - dble(i) * sinth(0,i) * costh(1,j) * costh(2,k) 
    		DDesignmat(1,ii) =   - dble(j) * costh(0,i) * sinth(1,j) * costh(2,k)
    		DDesignmat(2,ii) =   - dble(k) * costh(0,i) * costh(1,j) * sinth(2,k)
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
			DDesignmat(0,ii) =   - dble(i) * sinth(0,i) * (costh(1,j) * costh(2,k) + costh(1,k) * costh(2,j))
            DDesignmat(1,ii) =   costh(0,i) * (- dble(j) * sinth(1,j) * costh(2,k) - dble(k) * sinth(1,k) * costh(2,j))
            DDesignmat(2,ii) =   costh(0,i) * (- dble(k) * costh(1,j) * sinth(2,k) - dble(j) * sinth(2,j) * costh(1,k))

        elseif ((ii>=Ncoeff(1)).and.(ii <Ncoeff(2))) then
    		DDesignmat(0,ii) =	- dble(i) * sinth(0,i) * sinth(1,j) * sinth(2,k) 
    		DDesignmat(1,ii) =  + dble(j) * costh(0,i) * costh(1,j) * sinth(2,k)
    		DDesignmat(2,ii) =  + dble(k) * costh(0,i) * sinth(1,j) * costh(2,k)
        elseif ((ii>=Ncoeff(2)).and.(ii <Ncoeff(3))) then
    		DDesignmat(0,ii) =	- dble(i) * sinth(0,i) * (sinth(1,j) * sinth(2,k) +  sinth(2,j) * sinth(1,k))
    		DDesignmat(1,ii) =  costh(0,i) * (dble(j) * costh(1,j) * sinth(2,k) +  dble(k) * costh(1,k) * sinth(2,j))
    		DDesignmat(2,ii) =  costh(0,i) * (dble(k) * sinth(1,j) * costh(2,k) +  dble(j) * sinth(1,k) * costh(2,j))
        else
			DDesignmat(0,ii) =  dble(i)    * costh(0,i) * (sinth(1,k) * costh(2,j) + sinth(2,k) * costh(1,j))
			DDesignmat(1,ii) =  sinth(0,i) * (dble(k) * costh(1,k) * costh(2,j) - dble(j) * sinth(2,k) * sinth(1,j))
			DDesignmat(2,ii) =  sinth(0,i) * ( - dble(j) * sinth(1,k) * sinth(2,j) + dble(k) * costh(2,k) * costh(1,j))
		endif
	enddo
endsubroutine
!****************************************************************************************#
subroutine getDDesignMatEA(sinth,costh,fn_order,Ncoeff,stride_arr,DDesignMat)
    implicit none
	integer::ii,i,j,k,ctr
	integer,intent(in)::fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(out):: DDesignMat(0:n_fix_tors-1,0:fn_order-1) 
	DDesignMat  = 0.d0
	do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
    		DDesignmat(0,ii) =   - dble(i) * sinth(0,i) * (costh(1,j) * costh(2,k) - costh(1,k) * costh(2,j))
    		DDesignmat(1,ii) =   costh(0,i) * (-dble(j) * sinth(1,j) * costh(2,k) + dble(k) * sinth(1,k) * costh(2,j))
    		DDesignmat(2,ii) =   costh(0,i) * (-dble(k) * sinth(2,k) * costh(1,j) + dble(j) * sinth(2,j) * costh(1,k))
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
    		DDesignmat(0,ii) =    - dble(i) * sinth(0,i) * (sinth(1,j) * sinth(2,k) - sinth(2,j) * sinth(1,k))
    		DDesignmat(1,ii) =   costh(0,i) * (dble(j) * costh(1,j) * sinth(2,k) - dble(k) * costh(1,k) * sinth(2,j))
    		DDesignmat(2,ii) =   costh(0,i) * (dble(k) * costh(2,k) * sinth(1,j) - dble(j) * costh(2,j) * sinth(1,k))
        else
    		DDesignmat(0,ii) =   dble(i) * costh(0,i) * (sinth(1,k) * costh(2,j) - sinth(2,k) * costh(1,j))
    		DDesignmat(1,ii) =   sinth(0,i) * ( dble(k) * costh(1,k) * costh(2,j) +dble(j) *  sinth(2,k) * sinth(1,j))
    		DDesignmat(2,ii) =    sinth(0,i) * ( -dble(j) *  sinth(1,k) * sinth(2,j)  - dble(k) * costh(2,k) * costh(1,j))
		endif
	enddo
	
endsubroutine
!****************************************************************************************#
subroutine getDDesignMatOS(sinth,costh,fn_order,Ncoeff,stride_arr,DDesignMat)
    implicit none
	integer::ii,i,j,k,ctr
	integer,intent(in)::fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(out):: DDesignMat(0:n_fix_tors-1,0:fn_order-1) 
	DDesignMat  = 0.d0
	do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
    		DDesignmat(0,ii) =  dble(i) * costh(0,i) * costh(1,j) * costh(2,k) 
    		DDesignmat(1,ii) =  -dble(j) * sinth(0,i)  * sinth(1,j) * costh(2,k) 
    		DDesignmat(2,ii) =   -dble(k) *sinth(0,i)  * costh(1,j) *sinth(2,k)
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
    		DDesignmat(0,ii) =  dble(i) * costh(0,i)   * (costh(1,j) * costh(2,k) + costh(1,k) * costh(2,j))
    		DDesignmat(1,ii) =  sinth(0,i) * ( -dble(j) * sinth(1,j)  * costh(2,k) - dble(k) * sinth(1,k)  * costh(2,j)) 
    		DDesignmat(2,ii) =  sinth(0,i) * ( -dble(k) * sinth(2,k)  * costh(1,j) - dble(j) * sinth(2,j)  * costh(1,k)) 
        elseif ((ii>=Ncoeff(1)).and.(ii <Ncoeff(2))) then
    		DDesignmat(0,ii) =   dble(i) * costh(0,i) *   sinth(1,j) * sinth(2,k) 
    		DDesignmat(1,ii) =   dble(j) * sinth(0,i)  *  costh(1,j) * sinth(2,k)
    		DDesignmat(2,ii) =   dble(k) * sinth(0,i)  *  sinth(1,j) * costh(2,k)
        elseif ((ii>=Ncoeff(2)).and.(ii <Ncoeff(3))) then
    		DDesignmat(0,ii) =  dble(i) * costh(0,i) * ( sinth(1,j) * sinth(2,k) + sinth(1,k) * sinth(2,j))
    		DDesignmat(1,ii) =  sinth(0,i)  * (dble(j) * costh(1,j) * sinth(2,k) + dble(k) * costh(1,k) * sinth(2,j)) 
    		DDesignmat(2,ii) =  sinth(0,i)  * (dble(k) * costh(2,k) * sinth(1,j)  + dble(j) * costh(2,j) * sinth(1,k))
        else
    		DDesignmat(0,ii) =   -dble(i) * sinth(0,i) * (sinth(1,k) * costh(2,j) + costh(1,j) * sinth(2,k))
    		DDesignmat(1,ii) =   costh(0,i) * (dble(k) * costh(1,k) * costh(2,j) - dble(j) * sinth(1,j) * sinth(2,k))
    		DDesignmat(2,ii) =   costh(0,i) * (-dble(j) * sinth(1,k) * sinth(2,j) + dble(k) * costh(1,j) * costh(2,k))
		endif
	enddo
endsubroutine
!****************************************************************************************#
subroutine getDDesignMatOA(sinth,costh,fn_order,Ncoeff,stride_arr,DDesignMat)
    implicit none
	integer::ii,i,j,k,ctr
	integer,intent(in)::fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(out):: DDesignMat(0:n_fix_tors-1,0:fn_order-1) 
	DDesignMat  = 0.d0
	do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
    		DDesignmat(0,ii) =  dble(i) * costh(0,i) * ( costh(1,j) * costh(2,k) - costh(1,k) * costh(2,j)) 
    		DDesignmat(1,ii) =   sinth(0,i) * ( -dble(j) * sinth(1,j) * costh(2,k) + dble(k) * sinth(1,k) * costh(2,j))
    		DDesignmat(2,ii) =  sinth(0,i) * ( -dble(k) * costh(1,j) * sinth(2,k) + dble(j) * costh(1,k) * sinth(2,j))
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
    		DDesignmat(0,ii) =    dble(i) * costh(0,i) * (sinth(1,j) * sinth(2,k) - sinth(2,j) * sinth(1,k))
    		DDesignmat(1,ii) =   sinth(0,i) * (dble(j) * costh(1,j) * sinth(2,k) - dble(k) *  sinth(2,j) *costh(1,k)  )
    		DDesignmat(2,ii) =   sinth(0,i) * (dble(k) * sinth(1,j) * costh(2,k) - dble(j) * costh(2,j) * sinth(1,k))
        else
    		DDesignmat(0,ii) =   - dble(i) * sinth(0,i) * (sinth(1,k) * costh(2,j) - sinth(2,k) * costh(1,j))
    		DDesignmat(1,ii) =    costh(0,i) * ( dble(k) * costh(1,k) * costh(2,j) + dble(j) * sinth(2,k) * sinth(1,j))
    		DDesignmat(2,ii) =   costh(0,i) * (-dble(j) * sinth(1,k) * sinth(2,j) - dble(k) * costh(2,k) * costh(1,j))
		endif
	enddo
endsubroutine
!****************************************************************************************#
!****************************************************************************************#
subroutine getDesignMatES(sinth,costh,fn_order,Ncoeff,stride_arr,DesignMat)
    implicit none
	integer::ii,i,j,k,ctr
	integer,intent(in)::fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(out):: DesignMat(0:fn_order-1) 
	DesignMat  = 0.d0
	do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
    		DesignMat(ii) = costh(0,i) * costh(1,j) * costh(2,k) 
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
    		DesignMat(ii) = costh(0,i) * (costh(1,j) * costh(2,k) + costh(1,k) * costh(2,j))
        elseif ((ii>=Ncoeff(1)).and.(ii <Ncoeff(2))) then
    		DesignMat(ii) = costh(0,i) * sinth(1,j) * sinth(2,k)
        elseif ((ii>=Ncoeff(2)).and.(ii <Ncoeff(3))) then
    		DesignMat(ii) = costh(0,i) * (sinth(1,j) * sinth(2,k) + sinth(1,k) * sinth(2,j))
        else
    		DesignMat(ii) = sinth(0,i) * (sinth(1,k) * costh(2,j) + costh(1,j) * sinth(2,k))
		endif
	enddo
endsubroutine
!****************************************************************************************#
subroutine getDesignMatEA(sinth,costh,fn_order,Ncoeff,stride_arr,DesignMat)
    implicit none
	integer::ii,i,j,k,ctr
	integer,intent(in)::fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(out):: DesignMat(0:fn_order-1) 
	DesignMat  = 0.d0
	do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
			DesignMat(ii) = costh(0,i) * ( costh(1,j)*costh(2,k) - costh(1,k)*costh(2,j) ) 
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
			DesignMat(ii) = costh(0,i) * ( sinth(1,j)*sinth(2,k) - sinth(2,j)*sinth(1,k) )
        else
			DesignMat(ii) = sinth(0,i) * ( sinth(1,k)*costh(2,j) - sinth(2,k)*costh(1,j) )
		endif
	enddo
endsubroutine
!****************************************************************************************#
subroutine getDesignMatOS(sinth,costh,fn_order,Ncoeff,stride_arr,DesignMat)
    implicit none
	integer::ii,i,j,k,ctr
	integer,intent(in)::fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
	real(8),intent(out):: DesignMat(0:fn_order-1) 
	DesignMat  = 0.d0
	do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
			DesignMat(ii) = sinth(0,i) * costh(1,j) * costh(2,k) 
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
			DesignMat(ii) = sinth(0,i) * (costh(1,j) * costh(2,k) + costh(1,k) * costh(2,j))
        elseif ((ii>=Ncoeff(1)).and.(ii <Ncoeff(2))) then
			DesignMat(ii) = sinth(0,i) * sinth(1,j) * sinth(2,k)
        elseif ((ii>=Ncoeff(2)).and.(ii <Ncoeff(3))) then
			DesignMat(ii) = sinth(0,i) * ( sinth(1,j)*sinth(2,k) + sinth(1,k)*sinth(2,j) )
        else
			DesignMat(ii) = costh(0,i) * ( sinth(1,k)*costh(2,j) + costh(1,j)*sinth(2,k) )
		endif
	enddo
endsubroutine
!****************************************************************************************#
subroutine getDesignMatOA(sinth,costh,fn_order,Ncoeff,stride_arr,DesignMat)
    implicit none
	integer::ii,i,j,k,ctr
	integer,intent(in)::fn_order,Ncoeff(0:4),stride_arr(0:fn_order-1,0:2)
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(out):: DesignMat(0:fn_order-1) 
	DesignMat  = 0.d0
	do ii = 0,fn_order-1
        i = stride_arr(ii,0);j = stride_arr(ii,1);k = stride_arr(ii,2)
        if (ii<Ncoeff(0)) then
			DesignMat(ii)  = sinth(0,i) * ( costh(1,j)*costh(2,k) - costh(1,k)*costh(2,j) ) 
        elseif ((ii>=Ncoeff(0)).and.( ii<Ncoeff(1))) then
			DesignMat(ii)  = sinth(0,i) * ( sinth(1,j)*sinth(2,k) - sinth(2,j)*sinth(1,k) )
        else
			DesignMat(ii)  = costh(0,i) * ( sinth(1,k)*costh(2,j) - sinth(2,k)*costh(1,j) )
		endif
	enddo
endsubroutine
!****************************************************************************************#
subroutine getAllDM(sinth,costh,DMES_eqm,DMES_quad,DMES_gt5,DMES_qo,DMES_lt2,DMEA_eqm,DMEA_quad,DMEA_gt5,DMEA_qo,DMEA_lt2,DMOS_eqm,DMOS_quad,DMOS_gt5,DMOS_qo,DMOS_lt2,DMOA_eqm,DMOA_quad,DMOA_gt5,DMOA_qo,DMOA_lt2)
	implicit none
	integer::Nterms 
	real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(out):: DMES_eqm(0:fn_order_bonds_es-1),DMES_quad(0:fn_order_quad_es-1),DMES_gt5(0:fn_order_qubic_es_gt5-1),DMES_qo(0:fn_order_qo_es-1),DMES_lt2(0:fn_order_qubic_es_lt2-1) 
	real(8),intent(out):: DMEA_eqm(0:fn_order_bonds_ea-1),DMEA_quad(0:fn_order_quad_ea-1),DMEA_gt5(0:fn_order_qubic_ea_gt5-1),DMEA_qo(0:fn_order_qo_ea-1),DMEA_lt2(0:fn_order_qubic_ea_lt2-1)
	real(8),intent(out):: DMOS_eqm(0:fn_order_bonds_os-1),DMOS_quad(0:fn_order_quad_os-1),DMOS_gt5(0:fn_order_qubic_os_gt5-1),DMOS_qo(0:fn_order_qo_os-1),DMOS_lt2(0:fn_order_qubic_os_lt2-1)
	real(8),intent(out):: DMOA_eqm(0:fn_order_bonds_oa-1),DMOA_quad(0:fn_order_quad_oa-1),DMOA_gt5(0:fn_order_qubic_oa_gt5-1),DMOA_qo(0:fn_order_qo_oa-1),DMOA_lt2(0:fn_order_qubic_oa_lt2-1)
	DMES_eqm=0.d0;DMES_quad=0.d0;DMES_gt5=0.d0;DMES_qo=0.d0;DMES_lt2=0.d0 
	DMEA_eqm=0.d0;DMEA_quad=0.d0;DMEA_gt5=0.d0;DMEA_qo=0.d0;DMEA_lt2=0.d0
	DMOS_eqm=0.d0;DMOS_quad=0.d0;DMOS_gt5=0.d0;DMOS_qo=0.d0;DMOS_lt2=0.d0
	DMOA_eqm=0.d0;DMOA_quad=0.d0;DMOA_gt5=0.d0;DMOA_qo=0.d0;DMOA_lt2=0.d0

	Nterms = NtermsToNcoeffES(0,0)
	call getDesignMatES(sinth,costh,Nterms,Ncoeff_bonds_es,stride_arr_b_es,DMES_eqm)	
	Nterms = NtermsToNcoeffEA(0,0)
	call getDesignMatEA(sinth,costh,Nterms,Ncoeff_bonds_ea,stride_arr_b_ea,DMEA_eqm)	
	Nterms = NtermsToNcoeffOS(0,0)
	call getDesignMatOS(sinth,costh,Nterms,Ncoeff_bonds_os,stride_arr_b_os,DMOS_eqm)	
	Nterms = NtermsToNcoeffOA(0,0)
	call getDesignMatOA(sinth,costh,Nterms,Ncoeff_bonds_oa,stride_arr_b_oa,DMOA_eqm)	

	Nterms = NtermsToNcoeffES(1,0)
	call getDesignMatES(sinth,costh,Nterms,Ncoeff_quad_es,stride_arr_quad_es,DMES_quad)	
	Nterms = NtermsToNcoeffEA(1,0)
	call getDesignMatEA(sinth,costh,Nterms,Ncoeff_quad_ea,stride_arr_quad_ea,DMEA_quad)	
	Nterms = NtermsToNcoeffOS(1,0)
	call getDesignMatOS(sinth,costh,Nterms,Ncoeff_quad_os,stride_arr_quad_os,DMOS_quad)	
	Nterms = NtermsToNcoeffOA(1,0)
	call getDesignMatOA(sinth,costh,Nterms,Ncoeff_quad_oa,stride_arr_quad_oa,DMOA_quad)	

	Nterms = NtermsToNcoeffES(2,0)
	call getDesignMatES(sinth,costh,Nterms,Ncoeff_qubic_es_gt5,stride_arr_qubic_es_gt5,DMES_gt5)	
	Nterms = NtermsToNcoeffEA(2,0)
	call getDesignMatEA(sinth,costh,Nterms,Ncoeff_qubic_ea_gt5,stride_arr_qubic_ea_gt5,DMEA_gt5)	
	Nterms = NtermsToNcoeffOS(2,0)
	call getDesignMatOS(sinth,costh,Nterms,Ncoeff_qubic_os_gt5,stride_arr_qubic_os_gt5,DMOS_gt5)	
	Nterms = NtermsToNcoeffOA(2,0)
	call getDesignMatOA(sinth,costh,Nterms,Ncoeff_qubic_oa_gt5,stride_arr_qubic_oa_gt5,DMOA_gt5)	

	Nterms = NtermsToNcoeffES(5,0)
	call getDesignMatES(sinth,costh,Nterms,Ncoeff_qo_es,stride_arr_qo_es,DMES_qo)	
	Nterms = NtermsToNcoeffEA(5,0)                                                    
	call getDesignMatEA(sinth,costh,Nterms,Ncoeff_qo_ea,stride_arr_qo_ea,DMEA_qo)	
	Nterms = NtermsToNcoeffOS(5,0)                                                    
	call getDesignMatOS(sinth,costh,Nterms,Ncoeff_qo_os,stride_arr_qo_os,DMOS_qo)	
	Nterms = NtermsToNcoeffOA(5,0)                                                    
	call getDesignMatOA(sinth,costh,Nterms,Ncoeff_qo_oa,stride_arr_qo_oa,DMOA_qo)	

	DMES_lt2 = DMES_qo
	DMEA_lt2 = DMEA_qo
	DMOS_lt2 = DMOS_qo
	DMOA_lt2 = DMOA_qo
	
	!Nterms = NtermsToNcoeffES(4,0)
	!call getDesignMatES(sinth,costh,Nterms,Ncoeff_qubic_es_lt2,stride_arr_qubic_es_lt2,DMES_lt2)	
	!Nterms = NtermsToNcoeffEA(4,0)
	!call getDesignMatEA(sinth,costh,Nterms,Ncoeff_qubic_ea_lt2,stride_arr_qubic_ea_lt2,DMEA_lt2)	
	!Nterms = NtermsToNcoeffOS(4,0)
	!call getDesignMatOS(sinth,costh,Nterms,Ncoeff_qubic_os_lt2,stride_arr_qubic_os_lt2,DMOS_lt2)	
	!Nterms = NtermsToNcoeffOA(4,0)
	!call getDesignMatOA(sinth,costh,Nterms,Ncoeff_qubic_oa_lt2,stride_arr_qubic_oa_lt2,DMOA_lt2)	
endsubroutine 
!****************************************************************************************#
subroutine getAllDDM(sinth,costh,DDMES_eqm,DDMES_quad,DDMES_gt5,DDMES_qo,DDMES_lt2,DDMEA_eqm,DDMEA_quad,DDMEA_gt5,DDMEA_qo,DDMEA_lt2,DDMOS_eqm,DDMOS_quad,DDMOS_gt5,DDMOS_qo,DDMOS_lt2,DDMOA_eqm,DDMOA_quad,DDMOA_gt5,DDMOA_qo,DDMOA_lt2)
	implicit none
	integer::Nterms 
	real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
	real(8),intent(out):: DDMES_eqm(0:n_fix_tors-1,0:fn_order_bonds_es-1),DDMES_quad(0:n_fix_tors-1,0:fn_order_quad_es-1),DDMES_gt5(0:n_fix_tors-1,0:fn_order_qubic_es_gt5-1),DDMES_qo(0:n_fix_tors-1,0:fn_order_qo_es-1),DDMES_lt2(0:n_fix_tors-1,0:fn_order_qubic_es_lt2-1)
	real(8),intent(out):: DDMEA_eqm(0:n_fix_tors-1,0:fn_order_bonds_ea-1),DDMEA_quad(0:n_fix_tors-1,0:fn_order_quad_ea-1),DDMEA_gt5(0:n_fix_tors-1,0:fn_order_qubic_ea_gt5-1),DDMEA_qo(0:n_fix_tors-1,0:fn_order_qo_ea-1),DDMEA_lt2(0:n_fix_tors-1,0:fn_order_qubic_ea_lt2-1)
	real(8),intent(out):: DDMOS_eqm(0:n_fix_tors-1,0:fn_order_bonds_os-1),DDMOS_quad(0:n_fix_tors-1,0:fn_order_quad_os-1),DDMOS_gt5(0:n_fix_tors-1,0:fn_order_qubic_os_gt5-1),DDMOS_qo(0:n_fix_tors-1,0:fn_order_qo_os-1),DDMOS_lt2(0:n_fix_tors-1,0:fn_order_qubic_os_lt2-1)
	real(8),intent(out):: DDMOA_eqm(0:n_fix_tors-1,0:fn_order_bonds_oa-1),DDMOA_quad(0:n_fix_tors-1,0:fn_order_quad_oa-1),DDMOA_gt5(0:n_fix_tors-1,0:fn_order_qubic_oa_gt5-1),DDMOA_qo(0:n_fix_tors-1,0:fn_order_qo_oa-1),DDMOA_lt2(0:n_fix_tors-1,0:fn_order_qubic_oa_lt2-1)
	DDMES_eqm=0.d0;DDMES_quad=0.d0;DDMES_gt5=0.d0;DDMES_qo=0.d0;DDMES_lt2=0.d0 
	DDMEA_eqm=0.d0;DDMEA_quad=0.d0;DDMEA_gt5=0.d0;DDMEA_qo=0.d0;DDMEA_lt2=0.d0
	DDMOS_eqm=0.d0;DDMOS_quad=0.d0;DDMOS_gt5=0.d0;DDMOS_qo=0.d0;DDMOS_lt2=0.d0
	DDMOA_eqm=0.d0;DDMOA_quad=0.d0;DDMOA_gt5=0.d0;DDMOA_qo=0.d0;DDMOA_lt2=0.d0
	Nterms = NtermsToNcoeffES(0,0)
	call getDDesignMatES(sinth,costh,Nterms,Ncoeff_bonds_es,stride_arr_b_es,DDMES_eqm)	
	Nterms = NtermsToNcoeffEA(0,0)
	call getDDesignMatEA(sinth,costh,Nterms,Ncoeff_bonds_ea,stride_arr_b_ea,DDMEA_eqm)	
	Nterms = NtermsToNcoeffOS(0,0)
	call getDDesignMatOS(sinth,costh,Nterms,Ncoeff_bonds_os,stride_arr_b_os,DDMOS_eqm)	
	Nterms = NtermsToNcoeffOA(0,0)
	call getDDesignMatOA(sinth,costh,Nterms,Ncoeff_bonds_oa,stride_arr_b_oa,DDMOA_eqm)	

	Nterms = NtermsToNcoeffES(1,0)
	call getDDesignMatES(sinth,costh,Nterms,Ncoeff_quad_es,stride_arr_quad_es,DDMES_quad)	
	Nterms = NtermsToNcoeffEA(1,0)
	call getDDesignMatEA(sinth,costh,Nterms,Ncoeff_quad_ea,stride_arr_quad_ea,DDMEA_quad)	
	Nterms = NtermsToNcoeffOS(1,0)
	call getDDesignMatOS(sinth,costh,Nterms,Ncoeff_quad_os,stride_arr_quad_os,DDMOS_quad)	
	Nterms = NtermsToNcoeffOA(1,0)
	call getDDesignMatOA(sinth,costh,Nterms,Ncoeff_quad_oa,stride_arr_quad_oa,DDMOA_quad)	

	Nterms = NtermsToNcoeffES(2,0)
	call getDDesignMatES(sinth,costh,Nterms,Ncoeff_qubic_es_gt5,stride_arr_qubic_es_gt5,DDMES_gt5)	
	Nterms = NtermsToNcoeffEA(2,0)
	call getDDesignMatEA(sinth,costh,Nterms,Ncoeff_qubic_ea_gt5,stride_arr_qubic_ea_gt5,DDMEA_gt5)	
	Nterms = NtermsToNcoeffOS(2,0)
	call getDDesignMatOS(sinth,costh,Nterms,Ncoeff_qubic_os_gt5,stride_arr_qubic_os_gt5,DDMOS_gt5)	
	Nterms = NtermsToNcoeffOA(2,0)
	call getDDesignMatOA(sinth,costh,Nterms,Ncoeff_qubic_oa_gt5,stride_arr_qubic_oa_gt5,DDMOA_gt5)	

	Nterms = NtermsToNcoeffES(5,0)
	call getDDesignMatES(sinth,costh,Nterms,Ncoeff_qo_es,stride_arr_qo_es,DDMES_qo)	
	Nterms = NtermsToNcoeffEA(5,0)                                                    
	call getDDesignMatEA(sinth,costh,Nterms,Ncoeff_qo_ea,stride_arr_qo_ea,DDMEA_qo)	
	Nterms = NtermsToNcoeffOS(5,0)                                                    
	call getDDesignMatOS(sinth,costh,Nterms,Ncoeff_qo_os,stride_arr_qo_os,DDMOS_qo)	
	Nterms = NtermsToNcoeffOA(5,0)                                                    
	call getDDesignMatOA(sinth,costh,Nterms,Ncoeff_qo_oa,stride_arr_qo_oa,DDMOA_qo)	

	DDMES_lt2 = DDMES_qo
	DDMEA_lt2 = DDMEA_qo
	DDMOS_lt2 = DDMOS_qo
	DDMOA_lt2 = DDMOA_qo

	!Nterms = NtermsToNcoeffES(4,0)
	!call getDDesignMatES(sinth,costh,Nterms,Ncoeff_qubic_es_lt2,stride_arr_qubic_es_lt2,DDMES_lt2)	
	!Nterms = NtermsToNcoeffEA(4,0)
	!call getDDesignMatEA(sinth,costh,Nterms,Ncoeff_qubic_ea_lt2,stride_arr_qubic_ea_lt2,DDMEA_lt2)	
	!Nterms = NtermsToNcoeffOS(4,0)
	!call getDDesignMatOS(sinth,costh,Nterms,Ncoeff_qubic_os_lt2,stride_arr_qubic_os_lt2,DDMOS_lt2)	
	!Nterms = NtermsToNcoeffOA(4,0)
	!call getDDesignMatOA(sinth,costh,Nterms,Ncoeff_qubic_oa_lt2,stride_arr_qubic_oa_lt2,DDMOA_lt2)	
endsubroutine 
!****************************************************************************************#
!****************************************************************************************#
subroutine get_pot(X,egpot,VrsEg)
	implicit none 
	real(8),intent(in):: X(0:Ncarts-1)
	real(8) :: Xtmp(0:Ncarts-1)
	real(8),intent(out)::egpot,VrsEg
	real(8):: S(0:N_int-1) ,sym_S(0:N_int-1) ,S_eq(0:N_int-1) ,S_diff(0:N_int-1) ,sym_diff(0:N_int-1) ,sym_seq(0:N_int-1),dVdsym(0:N_int-1), dVds(0:N_int-1), dsdx(0:Ncarts-1,0:N_int-1) 
    real(8):: fix_tors(0:n_fix_tors-1),sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
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
	real(8):: fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1)
	real(8):: fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
	real(8):: fc_qo_es(0:n_fcs_qo_es-1),fc_qo_ea(0:n_fcs_qo_ea-1),fc_qo_os(0:n_fcs_qo_os-1),fc_qo_oa(0:n_fcs_qo_oa-1)
	real(8):: DMES_eqm(0:fn_order_bonds_es-1),DMES_quad(0:fn_order_quad_es-1),DMES_gt5(0:fn_order_qubic_es_gt5-1),DMES_qo(0:fn_order_qo_es-1),DMES_lt2(0:fn_order_qubic_es_lt2-1) 
	real(8):: DMEA_eqm(0:fn_order_bonds_ea-1),DMEA_quad(0:fn_order_quad_ea-1),DMEA_gt5(0:fn_order_qubic_ea_gt5-1),DMEA_qo(0:fn_order_qo_ea-1),DMEA_lt2(0:fn_order_qubic_ea_lt2-1)
	real(8):: DMOS_eqm(0:fn_order_bonds_os-1),DMOS_quad(0:fn_order_quad_os-1),DMOS_gt5(0:fn_order_qubic_os_gt5-1),DMOS_qo(0:fn_order_qo_os-1),DMOS_lt2(0:fn_order_qubic_os_lt2-1)
	real(8):: DMOA_eqm(0:fn_order_bonds_oa-1),DMOA_quad(0:fn_order_quad_oa-1),DMOA_gt5(0:fn_order_qubic_oa_gt5-1),DMOA_qo(0:fn_order_qo_oa-1),DMOA_lt2(0:fn_order_qubic_oa_lt2-1)
	integer::i
	DMES_eqm=0.d0;DMES_quad=0.d0;DMES_gt5=0.d0;DMES_qo=0.d0;DMES_lt2=0.d0 
	DMEA_eqm=0.d0;DMEA_quad=0.d0;DMEA_gt5=0.d0;DMEA_qo=0.d0;DMEA_lt2=0.d0
	DMOS_eqm=0.d0;DMOS_quad=0.d0;DMOS_gt5=0.d0;DMOS_qo=0.d0;DMOS_lt2=0.d0
	DMOA_eqm=0.d0;DMOA_quad=0.d0;DMOA_gt5=0.d0;DMOA_qo=0.d0;DMOA_lt2=0.d0
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
	fc_quartic_es_ltpt1 = 0.d0;fc_quartic_ea_ltpt1 = 0.d0;fc_quartic_os_ltpt1 = 0.d0;fc_quartic_oa_ltpt1 = 0.d0;	
	fc_qo_es = 0.d0;fc_qo_ea = 0.d0;fc_qo_os = 0.d0;fc_qo_oa = 0.d0 
	S =0.d0;S_eq = 0.d0; sym_Seq = 0.d0; sym_S = 0.d0
	egpot =0.d0;VrsEg=0.d0;sinth = 0.d0;costh=0.d0 
	!Convert Bohr to Angstrom
	Xtmp =  X * bohrToAng 
	!Convert cartesian to internal
	call cart_to_internal(Xtmp,S)
	!Fixed torsions array
	fix_tors = S(tors_idx) * degreeToRadian
	!Get sin and cosine arrays
	call get_csarr(fix_tors,max_id,sinth,costh)
	call getAllDM(sinth,costh,DMES_eqm,DMES_quad,DMES_gt5,DMES_qo,DMES_lt2,DMEA_eqm,DMEA_quad,DMEA_gt5,DMEA_qo,DMEA_lt2,DMOS_eqm,DMOS_quad,DMOS_gt5,DMOS_qo,DMOS_lt2,DMOA_eqm,DMOA_quad,DMOA_gt5,DMOA_qo,DMOA_lt2)
	!Get equlibrium symmetrized internal
	call get_se(sinth,costh,DMES_eqm,DMEA_eqm,DMOS_eqm,DMOA_eqm,Sym_Seq)
	!Convert eqilibrium symmetric internal to unsymmetrized internal
	sym_Seq(11) = sym_Seq(11) +360.0
	sym_Seq(18) = sym_Seq(18) +360.0
	S_eq  =  matmul(sym_mat_transpose,sym_Seq)
	do i = Nbonds+Nang,N_int-1
            !If S(i) and S_eq(i) have opposite signs
            if (S(i)*S_eq(i)<0.d0) then
               if (S_eq(i)>0.d0) then
                  S(i) = S(i) + 360.d0
               else
                  S(i) = S(i) - 360.d0
               endif
            endif
    enddo
    s_diff = S-S_eq
    Sym_diff = matmul(sym_mat,s_diff)
	!Convert to scaled coordinates
	call scale_coord(Sym_diff)
	!compute fitted force constants(fc)
	call get_all_fc(sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,DMES_quad,DMES_gt5,DMES_qo,DMES_lt2,DMEA_quad,DMEA_gt5,DMEA_qo,DMEA_lt2,DMOS_quad,DMOS_gt5,DMOS_qo,DMOS_lt2,DMOA_quad,DMOA_gt5,DMOA_qo,DMOA_lt2) 
	
	!Compute egpotential
	call  computepot(Sym_diff,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1 ,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,DMES_eqm,egpot,VrsEg) 
	egpot = egpot * CminvToHartee
    VrsEg = VrsEg * CminvToHartee
endsubroutine 
!****************************************************************************************#
subroutine get_force(X,dVdx,egpot,flag)
	!use allvars, only: currStep
	implicit none 
	real(8),intent(in) :: X(0:Ncarts-1)
	real(8) :: Xtmp(0:Ncarts-1)
	integer,intent(in) ::flag
	real(8),intent(out) :: dVdx(0:Ncarts-1),egpot
	integer::i
	real(8):: VrsEg,S(0:N_int-1) ,sym_S(0:N_int-1) ,sym_Seq(0:N_int-1) ,S_eq(0:N_int-1) ,S_diff(0:N_int-1) ,sym_diff(0:N_int-1) ,dVdsym(0:N_int-1), dVds(0:N_int-1), dsdx(0:Ncarts-1,0:N_int-1) 
    real(8):: fix_tors(0:n_fix_tors-1),sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
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
	real(8):: fc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1),fc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1)
	real(8):: fc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1),fc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1)
	real(8):: dfc_quartic_os_ltpt1(0:n_fcs_quartic_os_ltpt1-1,0:n_fix_tors-1),dfc_quartic_oa_ltpt1(0:n_fcs_quartic_oa_ltpt1-1,0:n_fix_tors-1)
	real(8):: dfc_quartic_es_ltpt1(0:n_fcs_quartic_es_ltpt1-1,0:n_fix_tors-1),dfc_quartic_ea_ltpt1(0:n_fcs_quartic_ea_ltpt1-1,0:n_fix_tors-1)
	real(8):: fc_qo_es(0:n_fcs_qo_es-1),fc_qo_ea(0:n_fcs_qo_ea-1),fc_qo_os(0:n_fcs_qo_os-1),fc_qo_oa(0:n_fcs_qo_oa-1)
	real(8):: dfc_qo_es(0:n_fcs_qo_es-1,0:n_fix_tors-1),dfc_qo_ea(0:n_fcs_qo_ea-1,0:n_fix_tors-1),dfc_qo_os(0:n_fcs_qo_os-1,0:n_fix_tors-1),dfc_qo_oa(0:n_fcs_qo_oa-1,0:n_fix_tors-1)
	real(8):: DMES_eqm(0:fn_order_bonds_es-1),DMES_quad(0:fn_order_quad_es-1),DMES_gt5(0:fn_order_qubic_es_gt5-1),DMES_qo(0:fn_order_qo_es-1),DMES_lt2(0:fn_order_qubic_es_lt2-1) 
	real(8):: DMEA_eqm(0:fn_order_bonds_ea-1),DMEA_quad(0:fn_order_quad_ea-1),DMEA_gt5(0:fn_order_qubic_ea_gt5-1),DMEA_qo(0:fn_order_qo_ea-1),DMEA_lt2(0:fn_order_qubic_ea_lt2-1)
	real(8):: DMOS_eqm(0:fn_order_bonds_os-1),DMOS_quad(0:fn_order_quad_os-1),DMOS_gt5(0:fn_order_qubic_os_gt5-1),DMOS_qo(0:fn_order_qo_os-1),DMOS_lt2(0:fn_order_qubic_os_lt2-1)
	real(8):: DMOA_eqm(0:fn_order_bonds_oa-1),DMOA_quad(0:fn_order_quad_oa-1),DMOA_gt5(0:fn_order_qubic_oa_gt5-1),DMOA_qo(0:fn_order_qo_oa-1),DMOA_lt2(0:fn_order_qubic_oa_lt2-1)
	real(8):: DDMES_eqm(0:n_fix_tors-1,0:fn_order_bonds_es-1),DDMES_quad(0:n_fix_tors-1,0:fn_order_quad_es-1),DDMES_gt5(0:n_fix_tors-1,0:fn_order_qubic_es_gt5-1),DDMES_qo(0:n_fix_tors-1,0:fn_order_qo_es-1),DDMES_lt2(0:n_fix_tors-1,0:fn_order_qubic_es_lt2-1)
	real(8):: DDMEA_eqm(0:n_fix_tors-1,0:fn_order_bonds_ea-1),DDMEA_quad(0:n_fix_tors-1,0:fn_order_quad_ea-1),DDMEA_gt5(0:n_fix_tors-1,0:fn_order_qubic_ea_gt5-1),DDMEA_qo(0:n_fix_tors-1,0:fn_order_qo_ea-1),DDMEA_lt2(0:n_fix_tors-1,0:fn_order_qubic_ea_lt2-1)
	real(8):: DDMOS_eqm(0:n_fix_tors-1,0:fn_order_bonds_os-1),DDMOS_quad(0:n_fix_tors-1,0:fn_order_quad_os-1),DDMOS_gt5(0:n_fix_tors-1,0:fn_order_qubic_os_gt5-1),DDMOS_qo(0:n_fix_tors-1,0:fn_order_qo_os-1),DDMOS_lt2(0:n_fix_tors-1,0:fn_order_qubic_os_lt2-1)
	real(8):: DDMOA_eqm(0:n_fix_tors-1,0:fn_order_bonds_oa-1),DDMOA_quad(0:n_fix_tors-1,0:fn_order_quad_oa-1),DDMOA_gt5(0:n_fix_tors-1,0:fn_order_qubic_oa_gt5-1),DDMOA_qo(0:n_fix_tors-1,0:fn_order_qo_oa-1),DDMOA_lt2(0:n_fix_tors-1,0:fn_order_qubic_oa_lt2-1)
	DMES_eqm=0.d0;DMES_quad=0.d0;DMES_gt5=0.d0;DMES_qo=0.d0;DMES_lt2=0.d0 
	DMEA_eqm=0.d0;DMEA_quad=0.d0;DMEA_gt5=0.d0;DMEA_qo=0.d0;DMEA_lt2=0.d0
	DMOS_eqm=0.d0;DMOS_quad=0.d0;DMOS_gt5=0.d0;DMOS_qo=0.d0;DMOS_lt2=0.d0
	DMOA_eqm=0.d0;DMOA_quad=0.d0;DMOA_gt5=0.d0;DMOA_qo=0.d0;DMOA_lt2=0.d0
	DDMES_eqm=0.d0;DDMES_quad=0.d0;DDMES_gt5=0.d0;DDMES_qo=0.d0;DDMES_lt2=0.d0 
	DDMEA_eqm=0.d0;DDMEA_quad=0.d0;DDMEA_gt5=0.d0;DDMEA_qo=0.d0;DDMEA_lt2=0.d0
	DDMOS_eqm=0.d0;DDMOS_quad=0.d0;DDMOS_gt5=0.d0;DDMOS_qo=0.d0;DDMOS_lt2=0.d0
	DDMOA_eqm=0.d0;DDMOA_quad=0.d0;DDMOA_gt5=0.d0;DDMOA_qo=0.d0;DDMOA_lt2=0.d0
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
	fc_qo_es = 0.d0;fc_qo_ea = 0.d0;fc_qo_os = 0.d0;fc_qo_oa = 0.d0 
	dfc_qo_es = 0.d0;dfc_qo_ea = 0.d0;dfc_qo_os = 0.d0;dfc_qo_oa = 0.d0 
	S =0.d0;S_eq = 0.d0; sym_Seq = 0.d0; sym_S = 0.d0
	VrsEg = 0.d0;dVdx = 0.d0;sinth = 0.d0;costh=0.d0;dVdsym=0.d0;dVds=0.d0; 
	!Convert Bohr to Angstrom
	Xtmp =  X * bohrToAng 
	!Convert cartesian to internal
	call cart_to_internal(Xtmp,S)
	!Fixed torsions array
	fix_tors = S(tors_idx) * degreeToRadian
	!Get sin and cosine arrays
	call get_csarr(fix_tors,max_id,sinth,costh)
	!Get equlibrium symmetrized internal
	call getAllDM(sinth,costh,DMES_eqm,DMES_quad,DMES_gt5,DMES_qo,DMES_lt2,DMEA_eqm,DMEA_quad,DMEA_gt5,DMEA_qo,DMEA_lt2,DMOS_eqm,DMOS_quad,DMOS_gt5,DMOS_qo,DMOS_lt2,DMOA_eqm,DMOA_quad,DMOA_gt5,DMOA_qo,DMOA_lt2)
	call getAllDDM(sinth,costh,DDMES_eqm,DDMES_quad,DDMES_gt5,DDMES_qo,DDMES_lt2,DDMEA_eqm,DDMEA_quad,DDMEA_gt5,DDMEA_qo,DDMEA_lt2,DDMOS_eqm,DDMOS_quad,DDMOS_gt5,DDMOS_qo,DDMOS_lt2,DDMOA_eqm,DDMOA_quad,DDMOA_gt5,DDMOA_qo,DDMOA_lt2)
	call get_se(sinth,costh,DMES_eqm,DMEA_eqm,DMOS_eqm,DMOA_eqm,Sym_Seq)
	!Convert eqilibrium symmetric internal to unsymmetrized internal
	sym_Seq(11) = sym_Seq(11) +360.0
	sym_Seq(18) = sym_Seq(18) +360.0
	S_eq  =  matmul(sym_mat_transpose,sym_Seq)
	do i = Nbonds+Nang,N_int-1
            !If S(i) and S_eq(i) have opposite signs
            if (S(i)*S_eq(i)<0.d0) then
               if (S_eq(i)>0.d0) then
                  S(i) = S(i) + 360.d0
               else
                  S(i) = S(i) - 360.d0
               endif
            endif
    enddo
    s_diff = S-S_eq
    Sym_diff = matmul(sym_mat,s_diff)
	!Convert to scaled coordinates
	call scale_coord(Sym_diff)
	!compute fitted force constants(fc)
	call get_all_fc(sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,DMES_quad,DMES_gt5,DMES_qo,DMES_lt2,DMEA_quad,DMEA_gt5,DMEA_qo,DMEA_lt2,DMOS_quad,DMOS_gt5,DMOS_qo,DMOS_lt2,DMOA_quad,DMOA_gt5,DMOA_qo,DMOA_lt2) 

	!Compute derivative of fitted fc (dfc)
	call  get_all_dfc(sinth,costh,dfc_quad_es,dfc_quad_os,dfc_quad_ea,dfc_quad_oa,dfc_qubic_es_gt5,dfc_qubic_es_lt5,dfc_qubic_es_lt2,dfc_qubic_ea_gt5,dfc_qubic_ea_lt5,dfc_qubic_ea_lt2,dfc_qubic_os_gt5,dfc_qubic_os_lt5,dfc_qubic_os_lt2,dfc_qubic_oa_gt5,dfc_qubic_oa_lt5,dfc_qubic_oa_lt2,dfc_quartic_es_gt5,dfc_quartic_es_lt5,dfc_quartic_es_lt2,dfc_quartic_ea_gt5,dfc_quartic_ea_lt5,dfc_quartic_ea_lt2,dfc_quartic_os_gt5,dfc_quartic_os_lt5,dfc_quartic_os_lt2,dfc_quartic_oa_gt5,dfc_quartic_oa_lt5,dfc_quartic_oa_lt2,dfc_qo_es,dfc_qo_ea,dfc_qo_os,dfc_qo_oa,DDMES_quad,DDMES_gt5,DDMES_qo,DDMES_lt2,DDMEA_quad,DDMEA_gt5,DDMEA_qo,DDMEA_lt2,DDMOS_quad,DDMOS_gt5,DDMOS_qo,DDMOS_lt2,DDMOA_quad,DDMOA_gt5,DDMOA_qo,DDMOA_lt2) 
	!Compute Dvdsym(cm-1)
	call getDvds(Sym_diff,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2 ,fc_quartic_ea_ltpt1, fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,dfc_quad_es,dfc_quad_ea,dfc_quad_os,dfc_quad_oa,dfc_qubic_es_gt5,dfc_qubic_es_lt5,dfc_qubic_es_lt2,dfc_qubic_ea_gt5,dfc_qubic_ea_lt5,dfc_qubic_ea_lt2,dfc_qubic_os_gt5,dfc_qubic_os_lt5,dfc_qubic_os_lt2,dfc_qubic_oa_gt5,dfc_qubic_oa_lt5,dfc_qubic_oa_lt2,dfc_quartic_es_gt5,dfc_quartic_es_lt5,dfc_quartic_es_lt2,dfc_quartic_es_ltpt1,dfc_quartic_ea_gt5,dfc_quartic_ea_lt5,dfc_quartic_ea_lt2,dfc_quartic_ea_ltpt1, dfc_quartic_os_gt5,dfc_quartic_os_lt5,dfc_quartic_os_ltpt1,dfc_quartic_os_lt2,dfc_quartic_oa_gt5,dfc_quartic_oa_lt5,dfc_quartic_oa_lt2,dfc_quartic_oa_ltpt1,dfc_qo_es,dfc_qo_ea,dfc_qo_os,dfc_qo_oa,DDMES_eqm,DDMEA_eqm,DDMOS_eqm,DDMOA_eqm,dVdsym) 
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
	!do i = 0,Ncarts-1
    !	if (abs(dVdx(i))>1.d0) then
    !    	print*,'Coord=',Xtmp
    !    	print*,'force=',dVdx
    !    	stop
    !	endif
    !enddo
	if (flag==4) then
		call  computepot(Sym_diff,sinth,costh,fc_quad_es,fc_quad_os,fc_quad_ea,fc_quad_oa,fc_qubic_es_gt5,fc_qubic_es_lt5,fc_qubic_es_lt2,fc_qubic_ea_gt5,fc_qubic_ea_lt5,fc_qubic_ea_lt2,fc_qubic_os_gt5,fc_qubic_os_lt5,fc_qubic_os_lt2,fc_qubic_oa_gt5,fc_qubic_oa_lt5,fc_qubic_oa_lt2,fc_quartic_es_gt5,fc_quartic_es_lt5,fc_quartic_es_lt2,fc_quartic_es_ltpt1,fc_quartic_ea_gt5,fc_quartic_ea_lt5,fc_quartic_ea_lt2,fc_quartic_ea_ltpt1 ,fc_quartic_os_gt5,fc_quartic_os_lt5,fc_quartic_os_lt2,fc_quartic_os_ltpt1,fc_quartic_oa_gt5,fc_quartic_oa_lt5,fc_quartic_oa_lt2,fc_quartic_oa_ltpt1,fc_qo_es,fc_qo_ea,fc_qo_os,fc_qo_oa,DMES_eqm,egpot,VrsEg) 
	endif
endsubroutine 
!****************************************************************************************#
subroutine numerical_dVdx(X,dVdx_num)
implicit none 
real(8),intent(inout)::X(0:Ncarts-1)
real(8),intent(out)::dVdx_num(0:Ncarts-1)
integer::i,istep,flag=1
real(8)::deltax,Xtmp(0:Ncarts-1),egpot_arr(0:4),egpot=0.d0,VrsEg,dum(0:Ncarts-1)
Xtmp = 0.d0
egpot = 0.d0
egpot_arr = 0.d0;dum = 0.d0
dVdx_num = 0.d0
!X(:)= X(:) * BohrToAng
deltax = AngToBohr * dx
!loop over the cartesian coordinates
do i = 0,Ncarts-1
	Xtmp(:) = X(:)
	egpot = 0.d0
	!displace five points and save the egpotential at each point
	do istep = 0,4
		Xtmp(i) = X(i)+ deltax*dble(istep-2)
		call get_pot(Xtmp,egpot,VrsEg)
		egpot_arr(istep) = egpot	
	enddo
	dVdx_num(i) = dot_product(five_pt_1st_deriv,egpot_arr) / deltax 
enddo
endsubroutine
end module 
