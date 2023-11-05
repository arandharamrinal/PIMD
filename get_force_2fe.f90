module fe_pot_force
use POTVARS
use constants 
use allvars, only : cartmass, total_mass_inv
implicit none 
contains
!****************************************************************************************#
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
                    ang = dacos(v_dot) * radian_to_degree
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
                    ang = dacos(v_dot) * radian_to_degree
                    ang = mod(ang,360.d0)
                    int_arr(Nbonds + a_ctr) = ang
                    !print*,int_arr(b_ctr);read(*,*)
                    a_ctr =a_ctr +1
                    a1(0) = v1(1)*v2(2) - v2(1)*v1(2)
                    a1(1) = v2(0)*v1(2) - v2(2)*v1(0)
                    a1(2) = v1(0)*v2(1) - v1(1)*v2(0)
                    a2(0) = v2(1)*v3(2) - v3(1)*v2(2)
                    a2(1) = v3(0)*v2(2) - v3(2)*v2(0)
                    a2(2) = v2(0)*v3(1) - v2(1)*v3(0)
                    v1_dot = norm2(v2) * dot_product(v1,a2)
                    v2_dot = dot_product(a1,a2)
                    dih = datan2(v1_dot,v2_dot) * radian_to_degree
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
!****************************************************************************************#
subroutine int_to_cart(S,cart_coord)
    implicit none
    real(8)::b2,a2,d2,r2,z,y
    real(8),allocatable::B(:),A(:),D(:),q(:),u1(:),u2(:)
    integer::i,a_ctr,b_ctr,d_ctr,Ndih,atm1,atm2,atm3,atm4
    real(8),intent(in),dimension(0:N_int-1)::S
    real(8),allocatable::n_v(:)
    real(8),intent(out),dimension(0:Ncarts-1)::cart_coord
    Ndih = N_int-Nbonds-Nang
    allocate(B(0:Nbonds-1),A(0:Nang-1),D(0:Ndih-1),q(0:2),n_v(0:2))
    allocate(u1(0:2),u2(0:2))
    cart_coord=0.0
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    B = S(0:Nbonds-1)
    A = S(Nbonds:Nbonds+Nang-1)
    D = S(Nbonds+Nang:N_int-1)
    do i = 0,Natoms-1
        if (atom_ind(i,1)==0) then
            !Place the first atom at the center.
            cart_coord(3*i:3*i+2) = 0.0

        else if (atom_ind(i,2)==0) then

            !Place the first atom at the center.
            cart_coord(3*i:3*i+1) = 0.0
            cart_coord(3*i+2)  = B(b_ctr)
            b_ctr = b_ctr + 1

        else if (atom_ind(i,3)==0) then

            !Place the first atom at the center.
            b2  = B(b_ctr)
            a2  = A(a_ctr)
            z = b2 * dcos( a2 * pi / 180.d0)
            y = dsqrt(b2*b2 - z*z)

            !Place the 2nd atom in the yz plane
            cart_coord(3*i) = 0.d0
            cart_coord(3*i+1) = y
            cart_coord(3*i+2) = z
            b_ctr =  b_ctr +1
            a_ctr =  a_ctr +1

        else
                atm4 =int( atom_ind(i,0))-1
                atm3 =int( atom_ind(i,1))-1
                atm2 =int( atom_ind(i,2))-1
                atm1 =int( atom_ind(i,3))-1

                !Get the bond length and angle involving atm4.
                r2   = B(b_ctr)
                a2   = A(a_ctr)

                !Get dihedral
                d2   = D(d_ctr)
                b_ctr = b_ctr +1
                a_ctr = a_ctr +1
                d_ctr = d_ctr +1

                !Compute the bond vectors
                u1   = cart_coord(3*atm2:3*atm2+2)-cart_coord(3*atm1:3*atm1+2)
                u2   = cart_coord(3*atm3:3*atm3+2)-cart_coord(3*atm2:3*atm2+2)
                !Get the cartesian coordinate with the 3rd atom as the origin.
                call internal_calculator(u1,u2,r2,a2,d2,q)

                !Compute the absolute cartesian position
                q   = q + cart_coord(3*atm3:3*atm3+2)
                cart_coord(3*i:3*i+2) = q(:)
        endif
    enddo
    deallocate(B,A,D,q,n_v,u1,u2)
end subroutine int_to_cart
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
        !write(*,*)ii,i,j,k,l;read(*,*)
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
!compute force constants for a given value of phi1 and phi2
!****************************************************************************************#
subroutine get_all_vsav(S,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even,Vsav)
        implicit none 
        real(8),intent(in):: S(0:N_int-1)
        real(8),intent(in):: fc_quad_even(0:n_fcs_quad_even-1),fc_quad_odd(0:n_fcs_quad_odd-1),fc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1),fc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1)
        real(8),intent(in):: fc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1),fc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1),fc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1),fc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1) 
        real(8),intent(in):: fc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1),fc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1),fc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1) 
        real(8),intent(in):: fc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1),fc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1),fc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1) 
        real(8),intent(in):: fc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1),fc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1) 
        real(8),intent(in):: fc_qo_ang_even(0:n_fcs_qo_ang_even-1),fc_qo_even(0:n_fcs_qo_even-1),fc_qo_odd(0:n_fcs_qo_odd-1) 
        real(8),intent(out) :: Vsav
        real(8) :: Vsav_quad_even,Vsav_quad_odd,Vsav_qubic_gt5_even,Vsav_qubic_lt5_even,Vsav_qubic_lt2_even,Vsav_qubic_gt5_odd
        real(8) :: Vsav_qubic_lt5_odd,Vsav_qubic_lt2_odd,Vsav_quartic_gt5_even,Vsav_quartic_lt5_even,Vsav_quartic_lt2_even
        real(8) ::  Vsav_quartic_gt5_odd,Vsav_quartic_lt5_odd,Vsav_quartic_lt2_odd
        real(8) :: Vsav_qo_even,Vsav_qo_ang_even,Vsav_qo_odd
        real(8) :: Vsav_quartic_ltpt1_even,Vsav_quartic_ltpt1_odd
        Vsav = 0.d0
        !quadratic 
        call getVsav(S,fc_idx_quad_even,fc_quad_even,n_fcs_quad_even,Vsav_quad_even)
        call getVsav(S,fc_idx_quad_odd,fc_quad_odd,n_fcs_quad_odd,Vsav_quad_odd)
        !Cubic even
        call getVsav(S,fc_idx_qubic_even_gt5,fc_qubic_even_gt5,n_fcs_qubic_even_gt5,Vsav_qubic_gt5_even)
        call getVsav(S,fc_idx_qubic_even_lt5,fc_qubic_even_lt5,n_fcs_qubic_even_lt5,Vsav_qubic_lt5_even)
        call getVsav(S,fc_idx_qubic_even_lt2,fc_qubic_even_lt2,n_fcs_qubic_even_lt2,Vsav_qubic_lt2_even)
        !Cubic Odd
        call getVsav(S,fc_idx_qubic_odd_gt5,fc_qubic_odd_gt5,n_fcs_qubic_odd_gt5,Vsav_qubic_gt5_odd)
        call getVsav(S,fc_idx_qubic_odd_lt5,fc_qubic_odd_lt5,n_fcs_qubic_odd_lt5,Vsav_qubic_lt5_odd)
        call getVsav(S,fc_idx_qubic_odd_lt2,fc_qubic_odd_lt2,n_fcs_qubic_odd_lt2,Vsav_qubic_lt2_odd)
        !Quartic even   
        call getVsav(S,fc_idx_quartic_even_gt5,fc_quartic_even_gt5,n_fcs_quartic_even_gt5,Vsav_quartic_gt5_even)
        call getVsav(S,fc_idx_quartic_even_lt5,fc_quartic_even_lt5,n_fcs_quartic_even_lt5,Vsav_quartic_lt5_even)
        call getVsav(S,fc_idx_quartic_even_lt2,fc_quartic_even_lt2,n_fcs_quartic_even_lt2,Vsav_quartic_lt2_even)
        !Quartic odd    
        call getVsav(S,fc_idx_quartic_odd_gt5,fc_quartic_odd_gt5,n_fcs_quartic_odd_gt5,Vsav_quartic_gt5_odd)
        !print*,fc_quartic_odd_gt5;read(*,*)
        call getVsav(S,fc_idx_quartic_odd_lt5,fc_quartic_odd_lt5,n_fcs_quartic_odd_lt5,Vsav_quartic_lt5_odd)
        call getVsav(S,fc_idx_quartic_odd_lt2,fc_quartic_odd_lt2,n_fcs_quartic_odd_lt2,Vsav_quartic_lt2_odd)
        call getVsav(S,fc_idx_quartic_even_ltpt1,fc_quartic_even_ltpt1,n_fcs_quartic_even_ltpt1,Vsav_quartic_ltpt1_even)
        call getVsav(S,fc_idx_quartic_odd_ltpt1,fc_quartic_odd_ltpt1,n_fcs_quartic_odd_ltpt1,Vsav_quartic_ltpt1_odd)
        call getVsav_qo(S,fc_idx_qo_even,fc_qo_even,n_fcs_qo_even,Vsav_qo_even)
        call getVsav_qo(S,fc_idx_qo_ang_even,fc_qo_ang_even,n_fcs_qo_ang_even,Vsav_qo_ang_even)
        call getVsav_qo(S,fc_idx_qo_odd,fc_qo_odd,n_fcs_qo_odd,Vsav_qo_odd)

        !print*,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2
        !Sum up all 
        Vsav = Vsav_quad_even + Vsav_quad_odd 
        !print*,'quad=',Vsav_quad_even , Vsav_quad_odd 
        Vsav = Vsav    + Vsav_qubic_gt5_even + Vsav_qubic_lt5_even + Vsav_qubic_lt2_even + Vsav_qubic_gt5_odd 
        !print*,'qubic1= ',Vsav_qubic_gt5_even , Vsav_qubic_lt5_even , Vsav_qubic_lt2_even , Vsav_qubic_gt5_odd 
        Vsav = Vsav + Vsav_qubic_lt5_odd + Vsav_qubic_lt2_odd + Vsav_quartic_gt5_even + Vsav_quartic_lt5_even + Vsav_quartic_lt2_even
        !print*, 'qubic2=',Vsav_qubic_lt5_odd , Vsav_qubic_lt2_odd  
        !print*, 'quartic1=',Vsav_quartic_gt5_even , Vsav_quartic_lt5_even , Vsav_quartic_lt2_even
        Vsav = Vsav + Vsav_quartic_gt5_odd + Vsav_quartic_lt5_odd + Vsav_quartic_lt2_odd 
        !print*, 'quartic1=',Vsav_quartic_gt5_odd , Vsav_quartic_lt5_odd , Vsav_quartic_lt2_odd 
        Vsav = Vsav + Vsav_quartic_ltpt1_even + Vsav_quartic_ltpt1_odd
        !print*,'quartic2=',Vsav_quartic_ltpt1_even , Vsav_quartic_ltpt1_odd
        Vsav = Vsav + Vsav_qo_even + Vsav_qo_ang_even + Vsav_qo_odd
        !print*,Vsav_qubic_gt5_even + Vsav_qubic_lt5_even + Vsav_qubic_lt2_even + Vsav_qubic_gt5_odd+Vsav_qubic_lt5_odd + Vsav_qubic_lt2_odd + Vsav_quartic_gt5_even + Vsav_quartic_lt5_even + Vsav_quartic_lt2_even+Vsav_quartic_gt5_odd + Vsav_quartic_lt5_odd + Vsav_quartic_lt2_odd+Vsav_qo_even + Vsav_qo_ang_even + Vsav_qo_odd
        !print*, 'qo=',Vsav_qo_even , Vsav_qo_ang_even , Vsav_qo_odd
        !print*,"Vsav=",Vsav;read(*,*)
endsubroutine  

!****************************************************************************************#
subroutine  compute_pot(S,sinth,costh,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even,pot,Vt_out)
    implicit none 
    real(8),intent(in)::S(0:N_int-1),sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
    real(8),intent(in):: fc_quad_even(0:n_fcs_quad_even-1),fc_quad_odd(0:n_fcs_quad_odd-1),fc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1),fc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1)
    real(8),intent(in):: fc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1),fc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1),fc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1),fc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1) 
    real(8),intent(in):: fc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1),fc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1),fc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1) 
    real(8),intent(in):: fc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1),fc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1),fc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1) 
    real(8),intent(in):: fc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1),fc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1) 
    real(8),intent(in):: fc_qo_ang_even(0:n_fcs_qo_ang_even-1),fc_qo_even(0:n_fcs_qo_even-1),fc_qo_odd(0:n_fcs_qo_odd-1) 
    real(8),intent(out)::pot,Vt_out
    real(8)::Vt=0.d0,Vsav=0.d0
    pot = 0.d0;Vt_out = 0.d0
    call getVt(sinth,costh,Vt_out)
    call get_all_vsav(S,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even,Vsav)
    !print*,Vsav
    pot  = Vt_out + Vsav 
    !print*,'pot=',pot,Vt_out,Vsav;read(*,*)
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
         do itors = 0,n_fix_tors-1
            dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + dfc(ii,itors) * S(i) * S(j) + fc(ii) * dS_dtors(i,itors) * S(j) + fc(ii) * dS_dtors(j,itors) * S(i)
        enddo
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
         do itors = 0,n_fix_tors-1
            dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + dfc(ii,itors) * S(i) * S(j) * S(k) + fc(ii) * dS_dtors(i,itors) * S(j) * S(k)  + fc(ii) * dS_dtors(j,itors) * S(i) * S(k) 
            dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(ii) * dS_dtors(k,itors) * S(i) * S(j)
        enddo
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
    !print*,n_params
    do ii =0,n_params-1
        i    = fc_idx(ii,0);j = fc_idx(ii,1);k = fc_idx(ii,2);l = fc_idx(ii,3);
        dVds(i) = dVds(i) + fc(ii) * S(j) * S(k) * S(l)
        dVds(j) = dVds(j) + fc(ii) * S(i) * S(k) * S(l)
        dVds(k) = dVds(k) + fc(ii) * S(i) * S(j) * S(l)
        dVds(l) = dVds(l) + fc(ii) * S(i) * S(j) * S(k)
         do itors = 0,n_fix_tors-1
            dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + dfc(ii,itors) * S(i) * S(j) * S(k) * S(l)
            dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(ii) * S(i) * dS_dtors(j,itors) * S(k) * S(l) + fc(ii) * S(j) * dS_dtors(i,itors) * S(k) * S(l) 
            dVds(tors_idx(itors)) = dVds(tors_idx(itors)) + fc(ii) * dS_dtors(k,itors) * S(i) * S(j) * S(l) + fc(ii) * dS_dtors(l,itors) * S(i) * S(j) * S(k)
        enddo
        !print*,'1 =',dVds(i),dVds(j),dVds(k),dVds(l),fc(ii)
        !print*,'2 =',dfc(ii,0),dS_dtors(j,0),dVds(tors_idx(0))
        !print*,'3 =',dfc(ii,1),dS_dtors(j,1),dVds(tors_idx(1))
    enddo
endsubroutine 
!****************************************************************************************#
subroutine  dvdsQO(fc_idx,fc,dS_dtors,dfc,S,icoord,n_params,dVds)
    implicit none 
    integer,intent(in)::icoord,n_params
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
subroutine getDvds(S,sinth,costh,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even,dfc_quad_even,dfc_quad_odd,dfc_qubic_even_gt5,dfc_qubic_even_lt5,dfc_qubic_even_lt2,dfc_qubic_odd_gt5,dfc_qubic_odd_lt5,dfc_qubic_odd_lt2,dfc_quartic_even_gt5,dfc_quartic_even_lt5,dfc_quartic_even_lt2,dfc_quartic_odd_gt5,dfc_quartic_odd_lt5,dfc_quartic_odd_lt2,dfc_quartic_even_ltpt1,dfc_quartic_odd_ltpt1,dfc_qo_even,dfc_qo_odd,dfc_qo_ang_even,dVds)
    implicit none 
    integer :: icoord
    real(8) :: dfcds(0:n_fix_tors-1),dVt(0:n_fix_tors-1),dS_dtors(0:N_int-1,0:n_fix_tors-1)
    integer :: t_id,i,idxs,b_ctr,a_ctr,d_ctr
    real(8) ::dVdstmp
    real(8),intent(in) :: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
    real(8),intent(out)::dVds(0:N_int-1)
    real(8),intent(in) :: S(0:N_int-1),fc_quad_even(0:n_fcs_quad_even-1),fc_quad_odd(0:n_fcs_quad_odd-1),fc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1),fc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1)
    real(8),intent(in) :: fc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1),fc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1),fc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1),fc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1) 
    real(8),intent(in) :: fc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1),fc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1),fc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1) 
    real(8),intent(in) :: fc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1),fc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1),fc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1) 
    real(8),intent(in) :: fc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1),fc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1) 
    real(8),intent(in) :: fc_qo_ang_even(0:n_fcs_qo_ang_even-1),fc_qo_even(0:n_fcs_qo_even-1),fc_qo_odd(0:n_fcs_qo_odd-1) 
    real(8),intent(in) :: dfc_quad_even(0:n_fcs_quad_even-1,0:n_fix_tors-1),dfc_quad_odd(0:n_fcs_quad_odd-1,0:n_fix_tors-1),dfc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1,0:n_fix_tors-1),dfc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1,0:n_fix_tors-1)
    real(8),intent(in) :: dfc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1,0:n_fix_tors-1),dfc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1,0:n_fix_tors-1),dfc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1,0:n_fix_tors-1),dfc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1,0:n_fix_tors-1) 
    real(8),intent(in) :: dfc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1,0:n_fix_tors-1),dfc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1,0:n_fix_tors-1),dfc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1,0:n_fix_tors-1) 
    real(8),intent(in) :: dfc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1,0:n_fix_tors-1),dfc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1,0:n_fix_tors-1),dfc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1,0:n_fix_tors-1) 
    real(8),intent(in) :: dfc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1,0:n_fix_tors-1),dfc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1,0:n_fix_tors-1) 
    real(8),intent(in) :: dfc_qo_ang_even(0:n_fcs_qo_ang_even-1,0:n_fix_tors-1),dfc_qo_even(0:n_fcs_qo_even-1,0:n_fix_tors-1),dfc_qo_odd(0:n_fcs_qo_odd-1,0:n_fix_tors-1) 
    !""" i'th element of dVds contains the derivatives of potential w.r.t to the i'th internal"""
    !Compute the derivatives of the i'th internal w.r.t the fixed coordinates
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    dS_dtors = 0.d0
    dVds = 0.d0
    do i = 0,size(odd_fn_idx)-1
        idxs = odd_fn_idx(i)
        if (idxs<Nbonds) then
            dfcds = 0.d0
            call dtrig_series_odd(sinth,costh,Ncoeff_bonds_odd,fitted_bonds_coeff_odd(:,b_ctr),stride_arr_b_odd,fn_order_bonds_odd,dfcds)
            dS_dtors(idxs,:) = -dfcds(:) 
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            dfcds = 0.d0
            call dtrig_series_odd(sinth,costh,Ncoeff_angs_odd,fitted_angs_coeff_odd(:,a_ctr),stride_arr_a_odd,fn_order_angs_odd,dfcds)
            dS_dtors(idxs,:) = -dfcds(:) 
            a_ctr = a_ctr + 1 
        else
            dfcds = 0.d0
            call dtrig_series_odd(sinth,costh,Ncoeff_dihs_odd,fitted_dihs_coeff_odd(:,d_ctr),stride_arr_d_odd,fn_order_dihs_odd,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:)
            d_ctr = d_ctr + 1 
        endif
    enddo
    b_ctr = 0
    a_ctr = 0
    d_ctr = 0
    do i = 0,size(even_fn_idx)-1
        idxs = even_fn_idx(i)
        if (idxs<Nbonds) then
            dfcds = 0.d0
            call dtrig_series_even(sinth,costh,Ncoeff_bonds_even,fitted_bonds_coeff_even(:,b_ctr),stride_arr_b_even,fn_order_bonds_even,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:)
            b_ctr = b_ctr + 1 
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            dfcds = 0.d0
            call dtrig_series_even(sinth,costh,Ncoeff_angs_even,fitted_angs_coeff_even(:,a_ctr),stride_arr_a_even,fn_order_angs_even,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:) 
            a_ctr = a_ctr + 1 
        else
            dfcds = 0.d0
            call dtrig_series_even(sinth,costh,Ncoeff_dihs_even,fitted_dihs_coeff_even(:,d_ctr),stride_arr_d_even,fn_order_dihs_even,dfcds)
            dS_dtors(idxs,:) =  -dfcds(:) 
            d_ctr = d_ctr + 1 
        endif
    enddo
    !Convert dS_dtors unit to 1/radian
    do i = 0,N_int-1 
        dS_dtors(i,:) = dS_dtors(i,:) * dim_scal_factor_inv(i)
        if (i>= Nbonds) then
            dS_dtors(i,:) = dS_dtors(i,:) * degree_to_radian
        endif
    enddo
    !Compute the derivative of the potential along the minimum path(zero'th order term in the expansion) w.r.t the fixed coordinates
    call getDVt(sinth,costh,dVt)
    dVds(tors_idx) =  dVds(tors_idx) +  dVt(:)
    !Compute derivative of the 2nd order, 3rd order and fourth order terms w.r.t to all internal coordinates
    !For quadratic even
    call dvdsQuad(fc_idx_quad_even,fc_quad_even,dS_dtors,dfc_quad_even,S,icoord,n_fcs_quad_even,dVds)
    !For quadratic Odd
    call dvdsQuad(fc_idx_quad_odd,fc_quad_odd,dS_dtors,dfc_quad_odd,S,icoord,n_fcs_quad_odd,dVds)
    !Cubic:
    !Cubic one body even terms
    call dvdsQubic(fc_idx_qubic_even_gt5,fc_qubic_even_gt5,dS_dtors,dfc_qubic_even_gt5,S,icoord,n_fcs_qubic_even_gt5,dVds)
    !Cubic two body even terms
    call dvdsQubic(fc_idx_qubic_even_lt5,fc_qubic_even_lt5,dS_dtors,dfc_qubic_even_lt5,S,icoord,n_fcs_qubic_even_lt5,dVds)
    !Cubic three body even terms
    call dvdsQubic(fc_idx_qubic_even_lt2,fc_qubic_even_lt2,dS_dtors,dfc_qubic_even_lt2,S,icoord,n_fcs_qubic_even_lt2,dVds)
    !Cubic one body Odd terms
    call dvdsQubic(fc_idx_qubic_odd_gt5,fc_qubic_odd_gt5,dS_dtors,dfc_qubic_odd_gt5,S,icoord,n_fcs_qubic_odd_gt5,dVds)
    !Cubic two body Odd terms
    call dvdsQubic(fc_idx_qubic_odd_lt5,fc_qubic_odd_lt5,dS_dtors,dfc_qubic_odd_lt5,S,icoord,n_fcs_qubic_odd_lt5,dVds)
    !Cubic three body Odd terms
    call dvdsQubic(fc_idx_qubic_odd_lt2,fc_qubic_odd_lt2,dS_dtors,dfc_qubic_odd_lt2,S,icoord,n_fcs_qubic_odd_lt2,dVds)
    !Quartic:    
    !Quartic one body even terms
    call dvdsQuartic(fc_idx_quartic_even_gt5,fc_quartic_even_gt5,dS_dtors,dfc_quartic_even_gt5,S,icoord,n_fcs_quartic_even_gt5,dVds)
    !Quartic two body even terms
    call dvdsQuartic(fc_idx_quartic_even_lt5,fc_quartic_even_lt5,dS_dtors,dfc_quartic_even_lt5,S,icoord,n_fcs_quartic_even_lt5,dVds)
    !Quartic three body even terms
    call dvdsQuartic(fc_idx_quartic_even_lt2,fc_quartic_even_lt2,dS_dtors,dfc_quartic_even_lt2,S,icoord,n_fcs_quartic_even_lt2,dVds)
    !Quartic two body Odd terms
    call dvdsQuartic(fc_idx_quartic_odd_gt5,fc_quartic_odd_gt5,dS_dtors,dfc_quartic_odd_gt5,S,icoord,n_fcs_quartic_odd_gt5,dVds)
    call dvdsQuartic(fc_idx_quartic_odd_lt5,fc_quartic_odd_lt5,dS_dtors,dfc_quartic_odd_lt5,S,icoord,n_fcs_quartic_odd_lt5,dVds)
    !Quartic three body Odd terms
    call dvdsQuartic(fc_idx_quartic_odd_lt2,fc_quartic_odd_lt2,dS_dtors,dfc_quartic_odd_lt2,S,icoord,n_fcs_quartic_odd_lt2,dVds)
    call dvdsQuartic(fc_idx_quartic_even_ltpt1,fc_quartic_even_ltpt1,dS_dtors,dfc_quartic_even_ltpt1,S,icoord,n_fcs_quartic_even_ltpt1,dVds)
    call dvdsQuartic(fc_idx_quartic_odd_ltpt1,fc_quartic_odd_ltpt1,dS_dtors,dfc_quartic_odd_ltpt1,S,icoord,n_fcs_quartic_odd_ltpt1,dVds)
    call dvdsQO(fc_idx_qo_even,fc_qo_even,dS_dtors,dfc_qo_even,S,icoord,n_fcs_qo_even,dVds)
    call dvdsQO(fc_idx_qo_odd,fc_qo_odd,dS_dtors,dfc_qo_odd,S,icoord,n_fcs_qo_odd,dVds)
    call dvdsQO(fc_idx_qo_ang_even,fc_qo_ang_even,dS_dtors,dfc_qo_ang_even,S,icoord,n_fcs_qo_ang_even,dVds)
    !print*,dVds
endsubroutine
!****************************************************************************************#
subroutine get_cc_deriv(i,j,sinth,costh,coeff,dfc_dtors)
    implicit none 
    integer,intent(in):: i,j
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:1)
    dfc_dtors = 0.d0
    dfc_dtors(0) =  - dble(i) * coeff * sinth(0,i) * costh(1,j)
    dfc_dtors(1) =  - dble(j) * coeff * costh(0,i) * sinth(1,j)
end subroutine
!****************************************************************************************#
subroutine  get_ss_deriv(i,j,sinth,costh,coeff,dfc_dtors)
    implicit none 
    integer,intent(in):: i,j
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:1) 
    dfc_dtors = 0.d0
    dfc_dtors(0) = + dble(i) * coeff * costh(0,i) * sinth(1,j)
    dfc_dtors(1) = + dble(j) * coeff * sinth(0,i) * costh(1,j)
    !print("sums=",th1,th2,coeff,sum_s[0],sum_s[1])
end subroutine
!****************************************************************************************#
subroutine get_cs_deriv(i,j,sinth,costh,coeff,dfc_dtors)
    implicit none 
    integer,intent(in):: i,j
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:1)
    dfc_dtors = 0.d0
    dfc_dtors(0) =  - dble(i) * coeff * sinth(0,i) * sinth(1,j)
    dfc_dtors(1) =  + dble(j) * coeff * costh(0,i) * costh(1,j)
end subroutine
!****************************************************************************************#
subroutine get_sc_deriv(i,j,sinth,costh,coeff,dfc_dtors)
    implicit none 
    integer,intent(in):: i,j
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff
    real(8),intent(out)::dfc_dtors(0:1) 
    dfc_dtors = 0.d0
    dfc_dtors(0) = + dble(i) * coeff * costh(0,i) * costh(1,j)
    dfc_dtors(1) = - dble(j) * coeff * sinth(0,i) * sinth(1,j)
    !print("sums=",th1,th2,coeff,sum_s[0],sum_s[1])
end subroutine
!****************************************************************************************#
subroutine dtrig_series_even(sinth,costh,Ncoeff,coeff_arr,stride_arr_even,fn_order_even,dfc_dtors)
    !"""
    !sum_s   : An array of length equal to the no. of fixed internals(len(fix_tors)). i'th element contains the 
    !theta   : An array of length equal to the no. of fixed internals. i'th element contains the value of the 
    !          i'th fixed coordinate.
    !idx     : An array of length equal to the no. of fixed internals. i'th element contains the value of the 
    !"""
    implicit none 
    integer,intent(in):: fn_order_even
    integer,intent(in)::Ncoeff(0:1),stride_arr_even(0:fn_order_even-1,0:1)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff_arr(0:fn_order_even-1)
    real(8),intent(out)::dfc_dtors(0:n_fix_tors-1)
    real(8)::dfdt(0:n_fix_tors-1)
    integer:: jj,i,j

    dfc_dtors = 0.d0 
    dfdt = 0.d0 
    do jj = 0,fn_order_even-1
        i = stride_arr_even(jj,0);j=stride_arr_even(jj,1)
        if (jj < Ncoeff(0)) then
            call get_cc_deriv(i,j,sinth,costh,coeff_arr(jj),dfdt) 
            dfc_dtors(:) = dfc_dtors(:) + dfdt(:)
        else
            call get_ss_deriv(i,j,sinth,costh,coeff_arr(jj),dfdt) 
            dfc_dtors(:) = dfc_dtors(:) + dfdt(:)
        endif
    enddo
end subroutine
!****************************************************************************************#
subroutine trig_series_even(sinth,costh,Ncoeff_even,coeff_arr,stride_arr_even,fn_order_even,fit_val)
    implicit none 
    integer,intent(in)::fn_order_even
    integer,intent(in)::Ncoeff_even(0:1),stride_arr_even(0:fn_order_even-1,0:1)
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),coeff_arr(0:fn_order_even-1)
    real(8),intent(out)::fit_val
    integer:: jj,i1,j1
    fit_val = 0.d0
    do jj = 0,fn_order_even-1
        i1 = stride_arr_even(jj,0); j1 = stride_arr_even(jj,1)
        if (jj<Ncoeff_even(0)) then 
            fit_val = fit_val + coeff_arr(jj) * costh(0,i1) * costh(1,j1)
        else
            fit_val = fit_val + coeff_arr(jj) * sinth(0,i1) * sinth(1,j1) 
        endif
    enddo
end subroutine
!****************************************************************************************#
subroutine  getVt(sinth,costh,Vt)
    implicit none 
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
    real(8),intent(out)::Vt
    Vt = 0.0
    call trig_series_even(sinth,costh,Ncoeff_pot,Vt_coeff,stride_arr_pot,fn_order_pot,Vt)
endsubroutine 
!****************************************************************************************#
subroutine  getDVt(sinth,costh,dVt)
    implicit none
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
    real(8),intent(out)::dVt(0:n_fix_tors-1)
    real(8) :: dfc_dtors(0:n_fix_tors-1)
    integer:: jj,i,j
    dVt = 0.d0;dfc_dtors = 0.d0;
    do jj = 0,fn_order_pot-1
        dfc_dtors = 0.d0;
        i = stride_arr_pot(jj,0);j=stride_arr_pot(jj,1)
        if (jj < Ncoeff_pot(0)) then
            call  get_cc_deriv(i,j,sinth,costh,Vt_coeff(jj),dfc_dtors)
            dVt = dVt + dfc_dtors
        else
            call  get_ss_deriv(i,j,sinth,costh,Vt_coeff(jj),dfc_dtors)
            dVt = dVt + dfc_dtors
        endif
    enddo
endsubroutine 
!Compute Vsav
!****************************************************************************************#
!given q6 and order of expansion compute the force constants:
subroutine get_fc(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,odd,fc)
    implicit none 
    integer::ii
    integer,intent(in)::Ncoeff(0:1),fn_order,n_params,odd,stride_arr(0:fn_order-1,0:1) 
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::fc(0:n_params-1) 
    fc = 0.d0
    !Even coefficients
    if (odd==1) then
        do ii=0,n_params-1
            call trig_series_odd(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,fc(ii)) 
         !print*,'odd=',fc(ii)
         !read(*,*)
        enddo 
    else
        do ii=0,n_params-1
            call trig_series_even(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,fc(ii)) 
         !print*,'even=',ii,fc(ii)
         !read(*,*)
        enddo 
    endif
endsubroutine 

!***************************************************************************************#
!given  function order compute the force constants:
subroutine get_dfc(sinth,costh,Ncoeff,fitted_fc_coeff,stride_arr,fn_order,n_params,odd,dfc)
    implicit none 
    integer::i1,j1,ii,jj
    real(8):: cs_deriv(0:1),sc_deriv(0:1),cc_deriv(0:1),ss_deriv(0:1),dfc_dtors(0:1)
    integer,intent(in):: Ncoeff(0:1),fn_order,n_params,odd,stride_arr(0:fn_order-1,0:1) 
    real(8),intent(in):: sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)),fitted_fc_coeff(0:fn_order-1,0:n_params-1)
    real(8),intent(out)::dfc(0:n_params-1,0:1)
    cs_deriv=0.d0; sc_deriv=0.d0;cc_deriv=0.d0;ss_deriv=0.d0;dfc=0.d0;dfc_dtors=0.d0
    !print*,n_params
    if (odd==1) then
       do ii = 0,n_params-1
           call dtrig_series_odd(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,dfc_dtors)
           dfc(ii,:) = dfc_dtors(:)
       enddo
    else
       do ii = 0,n_params-1
           call dtrig_series_even(sinth,costh,Ncoeff,fitted_fc_coeff(:,ii),stride_arr,fn_order,dfc_dtors)
           dfc(ii,:) = dfc_dtors(:)
        enddo
    endif
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
    do i = 0,size(odd_fn_idx)-1
        idxs = odd_fn_idx(i)
        !print(idxs,fitted_bonds_coeff_odd[i])
        if (idxs<Nbonds) then
            call trig_series_odd(sinth,costh,Ncoeff_bonds_odd,fitted_bonds_coeff_odd(:,i),stride_arr_b_odd,fn_order_bonds_odd,S_e(idxs))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            call trig_series_odd(sinth,costh,Ncoeff_angs_odd,fitted_angs_coeff_odd(:,i-b_ctr),stride_arr_a_odd,fn_order_angs_odd,S_e(idxs))  
            a_ctr = a_ctr + 1
        else
            call trig_series_odd(sinth,costh,Ncoeff_dihs_odd,fitted_dihs_coeff_odd(:,i-b_ctr-a_ctr),stride_arr_d_odd,fn_order_dihs_odd,S_e(idxs) )    
        endif
    enddo
    b_ctr =0;a_ctr = 0
    do i = 0,size(even_fn_idx)-1
        idxs = even_fn_idx(i)
        !print*,idxs,fitted_bonds_coeff_even(i,:)
        if (idxs<Nbonds) then
            call  trig_series_even(sinth,costh,Ncoeff_bonds_even,fitted_bonds_coeff_even(:,i),stride_arr_b_even,fn_order_bonds_even,S_e(idxs))
            b_ctr = b_ctr + 1
        else if ((idxs>=Nbonds).and.(idxs<Nbonds+Nang)) then
            call trig_series_even(sinth,costh,Ncoeff_angs_even,fitted_angs_coeff_even(:,i-b_ctr),stride_arr_a_even,fn_order_angs_even,S_e(idxs))  
            a_ctr = a_ctr + 1
        else
            call trig_series_even(sinth,costh,Ncoeff_dihs_even,fitted_dihs_coeff_even(:,i-b_ctr-a_ctr),stride_arr_d_even,fn_order_dihs_even,S_e(idxs))    
        endif
        !read(*,*)
    enddo
endsubroutine 
!****************************************************************************************************************************************************#
subroutine get_all_dfc(sinth,costh,dfc_quad_even,dfc_quad_odd,dfc_qubic_even_gt5,dfc_qubic_even_lt5,dfc_qubic_even_lt2,dfc_qubic_odd_gt5,dfc_qubic_odd_lt5,dfc_qubic_odd_lt2,dfc_quartic_even_gt5,dfc_quartic_even_lt5,dfc_quartic_even_lt2,dfc_quartic_odd_gt5,dfc_quartic_odd_lt5,dfc_quartic_odd_lt2,dfc_qo_even,dfc_qo_odd,dfc_qo_ang_even)
    implicit none
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
    real(8):: dfc_quad_even(0:n_fcs_quad_even-1,0:n_fix_tors-1),dfc_quad_odd(0:n_fcs_quad_odd-1,0:n_fix_tors-1)
    real(8):: dfc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1,0:n_fix_tors-1),dfc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1,0:n_fix_tors-1)
    real(8):: dfc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1,0:n_fix_tors-1),dfc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1,0:n_fix_tors-1)
    real(8):: dfc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1,0:n_fix_tors-1),dfc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1,0:n_fix_tors-1) 
    real(8):: dfc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1,0:n_fix_tors-1),dfc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1,0:n_fix_tors-1)
    real(8):: dfc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1,0:n_fix_tors-1) 
    real(8):: dfc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1,0:n_fix_tors-1)
    real(8):: dfc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1,0:n_fix_tors-1)
    real(8):: dfc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1,0:n_fix_tors-1) 
    real(8):: dfc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1,0:n_fix_tors-1),dfc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1,0:n_fix_tors-1)
    real(8):: dfc_qo_ang_even(0:n_fcs_qo_ang_even-1,0:n_fix_tors-1),dfc_qo_even(0:n_fcs_qo_even-1,0:n_fix_tors-1),dfc_qo_odd(0:n_fcs_qo_odd-1,0:n_fix_tors-1) 
    !Initialize to zero
    dfc_quad_even = 0.d0;dfc_quad_odd = 0.d0; dfc_qubic_even_gt5 = 0.d0;dfc_qubic_even_lt5 = 0.d0;dfc_qubic_even_lt2 = 0.d0;
    dfc_qubic_odd_gt5 = 0.d0;dfc_qubic_odd_lt5 = 0.d0;dfc_qubic_odd_lt2 = 0.d0;dfc_quartic_even_gt5 = 0.d0;
    dfc_quartic_even_lt5 = 0.d0;dfc_quartic_even_lt2 = 0.d0;dfc_quartic_odd_gt5 = 0.d0;dfc_quartic_odd_lt5 = 0.d0;dfc_quartic_odd_lt2 = 0.d0; 
    dfc_quartic_even_ltpt1 = 0.d0;dfc_quartic_odd_ltpt1=0.d0
    dfc_qo_ang_even=0.d0; dfc_qo_even=0.d0; dfc_qo_odd=0.d0
    !print*,size(dfc_quad_even),size(dfc_quad_odd)
    !print*,n_fcs_quad_even,n_fcs_quad_odd
    !print*,size(dfc_qubic_even_gt5),size(dfc_qubic_even_lt5),size(dfc_qubic_even_lt2)
    !print*,n_fcs_qubic_even_gt5,n_fcs_qubic_even_lt5,n_fcs_qubic_even_lt2
    !print*,size(dfc_qubic_odd_gt5),size(dfc_qubic_odd_lt5),size(dfc_qubic_odd_lt2)
    !print*,n_fcs_qubic_odd_gt5,n_fcs_qubic_odd_lt5,n_fcs_qubic_odd_lt2
    !print*,size(dfc_quartic_even_gt5),size(dfc_quartic_even_lt5),size(dfc_quartic_even_lt2),size(dfc_quartic_even_ltpt1)
    !print*,n_fcs_quartic_odd_gt5,n_fcs_quartic_odd_lt5,n_fcs_quartic_odd_lt2,n_fcs_quartic_odd_ltpt1
    !print*,size(dfc_quartic_odd_gt5),size(dfc_quartic_odd_lt5),size(dfc_quartic_odd_lt2),size(dfc_quartic_odd_ltpt1)
    !print*,n_fcs_quartic_even_gt5,n_fcs_quartic_even_lt5,n_fcs_quartic_even_lt2,n_fcs_quartic_even_ltpt1
    !print*,size(dfc_qo_even),size(dfc_qo_odd),size(dfc_qo_ang_even)
    !print*,n_fcs_qo_ang_even,n_fcs_qo_even,n_fcs_qo_odd
    !read(*,*)
    !print*,size(dfc_quartic_even_lt5)
    !For quadratic
    !Even
    call get_dfc(sinth,costh,Ncoeff_fij_even,fitted_fc_coeff_ij_even,stride_arr_quad_even,fn_order_fij_even,n_fcs_quad_even,0,dfc_quad_even)
    !Odd
    call get_dfc(sinth,costh,Ncoeff_fij_odd,fitted_fc_coeff_ij_odd,stride_arr_quad_odd,fn_order_fij_odd,n_fcs_quad_odd,1,dfc_quad_odd)
        
    !For qubic
    !Even
    !One body terms
    call get_dfc(sinth,costh,Ncoeff_fijk_even_gt5,fitted_fc_coeff_ijk_even_gt5,stride_arr_qubic_even_gt5,fn_order_fijk_even_gt5,n_fcs_qubic_even_gt5,0,dfc_qubic_even_gt5)
    
    !Two body terms
    call get_dfc(sinth,costh,Ncoeff_fijk_even_lt5,fitted_fc_coeff_ijk_even_lt5,stride_arr_qubic_even_lt5,fn_order_fijk_even_lt5,n_fcs_qubic_even_lt5,0,dfc_qubic_even_lt5)
    
    !Three body terms
    call get_dfc(sinth,costh,Ncoeff_fijk_even_lt2,fitted_fc_coeff_ijk_even_lt2,stride_arr_qubic_even_lt2,fn_order_fijk_even_lt2,n_fcs_qubic_even_lt2,0,dfc_qubic_even_lt2)
    
    !Odd:
    !One body terms
    call get_dfc(sinth,costh,Ncoeff_fijk_odd_gt5,fitted_fc_coeff_ijk_odd_gt5,stride_arr_qubic_odd_gt5,fn_order_fijk_odd_gt5,n_fcs_qubic_odd_gt5,1,dfc_qubic_odd_gt5)
    
    !Two body terms
    call get_dfc(sinth,costh,Ncoeff_fijk_odd_lt5,fitted_fc_coeff_ijk_odd_lt5,stride_arr_qubic_odd_lt5,fn_order_fijk_odd_lt5,n_fcs_qubic_odd_lt5,1,dfc_qubic_odd_lt5)
    
    !Three body terms
    call get_dfc(sinth,costh,Ncoeff_fijk_odd_lt2,fitted_fc_coeff_ijk_odd_lt2,stride_arr_qubic_odd_lt2,fn_order_fijk_odd_lt2,n_fcs_qubic_odd_lt2,1,dfc_qubic_odd_lt2)
    
    !For quartic
    !Even:
    !One body terms
    call get_dfc(sinth,costh,Ncoeff_fijkl_even_gt5,fitted_fc_coeff_ijkl_even_gt5,stride_arr_quartic_even_gt5,fn_order_fijkl_even_gt5,n_fcs_quartic_even_gt5,0,dfc_quartic_even_gt5)
       !print*,fitted_fc_coeff_ijkl_even_gt5 
    !print*,n_fcs_quartic_even_gt5;read(*,*)
    !print*,size(dfc_quartic_even_gt5),n_fcs_quartic_even_gt5
    !Two body terms
    call get_dfc(sinth,costh,Ncoeff_fijkl_even_lt5,fitted_fc_coeff_ijkl_even_lt5,stride_arr_quartic_even_lt5,fn_order_fijkl_even_lt5,n_fcs_quartic_even_lt5,0,dfc_quartic_even_lt5)
    
    !Three body terms
    call get_dfc(sinth,costh,Ncoeff_fijkl_even_lt2,fitted_fc_coeff_ijkl_even_lt2,stride_arr_quartic_even_lt2,fn_order_fijkl_even_lt2,n_fcs_quartic_even_lt2,0,dfc_quartic_even_lt2)
    
    !Odd:
    !Two body terms
    call get_dfc(sinth,costh,Ncoeff_fijkl_odd_gt5,fitted_fc_coeff_ijkl_odd_gt5,stride_arr_quartic_odd_gt5,fn_order_fijkl_odd_gt5,n_fcs_quartic_odd_gt5,1,dfc_quartic_odd_gt5)
    call get_dfc(sinth,costh,Ncoeff_fijkl_odd_lt5,fitted_fc_coeff_ijkl_odd_lt5,stride_arr_quartic_odd_lt5,fn_order_fijkl_odd_lt5,n_fcs_quartic_odd_lt5,1,dfc_quartic_odd_lt5)
    call get_dfc(sinth,costh,Ncoeff_fijkl_odd_lt2,fitted_fc_coeff_ijkl_odd_lt2,stride_arr_quartic_odd_lt2,fn_order_fijkl_odd_lt2,n_fcs_quartic_odd_lt2,1,dfc_quartic_odd_lt2)
    
    !Three body terms
    call get_dfc(sinth,costh,Ncoeff_qo_even,fitted_fc_coeff_qo_even,stride_arr_qo_even,fn_order_qo_even,n_fcs_qo_even,0,dfc_qo_even)
    call get_dfc(sinth,costh,Ncoeff_qo_ang_even,fitted_fc_coeff_qo_ang_even,stride_arr_qo_ang_even,fn_order_qo_ang_even,n_fcs_qo_ang_even,0,dfc_qo_ang_even)
    call get_dfc(sinth,costh,Ncoeff_qo_odd,fitted_fc_coeff_qo_odd,stride_arr_qo_odd,fn_order_qo_odd,n_fcs_qo_odd,1,dfc_qo_odd)
endsubroutine
!****************************************************************************************************************************************************#
subroutine get_all_fc(sinth,costh,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even)
    implicit none 
    real(8),intent(in)::sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0)) 
    real(8),intent(out):: fc_quad_even(0:n_fcs_quad_even-1),fc_quad_odd(0:n_fcs_quad_odd-1),fc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1),fc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1)
    real(8),intent(out):: fc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1),fc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1),fc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1),fc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1) 
    real(8),intent(out):: fc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1),fc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1),fc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1) 
    real(8),intent(out):: fc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1),fc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1),fc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1) 
    real(8),intent(out):: fc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1),fc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1) 
    real(8),intent(out):: fc_qo_ang_even(0:n_fcs_qo_ang_even-1),fc_qo_even(0:n_fcs_qo_even-1),fc_qo_odd(0:n_fcs_qo_odd-1) 
    !Initialize to zero
    fc_quad_even = 0.d0;fc_quad_odd = 0.d0; fc_qubic_even_gt5 = 0.d0;fc_qubic_even_lt5 = 0.d0;fc_qubic_even_lt2 = 0.d0;
    fc_qubic_odd_gt5 = 0.d0;fc_qubic_odd_lt5 = 0.d0;fc_qubic_odd_lt2 = 0.d0;fc_quartic_even_gt5 = 0.d0;
    fc_quartic_even_lt5 = 0.d0;fc_quartic_even_lt2 = 0.d0;fc_quartic_odd_gt5 = 0.d0;fc_quartic_odd_lt5 = 0.d0;fc_quartic_odd_lt2 = 0.d0; 
    fc_quartic_even_ltpt1 =0.d0;fc_quartic_odd_ltpt1 =0.d0;fc_qo_ang_even=0.d0;fc_qo_even=0.d0;fc_qo_odd=0.d0
    !For quadratic
    !Even
    call get_fc(sinth,costh,Ncoeff_fij_even,fitted_fc_coeff_ij_even,stride_arr_quad_even,fn_order_fij_even,n_fcs_quad_even,0,fc_quad_even)
    !Odd
    call get_fc(sinth,costh,Ncoeff_fij_odd,fitted_fc_coeff_ij_odd,stride_arr_quad_odd,fn_order_fij_odd,n_fcs_quad_odd,1,fc_quad_odd)
    !For qubic
    !Even
    !One body terms
    call get_fc(sinth,costh,Ncoeff_fijk_even_gt5,fitted_fc_coeff_ijk_even_gt5,stride_arr_qubic_even_gt5,fn_order_fijk_even_gt5,n_fcs_qubic_even_gt5,0,fc_qubic_even_gt5)
    
    !Two body terms
       call get_fc(sinth,costh,Ncoeff_fijk_even_lt5,fitted_fc_coeff_ijk_even_lt5,stride_arr_qubic_even_lt5,fn_order_fijk_even_lt5,n_fcs_qubic_even_lt5,0,fc_qubic_even_lt5)
    
    !Three body terms
    call get_fc(sinth,costh,Ncoeff_fijk_even_lt2,fitted_fc_coeff_ijk_even_lt2,stride_arr_qubic_even_lt2,fn_order_fijk_even_lt2,n_fcs_qubic_even_lt2,0,fc_qubic_even_lt2)
    
    !Odd:
    !One body terms
    call get_fc(sinth,costh,Ncoeff_fijk_odd_gt5,fitted_fc_coeff_ijk_odd_gt5,stride_arr_qubic_odd_gt5,fn_order_fijk_odd_gt5,n_fcs_qubic_odd_gt5,1,fc_qubic_odd_gt5)
    
    !Two body terms
    call get_fc(sinth,costh,Ncoeff_fijk_odd_lt5,fitted_fc_coeff_ijk_odd_lt5,stride_arr_qubic_odd_lt5,fn_order_fijk_odd_lt5,n_fcs_qubic_odd_lt5,1,fc_qubic_odd_lt5)
    
    !Three body terms
    call get_fc(sinth,costh,Ncoeff_fijk_odd_lt2,fitted_fc_coeff_ijk_odd_lt2,stride_arr_qubic_odd_lt2,fn_order_fijk_odd_lt2,n_fcs_qubic_odd_lt2,1,fc_qubic_odd_lt2)
    
    !For quartic
    !Even:
    !One body terms
    call get_fc(sinth,costh,Ncoeff_fijkl_even_gt5,fitted_fc_coeff_ijkl_even_gt5,stride_arr_quartic_even_gt5,fn_order_fijkl_even_gt5,n_fcs_quartic_even_gt5,0,fc_quartic_even_gt5)
    
    !Two body terms
    call get_fc(sinth,costh,Ncoeff_fijkl_even_lt5,fitted_fc_coeff_ijkl_even_lt5,stride_arr_quartic_even_lt5,fn_order_fijkl_even_lt5,n_fcs_quartic_even_lt5,0,fc_quartic_even_lt5)
    
    !Three body terms
    call get_fc(sinth,costh,Ncoeff_fijkl_even_lt2,fitted_fc_coeff_ijkl_even_lt2,stride_arr_quartic_even_lt2,fn_order_fijkl_even_lt2,n_fcs_quartic_even_lt2,0,fc_quartic_even_lt2)
    
    !Odd:
    !Two body terms
    !read(*,*)
    call get_fc(sinth,costh,Ncoeff_fijkl_odd_gt5,fitted_fc_coeff_ijkl_odd_gt5,stride_arr_quartic_odd_gt5,fn_order_fijkl_odd_gt5,n_fcs_quartic_odd_gt5,1,fc_quartic_odd_gt5)
    !!print*,'odd_gt5=',Ncoeff_fijkl_odd_gt5,fn_order_fijkl_odd_gt5,n_fcs_quartic_odd_gt5,stride_arr_quartic_odd_gt5;read(*,*)
    call get_fc(sinth,costh,Ncoeff_fijkl_odd_lt5,fitted_fc_coeff_ijkl_odd_lt5,stride_arr_quartic_odd_lt5,fn_order_fijkl_odd_lt5,n_fcs_quartic_odd_lt5,1,fc_quartic_odd_lt5)
    
    !Three body terms
    call get_fc(sinth,costh,Ncoeff_fijkl_odd_lt2,fitted_fc_coeff_ijkl_odd_lt2,stride_arr_quartic_odd_lt2,fn_order_fijkl_odd_lt2,n_fcs_quartic_odd_lt2,1,fc_quartic_odd_lt2)
    fc_quartic_odd_ltpt1 = fitted_fc_coeff_ijkl_odd_ltpt1
    fc_quartic_even_ltpt1 = fitted_fc_coeff_ijkl_even_ltpt1
 
    call get_fc(sinth,costh,Ncoeff_qo_odd,fitted_fc_coeff_qo_odd,stride_arr_qo_odd,fn_order_qo_odd,n_fcs_qo_odd,1,fc_qo_odd)
    !print*,'qo=',fitted_fc_coeff_qo_odd;read(*,*)
    call get_fc(sinth,costh,Ncoeff_qo_even,fitted_fc_coeff_qo_even,stride_arr_qo_even,fn_order_qo_even,n_fcs_qo_even,0,fc_qo_even)
    !print*,'qo=',fitted_fc_coeff_qo_even;read(*,*)
    call get_fc(sinth,costh,Ncoeff_qo_ang_even,fitted_fc_coeff_qo_ang_even,stride_arr_qo_ang_even,fn_order_qo_ang_even,n_fcs_qo_ang_even,0,fc_qo_ang_even)
endsubroutine

!!****************************************************************************************#
!subroutine get_dsymdx(X,dsymdx):
!    implicit none
!    real(8),intent(in)::X(0:Ncarts-1)
!    real(8),intent(out)::dsymdx(0:Ncarts-1,0:N_int-1)
!    integer::i,j,istep
!    real(8)::Xtmp(0:Ncarts-1),dsym_ar(0:4,0:N_int-1),stmp(0:N_int-1),symtmp(0:N_int-1),S(0:N_int-1)
!    dsymdx = 0.d0;Xtmp=0.d0;,dsym_ar=0.d0;stmp=0.d0;symtmp=0.d0;S=0.d0
!    !Compute internals
!    cart_to_internal(X,S)
!    !SYmmetrized coordinate
!    sym_S = matmul(sym_mat,S)
!    
!    do i = 0,Ncarts-1
!        dsym_ar = 0.d0
!
!        do istep =-2,2
!            Xtmp(:) = X(:)
!            Xtmp(i) = X(i) + dble(istep) * dx
!            Stmp    = cart_to_internal(Xtmp,Stmp)
!            symtmp = matmul(sym_mat,Stmp)
!            !print("sym_tmp=",i,istep,sym_tmp)
!            do  j = Nbonds+Nang,N_int-1
!                if ((abs(sym_S(j)-180.0)<1.e-2) .or. (abs(sym_S[j])<1.e-2))
!                    sym_tmp(j) = sym_tmp(j) + 0.0
!                    !Stmp[j]   = mod(Stmp[j],360.0)
!                endif
!            enddo
!            dsym_ar(istep+2)(:)   = sym_tmp(:)
!        enddo
!        dsym_dx(i,:) = matmul(five_pt_first_deriv(:),dsym_ar(:,:)) * dx_inv
!    
!    enddo
!endsubroutine

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
    do j = 1,max_id(0)
        do i = 0,n_fix_tors-1
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
        S(i) = S(i) * degree_to_radian
    endif
enddo
endsubroutine 
!****************************************************************************************#
subroutine get_pot(X,pot,Vt_out)
    implicit none 
   integer::i
    real(8),intent(in):: X(0:Ncarts-1)
    real(8) :: Xtmp(0:Ncarts-1)
    real(8),intent(out)::pot,Vt_out
    real(8):: S(0:N_int-1),sym_S(0:N_int-1) ,sym_Seq(0:N_int-1) ,S_eq(0:N_int-1) ,Sym_diff(0:N_int-1) ,S_diff(0:N_int-1) ,dVdsym(0:N_int-1), dVds(0:N_int-1), dsdx(0:Ncarts-1,0:N_int-1) 
    real(8) :: fix_tors(0:n_fix_tors-1),sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
    real(8):: dfc_quad_even(0:n_fcs_quad_even-1,0:n_fix_tors-1),dfc_quad_odd(0:n_fcs_quad_odd-1,0:n_fix_tors-1),dfc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1,0:n_fix_tors-1),dfc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1,0:n_fix_tors-1)
    real(8):: dfc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1,0:n_fix_tors-1),dfc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1,0:n_fix_tors-1),dfc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1,0:n_fix_tors-1),dfc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1,0:n_fix_tors-1) 
    real(8):: dfc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1,0:n_fix_tors-1),dfc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1,0:n_fix_tors-1),dfc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1,0:n_fix_tors-1) 
    real(8):: dfc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1,0:n_fix_tors-1),dfc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1,0:n_fix_tors-1),dfc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1,0:n_fix_tors-1) 
    real(8):: fc_quad_even(0:n_fcs_quad_even-1),fc_quad_odd(0:n_fcs_quad_odd-1),fc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1),fc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1)
    real(8):: fc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1),fc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1),fc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1),fc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1) 
    real(8):: fc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1),fc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1),fc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1) 
    real(8):: fc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1),fc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1),fc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1) 
    real(8):: dfc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1,0:n_fix_tors-1),dfc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1,0:n_fix_tors-1),dfc_qo_even(0:n_fcs_qo_even-1,0:n_fix_tors-1),dfc_qo_odd(0:n_fcs_qo_odd-1,0:n_fix_tors-1),dfc_qo_ang_even(0:n_fcs_qo_ang_even-1,0:n_fix_tors-1) 
    real(8):: fc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1),fc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1) 
    real(8):: fc_qo_ang_even(0:n_fcs_qo_ang_even-1),fc_qo_even(0:n_fcs_qo_even-1),fc_qo_odd(0:n_fcs_qo_odd-1) 
    fc_quad_even = 0.d0;fc_quad_odd = 0.d0; fc_qubic_even_gt5 = 0.d0;fc_qubic_even_lt5 = 0.d0;fc_qubic_even_lt2 = 0.d0;
    fc_qubic_odd_gt5 = 0.d0;fc_qubic_odd_lt5 = 0.d0;fc_qubic_odd_lt2 = 0.d0;fc_quartic_even_gt5 = 0.d0;
    fc_quartic_even_lt5 = 0.d0;fc_quartic_even_lt2 = 0.d0;fc_quartic_odd_gt5 = 0.d0;fc_quartic_odd_lt5 = 0.d0;fc_quartic_odd_lt2 = 0.d0; 
    dfc_quad_even = 0.d0;dfc_quad_odd = 0.d0; dfc_qubic_even_gt5 = 0.d0;dfc_qubic_even_lt5 = 0.d0;dfc_qubic_even_lt2 = 0.d0;
    dfc_qubic_odd_gt5 = 0.d0;dfc_qubic_odd_lt5 = 0.d0;dfc_qubic_odd_lt2 = 0.d0;dfc_quartic_even_gt5 = 0.d0;
    dfc_quartic_even_lt5 = 0.d0;dfc_quartic_even_lt2 = 0.d0;dfc_quartic_odd_gt5 = 0.d0;dfc_quartic_odd_lt5 = 0.d0;dfc_quartic_odd_lt2 = 0.d0; 
    dfc_qo_ang_even=0.d0;dfc_qo_even=0.d0;dfc_qo_odd=0.d0
    dfc_quartic_even_ltpt1 =0.d0;dfc_quartic_odd_ltpt1 =0.d0
    fc_quartic_even_ltpt1 =0.d0;fc_quartic_odd_ltpt1 =0.d0;fc_qo_ang_even=0.d0;fc_qo_even=0.d0;fc_qo_odd=0.d0
    pot =0.d0;Vt_out=0.d0;sinth = 0.d0;costh=0.d0 
    Xtmp = 0.d0;S=0.d0;sym_Seq = 0.d0 ;sym_diff = 0.d0
    !Convert Bohr to Angstrom
    Xtmp =  X * BohrToAng
    !print*,'X=',Xtmp;read(*,*)
    !Convert cartesian to internal
    !print*,'S=',S;read(*,*)
    call cart_to_internal(Xtmp,S)
    !print*,'S=',S;read(*,*)
    !Fixed torsions array
    fix_tors = S(tors_idx) * degree_to_radian
    !fix_tors = fix_tors * 
    !Get sin and cosine arrays
    call get_csarr(fix_tors,max_id,sinth,costh)

    !Get equlibrium S
    call get_se(sinth,costh,S_eq)
    sym_Seq  =  matmul(sym_mat,S_eq)
    do i = Nbonds+Nang,N_int-1
        !Check if S(i) is close to 180 degree
        if (abs(abs(S(i))-180.d0)<90.d0) then
            !If S(i) and S_eq(i) have opposite signs
            if (S(i)*sym_Seq(i)<0.d0) then
               if (sym_Seq(i)>0.d0) then
                  S(i) = S(i) + 360.d0
               else
                  S(i) = S(i) - 360.d0
               endif                                       
            endif
        endif
    enddo
    !Convert internal to symmetric internal
    !sym_S = matmul(sym_mat,S)
    sym_diff = S-sym_Seq
    S_diff = matmul(sym_mat,sym_diff)
    !S_diff = sym_S-S_eq
    !Convert to scaled coordinates
    call scale_coord(S_diff)

    !compute fitted force constants(fc)
    call get_all_fc(sinth,costh,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even)
    call  compute_pot(S_diff,sinth,costh,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even,pot,Vt_out)
    pot = pot * CminvToHartee
    Vt_out = Vt_out * CminvToHartee
endsubroutine 
!****************************************************************************************#
subroutine get_force(X,dVdx,pot,flag)
    implicit none 
    integer::i
    real(8),intent(in):: X(0:Ncarts-1)
    real(8) :: Xtmp(0:Ncarts-1),Vt_out 
    integer,intent(in)::flag
    real(8),intent(out):: dVdx(0:Ncarts-1),pot
    real(8):: S(0:N_int-1) ,sym_S(0:N_int-1) ,sym_Seq(0:N_int-1) ,S_eq(0:N_int-1) ,Sym_diff(0:N_int-1) 
    real(8):: S_diff(0:N_int-1) ,dVdsym(0:N_int-1), dVds(0:N_int-1), dsdx(0:Ncarts-1,0:N_int-1) 
    real(8):: fix_tors(0:n_fix_tors-1),sinth(0:n_fix_tors-1,0:max_id(0)),costh(0:n_fix_tors-1,0:max_id(0))
    real(8):: dfc_quad_even(0:n_fcs_quad_even-1,0:n_fix_tors-1),dfc_quad_odd(0:n_fcs_quad_odd-1,0:n_fix_tors-1)
    real(8):: dfc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1,0:n_fix_tors-1),dfc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1,0:n_fix_tors-1)
    real(8):: dfc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1,0:n_fix_tors-1),dfc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1,0:n_fix_tors-1)
    real(8):: dfc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1,0:n_fix_tors-1),dfc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1,0:n_fix_tors-1) 
    real(8):: dfc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1,0:n_fix_tors-1),dfc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1,0:n_fix_tors-1)
    real(8):: dfc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1,0:n_fix_tors-1) 
    real(8):: dfc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1,0:n_fix_tors-1)
    real(8):: dfc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1,0:n_fix_tors-1)
    real(8):: dfc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1,0:n_fix_tors-1) 
    real(8):: dfc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1,0:n_fix_tors-1),dfc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1,0:n_fix_tors-1)
    real(8):: dfc_qo_ang_even(0:n_fcs_qo_ang_even-1,0:n_fix_tors-1),dfc_qo_even(0:n_fcs_qo_even-1,0:n_fix_tors-1),dfc_qo_odd(0:n_fcs_qo_odd-1,0:n_fix_tors-1) 
    real(8):: fc_quad_even(0:n_fcs_quad_even-1),fc_quad_odd(0:n_fcs_quad_odd-1),fc_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1)
    real(8):: fc_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1),fc_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1),fc_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1) 
    real(8):: fc_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1),fc_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1) 
    real(8):: fc_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1),fc_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1),fc_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1) 
    real(8):: fc_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1),fc_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1),fc_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1) 
    real(8):: fc_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1),fc_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1) 
    real(8):: fc_qo_ang_even(0:n_fcs_qo_ang_even-1),fc_qo_even(0:n_fcs_qo_even-1),fc_qo_odd(0:n_fcs_qo_odd-1) 
    fc_quad_even = 0.d0;fc_quad_odd = 0.d0; fc_qubic_even_gt5 = 0.d0;fc_qubic_even_lt5 = 0.d0;fc_qubic_even_lt2 = 0.d0;
    fc_qubic_odd_gt5 = 0.d0;fc_qubic_odd_lt5 = 0.d0;fc_qubic_odd_lt2 = 0.d0;fc_quartic_even_gt5 = 0.d0;
    fc_quartic_even_lt5 = 0.d0;fc_quartic_even_lt2 = 0.d0;fc_quartic_odd_gt5 = 0.d0;fc_quartic_odd_lt5 = 0.d0;fc_quartic_odd_lt2 = 0.d0; 
    fc_quartic_even_ltpt1 =0.d0;fc_quartic_odd_ltpt1 =0.d0;fc_qo_ang_even=0.d0;fc_qo_even=0.d0;fc_qo_odd=0.d0
    dfc_quad_even = 0.d0;dfc_quad_odd = 0.d0; dfc_qubic_even_gt5 = 0.d0;dfc_qubic_even_lt5 = 0.d0;dfc_qubic_even_lt2 = 0.d0;
    dfc_qubic_odd_gt5 = 0.d0;dfc_qubic_odd_lt5 = 0.d0;dfc_qubic_odd_lt2 = 0.d0;dfc_quartic_even_gt5 = 0.d0;
    dfc_quartic_even_lt5 = 0.d0;dfc_quartic_even_lt2 = 0.d0;dfc_quartic_odd_gt5 = 0.d0;dfc_quartic_odd_lt5 = 0.d0;dfc_quartic_odd_lt2 = 0.d0; 
    dfc_qo_ang_even=0.d0;dfc_qo_even=0.d0;dfc_qo_odd=0.d0
    dfc_quartic_even_ltpt1 =0.d0;dfc_quartic_odd_ltpt1 =0.d0

    dVdx = 0.d0;sinth = 0.d0;costh=0.d0;dVdsym=0.d0;dVds=0.d0;
    Vt_out =0.d0;    sym_Seq = 0.d0 ;sym_diff = 0.d0
    !print*,size(dfc_quartic_even_gt5),n_fcs_quartic_even_gt5
    !print*,n_fcs_quartic_even_gt5;read(*,*)
    !Convert Bohr to Angstrom
    Xtmp = BohrToAng * X 
    !Convert cartesian to internal (Angstrom and degree)
    call cart_to_internal(Xtmp,S)
    !Fixed torsions array (in radian)
    fix_tors = S(tors_idx) * degree_to_radian
    !Get sin and cosine arrays
    call get_csarr(fix_tors,max_id,sinth,costh)
    !print*,"sinth=",sinth
    !print*,"costh=",costh
    !Get equlibrium S(Angstrom and degree)
    call get_se(sinth,costh,S_eq)
    sym_Seq    =  matmul(sym_mat,S_eq)
    do i = Nbonds+Nang,N_int-1
        !Check if S(i) is close to 180 degree
        if (abs(abs(S(i))-180.d0)<90.d0) then
            !If S(i) and S_eq(i) have opposite signs
            if (S(i)*sym_Seq(i)<0.d0) then
               if (sym_Seq(i)>0.d0) then
                  S(i) = S(i) + 360.d0
               else
                  S(i) = S(i) - 360.d0
               endif                                       
            endif
        endif
    enddo
    !print*,'S=',S
    !print*,'S_eq=',S_eq;
    !print*,'Sym_eq=',matmul(sym_mat,S_eq)
    !Convert internal to symmetric internal (Angstrom and degree)
    !sym_S = matmul(sym_mat,S)
    !read(*,*)
    !print*,'sym_S=',sym_S
    !Compute deviation from the equilibrium (Angstrom and degree)
    !S_diff =sym_S-S_eq
    sym_diff = S-sym_Seq
    S_diff = matmul(sym_mat,sym_diff)
    !print*,'Sym_diff =',S_diff;read(*,*)
    !Convert to scaled coordinates (dimensionless)
    call scale_coord(S_diff)
    !compute fitted force constants(fc)
    call get_all_fc(sinth,costh,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even)
    !print*,'fc'
    !Compute derivative of fitted fc (dfc)
    call get_all_dfc(sinth,costh,dfc_quad_even,dfc_quad_odd,dfc_qubic_even_gt5,dfc_qubic_even_lt5,dfc_qubic_even_lt2,dfc_qubic_odd_gt5,dfc_qubic_odd_lt5,dfc_qubic_odd_lt2,dfc_quartic_even_gt5,dfc_quartic_even_lt5,dfc_quartic_even_lt2,dfc_quartic_odd_gt5,dfc_quartic_odd_lt5,dfc_quartic_odd_lt2,dfc_qo_even,dfc_qo_odd,dfc_qo_ang_even)
    !print*,'dfc'
    !Compute Dvdsym(cm-1)
    call getDvds(S_diff,sinth,costh,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even,dfc_quad_even,dfc_quad_odd,dfc_qubic_even_gt5,dfc_qubic_even_lt5,dfc_qubic_even_lt2,dfc_qubic_odd_gt5,dfc_qubic_odd_lt5,dfc_qubic_odd_lt2,dfc_quartic_even_gt5,dfc_quartic_even_lt5,dfc_quartic_even_lt2,dfc_quartic_odd_gt5,dfc_quartic_odd_lt5,dfc_quartic_odd_lt2,dfc_quartic_even_ltpt1,dfc_quartic_odd_ltpt1,dfc_qo_even,dfc_qo_odd,dfc_qo_ang_even,dVdsym)
    !Convert to cm^-1/degree or cm^-1/angstrom 
    do i = 0,N_int-1
        dVdsym(i) = dVdsym(i) * dim_scal_factor_inv(i)
        !Radian to degree
        if (i>= Nbonds) then
                dVdsym(i) = dVdsym(i) * degree_to_radian
        endif
    enddo
    !Compute dVds (cm^-1)
    dVds = matmul(sym_mat_transpose,dVdsym)
    !Get dsdx (degree/angstrom)
    call get_dsdx(Xtmp,dsdx)
    !do i = 0,Ncarts-1
    !   print*,'dsdx=',i,dsdx(i,:)
    !enddo
    !get cart_grad
    call int_to_cart_grad(dsdx,dVds,dVdx)
    dVdx = dVdx * CminvToHartee * BohrToAng 
    if (flag==1) then
        call  compute_pot(S_diff,sinth,costh,fc_quad_even,fc_quad_odd,fc_qubic_even_gt5,fc_qubic_even_lt5,fc_qubic_even_lt2,fc_qubic_odd_gt5,fc_qubic_odd_lt5,fc_qubic_odd_lt2,fc_quartic_even_gt5,fc_quartic_even_lt5,fc_quartic_even_lt2,fc_quartic_odd_gt5,fc_quartic_odd_lt5,fc_quartic_odd_lt2,fc_quartic_even_ltpt1,fc_quartic_odd_ltpt1,fc_qo_even,fc_qo_odd,fc_qo_ang_even,pot,Vt_out)
    endif
endsubroutine 
!****************************************************************************************#
subroutine numerical_dVdx(X,dVdx_num)
implicit none 
real(8),intent(in)::X(0:Ncarts-1)
real(8),intent(out)::dVdx_num(0:Ncarts-1)
integer::i,istep,flag=1
real(8)::Vt=0.d0,Xinp(0:Ncarts-1),Xtmp(0:Ncarts-1),pot_arr(0:4),pot=0.d0,dum(0:Ncarts-1)
Xtmp = 0.d0;Xinp=0.d0
pot = 0.d0
pot_arr = 0.d0;dum = 0.d0
dVdx_num = 0.d0
Xinp(:)= X(:) * BohrToAng
!loop over the cartesian coordinates
do i = 0,Ncarts-1
    Xtmp(:) = Xinp(:)
    pot = 0.d0
    !displace five points and save the potential at each point
    do istep = 0,4
        Xtmp(i) = Xinp(i)+ dx*dble(istep-2)
        call get_pot(Xtmp,pot,vt)
        pot_arr(istep) = pot    
    enddo
    dVdx_num(i) = dx_inv * dot_product(five_pt_1st_deriv,pot_arr)  
enddo
dVdx_num = dVdx_num * BohrToAng
endsubroutine
end module 
