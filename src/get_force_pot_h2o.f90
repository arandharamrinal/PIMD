!************************************************************************************************!
!************************************************************************************************!
subroutine get_H2O_pot()
!Given the cartesian coordinates,force constants, this routine computes
!the cartesian gradient and potential
!for H2O molecule
use allvars
use potvars
use constants 
implicit none
real(8),allocatable ::tmp_r(:),bond_vec(:),bond_uvec(:),bond_vec1(:),bond_vec2(:),bond_uvec1(:),bond_uvec2(:)
real(8) :: bond_len,bond_dr,bond_dot_prod,cos_theta,theta
real(8)::rtmp1(0:2),rtmp2(0:2),del_theta,sinth_inv,bond_dr1,bond_dr2,bond_len1,bond_len2
integer  ::i_tba,a0,a1,a2,b0,b1,i_ang,i_bond,i_bead
allocate(tmp_r(0:Ndim-1),bond_vec(0:Ndim-1),bond_uvec(0:Ndim-1),bond_vec1(0:Ndim-1),bond_vec2(0:Ndim-1),bond_uvec1(0:Ndim-1),bond_uvec2(0:Ndim-1))
dVdq = 0.d0
ext_V= 0.d0
do i_tba = 0,tot_bond_angles-1
    !Check whether it is a bond or an angle and then compute the
    !gradient and potential
        if (bond_pair(i_tba,2)==0) then
        !For bonds
           do i_bead = 0,Nbeads-1
              b0 = bond_pair(i_tba,0)
              b1 = bond_pair(i_tba,1)
              bond_vec            = q(3*b0:3*b0+2,i_bead) - q(3*b1:3*b1+2,i_bead)
              bond_len            = dsqrt(dot_product(bond_vec,bond_vec))
              bond_uvec           = bond_vec/bond_len
              bond_dr             = bond_len - bond_eq
              ext_V(i_bead)       = ext_V(i_bead)  + 0.5d0 * k_bond * bond_dr * bond_dr
              dVdq(3*b0:3*b0+2,i_bead)   = dVdq(3*b0:3*b0+2,i_bead) + bond_uvec(:) * k_bond * bond_dr
              dVdq(3*b1:3*b1+2,i_bead)   = dVdq(3*b1:3*b1+2,i_bead) - bond_uvec(:) * k_bond * bond_dr
           enddo
        else if (bond_pair(i_tba,3)==0) then
        !For angles
          do i_bead  = 0, Nbeads-1
              a0 = bond_pair(i_tba,0)
              a1 = bond_pair(i_tba,1)
              a2 = bond_pair(i_tba,2)
              bond_vec1    = q(3*a1:3*a1+2,i_bead) - q(3*a0:3*a0+2,i_bead)
              bond_vec2    = q(3*a1:3*a1+2,i_bead) - q(3*a2:3*a2+2,i_bead)
              bond_len1       = dsqrt(dot_product(bond_vec1,bond_vec1))
              bond_len2       = dsqrt(dot_product(bond_vec2,bond_vec2))
              bond_uvec1      = bond_vec1/bond_len1
              bond_uvec2      = bond_vec2/bond_len2
              !bond_dr1        = bond_len1 - bond_eq
              !bond_dr2        = bond_len2 - bond_eq
              bond_dot_prod   = dot_product(bond_vec1,bond_vec2)
              cos_theta       = bond_dot_prod/(bond_len1*bond_len2)
              theta           = dacos(cos_theta)
              del_theta       = theta - ang_eq
              ext_V(i_bead)   = ext_V(i_bead) + 0.5d0 * k_ang * del_theta * del_theta
              sinth_inv = 1.d0/dsin(theta)
              rtmp1 = - k_ang*del_theta*( bond_dot_prod * bond_vec1 / (bond_len1**3 * bond_len2 ) - bond_vec2 / (bond_len1*bond_len2) ) * sinth_inv
              rtmp2 = - k_ang*del_theta*( bond_dot_prod * bond_vec2 / (bond_len1 * bond_len2**3 ) - bond_vec1 / (bond_len1*bond_len2) ) * sinth_inv

              dVdq(3*a0:3*a0+2,i_bead) = dVdq(3*a0:3*a0+2,i_bead) + rtmp1
              dVdq(3*a2:3*a2+2,i_bead) = dVdq(3*a2:3*a2+2,i_bead) + rtmp2
              dVdq(3*a1:3*a1+2,i_bead) = dVdq(3*a1:3*a1+2,i_bead) - rtmp1 - rtmp2
           enddo
        endif
enddo
!pot = sum(ext_V)
end subroutine
!************************************************************************************************!

