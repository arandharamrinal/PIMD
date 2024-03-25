!****************************************************************************************************#
!******************Unitary matrix for cartesian to Normaal mode Transformation*******************!
!************************************************************************************************!
subroutine getCmat()
use allvars, only : Nbeads,Cmat, NbeadsInv, CmatInv
use constants,only : pi
implicit none 
integer::i,j,k
real(8)::dum,di,dj,dk
do j = 0,Nbeads-1
   dj = dble(j)
   do k = 0,Nbeads-1
        dk = dble(k)
        if (k==0) then
            CmatInv(j,k)= dsqrt(1.d0*NbeadsInv)
        else if ((k>0).and.(k < anint(0.5d0*Nbeads))) then
                CmatInv(j,k) =dsqrt(2.d0*NbeadsInv)*dcos(2.d0*pi*dj*dk*NbeadsInv)
        else if (k==anint(0.5d0*Nbeads))  then
             CmatInv(j,k)= dsqrt(1.d0*NbeadsInv)*(-1.d0)**j
        else if ((k > anint(0.5d0*Nbeads).and.(k < Nbeads))) then
             CmatInv(j,k) = dsqrt(2.d0*NbeadsInv)*dsin(2.d0*pi*dj*dk*NbeadsInv)
        endif
   enddo
enddo
!Cmat is the orthogonal matrix of dimension (Nbeads * Nbeads)  for generating the Normal modes
!CmatInv the inverse of Cmat and since Cmat is orthogonal and real, CmatInv is same as the 
!transpose of Cmat
Cmat  = transpose(CmatInv)
end subroutine
!************************************************************************************************!
!***************Compute the eigen values of the transformation matrix****************************!
!************************************************************************************************!
subroutine eigVals()
use allvars, only: Nbeads, lamda, Nbeads, piN
implicit none 
integer :: i,k
real(8) :: di
do k =0,Nbeads-1
   lamda(k) = 4.d0 * Nbeads * dsin ( dble(k) * piN ) * dsin( dble(k) * piN )   
enddo
end subroutine
!****************************************************************************************************#
