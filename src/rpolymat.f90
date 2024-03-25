!Matrix for free ring polymer dynamics 
!************************************************************************************************!
subroutine polyMat()
use allvars, only:Natoms,poly,polyHalf,massInv,dt,Nbeads,omegaK,mass
implicit none
double precision ::  twown,  wk, wt, wm, cos_wt, sinwt
integer :: i, j, k
real(8) :: dk
poly = 0.d0;polyHalf = 0.d0
do j = 0, Natoms-1
        poly(j,0,0) = 1.0d0
        poly(j,1,0) = 0.d0
        poly(j,2,0) = dt * massInv(j)
        poly(j,3,0) = 1.0d0
        polyHalf(j,0,0) = 1.0d0
        polyHalf(j,1,0) = 0.d0
        polyHalf(j,2,0) = 0.5d0 * dt * massInv(j)
        polyHalf(j,3,0) = 1.0d0
end do
do k = 1,Nbeads-1
   do j = 0, Natoms-1
        poly(j,0,k) = dcos(omegaK(k)*dt)
        poly(j,1,k) = -mass(j) * omegaK(k) * dsin(omegaK(k)*dt)
        poly(j,2,k) = massInv(j)/omegaK(k)*dsin(omegaK(k)*dt)
        poly(j,3,k) = dcos(omegaK(k)*dt)
        polyHalf(j,0,k) = dcos(omegaK(k)*0.5d0 * dt)
        polyHalf(j,1,k) = -mass(j) * omegaK(k) * dsin(omegaK(k)*0.5d0 * dt)
        polyHalf(j,2,k) = massInv(j)/omegaK(k)*dsin(omegaK(k)*0.5d0 * dt)
        polyHalf(j,3,k) = dcos(omegaK(k)*0.5d0 * dt)
   end do
end do
end subroutine
!****************************************************************************************************#  

!Matrix for free ring polymer dynamics 
!************************************************************************************************!
subroutine CayleyPolyMat()
use allvars, only:Natoms,CayleyPoly,sqrtCayleyPoly,massInv,dt,Nbeads,omegaK,mass
implicit none
real(8) :: detA,sqrtDetA, fact, sqrt_fact, a3,b3,c3,d3,a1,b1,c1,d1,a2,b2,c2,d2,twown,  wk, wt, wm, cos_wt, sinwt
integer :: i, j, k
real(8) :: dk
CayleyPoly = 0.d0
sqrtCayleyPoly = 0.d0
sqrt_fact = 0.d0
fact = 0.d0
do j = 0, Natoms-1
   do k = 0,Nbeads-1
        !I + 0.5*dt*A
        a1 = dt  * omegaK(k) * omegaK(k)
        b1 = dt  * dt * omegaK(k) * omegaK(k)
        c1 = dt  
        detA = 1.d0 / (4.d0 + b1)
        sqrtDetA = dsqrt(detA) 
        CayleyPoly(j,0,k) =  (4.d0 - b1)*detA
        CayleyPoly(j,1,k) =  (4.d0*c1)*detA
        CayleyPoly(j,2,k) =  -4.d0*a1*detA
        CayleyPoly(j,3,k) =  (4.d0 - b1)*detA
        sqrtCayleyPoly(j,0,k) =  2.d0*sqrtDetA
        sqrtCayleyPoly(j,1,k) =  c1*sqrtDetA 
        sqrtCayleyPoly(j,2,k) =  -a1*sqrtDetA 
        sqrtCayleyPoly(j,3,k) =  2.d0*sqrtDetA
        CayleyPoly(j,1,k) =  CayleyPoly(j,1,k)/mass(j)
        CayleyPoly(j,2,k) =  CayleyPoly(j,2,k)*mass(j)
        sqrtCayleyPoly(j,1,k) =  sqrtCayleyPoly(j,1,k)/mass(j)
        sqrtCayleyPoly(j,2,k) =  sqrtCayleyPoly(j,2,k)*mass(j)
   end do
end do

end subroutine
!****************************************************************************************************#  

