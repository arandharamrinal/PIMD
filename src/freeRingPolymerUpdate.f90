!************************************************************************************************!
subroutine freeRPFullStep()
    use allvars, only: Ndim,Natoms,Nbeads,poly,u,up
    !updates the normal mode coordinate and momenta of free ringpolymer i.e evolution is due to only 
    !the harmonic forces between the beads 
    implicit none
    real(8) ::  pNew
    integer :: i, j, k
    ! Transform cartesian coordinates and momenta to normal mode space
    call CNtrans()
    call CNtransP()
    do k = 0, Nbeads-1
       do j = 0, Natoms-1
            do i = 0,Ndim-1
                pNew = up(3*j+i,k) * poly(j,0,k) + u(3*j+i,k) * poly(j,1,k)
                u(3*j+i,k) = up(3*j+i,k) * poly(j,2,k) + u(3*j+i,k) * poly(j,3,k)
                up(3*j+i,k) = pNew
            end do
        end do
    end do
    call NCtrans()
    call NCtransP()
end subroutine freeRPFullStep
!Matrix for free ring polymer dynamics 
!************************************************************************************************!
subroutine freeRPHalfStep()
    use allvars, only: Ndim,Natoms,Nbeads,polyHalf,u,up
    !updates the normal mode coordinate and momenta of free ringpolymer i.e evolution is due to only 
    !the harmonic forces between the beads 
    implicit none
    real(8) ::  pNew
    integer :: i, j, k
    ! Transform cartesian coordinates and momenta to normal mode space
    call CNtrans()
    call CNtransP()
    do k = 0, Nbeads-1
       do j = 0, Natoms-1
            do i = 0,Ndim-1
                pNew = up(3*j+i,k) * polyHalf(j,0,k) + u(3*j+i,k) * polyHalf(j,1,k)
                u(3*j+i,k) = up(3*j+i,k) * polyHalf(j,2,k) + u(3*j+i,k) * polyHalf(j,3,k)
                up(3*j+i,k) = pNew
            end do
        end do
    end do
    call NCtrans()
    call NCtransP()
end subroutine freeRPHalfStep
!Matrix for free ring polymer dynamics 
!************************************************************************************************!
subroutine freeRPCayleyFullStep()
    use allvars, only: Ndim,Natoms,Nbeads,CayleyPoly,u,up
    !updates the normal mode coordinate and momenta of free ringpolymer i.e evolution is due to only 
    !the harmonic forces between the beads 
    implicit none
    real(8) ::  pNew
    integer :: i, j, k
    ! Transform cartesian coordinates and momenta to normal mode space
    call CNtrans()
    call CNtransP()
    do k = 0, Nbeads-1
       do j = 0, Natoms-1
            do i = 0,Ndim-1
                !print*,j,k, CayleyPoly(j,:,k);read(*,*)
                pNew       =  u(3*j+i,k) * CayleyPoly(j,2,k) + up(3*j+i,k) * CayleyPoly(j,3,k)
                u(3*j+i,k)  =  u(3*j+i,k) * CayleyPoly(j,0,k) + up(3*j+i,k) * CayleyPoly(j,1,k)
                up(3*j+i,k) = pNew
            end do
        end do
    end do
    call NCtrans()
    call NCtransP()
end subroutine freeRPCayleyFullStep
!Matrix for free ring polymer dynamics 
!************************************************************************************************!
subroutine freeRPCayleyHalfStep()
    use allvars, only: Ndim,Natoms,Nbeads,sqrtCayleyPoly,u,up
    !updates the normal mode coordinate and momenta of free ringpolymer i.e evolution is due to only 
    !the harmonic forces between the beads 
    implicit none
    real(8) ::  pNew
    integer :: i, j, k
    ! Transform cartesian coordinates and momenta to normal mode space
    call CNtrans()
    call CNtransP()
    do k = 0, Nbeads-1
       do j = 0, Natoms-1
            do i = 0,Ndim-1
                pNew       =  u(3*j+i,k) * sqrtCayleyPoly(j,2,k) + up(3*j+i,k) * sqrtCayleyPoly(j,3,k)
                u(3*j+i,k)  =  u(3*j+i,k) * sqrtCayleyPoly(j,0,k) + up(3*j+i,k) * sqrtCayleyPoly(j,1,k)
                up(3*j+i,k) = pNew
            end do
        end do
    end do
    call NCtrans()
    call NCtransP()
end subroutine freeRPCayleyHalfStep
!Matrix for free ring polymer dynamics 
!************************************************************************************************!
