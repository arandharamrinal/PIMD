!************************************************************************************************!
!Properties
!************************************************************************************************!
subroutine getCentroid()
        use allvars, only: centroid, PreToNonpreNM,u
        !Computes centroid postion of beads for each atom.
        implicit none
        integer :: j
        centroid = 0.0d0
        centroid(:) = PreToNonpreNM * u(:,0) 
end subroutine

!************************************************************************************************!
subroutine getCentP()
        use allvars,only: centP, PreToNonpreNM, up
        !Computes centroid momentum of beads for each atom. 
        implicit none
        integer :: j
        centP = 0.0d0
        centP(:) = PreToNonpreNM * up(:,0)
end subroutine
!************************************************************************************************!
subroutine getCentVel()
    use allvars, only : Ncarts, centVel, PreToNonpreNM, up, massKprimeInv
    !Computes centroid velocity of beads for each atom.
    implicit none
    integer ::  j
    centVel = 0.0d0
    do j = 0, Ncarts-1
        centVel(j) = PreToNonpreNM * up(j,0) * massKprimeInv(j,0)
    end do
    end subroutine
!************************************************************************************************!
!*************************  Compute center of mass in the cartesian  space  *********************!
!************************************************************************************************!
!Which mass ?
subroutine getCoM()
    use allvars, only : Nbeads, Natoms, Ndim, u, massK, CoM
    !Computes center of mass in the cartesian space.
    implicit none 
    integer :: i, j, k
    real(8) :: tmass = 0.d0
    CoM = 0.d0
    do k = 0, Nbeads-1
        tmass = 0.d0
        do j = 0, Natoms-1
            do i = 0,Ndim-1
                CoM(i,k) = CoM(i,k) + u(3*j+i,k) * massK(3*j,k)
            end do
            tmass = tmass + massK(3*j,k)
        end do
        !print*,tot_mass
        !CoM(:,k) = CoM(:,k) * totalMassInv
        CoM(:,k) = CoM(:,k) / tmass
    end do
end subroutine getCoM
!************************************************************************************************!
!Which mass ?
subroutine getCoMCart()
    use allvars, only : Nbeads, Natoms, Ndim, CoMCart, q, mass, CoMCart, totalMassInv, NbeadsInv
    !Computes center of mass in the cartesian space.
    implicit none
    integer :: i, j, k
    CoMCart = 0.d0
    do k = 0, Nbeads-1
        do j = 0, Natoms-1
            do i = 0,Ndim-1
                CoMCart(i) = CoMCart(i) + q(3*j+i,k) * mass(j)
            end do
        end do
    end do
    CoMCart(:) = CoMCart(:) * totalMassInv * NbeadsInv
end subroutine getCoMCart
!************************************************************************************************!
subroutine getCoMVel()
use allvars, only: Nbeads, Natoms, Ndim, CoMVel, up, totalMassInv
!Compute velocity of CoM.
implicit none
!real(8) :: sys_mass
integer :: i, j, k
!sys_mass = total_mass * dble(Nbeads)
CoMVel = 0.d0
do k = 0, Nbeads-1
do j = 0, Natoms-1
	do i = 0,Ndim-1
		CoMVel(i,k) = CoMVel(i,k) + up(3*j+i,k)
	end do
end do
CoMVel(:,k) = CoMVel(:,k) * totalMassInv
    end do
end subroutine getCoMVel
!************************************************************************************************!
!************************************************************************************************!
subroutine getCoMP()
    use allvars, only : Nbeads, Natoms, Ndim, CoMP, massKprime, up, totalMassInv
    !Compute momentum of CoM.
    implicit none
    integer :: i, j, k
    CoMP = 0.d0
    do k = 0, Nbeads-1
        do j = 0, Natoms-1
            do i = 0,Ndim-1
                CoMP(i,k) = CoMP(i,k) + massKprime(3*j+i,k) * up(3*j+i,k)
            end do
        end do
        CoMP(:,k) = CoMP(:,k) * totalMassInv
    end do
end subroutine getCoMP
!************************************************************************************************!
subroutine getGyrationRadii()
     use allvars, only : Natoms, Ndim, Nbeads, NbeadsInv, q, centroid, Rg
    !Computes the radius of gyration for each atom
    implicit none
    real(8) ::  rDist(0:Ndim-1)
    integer :: i, j, k

    call getCentroid()
    Rg = 0.0d0

    do j = 0, Natoms-1
        do i = 0,Ndim-1
            do k = 0, Nbeads-1
                rDist = q(3*j:3*j+2,k) - centroid(3*j:3*j+2)
                Rg(j) = Rg(j) + dot_product(rDist,rDist)
            end do
        end do
        Rg(j) = sqrt(Rg(j) * NbeadsInv)
    end do

end subroutine getGyrationRadii
!************************************************************************************************!
