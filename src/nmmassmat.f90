subroutine getMass()
    use allvars, only: Natoms,Ndim, Nbeads, method, massK, massKprime, massKprimeInv, lamda, & 
                   pimdTorpmdMass, mass, massInv, sysNMmass, totalMass, totalMassInv
    implicit none
    integer:: jdim,iatom,ibead
    sysNMmass = 0.d0
    totalMass  = 0.d0
    !print*,mass
    do iatom = 0,Natoms-1
       do ibead = 0,Nbeads-1
          do jdim = 0,Ndim-1
                if ((method==0).or.(((method==1).or.(method==3)).or.(method==4))) then
                    massK(Ndim*iatom+jdim,ibead)=mass(iatom)
                    massKprime(Ndim*iatom+jdim,ibead)=mass(iatom)
                else if (method == 2) then
                    massK(Ndim*iatom+jdim,ibead) = lamda(ibead) * mass(iatom)
                    if (ibead == 0) then
                         massKprime(Ndim*iatom+jdim,ibead) = mass(iatom)
                    else
                         massKprime(Ndim*iatom+jdim,ibead) = massK(Ndim*iatom+jdim,ibead)
                    endif
                endif
                massKprimeInv(Ndim*iatom+jdim,ibead) =  1.d0 / massKprime(Ndim*iatom+jdim,ibead)
                sysNMmass = sysNMmass + massKprime(Ndim*iatom,ibead)
                if (ibead == 0) then
                         pimdTorpmdMass(Ndim*iatom+jdim,ibead) = 1.d0 
                else
                         pimdTorpmdMass(Ndim*iatom+jdim,ibead) = 1.d0/dsqrt(lamda(ibead))
                endif
          enddo
       enddo
       totalMass = totalMass + massKprime(3*iatom,0)
    enddo
    totalMassInv = 1.d0/totalMass
end subroutine
