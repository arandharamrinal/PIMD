!Global translation and rotation removal 
!For reference check 10.1063/1.3125009 by Dominik Marx 
!Everything is in the Normal Mode represantation. 
!Following routines are for correcting the centroid velocity(Angular  momentum) and  centroid force(Torque).
!It involves the following steps:
!        Transform to the center of mass frame of the whole system
!        For Velocity correction
!              1) Computes the center of mass velocity(qDot) and removes the same from the centroid velocity of each atom
!              2) Computes the centroid angular velocity (by the relation (inverse of MoI matrix multiplied by  centroid angular 
!                 momentum)
!              3) Computes the linear velocity due to rotation(qDotDoublePrime) by taking the cross product of centroid velocity 
!                 and centroid position
!              4) Subtract qDotDoublePrime from the centroid velocity of each atom
!        For force correction:
!             1) Compute the centroid force of the whole system and subtract it from the centroid force on each atom
!              2) Compute the time derivative of centroid angular velocity (by the relation (inverse of MoI matrix multiplied by  centroid torque))
!              3) Compute Force correction due to rotation which is forceDoublePrime and subtract this from the centroid force on each atom  
!************************************************************************************************!
subroutine getOuterProd(a,b,Ndim)
        !Computes the outer product of two vectors. Here both are a. 
        implicit none
        integer :: i,j,k
		integer,intent(in) ::Ndim
        real(8),intent(in) :: a(0:Ndim-1)
        real(8),intent(out) :: b(0:Ndim-1,0:Ndim-1)
        b = 0.d0
        do i = 0 , Ndim-1
           do j = 0 , Ndim-1
              b(i,j) = a(i)*a(j)
           enddo
        enddo
endsubroutine
!************************************************************************************************!
subroutine  removeCoM()
        use allvars, only : Natoms,u,CoM
        !Shift CoM to the origin
        implicit none 
        integer :: j,k
        call getCoM()
        !do k = 0,Nbeads-1
        !    do j =0,Natoms-1
        !            u(3*j:3*j+2,k) = u(3*j:3*j+2,k) - CoM(:,k)
        !    enddo
        !enddo
        do j =0,Natoms-1
                u(3*j:3*j+2,0) = u(3*j:3*j+2,0) - CoM(:,0)
        enddo
end subroutine 
!************************************************************************************************!
subroutine  removeCoMCart()
        use allvars, only : Nbeads, Natoms, q, CoMCart
        !Shift CoM to the origin
        integer :: j,k
        call getCoMCart()
        do k = 0,Nbeads-1
            do j =0,Natoms-1
                    q(3*j:3*j+2,k) = q(3*j:3*j+2,k) - CoMCart(:)
            enddo
        enddo
end subroutine
!************************************************************************************************!
subroutine  removeCoMVel()
        use allvars, only: Nbeads,Natoms,Ndim,Ncarts,up,massKprime
        !Shift CoM to the origin
        implicit none
        integer :: i,j,k
        real(8) :: CoMv(0:Ndim-1),up_tmp(0:Ncarts-1,0:Nbeads-1),tmass=0.d0
        do k = 0, Nbeads-1
            CoMv = 0.d0
            tmass=0.d0
            do j = 0, Natoms-1
                do i = 0,Ndim-1
                    CoMv(i) = CoMv(i) +  up(3*j+i,k)
                end do
                tmass=tmass + massKprime(3*j,k) 
            end do
            CoMv= CoMv / tmass
            do j =0,Natoms-1
               up(3*j:3*j+2,k) = up(3*j:3*j+2,k) - CoMv(:) * massKprime(3*j,k)
            enddo
       enddo
end subroutine
!************************************************************************************************!
subroutine getqDotPrime()
        use allvars, only : qDotPrime, Natoms, PreToNonpreNM, up, totalMassInv
        !Computes centroid velocity of the whole system
        implicit none
        integer :: i,j,k
        !total mass = mass of all atoms in the system
        qDotPrime = 0.d0
        do j = 0,Natoms-1
           qDotPrime(:) = qDotPrime(:) + PreToNonpreNM * up(3*j:3*j+2,0)
        enddo
        qDotPrime = qDotPrime * totalMassInv
end subroutine
!************************************************************************************************!
subroutine getFcent()
        use allvars, only : dVdu, PreToNonpreNM, Fcent
        !Computes force on the centroid of each atom, which is just the force on the zeroth normal mode 
        implicit none
        integer :: i, j, k
        Fcent = 0.0d0
        Fcent(:) = PreToNonpreNM * dVdu(:,0)              

end subroutine
!************************************************************************************************!
subroutine getforcePrime()
        use allvars, only : forcePrime, Natoms, PreToNonpreNM, dVdu, NatomsInv
        !Computes force on the centroid of the whole system
        implicit none
        integer:: j,i
        forcePrime = 0.d0
        do j = 0,Natoms-1 
               forcePrime(:) = forcePrime(:) +  PreToNonpreNM * dVdu(3*j:3*j+2,0)
        enddo
        forcePrime = forcePrime*NatomsInv
end subroutine
!************************************************************************************************!
subroutine removeFcent()
        use allvars, only: Natoms, NatomsInv, PreToNonpreNM, dVdu, PreToNonpreCart
    !Computes total centroid force and removes it from the centroid force of each atom
       implicit none
       integer :: i, j, k
       real(8) :: Fcent_u(0:2)

       Fcent_u = 0.0d0
       do j = 0, Natoms-1
          Fcent_u(:) = Fcent_u(:) + PreToNonpreNM * dVdu(3*j:3*j+2,0)
       end do
       Fcent_u = Fcent_u * NatomsInv
       do j = 0,Natoms-1
          dVdu(3*j:3*j+2,0) = dVdu(3*j:3*j+2,0) - PreToNonpreCart * Fcent_u(:)
       end do
       call NCpotTrans()
end subroutine removeFcent 

!************************************************************************************************!
subroutine removeRotation()
        use allvars, only: Natoms, Ndim, up, Cartmass, MoI, mass, identityMat,  Jcent, & 
                 delqDot, centVel, qDotPrime, centroid, delqDot, eigMoI, JinPA, eigMoIinv, AngVelcent, &
                 qDotDoublePrime, PreToNonpreNM, PreToNonpreCart 
        implicit none
        real(8) :: c(0:Ndim-1),outProd(0:Ndim-1,0:Ndim-1)
        integer :: i, j ,k
        !set center of mass of the system to be the origin
        !call CNtrans()
        !call CNtransP()
        !call removeCoMCart()
        !call CNtrans()
        call removeCoM()
        call NCtrans()
        !call removeCoMVel()
        !compute moment of inertia of beads for each atom
        call getcentroid()
        call getCentP()
        call getcentVel()
        MoI = 0.d0
        do j = 0,Natoms-1
           c = centroid(3*j:3*j+2)
           call getOuterProd(c,Ndim)
           MoI = MoI + mass(j)*(dot_product(c,c)*identityMat - outProd)
        enddo
        !compute the translational corretions to normal mode velocity
        call getqDotPrime()
        !calculate the angular momentum about the CoM of the whole system
        Jcent  = 0.d0
        do j = 0,Natoms-1
           delqDot(:)  = centVel(3*j:3*j+2) - qDotPrime(:)
           Jcent(0) = Jcent(0) +  mass(j) * (  centroid(3*j+1) * delqDot(2) - centroid(3*j+2) * delqDot(1)   )
           Jcent(1) = Jcent(1) +  mass(j) * (  centroid(3*j+2) * delqDot(0) - centroid(3*j)   * delqDot(2)   )
           Jcent(2) = Jcent(2) +  mass(j) * (  centroid(3*j)   * delqDot(1) - centroid(3*j+1) * delqDot(0)   )
        enddo
        !compute inverse of moment of inertia matrix
        !Diagonalise the moment of inertia matrix toget the principle axis 
        call HOUSEDIAG(MoI,Ndim,Ndim,eigMoI)
        forall(i=0:2) eigMoIinv(i) = 1.d0 /eigMoI(i)
        !Transform angular momentum from lab frame to principle axis frame 
        JinPA(0) = MoI(0,0)*Jcent(0) + MoI(1,0)*Jcent(1) + MoI(2,0)*Jcent(2)
        JinPA(1) = MoI(0,1)*Jcent(0) + MoI(1,1)*Jcent(1) + MoI(2,1)*Jcent(2)
        JinPA(2) = MoI(0,2)*Jcent(0) + MoI(1,2)*Jcent(1) + MoI(2,2)*Jcent(2)
        !compute angular momentum / moment of inertia to get the centroid angular velocity
        if ( Natoms == 1 ) then
           JinPA(0) = 0.d0
           JinPA(1) = 0.d0
           JinPA(2) = 0.d0
        else if ( Natoms == 2 ) then
           JinPA(0) = 0.d0
           JinPA(1) = JinPA(1)*eigMoIinv(1)
           JinPA(2) = JinPA(2)*eigMoIinv(2)
        else
           JinPA(0) = JinPA(0)*eigMoIinv(0)
           JinPA(1) = JinPA(1)*eigMoIinv(1)
           JinPA(2) = JinPA(2)*eigMoIinv(2)
        end if                                 
        !Transform centroid angular velocities from PA frame to lab frame.
        AngVelcent(0) = MoI(0,0)*JinPA(0) + MoI(0,1)*JinPA(1) + MoI(0,2)*JinPA(2)
        AngVelcent(1) = MoI(1,0)*JinPA(0) + MoI(1,1)*JinPA(1) + MoI(1,2)*JinPA(2)
        AngVelcent(2) = MoI(2,0)*JinPA(0) + MoI(2,1)*JinPA(1) + MoI(2,2)*JinPA(2)
        !Compute the rotational correction to velocity
        qDotDoublePrime =0.d0
        do j = 0,Natoms-1
           qDotDoublePrime(3*j)   = (  AngVelcent(1) * centroid(3*j+2) - AngVelcent(2) * centroid(3*j+1)   )
           qDotDoublePrime(3*j+1) = (  AngVelcent(2) * centroid(3*j)   - AngVelcent(0) * centroid(3*j+2)   )
           qDotDoublePrime(3*j+2) = (  AngVelcent(0) * centroid(3*j+1) - AngVelcent(1) * centroid(3*j)     )
        enddo
        do j = 0,Natoms-1
           up(3*j:3*j+2,0) = up(3*j:3*j+2,0) - PreToNonpreCart * mass(j) * qDotPrime(:) - PreToNonpreCart * mass(j) * qDotDoublePrime(3*j:3*j+2)  
        enddo
endsubroutine 
!************************************************************************************************!
subroutine removeForce()
        use allvars, only: Natoms, Ndim, dVdu, MoI, mass, identityMat,  torqueCent, & 
                 delqDot, centVel, delFcent, centroid, delqDot, eigMoI, torquePA, eigMoIinv, AngACent, & 
                 forceDoublePrime, forcePrime,  PreToNonpreNM, PreToNonpreCart 
        implicit none
        real(8) :: b(0:Ndim-1)
        integer :: i, j ,k
        real(8) :: outProd(0:Ndim-1,0:Ndim-1)
        !set center of mass of the system to be the origin
        !call removeCoMCart()
        !call CNtrans()
        call removeCoM()
        call NCtrans()
        call getcentroid()
        call getCentP()
        !compute moment of inertia of beads for each atom
        MoI = 0.d0
        do j = 0,Natoms-1
           b = centroid(3*j:3*j+2)
           call getOuterProd(b,Ndim)
           MoI = MoI + mass(j)*(dot_product(b,b)*identityMat - outProd)
        enddo
        !compute the translational corrections to normal mode force 
        call getforcePrime()
        !compute torque about the CoM of the whole system
        torqueCent = 0.d0
        do j = 0,Natoms-1
          delFcent(:)  = PreToNonpreNM*dVdu(3*j:3*j+2,0) - forcePrime(:)
          torqueCent(0) = torqueCent(0) +   (  centroid(3*j+1) * delFcent(2) - centroid(3*j+2) * delFcent(1) )
          torqueCent(1) = torqueCent(1) +   (  centroid(3*j+2) * delFcent(0) - centroid(3*j) * delFcent(2) )
          torqueCent(2) = torqueCent(2) +   (  centroid(3*j)   * delFcent(1) - centroid(3*j+1)   * delFcent(0) )
        enddo
        !Diagonalise the moment of inertia matrix to get the principle axis 
        call HOUSEDIAG(MoI,Ndim,Ndim,eigMoI)
        forall(i=0:2) eigMoIinv(i) = 1.d0 /eigMoI(i)
        !Transform centroid torque from lab frame to principle axis frame 
        torquePA(0) = MoI(0,0)*torqueCent(0) + MoI(1,0)*torqueCent(1) + MoI(2,0)*torqueCent(2)
        torquePA(1) = MoI(0,1)*torqueCent(0) + MoI(1,1)*torqueCent(1) + MoI(2,1)*torqueCent(2)
        torquePA(2) = MoI(0,2)*torqueCent(0) + MoI(1,2)*torqueCent(1) + MoI(2,2)*torqueCent(2)
        !compute angular momentum / moment of inertia to get the centroid angular velocity
        if ( Natoms == 1 ) then
           torquePA(0) = 0.d0
           torquePA(1) = 0.d0
           torquePA(2) = 0.d0
        else if ( Natoms == 2 ) then
           torquePA(0) = 0.d0
           torquePA(1) = torquePA(1)*eigMoIinv(1)
           torquePA(2) = torquePA(2)*eigMoIinv(2)
        else
           torquePA(0) = torquePA(0)*eigMoIinv(0)
           torquePA(1) = torquePA(1)*eigMoIinv(1)
           torquePA(2) = torquePA(2)*eigMoIinv(2)
        end if                                 
        !Transform centroid angular velocities from PA frame to lab frame.
        AngACent(0) = MoI(0,0)*torquePA(0) + MoI(0,1)*torquePA(1) + MoI(0,2)*torquePA(2)
        AngACent(1) = MoI(1,0)*torquePA(0) + MoI(1,1)*torquePA(1) + MoI(1,2)*torquePA(2)
        AngACent(2) = MoI(2,0)*torquePA(0) + MoI(2,1)*torquePA(1) + MoI(2,2)*torquePA(2)
        !Compute rotation correction to force
        forceDoublePrime = 0.d0
        do j = 0,Natoms-1
          forceDoublePrime(3*j)   = mass(j) * (  AngACent(1) * centroid(3*j+2) - AngACent(2) * centroid(3*j+1)   )
          forceDoublePrime(3*j+1) = mass(j) * (  AngACent(2) * centroid(3*j)   - AngACent(0) * centroid(3*j+2)   )
          forceDoublePrime(3*j+2) = mass(j) * (  AngACent(0) * centroid(3*j+1) - AngACent(1) * centroid(3*j)     )
        enddo
        do j = 0,Natoms-1
           dVdu(3*j:3*j+2,0) = dVdu(3*j:3*j+2,0) - PreToNonpreCart * forcePrime(:) - PreToNonpreCart * forceDoublePrime(3*j:3*j+2)  
        enddo
endsubroutine 
!*****************************************************************************************************!
!************************************************************************************************!
subroutine removeForceClassical()
        use allvars, only: Natoms, Ndim, dVdq, MoI, mass, cartmass, identityMat,  torqueCent, & 
                 delqDot, centVel, delFcent,q, delqDot, eigMoI, torquePA, eigMoIinv, AngACent, & 
                 forceDoublePrime, forceDoublePrime, forcePrime,  PreToNonpreNM, PreToNonpreCart 
        implicit none
        real(8) :: b(0:Ndim-1)
        real(8) :: outProd(0:Ndim-1,0:Ndim-1)
        integer :: i, j ,k
        !set center of mass of the system to be the origin
        call removeCoMCart()
        !compute moment of inertia of beads for each atom
        MoI = 0.d0
        do j = 0,Natoms-1
           b = q(3*j:3*j+2,0)
           call getOuterProd(b,Ndim)
           MoI = MoI + mass(j)*(dot_product(b,b)*identityMat - outProd)
        enddo
        !compute the translational corrections to normal mode force 
        forcePrime = 0.d0
        do j = 0,Natoms-1 
               forcePrime(:) = forcePrime(:) +  dVdq(3*j:3*j+2,0)
        enddo
        !compute torque about the CoM of the whole system
        torqueCent = 0.d0
        do j = 0,Natoms-1
          delFcent(:)  = dVdq(3*j:3*j+2,0) - forcePrime(:)
          torqueCent(0) = torqueCent(0) +   (  q(3*j+1,0) * delFcent(2) - q(3*j+2,0) * delFcent(1) )
          torqueCent(1) = torqueCent(1) +   (  q(3*j+2,0) * delFcent(0) - q(3*j,0) * delFcent(2) )
          torqueCent(2) = torqueCent(2) +   (  q(3*j,0)   * delFcent(1) - q(3*j+1,0)   * delFcent(0) )
        enddo
        !Diagonalise the moment of inertia matrix to get the principle axis 
        call HOUSEDIAG(MoI,Ndim,Ndim,eigMoI)
        forall(i=0:2) eigMoIinv(i) = 1.d0 /eigMoI(i)
        !Transform centroid torque from lab frame to principle axis frame 
        torquePA(0) = MoI(0,0)*torqueCent(0) + MoI(1,0)*torqueCent(1) + MoI(2,0)*torqueCent(2)
        torquePA(1) = MoI(0,1)*torqueCent(0) + MoI(1,1)*torqueCent(1) + MoI(2,1)*torqueCent(2)
        torquePA(2) = MoI(0,2)*torqueCent(0) + MoI(1,2)*torqueCent(1) + MoI(2,2)*torqueCent(2)
        !compute angular momentum / moment of inertia to get the centroid angular velocity
        if ( Natoms == 1 ) then
           torquePA(0) = 0.d0
           torquePA(1) = 0.d0
           torquePA(2) = 0.d0
        else if ( Natoms == 2 ) then
           torquePA(0) = 0.d0
           torquePA(1) = torquePA(1)*eigMoIinv(1)
           torquePA(2) = torquePA(2)*eigMoIinv(2)
        else
           torquePA(0) = torquePA(0)*eigMoIinv(0)
           torquePA(1) = torquePA(1)*eigMoIinv(1)
           torquePA(2) = torquePA(2)*eigMoIinv(2)
        end if                                 
        !Transform centroid angular velocities from PA frame to lab frame.
        AngACent(0) = MoI(0,0)*torquePA(0) + MoI(0,1)*torquePA(1) + MoI(0,2)*torquePA(2)
        AngACent(1) = MoI(1,0)*torquePA(0) + MoI(1,1)*torquePA(1) + MoI(1,2)*torquePA(2)
        AngACent(2) = MoI(2,0)*torquePA(0) + MoI(2,1)*torquePA(1) + MoI(2,2)*torquePA(2)
        !Compute rotation correction to force
        forceDoublePrime = 0.d0
        do j = 0,Natoms-1
          forceDoublePrime(3*j)   = mass(j) * (  AngACent(1) * q(3*j+2,0) - AngACent(2) * q(3*j+1,0)   )
          forceDoublePrime(3*j+1) = mass(j) * (  AngACent(2) * q(3*j,0)   - AngACent(0) * q(3*j+2,0)   )
          forceDoublePrime(3*j+2) = mass(j) * (  AngACent(0) * q(3*j+1,0) - AngACent(1) * q(3*j,0)     )
        enddo
        do j = 0,Natoms-1
           dVdq(3*j:3*j+2,0) = dVdq(3*j:3*j+2,0) - forcePrime(:) -  forceDoublePrime(3*j:3*j+2)  
        enddo
endsubroutine 
!*****************************************************************************************************!
!************************************************************************************************!
subroutine removeRotationClassical()
        use allvars, only: Natoms, Ndim, p, MoI, mass,cartmass, identityMat,  Jcent, & 
                 delqDot, centVel, qDotPrime, delqDot, eigMoI,q, JinPA, eigMoIinv, AngVelcent, &
                 qDotDoublePrime, PreToNonpreNM, PreToNonpreCart,totalMassInv 
        implicit none
        real(8) :: c(0:Ndim-1)
        integer :: i, j ,k
        real(8) :: outProd(0:Ndim-1,0:Ndim-1)
        !set center of mass of the system to be the origin
        call removeCoMCart()
        !compute moment of inertia of beads for each atom
        MoI = 0.d0
        do j = 0,Natoms-1
           c = q(3*j:3*j+2,0)
           call getOuterProd(c,Ndim)
           MoI = MoI + mass(j)*(dot_product(c,c)*identityMat - outProd)
        enddo
        !compute the translational corretions to normal mode velocity
        qDotPrime = 0.d0
        do j = 0,Natoms-1
           qDotPrime(:) = qDotPrime(:) + p(3*j:3*j+2,0)
        enddo
        qDotPrime = qDotPrime * totalMassInv
        !calculate the angular momentum about the CoM of the whole system
        Jcent  = 0.d0
        centVel = p(:,0) / Cartmass
        do j = 0,Natoms-1
           delqDot(:)  = centVel(3*j:3*j+2) - qDotPrime(:)
           Jcent(0) = Jcent(0) +  mass(j) * (  q(3*j+1,0) * delqDot(2) - q(3*j+2,0) * delqDot(1)   )
           Jcent(1) = Jcent(1) +  mass(j) * (  q(3*j+2,0) * delqDot(0) - q(3*j,0)   * delqDot(2)   )
           Jcent(2) = Jcent(2) +  mass(j) * (  q(3*j,0)   * delqDot(1) - q(3*j+1,0) * delqDot(0)   )
        enddo
        !compute inverse of moment of inertia matrix
        !Diagonalise the moment of inertia matrix toget the principle axis 
        call HOUSEDIAG(MoI,Ndim,Ndim,eigMoI)
        forall(i=0:2) eigMoIinv(i) = 1.d0 /eigMoI(i)
        !Transform angular momentum from lab frame to principle axis frame 
        JinPA(0) = MoI(0,0)*Jcent(0) + MoI(1,0)*Jcent(1) + MoI(2,0)*Jcent(2)
        JinPA(1) = MoI(0,1)*Jcent(0) + MoI(1,1)*Jcent(1) + MoI(2,1)*Jcent(2)
        JinPA(2) = MoI(0,2)*Jcent(0) + MoI(1,2)*Jcent(1) + MoI(2,2)*Jcent(2)
        !compute angular momentum / moment of inertia to get the centroid angular velocity
        if ( Natoms == 1 ) then
           JinPA(0) = 0.d0
           JinPA(1) = 0.d0
           JinPA(2) = 0.d0
        else if ( Natoms == 2 ) then
           JinPA(0) = 0.d0
           JinPA(1) = JinPA(1)*eigMoIinv(1)
           JinPA(2) = JinPA(2)*eigMoIinv(2)
        else
           JinPA(0) = JinPA(0)*eigMoIinv(0)
           JinPA(1) = JinPA(1)*eigMoIinv(1)
           JinPA(2) = JinPA(2)*eigMoIinv(2)
        end if                                 
        !Transform centroid angular velocities from PA frame to lab frame.
        AngVelcent(0) = MoI(0,0)*JinPA(0) + MoI(0,1)*JinPA(1) + MoI(0,2)*JinPA(2)
        AngVelcent(1) = MoI(1,0)*JinPA(0) + MoI(1,1)*JinPA(1) + MoI(1,2)*JinPA(2)
        AngVelcent(2) = MoI(2,0)*JinPA(0) + MoI(2,1)*JinPA(1) + MoI(2,2)*JinPA(2)
        !Compute the rotational correction to velocity
        qDotDoublePrime =0.d0
        do j = 0,Natoms-1
           qDotDoublePrime(3*j)   = (  AngVelcent(1) * q(3*j+2,0) - AngVelcent(2) * q(3*j+1,0)   )
           qDotDoublePrime(3*j+1) = (  AngVelcent(2) * q(3*j,0)   - AngVelcent(0) * q(3*j+2,0)   )
           qDotDoublePrime(3*j+2) = (  AngVelcent(0) * q(3*j+1,0) - AngVelcent(1) * q(3*j,0)     )
        enddo
        do j = 0,Natoms-1
           p(3*j:3*j+2,0) = p(3*j:3*j+2,0) -  mass(j) * qDotPrime(:) -  mass(j) * qDotDoublePrime(3*j:3*j+2)  
        enddo
endsubroutine 
